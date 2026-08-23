#!/usr/bin/env python3
"""
PPPR parameter sweep driver.

Runs a matrix of OpenFAST cases with different PPPR settings, writing each case
into its own self-contained directory under SweepRuns/cases/ and collecting
per-case metrics into SweepRuns/results.csv.

Usage
-----
    conda activate openfast_env
    python run_pppr_sweep.py                 # run everything
    python run_pppr_sweep.py --list          # show the matrix, run nothing
    python run_pppr_sweep.py --only A_,C_    # run only ids containing A_ or C_
    python run_pppr_sweep.py --analyze-only  # recompute metrics from existing .out files
    python run_pppr_sweep.py --results r.csv # write to SweepRuns/r.csv instead
    python run_pppr_sweep.py -j 14           # 14 OpenFAST instances at once
    python run_pppr_sweep.py -j 0            # all cores minus two

A partial run needs the baseline in the filter (--only BASE,A_), because the
baseline is what measures the equilibrium tilt the PPPR cases offset from and
what the delta-power column is referenced against. Without it the script falls
back to --default-tilt and leaves pwr_vs_base_pct empty.

Design notes
------------
* Baselines (PPPR off, stock ROSCO) run FIRST. They serve two purposes: they
  give the power reference for the delta-power column, and they measure the
  equilibrium platform tilt at each wind speed, which the PPPR cases need for
  PPPR_offset_phi. You do not have to hand-measure the tilt.

* The reference phi is a DEVIATION from the equilibrium tilt (the linearised
  model has phi = 0 at equilibrium), while ROSCO's PtfmRDY is absolute. Cases
  therefore specify `offset_dev_deg` and the script adds the measured tilt.

* Gains follow the reference gain formula and depend on BOTH wind speed and forcing
  frequency, so they are recomputed per case. Set `gains=` explicitly on a case
  to override and hold them fixed while sweeping frequency (useful for
  separating "resonance is the problem" from "the gain formula overshoots").

* -j/--jobs runs several cases concurrently. OpenFAST is single-threaded here
  (no OpenMP) and one instance peaks near 50 MB, so cores are the only limit --
  memory is not. Each case owns its directory and each run is a separate
  process, so ROSCO's module-level state cannot leak between concurrent runs.
  Baselines always complete before any PPPR case starts, because they measure
  the equilibrium tilt that gets baked into PPPR_offset_phi.

* Only diagnostics/performance channels are written. Per-blade-node output is
  disabled, which is what drops the .out files from ~350 MB to a few MB.
"""

import argparse
import concurrent.futures as futures
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.join(HERE, "SweepRuns")
CASEDIR = os.path.join(SWEEP, "cases")
SIBLING = os.path.abspath(os.path.join(HERE, "..", "IEA-15-240-RWT"))
RESULTS = os.path.join(SWEEP, "results.csv")

FST = "IEA-15-240-RWT-UMaineSemi.fst"
DISCON = "IEA-15-240-RWT-UMaineSemi_DISCON.IN"
ED = "IEA-15-240-RWT-UMaineSemi_ElastoDyn.dat"
AD = "IEA-15-240-RWT-UMaineSemi_AeroDyn15.dat"
SD = "IEA-15-240-RWT-UMaineSemi_ServoDyn.dat"
SEA = "IEA-15-240-RWT-UMaineSemi_SeaState.dat"
INFLOW = "IEA-15-240-RWT_InflowFile.dat"

# Case-directory files copied verbatim from HERE. Blade/airfoil data stays in
# the sibling folder and is reached through a symlink (see ensure_layout).
COPY_FILES = [
    FST, DISCON, ED, AD, SD, SEA,
    "IEA-15-240-RWT-UMaineSemi_ElastoDyn_tower.dat",
    "IEA-15-240-RWT-UMaineSemi_HydroDyn.dat",
    "IEA-15-240-RWT-UMaineSemi_MoorDyn.dat",
    "IEA-15-240-RWT-UMaineSemi_MAP.dat",
]

# ----------------------------------------------------------------------------
# Output channels. Trimmed to diagnostics + performance + loads.
# ----------------------------------------------------------------------------
ED_CHANNELS = [
    "BldPitch1", "GenSpeed", "RotSpeed",
    "PtfmPitch", "PtfmSurge", "PtfmHeave",
    "RotThrust", "RotTorq",
    "TwrBsMyt", "TwrBsFxt", "RootMyb1",
]
AD_CHANNELS = ["RtFldCp", "RtFldCt", "RtTSR", "RtVAvgxh", "RtSpeed"]
SD_CHANNELS = ["GenPwr", "GenTq"]

# ----------------------------------------------------------------------------
# Turbine constants and the aerodynamic sensitivity surface.
#
# calculateGains.m interpolates dCp/dTSR, dCp/dbeta and Cp on the IEA-15MW
# steady surface at the case's (TSR0, beta0) rather than assuming a fixed
# operating point. Doing the same here is what lets TSR and wind speed become
# swept axes instead of baked-in constants. Verified to reproduce
# calculateGains() to six digits at its defaults: kp = 0.463407,
# kp_Tg = 2.63639e7 at U = 10, freq = 0.1, zeta = 0.7.
# ----------------------------------------------------------------------------
RHO = 1.2
RROT = 120.0
JR = 3.525e8
NG = 1.0
CPFILE = os.path.join(HERE, "MatlabSimulations", "IEA15MW_Cp_Ct_Cq.mat")

# Fallback values at TSR0 = 8.5, beta0 = 0, used only if the .mat or scipy is
# unavailable. Valid at that operating point ONLY.
TSR0 = 8.5
CP0 = 0.46969
DCP_DTSR = 0.005270
DCP_DBETA = -0.004124          # per DEGREE

_CP = {}


def _cp_sensitivities(tsr0, beta0):
    """(dCp_dTSR, dCp_dbeta_per_deg, Cp0) at an operating point."""
    if "interp" not in _CP:
        try:
            import numpy as np
            from scipy.io import loadmat
            from scipy.interpolate import RegularGridInterpolator
            d = loadmat(CPFILE)
            tsrs, angs, cp = d["TSRs"].ravel(), d["angles"].ravel(), d["Cp"]
            # MATLAB gradient(Cp, angles, TSRs) -> d/d(angles) on dim 2,
            # d/d(TSRs) on dim 1. numpy returns them in axis order.
            g_tsr, g_beta = np.gradient(cp, tsrs, angs)
            mk = lambda Z: RegularGridInterpolator((tsrs, angs), Z, method="linear",
                                                   bounds_error=False, fill_value=None)
            _CP["interp"] = (mk(g_tsr), mk(g_beta), mk(cp))
        except Exception as e:                      # noqa: BLE001
            _CP["interp"] = None
            _CP["why"] = str(e)
    it = _CP["interp"]
    if it is None:
        if abs(tsr0 - TSR0) > 1e-9 or abs(beta0) > 1e-9:
            print("  WARNING: Cp surface unavailable ({}); gains at TSR0={} beta0={} "
                  "fall back to constants valid only at TSR0=8.5, beta0=0."
                  .format(_CP.get("why", "?"), tsr0, beta0))
        return DCP_DTSR, DCP_DBETA, CP0
    pt = [[tsr0, beta0]]
    return (float(it[0](pt)[0]), float(it[1](pt)[0]), float(it[2](pt)[0]))


def pir_gains(U, omega, zeta=0.7, ratio=10.0, tsr0=TSR0, beta0=0.0):
    """PIR gains at an operating point -- the calculateGains.m formulae.

    Returns ROSCO-native units: kp/kr in rad-pitch per rad-error, kp_Tg/kr_Tg
    in N-m per rad/s.
    """
    dCp_dTSR, dCp_dbeta_deg, Cp0 = _cp_sensitivities(tsr0, beta0)
    dCp_dbeta = dCp_dbeta_deg * 180.0 / math.pi          # per-deg -> per-rad
    A = 0.5 * RHO * math.pi * RROT**4 * U / (JR * tsr0**2) * (dCp_dTSR * tsr0 - Cp0)
    B_beta = NG / (2 * JR * tsr0**2) * RHO * math.pi * RROT**3 * U**2 * (dCp_dbeta * tsr0)
    B_Tg = -NG**2 / JR
    root = 2 * zeta * omega + A
    kp = -1.0 / (2 * math.pi * B_beta) * root
    kp_tg = -1.0 / B_Tg * root
    return dict(kp=kp, kr=kp / ratio, kp_tg=kp_tg, kr_tg=kp_tg / ratio)


# ============================================================================
# CASE MATRIX -- edit here
# ============================================================================
DEFAULTS = dict(
    U=10.0,
    tsr0=8.5,               # operating tip-speed ratio
    beta0=0.0,              # operating blade pitch [deg]
    omega=0.20,             # forcing frequency [rad/s]
    amp_phi_deg=1.0,
    amp_omega=0.01,         # [rad/s]
    offset_dev_deg=0.0,     # phi reference mean, as DEVIATION from equilibrium tilt
    phi_phase_deg=0.0,
    omega_phase_deg=-90.0,  # reference is sin(wt - 90deg); ROSCO adds the offset
    zeta=0.7,
    ratio=10.0,             # kp/kr
    gains=None,             # None -> compute; or dict(kp=,kr=,kp_tg=,kr_tg=)
    waves=False,            # False -> WaveMod 0 (calm), isolates the controller
    tmax=600.0,
    dt=0.025,
    dt_out=0.05,
    pppr=True,
    # Blade-pitch floor. This is the single most consequential setting in the file:
    # the reference tree shipped PC_MinPit = -1.57 rad / PC_FinePit = -0.35 rad,
    # which lets the controller command NEGATIVE pitch, while this tree shipped 0/0.
    # PC_FinePit is the operative floor (PS_Mode = 0 sets LocalVar%PC_MinPit from
    # it); PC_MinPit is only the hardware backstop applied to PitComAct. Holding
    # the platform near equilibrium tilt needs ~0 deg of steady pitch, so with a
    # 0 floor the command clips and the controller diverges -- which is why this
    # tree needed offset_dev_deg < 0 and the reference tree did not. Defaults now
    # track the merged DISCON.IN; sweep D varies them deliberately. Leaving these
    # at 0 would have silently re-imposed the old floor on every case and
    # reproduced the old divergences.
    pc_minpit_rad=-1.57,
    pc_finepit_rad=-0.35,
)


# Operating-point deviation used by every sweep that is NOT deliberately varying
# it. The first campaign ran A and C at dev = 0, which is unusable: holding the
# platform at its *equilibrium* tilt needs ~0 deg of steady blade pitch, the
# turbine already sits at 0.5-1.25 deg, so the command clips on the PC_MinPit = 0
# floor, the clipped value is fed back into the controller's own stored state,
# and it diverges. Every dev >= 0 case in results_dev0.csv failed and every
# dev < 0 case survived, with nothing else in the matrix predicting the outcome.
# -2 deg is the deepest deviation that still holds a workable amount of power
# (-20.7% vs -59.1% at -4 deg), so it is the operating point to sweep around.
DEV = -2.0


def build_cases():
    cases = []

    def add(cid, group, **kw):
        c = dict(DEFAULTS)
        c.update(kw)
        c["id"] = cid
        c["group"] = group
        cases.append(c)

    # --- Baselines: PPPR off, stock ROSCO. Run first; measure equilibrium tilt.
    for U in (10.0,):
        add(f"BASE_U{U:g}", "baseline", U=U, pppr=False)

    # --- Sweep A: amplitude (AC drive strength) about the DEV operating point.
    #     Expect a ceiling somewhere near 4 deg: at dev = -2 the baseline holds
    #     ~5.5 deg of pitch, and the observed exchange rate is ~1.2 deg of pitch
    #     swing per deg of platform amplitude, so amp ~4.5 puts the trough back
    #     on the 0 deg floor. That ceiling is now a real amplitude limit rather
    #     than the dev = 0 artefact.
    for a in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0):
        add(f"A_amp{a:g}_d{DEV:+g}", "A_amplitude", amp_phi_deg=a,
            offset_dev_deg=DEV, omega=0.20)

    # --- Sweep B: DC offset. Directly tests the PIR integrator, which is the
    #     entire reason for the I term -- a pure PR cannot hold any of these.
    #     This is the one sweep that varies dev, so it does not take DEV.
    for d in (0.0, -1.0, -2.0, -4.0, 2.0):
        add(f"B_dev{d:+g}", "B_dc_offset", offset_dev_deg=d, amp_phi_deg=1.0, omega=0.20)

    # --- Sweep C: frequency. 0.213 rad/s is the UMaineSemi platform mode and the
    #     abstract's actual target. Gains move with frequency (see pir_gains).
    for w in (0.15, 0.18, 0.20, 0.213, 0.24, 0.28):
        add(f"C_w{w:g}_d{DEV:+g}", "C_frequency", omega=w, amp_phi_deg=1.0,
            offset_dev_deg=DEV)

    # --- Sweep C2: frequency with gains PINNED at the 0.20 rad/s values, so that
    #     frequency is the only thing moving. Only meaningful if C shows a break.
    g20 = pir_gains(10.0, 0.20)
    for w in (0.213, 0.24):
        add(f"C2_w{w:g}_fixedgain_d{DEV:+g}", "C2_frequency_fixedgain", omega=w,
            gains=g20, offset_dev_deg=DEV)

    # --- Sweep D: the pitch floor. This is the direct test of why the reference
    #     tree works near equilibrium tilt and this one does not. Each row pairs a
    #     floor with dev = 0, the operating point that diverged at t=123 s under a
    #     0 floor. If the -0.35 rad rows survive, the floor is the whole story.
    for name, mn, fn in (("0", 0.0, 0.0),            # this tree's shipped floor
                         ("fine-0.35", -1.57, -0.35),  # reference tree's shipped pair
                         ("fine-0.09", -1.57, -0.0873)):  # -5 deg, a defensible limit
        add(f"D_floor{name}_dev0", "D_pitch_floor", offset_dev_deg=0.0,
            amp_phi_deg=1.0, omega=0.20, pc_minpit_rad=mn, pc_finepit_rad=fn)
        add(f"D_floor{name}_dev-2", "D_pitch_floor", offset_dev_deg=-2.0,
            amp_phi_deg=1.0, omega=0.20, pc_minpit_rad=mn, pc_finepit_rad=fn)

    # --- Sweep E: the amplitude x frequency grid for the contour/regime map.
    #     Run at dev = 0 -- the equilibrium operating point -- which is viable now
    #     that the merged tree ships a negative pitch floor. 20 cases; at -j 14
    #     that is two batches.
    for w in (0.10, 0.15, 0.20, 0.24, 0.28):
        for a in (0.5, 1.0, 2.0, 3.0):
            add(f"E_w{w:g}_a{a:g}", "E_grid", omega=w, amp_phi_deg=a,
                offset_dev_deg=0.0)

    # de-duplicate ids (sweeps overlap at the nominal point)
    seen, out = set(), []
    for c in cases:
        if c["id"] not in seen:
            seen.add(c["id"])
            out.append(c)
    return out


# ============================================================================
# File patching helpers
# ============================================================================
def set_fast(lines, label, value, fname=""):
    """Set a value in OpenFAST '<value>  <Label>  - description' format."""
    # The label must be matched in the LABEL POSITION only -- i.e. the last token
    # before the ' - ' description separator. Matching anywhere on the line is
    # unsafe: e.g. the MSL2SWL line's description contains "unused when WaveMod
    # = 6", which a looser pattern happily clobbers.
    for i, ln in enumerate(lines):
        m = re.search(r"\s+-\s", ln)
        head = ln[:m.start()] if m else ln.rstrip("\n")
        tail = ln[m.start():] if m else "\n"
        toks = head.split()
        if len(toks) >= 2 and toks[-1] == label:
            indent = " " * (len(head) - len(head.lstrip()))
            val = str(value).rstrip()
            lines[i] = indent + val + " " * max(2, 22 - len(val)) + label + tail
            if not lines[i].endswith("\n"):
                lines[i] += "\n"
            return True
    print("  WARNING: '{}' not found in {}".format(label, fname))
    return False


def set_discon(lines, label, value, fname=DISCON):
    """Set a value in DISCON.IN '<value>   ! <Label>  - description' format."""
    pat = re.compile(r"^\s*.*?!\s*" + re.escape(label) + r"\b")
    for i, ln in enumerate(lines):
        if pat.match(ln) and "!" in ln:
            comment = ln[ln.index("!"):]
            # always keep >=2 spaces before '!' -- ROSCO tokenises on whitespace,
            # so a long value butting against the comment marker breaks parsing
            val = str(value)
            lines[i] = val + " " * max(2, 20 - len(val)) + comment
            return True
    print("  WARNING: '{}' not found in {}".format(label, fname))
    return False


def replace_outlist(lines, channels, occurrence=0, quote=True):
    """Replace the Nth OutList block (header .. END) with `channels`."""
    idx = [i for i, ln in enumerate(lines) if re.search(r"\bOutList\b", ln)]
    if occurrence >= len(idx):
        return False
    start = idx[occurrence]
    end = start + 1
    while end < len(lines) and not lines[end].lstrip().upper().startswith("END"):
        end += 1
    fmt = '"{}"' if quote else "{}"
    block = [fmt.format(c) + "\n" for c in channels]
    lines[start + 1:end] = block
    return True


def read(p):
    with open(p, "r", errors="replace") as f:
        return f.readlines()


def write(p, lines):
    with open(p, "w") as f:
        f.writelines(lines)


# ============================================================================
# Case setup
# ============================================================================
def ensure_layout():
    os.makedirs(CASEDIR, exist_ok=True)
    # Case dirs sit one level below HERE, so '../IEA-15-240-RWT/...' references
    # inside AeroDyn/ElastoDyn resolve through this symlink.
    link = os.path.join(CASEDIR, "IEA-15-240-RWT")
    if not os.path.exists(link):
        os.symlink(SIBLING, link)


def setup_case(c):
    d = os.path.join(CASEDIR, c["id"])
    os.makedirs(d, exist_ok=True)
    for fn in COPY_FILES:
        src = os.path.join(HERE, fn)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(d, fn))
    shutil.copy2(os.path.join(SIBLING, INFLOW), os.path.join(d, INFLOW))
    # HydroDyn reads potential-flow data via './HydroData/...' relative to the
    # case directory, so it has to be reachable from here. Symlink, don't copy.
    for sub in ("HydroData",):
        src, dst = os.path.join(HERE, sub), os.path.join(d, sub)
        if os.path.isdir(src) and not os.path.exists(dst):
            os.symlink(src, dst)

    freq_hz = c["omega"] / (2 * math.pi)
    g = c["gains"] or pir_gains(c["U"], c["omega"], c["zeta"], c["ratio"],
                                c["tsr0"], c["beta0"])
    tilt = c["tilt_deg"]
    # (offset is written in degrees now -- see the DISCON block below)
    omega_ref = c["tsr0"] * c["U"] / RROT

    # ---- .fst
    p = os.path.join(d, FST)
    L = read(p)
    set_fast(L, "TMax", "{:<22.1f}".format(c["tmax"]), FST)
    set_fast(L, "DT", "{:<22.4f}".format(c["dt"]), FST)
    set_fast(L, "DT_Out", "{:<22.4f}".format(c["dt_out"]), FST)
    set_fast(L, "InflowFile", '"{}"'.format(INFLOW), FST)
    write(p, L)

    # ---- InflowWind
    p = os.path.join(d, INFLOW)
    L = read(p)
    set_fast(L, "HWindSpeed", "{:<22.2f}".format(c["U"]), INFLOW)
    write(p, L)

    # ---- ElastoDyn: initial conditions matched to the references, trimmed output
    p = os.path.join(d, ED)
    L = read(p)
    set_fast(L, "RotSpeed", "{:<11.4f}".format(omega_ref * 60 / (2 * math.pi)), ED)
    set_fast(L, "PtfmPitch", "{:<11.4f}".format(tilt if c["pppr"] else 0.0), ED)
    set_fast(L, "BlPitch(1)", "{:<11.2f}".format(0.0), ED)
    set_fast(L, "BlPitch(2)", "{:<11.2f}".format(0.0), ED)
    set_fast(L, "BlPitch(3)", "{:<11.2f}".format(0.0), ED)
    set_fast(L, "BldNd_BladesOut", "{:<11d}".format(0), ED)   # kill per-node output
    replace_outlist(L, ED_CHANNELS, occurrence=0)
    write(p, L)

    # ---- AeroDyn: kill the per-node block (this is the file-size hog)
    p = os.path.join(d, AD)
    L = read(p)
    set_fast(L, "BldNd_BladesOut", "{:<15d}".format(0), AD)
    replace_outlist(L, AD_CHANNELS, occurrence=0)
    write(p, L)

    # ---- ServoDyn
    p = os.path.join(d, SD)
    L = read(p)
    replace_outlist(L, SD_CHANNELS, occurrence=0)
    write(p, L)

    # ---- SeaState
    p = os.path.join(d, SEA)
    L = read(p)
    set_fast(L, "WaveMod", "{:<10d}".format(2 if c["waves"] else 0), SEA)
    # Pin seed 1 so wave-on runs are comparable. Seed 2 is left alone -- it holds
    # "RANLUX" (an RNG selector, not a number), and seed 1 fully determines the
    # realisation in that mode.
    set_fast(L, "WaveSeed(1)", "{:<10d}".format(123456789), SEA)
    write(p, L)

    # ---- DISCON
    p = os.path.join(d, DISCON)
    L = read(p)
    if c["pppr"]:
        # Unit conventions below follow the merged (reference-tree) controller:
        #   PPPR_Mode        1, not 2 -- the additive/delayed mode no longer exists
        #   PPPR_amp_phi     DEGREES  (was radians)
        #   PPPR_offset_phi  DEGREES  (was radians)
        #   PPPR_freq_*      rad/s    (was Hz)
        #   PPPR_CntrGains_* 3 values [Kp, Kr, omega_z]; PPPR_fz_* no longer exist
        #   phase enters as sin(wt - phase*D2R), so the sign is negated here to
        #   keep a case's phi_phase_deg/omega_phase_deg meaning what it did before.
        omega_z = c["omega"] / 10.0                     # PI zero, 0.1x forcing freq
        set_discon(L, "PPPR_Mode", 1)
        set_discon(L, "PC_ControlMode", 0)
        set_discon(L, "VS_ControlMode", 0)
        set_discon(L, "SS_Mode", 0)
        set_discon(L, "PS_Mode", 0)
        set_discon(L, "Fl_Mode", 0)
        set_discon(L, "PC_MinPit", "{:.6f}".format(c["pc_minpit_rad"]))
        set_discon(L, "PC_FinePit", "{:.6f}".format(c["pc_finepit_rad"]))
        set_discon(L, "PPPR_amp_phi", "{:.6f}".format(c["amp_phi_deg"]))
        set_discon(L, "PPPR_freq_phi", "{:.6f}".format(c["omega"]))
        set_discon(L, "PPPR_amp_omega", "{:.6f}".format(c["amp_omega"]))
        set_discon(L, "PPPR_freq_omega", "{:.6f}".format(c["omega"]))
        set_discon(L, "Phi_phaseoffset", "{:.2f}".format(-c["phi_phase_deg"]))
        set_discon(L, "Omega_phaseoffset", "{:.2f}".format(-c["omega_phase_deg"]))
        set_discon(L, "PPPR_offset_phi", "{:.6f}".format(tilt + c["offset_dev_deg"]))
        set_discon(L, "PPPR_offset_omega", "{:.6f}".format(omega_ref))
        set_discon(L, "PPPR_CntrGains_phi",
                   "{:.6f}  {:.6f}  {:.6f}".format(g["kp"], g["kr"], omega_z))
        set_discon(L, "PPPR_CntrGains_omega",
                   "{:.5e}  {:.5e}  {:.6f}".format(g["kp_tg"], g["kr_tg"], omega_z))
    else:
        set_discon(L, "PPPR_Mode", 0)
        set_discon(L, "PC_ControlMode", 1)
        set_discon(L, "VS_ControlMode", 1)
        set_discon(L, "SS_Mode", 1)
        set_discon(L, "PS_Mode", 1)
        set_discon(L, "Fl_Mode", 0)
        # Pin the baseline to the stock physical floor regardless of what the
        # merged DISCON.IN ships, or the power reference moves with the PPPR
        # cases and every pwr_vs_base_pct becomes meaningless.
        set_discon(L, "PC_MinPit", "{:.6f}".format(0.0))
        set_discon(L, "PC_FinePit", "{:.6f}".format(0.0))
    write(p, L)

    c["_dir"] = d
    c["_gains"] = g
    c["_freq_hz"] = freq_hz
    c["_omega_ref"] = omega_ref
    c["_offset_phi_deg"] = tilt + c["offset_dev_deg"]
    return d


# ============================================================================
# Run + analyse
# ============================================================================
def run_case(c, timeout=7200):
    d = c["_dir"]
    log = os.path.join(d, "openfast.log")
    t0 = time.time()
    with open(log, "w") as lf:
        try:
            r = subprocess.run(["openfast", FST], cwd=d, stdout=lf,
                               stderr=subprocess.STDOUT, timeout=timeout)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = -1
            lf.write("\n*** TIMEOUT after {} s ***\n".format(timeout))
    return rc, time.time() - t0


def load_out(path):
    """Parse an OpenFAST .out file into {channel: [values]} using stdlib only."""
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()
    hdr = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("Time"):
            hdr = i
            break
    if hdr is None:
        return {}
    names = lines[hdr].split()
    data = {n: [] for n in names}
    for ln in lines[hdr + 2:]:            # +2 skips the units row
        p = ln.split()
        if len(p) != len(names):
            continue
        try:
            vals = [float(x) for x in p]
        except ValueError:
            continue
        for n, v in zip(names, vals):
            data[n].append(v)
    return data


def fit_sinusoid(t, x, f):
    """Single-harmonic projection: x ~ mean + A*sin(2*pi*f*t + psi).

    Returns (mean, A, psi_deg). psi uses the same convention as the reference
    waveform, which ROSCO builds as offset + amp*sin(w*t - Phi_phaseoffset*D2R);
    setup_case writes Phi_phaseoffset = -phi_phase_deg, so a case's
    phi_phase_deg IS the commanded psi and the two are directly comparable.
    """
    n = len(x)
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    m = sum(x) / n
    cc = ss = 0.0
    for ti, xi in zip(t, x):
        a = 2 * math.pi * f * ti
        cc += (xi - m) * math.cos(a)
        ss += (xi - m) * math.sin(a)
    return m, 2.0 * math.hypot(cc, ss) / n, math.degrees(math.atan2(cc, ss))


def amp_at(t, x, f):
    """Amplitude only -- kept for callers that do not need the phase."""
    return fit_sinusoid(t, x, f)[1]


def wrap180(a):
    """Fold an angle into (-180, 180] so phase errors do not read as ~360."""
    return (a + 180.0) % 360.0 - 180.0


def rmse(a, b):
    n = len(a)
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / n) if n else float("nan")


def phase_average(t, x, f, nbins=72):
    """Fold x onto a single forcing cycle.

    Returns (phase_deg, mean_in_bin, count_in_bin). This is the phase-averaged
    quantity: it suppresses wave and turbulence content that is incoherent with
    the forcing, leaving the cycle the controller is actually driving.
    """
    T = 1.0 / f
    acc = [0.0] * nbins
    cnt = [0] * nbins
    for ti, xi in zip(t, x):
        k = int(((ti % T) / T) * nbins) % nbins
        acc[k] += xi
        cnt[k] += 1
    return ([(k + 0.5) * 360.0 / nbins for k in range(nbins)],
            [acc[k] / cnt[k] if cnt[k] else float("nan") for k in range(nbins)],
            cnt)


def stat(d, key, sl):
    v = d.get(key)
    if not v:
        return None
    return v[sl]


def analyse(c, outpath):
    d = load_out(outpath)
    if not d or "Time" not in d or len(d["Time"]) < 10:
        return dict(completed=0, note="no/short output")
    t = d["Time"]
    tmax = t[-1]
    i0 = next((i for i, ti in enumerate(t) if ti >= tmax / 3.0), 0)   # drop first third
    sl = slice(i0, None)
    ts = t[sl]
    r = dict(completed=1 if tmax >= 0.98 * c["tmax"] else 0, t_end=round(tmax, 1))

    def col(k):
        return d[k][sl] if k in d else None

    bp = col("BldPitch1")
    if bp:
        n = len(bp)
        r["pitch_lo_pct"] = round(100.0 * sum(1 for v in bp if v < 0.01) / n, 1)
        r["pitch_hi_pct"] = round(100.0 * sum(1 for v in bp if v > 89.9) / n, 1)
        r["pitch_mean"] = round(sum(bp) / n, 3)
        r["pitch_min"] = round(min(bp), 3)
        r["pitch_max"] = round(max(bp), 3)

    pp = col("PtfmPitch")
    if pp:
        r["ptfm_mean"] = round(sum(pp) / len(pp), 4)
        if c["pppr"]:
            r["ptfm_ref_mean"] = round(c["_offset_phi_deg"], 4)
            r["ptfm_mean_err"] = round(r["ptfm_mean"] - c["_offset_phi_deg"], 4)
            _, amp, psi = fit_sinusoid(ts, pp, c["_freq_hz"])
            r["ptfm_amp"] = round(amp, 4)
            r["ptfm_ref_amp"] = round(c["amp_phi_deg"], 4)
            r["ptfm_phase"] = round(psi, 2)
            r["ptfm_phase_err"] = round(wrap180(psi - c["phi_phase_deg"]), 2)
            if c["amp_phi_deg"] > 0:
                r["ptfm_amp_ratio"] = round(amp / c["amp_phi_deg"], 3)
            w = 2 * math.pi * c["_freq_hz"]
            ref = [c["_offset_phi_deg"] + c["amp_phi_deg"]
                   * math.sin(w * ti + math.radians(c["phi_phase_deg"])) for ti in ts]
            r["ptfm_rmse"] = round(rmse(pp, ref), 4)

    gs = col("GenSpeed")
    if gs:
        r["gen_mean_rpm"] = round(sum(gs) / len(gs), 4)
        if c["pppr"]:
            ref = c["_omega_ref"] * 60 / (2 * math.pi)
            r["gen_ref_rpm"] = round(ref, 4)
            r["gen_mean_err"] = round(r["gen_mean_rpm"] - ref, 4)
            _, gamp, gpsi = fit_sinusoid(ts, gs, c["_freq_hz"])
            r["gen_amp"] = round(gamp, 5)
            r["gen_ref_amp"] = round(c["amp_omega"] * 60 / (2 * math.pi), 5)
            r["gen_phase"] = round(gpsi, 2)
            r["gen_phase_err"] = round(wrap180(gpsi - c["omega_phase_deg"]), 2)
            if c["amp_omega"] > 0:
                r["gen_amp_ratio"] = round(gamp / (c["amp_omega"] * 60 / (2 * math.pi)), 3)
            w = 2 * math.pi * c["_freq_hz"]
            gref = [ref + c["amp_omega"] * 60 / (2 * math.pi)
                    * math.sin(w * ti + math.radians(c["omega_phase_deg"])) for ti in ts]
            r["gen_rmse"] = round(rmse(gs, gref), 5)

    gt = col("GenTq")
    if gt and c["pppr"]:
        _, tamp, tpsi = fit_sinusoid(ts, gt, c["_freq_hz"])
        r["gentq_mean"] = round(sum(gt) / len(gt), 1)
        r["gentq_amp"] = round(tamp, 1)
        r["gentq_phase"] = round(tpsi, 2)

    pw = col("GenPwr")
    if pw:
        r["pwr_mean_kW"] = round(sum(pw) / len(pw), 1)
        r["pwr_min_kW"] = round(min(pw), 1)
        r["pwr_max_kW"] = round(max(pw), 1)

    ts_ = col("RtTSR")
    if ts_:
        r["tsr_min"] = round(min(ts_), 2)
        r["tsr_mean"] = round(sum(ts_) / len(ts_), 2)

    for k, lbl in (("TwrBsMyt", "twrbs_max"), ("RootMyb1", "rootmy_max"),
                   ("RotThrust", "thrust_mean")):
        v = col(k)
        if v:
            r[lbl] = round(max(abs(x) for x in v) if "max" in lbl
                           else sum(v) / len(v), 1)

    # Phase-averaged cycle, written next to the .out. Folding the settled window
    # onto one forcing period averages away wave and turbulence content that is
    # incoherent with the forcing, leaving the cycle the controller is driving.
    if c["pppr"] and c.get("_freq_hz", 0) > 0:
        chans = [(k, v) for k, v in (("BldPitch1", bp), ("PtfmPitch", pp),
                                     ("GenSpeed", gs), ("GenTq", gt)) if v]
        if chans:
            ph, cols = None, {}
            for k, v in chans:
                ph, cols[k], _ = phase_average(ts, v, c["_freq_hz"])
            with open(os.path.join(os.path.dirname(outpath), "phase_avg.csv"),
                      "w", newline="") as f:
                wr = csv.writer(f)
                wr.writerow(["phase_deg"] + list(cols))
                for j in range(len(ph)):
                    wr.writerow(["{:.2f}".format(ph[j])]
                                + ["{:.6g}".format(cols[k][j]) for k in cols])
    return r


def measure_tilt(outpath):
    """Mean PtfmPitch over the settled portion of a baseline run."""
    d = load_out(outpath)
    if not d or "PtfmPitch" not in d or not d["PtfmPitch"]:
        return None
    t, p = d["Time"], d["PtfmPitch"]
    i0 = next((i for i, ti in enumerate(t) if ti >= t[-1] / 3.0), 0)
    seg = p[i0:]
    return sum(seg) / len(seg)


# ============================================================================
def write_results(path, rows):
    """Rewrite the whole CSV. Called under `lock` whenever rows changes."""
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def process_case(c, args, tilts, rows, lock, results, idx, total):
    """Set up, run and analyse one case, appending its row.

    Safe to call from several threads at once: every case owns its directory,
    and each OpenFAST instance is a separate process, so ROSCO's module-level
    state (the instRes_* counters in particular) is private to each run. Only
    the shared `rows`/`tilts`/stdout touches are serialised, via `lock`.
    """
    c["tilt_deg"] = tilts.get(c["U"], args.default_tilt)
    d = setup_case(c)
    outpath = os.path.join(d, FST.replace(".fst", ".out"))
    tag = "[{}/{}] {:<24}".format(idx, total, c["id"])

    if args.analyze_only:
        if not os.path.exists(outpath):
            with lock:
                print(tag + " no .out, skipped", flush=True)
            return
        rc, dur = 0, 0.0
    else:
        rc, dur = run_case(c, args.timeout)

    if not os.path.exists(outpath):
        with lock:
            print(tag + " FAILED (rc={}, {:.0f}s) - see openfast.log".format(rc, dur),
                  flush=True)
            rows.append(dict(id=c["id"], group=c["group"], completed=0,
                             note="no output, rc={}".format(rc)))
            write_results(results, rows)
        return

    m = analyse(c, outpath)
    note = ""
    if not c["pppr"]:
        tv = measure_tilt(outpath)
        if tv is not None:
            tilts[c["U"]] = tv
            note = "[tilt {:.2f} deg] ".format(tv)

    row = dict(id=c["id"], group=c["group"], U=c["U"], omega=c["omega"],
               freq_hz=round(c["_freq_hz"], 6), amp_phi_deg=c["amp_phi_deg"],
               amp_omega=c["amp_omega"], offset_dev_deg=c["offset_dev_deg"],
               tilt_deg=round(c["tilt_deg"], 3), waves=int(c["waves"]),
               pc_finepit=c["pc_finepit_rad"],
               kp=round(c["_gains"]["kp"], 6), kr=round(c["_gains"]["kr"], 6),
               kp_tg="{:.4e}".format(c["_gains"]["kp_tg"]),
               kr_tg="{:.4e}".format(c["_gains"]["kr_tg"]),
               runtime_s=round(dur, 1))
    row.update(m)

    with lock:
        # power delta vs the baseline at the same wind speed
        base = next((r for r in rows if r.get("group") == "baseline"
                     and r.get("U") == c["U"]), None)
        if base and base.get("pwr_mean_kW") and m.get("pwr_mean_kW"):
            row["pwr_vs_base_pct"] = round(
                100.0 * (m["pwr_mean_kW"] / base["pwr_mean_kW"] - 1.0), 2)
        rows.append(row)
        print(tag + " {}{} t={}s P={} kW {}".format(
            note, "OK " if m.get("completed") else "DIVERGED",
            m.get("t_end"), m.get("pwr_mean_kW"),
            "amp_ratio={}".format(m.get("ptfm_amp_ratio")) if c["pppr"] else ""),
            flush=True)
        write_results(results, rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None,
                    help="comma-separated substring filters on case id, e.g. 'A_,C_'")
    ap.add_argument("--list", action="store_true", help="print matrix and exit")
    ap.add_argument("--analyze-only", action="store_true",
                    help="recompute metrics from existing .out files")
    ap.add_argument("--default-tilt", type=float, default=5.3,
                    help="fallback equilibrium tilt [deg] if no baseline was run")
    ap.add_argument("--results", default=RESULTS,
                    help="output CSV (relative names land in SweepRuns/)")
    ap.add_argument("--jobs", "-j", type=int, default=1,
                    help="OpenFAST instances to run at once. Each is single-threaded "
                         "and peaks near 50 MB, so the limit is physical cores; -j0 "
                         "uses all of them minus two.")
    ap.add_argument("--timeout", type=int, default=7200)
    args = ap.parse_args()

    results = args.results
    if not os.path.isabs(results):
        results = os.path.join(SWEEP, results)

    jobs = args.jobs
    if jobs <= 0:
        jobs = max(1, (os.cpu_count() or 2) - 2)

    cases = build_cases()
    if args.only:
        pats = [p.strip() for p in args.only.split(",") if p.strip()]
        cases = [c for c in cases if any(p in c["id"] for p in pats)]
    # baselines first: they supply the tilt and the power reference
    cases.sort(key=lambda c: (c["pppr"], c["id"]))

    if args.list:
        print("{:<24} {:<22} {:>6} {:>8} {:>8} {:>8} {:>7}".format(
            "id", "group", "U", "w", "amp", "dev", "waves"))
        for c in cases:
            print("{:<24} {:<22} {:>6.1f} {:>8.3f} {:>8.2f} {:>8.1f} {:>7}".format(
                c["id"], c["group"], c["U"], c["omega"],
                c["amp_phi_deg"], c["offset_dev_deg"], str(c["waves"])))
        print("\n{} cases".format(len(cases)))
        return

    ensure_layout()
    tilts, rows = {}, []
    lock = threading.Lock()
    total = len(cases)
    order = {c["id"]: i for i, c in enumerate(cases)}
    print("{} cases -> {}   ({} at a time)\n".format(total, SWEEP, jobs))

    def run_group(group, start):
        if not group:
            return
        if jobs == 1:
            for k, c in enumerate(group):
                process_case(c, args, tilts, rows, lock, results, start + k, total)
            return
        with futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            fs = [ex.submit(process_case, c, args, tilts, rows, lock, results,
                            start + k, total) for k, c in enumerate(group)]
            for f in fs:
                f.result()          # re-raise anything that blew up in a worker

    # Baselines must ALL finish before any PPPR case is set up: they measure the
    # equilibrium tilt that every PPPR reference is built on, and setup_case
    # bakes that tilt into PPPR_offset_phi. Running the two groups concurrently
    # would silently fall back to --default-tilt.
    baselines = [c for c in cases if not c["pppr"]]
    run_group(baselines, 1)
    run_group([c for c in cases if c["pppr"]], len(baselines) + 1)

    # Completion order is nondeterministic under -j; restore matrix order.
    with lock:
        rows.sort(key=lambda r: order.get(r.get("id"), 10 ** 6))
        write_results(results, rows)
    print("\nresults -> {}".format(results))


if __name__ == "__main__":
    sys.exit(main())
