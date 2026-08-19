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

* Only diagnostics/performance channels are written. Per-blade-node output is
  disabled, which is what drops the .out files from ~350 MB to a few MB.
"""

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
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
# Turbine constants for the gain formula (evaluated at TSR0=8.5, beta0=0
# from IEA15MW_Cp_Ct_Cq.mat -- see MatlabSimulations/PR_control_nonlinear_sim.m).
# ----------------------------------------------------------------------------
TSR0 = 8.5
RHO = 1.2
RROT = 120.0
JR = 3.525e8
NG = 1.0
CP0 = 0.46969
DCP_DTSR = 0.005270
DCP_DBETA = -0.004124          # per DEGREE (data.angles is in degrees)


def pir_gains(U, omega, zeta=0.7, ratio=10.0):
    """Modified-Abbas PIR gains at wind speed U and forcing freq omega [rad/s].

    Returns ROSCO-native units: kp/kr in rad-pitch per rad-error (the reference
    formula is deg-pitch per rad-error, hence the pi/180), kp_Tg/kr_Tg in N-m
    per rad/s.
    """
    A = 0.5 * RHO * math.pi * RROT**4 * U / (JR * TSR0**2) * (DCP_DTSR * TSR0 - CP0)
    B = NG / (2 * JR * TSR0**2) * RHO * math.pi * RROT**3 * U**2 * (DCP_DBETA * TSR0)
    root = 2 * zeta * omega + A
    kp = -1.0 / (2 * math.pi * B) * root * math.pi / 180.0
    kp_tg = JR / NG**2 * root
    return dict(kp=kp, kr=kp / ratio, kp_tg=kp_tg, kr_tg=kp_tg / ratio)


# ============================================================================
# CASE MATRIX -- edit here
# ============================================================================
DEFAULTS = dict(
    U=10.0,
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
    g = c["gains"] or pir_gains(c["U"], c["omega"], c["zeta"], c["ratio"])
    tilt = c["tilt_deg"]
    offset_phi_rad = math.radians(tilt + c["offset_dev_deg"])
    omega_ref = TSR0 * c["U"] / RROT

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
        set_discon(L, "PPPR_Mode", 2)
        set_discon(L, "PC_ControlMode", 0)
        set_discon(L, "VS_ControlMode", 0)
        set_discon(L, "SS_Mode", 0)
        set_discon(L, "PS_Mode", 0)
        set_discon(L, "Fl_Mode", 0)
        set_discon(L, "PPPR_amp_phi", "{:.7f}".format(math.radians(c["amp_phi_deg"])))
        set_discon(L, "PPPR_freq_phi", "{:.7f}".format(freq_hz))
        set_discon(L, "PPPR_amp_omega", "{:.6f}".format(c["amp_omega"]))
        set_discon(L, "PPPR_freq_omega", "{:.7f}".format(freq_hz))
        set_discon(L, "Phi_phaseoffset", "{:.2f}".format(c["phi_phase_deg"]))
        set_discon(L, "Omega_phaseoffset", "{:.2f}".format(c["omega_phase_deg"]))
        set_discon(L, "PPPR_fz_phi", "{:.8f}".format(freq_hz / 10.0))
        set_discon(L, "PPPR_fz_omega", "{:.8f}".format(freq_hz / 10.0))
        set_discon(L, "PPPR_offset_phi", "{:.7f}".format(offset_phi_rad))
        set_discon(L, "PPPR_offset_omega", "{:.6f}".format(omega_ref))
        set_discon(L, "PPPR_CntrGains_phi", "{:.6f}  {:.6f}".format(g["kp"], g["kr"]))
        set_discon(L, "PPPR_CntrGains_omega", "{:.5e}  {:.5e}".format(g["kp_tg"], g["kr_tg"]))
    else:
        set_discon(L, "PPPR_Mode", 0)
        set_discon(L, "PC_ControlMode", 1)
        set_discon(L, "VS_ControlMode", 1)
        set_discon(L, "SS_Mode", 1)
        set_discon(L, "PS_Mode", 1)
        set_discon(L, "Fl_Mode", 0)
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


def amp_at(t, x, f):
    """Amplitude of the component of x at frequency f [Hz] (single-bin projection)."""
    n = len(x)
    if n == 0:
        return float("nan")
    m = sum(x) / n
    re_ = im = 0.0
    for ti, xi in zip(t, x):
        a = 2 * math.pi * f * ti
        re_ += (xi - m) * math.cos(a)
        im -= (xi - m) * math.sin(a)
    return 2.0 * math.hypot(re_, im) / n


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
            r["ptfm_amp"] = round(amp_at(ts, pp, c["_freq_hz"]), 4)
            r["ptfm_ref_amp"] = round(c["amp_phi_deg"], 4)
            if c["amp_phi_deg"] > 0:
                r["ptfm_amp_ratio"] = round(r["ptfm_amp"] / c["amp_phi_deg"], 3)

    gs = col("GenSpeed")
    if gs:
        r["gen_mean_rpm"] = round(sum(gs) / len(gs), 4)
        if c["pppr"]:
            ref = c["_omega_ref"] * 60 / (2 * math.pi)
            r["gen_ref_rpm"] = round(ref, 4)
            r["gen_mean_err"] = round(r["gen_mean_rpm"] - ref, 4)
            r["gen_amp"] = round(amp_at(ts, gs, c["_freq_hz"]), 5)
            r["gen_ref_amp"] = round(c["amp_omega"] * 60 / (2 * math.pi), 5)

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
    ap.add_argument("--timeout", type=int, default=7200)
    args = ap.parse_args()

    results = args.results
    if not os.path.isabs(results):
        results = os.path.join(SWEEP, results)

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
    print("{} cases -> {}\n".format(len(cases), SWEEP))

    for i, c in enumerate(cases, 1):
        c["tilt_deg"] = tilts.get(c["U"], args.default_tilt)
        d = setup_case(c)
        outpath = os.path.join(d, FST.replace(".fst", ".out"))

        print("[{}/{}] {:<24}".format(i, len(cases), c["id"]), end=" ", flush=True)
        if args.analyze_only:
            if not os.path.exists(outpath):
                print("no .out, skipped")
                continue
            rc, dur = 0, 0.0
        else:
            rc, dur = run_case(c, args.timeout)

        if not os.path.exists(outpath):
            print("FAILED (rc={}, {:.0f}s) - see openfast.log".format(rc, dur))
            rows.append(dict(id=c["id"], group=c["group"], completed=0,
                             note="no output, rc={}".format(rc)))
            continue

        m = analyse(c, outpath)
        if not c["pppr"]:
            tv = measure_tilt(outpath)
            if tv is not None:
                tilts[c["U"]] = tv
                print("[tilt {:.2f} deg] ".format(tv), end="")

        row = dict(id=c["id"], group=c["group"], U=c["U"], omega=c["omega"],
                   freq_hz=round(c["_freq_hz"], 6), amp_phi_deg=c["amp_phi_deg"],
                   amp_omega=c["amp_omega"], offset_dev_deg=c["offset_dev_deg"],
                   tilt_deg=round(c["tilt_deg"], 3), waves=int(c["waves"]),
                   kp=round(c["_gains"]["kp"], 6), kr=round(c["_gains"]["kr"], 6),
                   kp_tg="{:.4e}".format(c["_gains"]["kp_tg"]),
                   kr_tg="{:.4e}".format(c["_gains"]["kr_tg"]),
                   runtime_s=round(dur, 1))
        row.update(m)
        # power delta vs the baseline at the same wind speed
        base = next((r for r in rows if r.get("group") == "baseline"
                     and r.get("U") == c["U"]), None)
        if base and base.get("pwr_mean_kW") and m.get("pwr_mean_kW"):
            row["pwr_vs_base_pct"] = round(
                100.0 * (m["pwr_mean_kW"] / base["pwr_mean_kW"] - 1.0), 2)
        rows.append(row)

        print("{} t={}s P={} kW {}".format(
            "OK " if m.get("completed") else "DIVERGED",
            m.get("t_end"), m.get("pwr_mean_kW"),
            "amp_ratio={}".format(m.get("ptfm_amp_ratio")) if c["pppr"] else ""))

        # write incrementally so an interrupted sweep still leaves usable results
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(results, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)

    print("\nresults -> {}".format(results))


if __name__ == "__main__":
    sys.exit(main())
