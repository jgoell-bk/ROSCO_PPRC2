#!/usr/bin/env python3
"""
Reproduce the reference 4-panel PPPR figure: blade pitch, platform tilt,
rotor speed, power -- with the demand waveforms overlaid as dashed lines.

The reference curves are rebuilt from the DISCON.IN that the run used, so the
plot is self-configuring: point it at a .out and its DISCON and the dashed
lines are whatever that case actually commanded.

Usage
-----
    python plot_pppr_compare.py                          # defaults, this folder
    python plot_pppr_compare.py --out RUN.out --discon RUN_DISCON.IN
    python plot_pppr_compare.py --tmin 0               # crop the transient
    python plot_pppr_compare.py --save fig.png --no-show

Conventions follow the merged controller: PPPR_freq_* in rad/s, PPPR_amp_phi
and PPPR_offset_phi in degrees, phase entering as sin(w*t - phase*D2R).
"""

import argparse
import math
import os
import re

import matplotlib
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_OUT = os.path.join(HERE, "IEA-15-240-RWT-UMaineSemi.out")
DEF_DISCON = os.path.join(HERE, "IEA-15-240-RWT-UMaineSemi_DISCON.IN")


def read_out(path):
    """Parse an OpenFAST .out into {channel: [floats]} (stdlib only)."""
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()
    hdr = next((i for i, ln in enumerate(lines) if ln.strip().startswith("Time")), None)
    if hdr is None:
        raise SystemExit("no channel header found in {}".format(path))
    names = lines[hdr].split()
    data = {n: [] for n in names}
    for ln in lines[hdr + 2:]:                     # +2 skips the units row
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


def read_discon(path):
    """Pull '<value>  ! <Label>' records out of a DISCON.IN."""
    out = {}
    with open(path, "r", errors="replace") as f:
        for ln in f:
            if "!" not in ln:
                continue
            val, rest = ln.split("!", 1)
            m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)", rest)
            if not m:
                continue
            nums = []
            for tok in val.split():
                try:
                    nums.append(float(tok))
                except ValueError:
                    break
            if nums:
                out[m.group(1)] = nums
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEF_OUT)
    ap.add_argument("--discon", default=DEF_DISCON)
    ap.add_argument("--tmin", type=float, default=None, help="crop start [s]")
    ap.add_argument("--tmax", type=float, default=None, help="crop end [s]")
    ap.add_argument("--save", default=None, help="write PNG here")
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--autoscale", action="store_true",
                    help="let matplotlib pick y-limits instead of the reference figure's")
    a = ap.parse_args()

    if a.no_show:
        matplotlib.use("Agg")

    d = read_out(a.out)
    p = read_discon(a.discon)

    def par(name, i=0, default=0.0):
        v = p.get(name)
        return v[i] if v and len(v) > i else default

    t = d["Time"]
    lo = a.tmin if a.tmin is not None else t[0]
    hi = a.tmax if a.tmax is not None else t[-1]
    keep = [i for i, ti in enumerate(t) if lo <= ti <= hi]
    ts = [t[i] for i in keep]

    def ch(name, scale=1.0):
        if name not in d:
            raise SystemExit("channel '{}' not in {} -- check the OutList".format(name, a.out))
        return [d[name][i] * scale for i in keep]

    beta = ch("BldPitch1")                       # deg
    tilt = ch("PtfmPitch")                       # deg
    omega = ch("RotSpeed", 2 * math.pi / 60.0)   # rpm -> rad/s
    power = ch("GenPwr", 1e-3)                   # kW  -> MW

    D2R = math.pi / 180.0
    tilt_ref = [par("PPPR_offset_phi") +
                par("PPPR_amp_phi") * math.sin(ti * par("PPPR_freq_phi")
                                               - par("Phi_phaseoffset") * D2R) for ti in ts]
    omega_ref = [par("PPPR_offset_omega") +
                 par("PPPR_amp_omega") * math.sin(ti * par("PPPR_freq_omega")
                                                  - par("Omega_phaseoffset") * D2R) for ti in ts]

    fig, ax = plt.subplots(4, 1, figsize=(13, 9), sharex=True)
    ax[0].plot(ts, beta, lw=1.0)
    ax[0].set_ylabel(r"$\beta$")
    ax[0].axhline(0.0, color="0.7", lw=0.8, ls=":")     # the pitch floor that used to bind

    ax[1].plot(ts, tilt, lw=1.0)
    ax[1].plot(ts, tilt_ref, "r--", lw=1.2)
    ax[1].set_ylabel("Tilt")

    ax[2].plot(ts, omega, lw=1.0)
    ax[2].plot(ts, omega_ref, "r--", lw=1.2)
    ax[2].set_ylabel(r"$\omega$")

    ax[3].plot(ts, power, lw=1.0)
    ax[3].set_ylabel("Power (MW)")
    ax[3].set_xlabel("Time (s)")

    for x in ax:
        x.grid(True, alpha=0.3)
        x.set_xlim(ts[0], ts[-1])

    # Y-scales read off the reference figure, so the two plots can be laid side
    # by side. Clamping can push data off-panel -- notably omega, which
    # overshoots to ~1.03 rad/s here against a 0.855 ceiling there -- so every
    # series that leaves its box is reported rather than silently cropped.
    REF_YLIM = [(-4.5, 2.9,  [-4, -2, 0, 2]),
                (0.0, 10.0,  [0, 5, 10]),
                (0.695, 0.855, [0.70, 0.75, 0.80, 0.85]),
                (0.0, 15.0,  [0, 5, 10, 15])]
    if not a.autoscale:
        clipped = []
        for x, series, label, (ylo, yhi, ticks) in zip(
                ax, (beta, tilt, omega, power),
                ("beta", "tilt", "omega", "power"), REF_YLIM):
            dlo, dhi = min(series), max(series)
            if dlo < ylo or dhi > yhi:
                clipped.append("{} spans {:.3f}..{:.3f}, panel shows {:.3f}..{:.3f}"
                               .format(label, dlo, dhi, ylo, yhi))
            x.set_ylim(ylo, yhi)
            x.set_yticks(ticks)
        if clipped:
            print("NOTE: data outside the reference figure's y-limits "
                  "(re-run with --autoscale to see it all):")
            for c in clipped:
                print("  " + c)

    fig.tight_layout()

    # Settled-window summary -- the numbers you actually compare against a
    # collaborator's run. Uses the last two thirds unless --tmin was given.
    i0 = 0 if a.tmin is not None else next(
        (k for k, ti in enumerate(ts) if ti >= ts[0] + (ts[-1] - ts[0]) / 3.0), 0)
    n = len(ts) - i0

    def mean(v):
        return sum(v[i0:]) / n

    print("settled window: {:.0f}-{:.0f} s".format(ts[i0], ts[-1]))
    print("  beta   mean {:8.3f}  min {:8.3f}  max {:8.3f} deg".format(
        mean(beta), min(beta[i0:]), max(beta[i0:])))
    print("  below 0 deg: {:.1f}% of the time".format(
        100.0 * sum(1 for v in beta[i0:] if v < 0) / n))
    print("  tilt   mean {:8.3f} deg  (ref {:8.3f})".format(mean(tilt), mean(tilt_ref)))
    print("  omega  mean {:8.4f} rad/s (ref {:8.4f})".format(mean(omega), mean(omega_ref)))
    print("  power  mean {:8.3f} MW   min {:7.3f}  max {:7.3f}".format(
        mean(power), min(power[i0:]), max(power[i0:])))

    if a.save:
        fig.savefig(a.save, dpi=130)
        print("wrote {}".format(a.save))
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
