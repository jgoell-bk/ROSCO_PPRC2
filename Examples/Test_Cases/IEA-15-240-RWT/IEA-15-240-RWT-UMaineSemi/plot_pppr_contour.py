#!/usr/bin/env python3
"""
Regime map: controller tracking error as a contour over the forcing parameter
space, which is what says where the approach works and where it does not.

Reads a results CSV written by run_pppr_sweep.py and contours the demand-vs-
simulation errors against two swept inputs (tilt amplitude and frequency by
default). Cases that diverged are drawn as hollow markers rather than
interpolated through -- a diverged run has no meaningful tracking error, and
smoothing over it would paint a stable-looking region that is not.

Usage
-----
    python plot_pppr_contour.py                       # E_grid from results.csv
    python plot_pppr_contour.py --results r.csv --group E_grid
    python plot_pppr_contour.py --x amp_phi_deg --y omega
    python plot_pppr_contour.py --save regime.png --no-show
"""

import argparse
import csv
import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_RESULTS = os.path.join(HERE, "SweepRuns", "results.csv")

# (column, panel title, colormap, symmetric-about-this or None)
PANELS = [
    ("ptfm_rmse",      "Platform tilt RMSE [deg]",      "viridis_r", None),
    ("ptfm_amp_ratio", "Tilt amplitude ratio [-]",      "RdBu_r",    1.0),
    ("ptfm_phase_err", "Tilt phase error [deg]",        "RdBu_r",    0.0),
    ("pwr_vs_base_pct", "Power vs baseline [%]",        "RdBu_r",    0.0),
]


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=DEF_RESULTS)
    ap.add_argument("--group", default="E_grid", help="group column to select")
    ap.add_argument("--x", default="amp_phi_deg")
    ap.add_argument("--y", default="omega")
    ap.add_argument("--save", default=None)
    ap.add_argument("--no-show", action="store_true")
    a = ap.parse_args()

    if a.no_show:
        matplotlib.use("Agg")

    path = a.results if os.path.isabs(a.results) else os.path.join(HERE, "SweepRuns", a.results)
    rows = [r for r in csv.DictReader(open(path)) if r.get("group") == a.group]
    if not rows:
        raise SystemExit("no rows with group == '{}' in {}".format(a.group, path))

    ok = [r for r in rows if r.get("completed") == "1"]
    bad = [r for r in rows if r.get("completed") != "1"]
    print("{}: {} cases, {} completed, {} diverged".format(a.group, len(rows), len(ok), len(bad)))
    if len(ok) < 3:
        raise SystemExit("need at least 3 completed cases to contour; got {}".format(len(ok)))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (col, title, cmap, sym) in zip(axes.ravel(), PANELS):
        pts = [(fnum(r[a.x]), fnum(r[a.y]), fnum(r.get(col))) for r in ok]
        pts = [p for p in pts if None not in p]
        if len(pts) < 3:
            ax.set_title(title + "\n(insufficient data)")
            ax.set_axis_off()
            continue
        xs, ys, zs = zip(*pts)
        tri = mtri.Triangulation(xs, ys)
        kw = {}
        if sym is not None:
            half = max(abs(min(zs) - sym), abs(max(zs) - sym)) or 1.0
            kw = dict(vmin=sym - half, vmax=sym + half)
        cf = ax.tricontourf(tri, zs, levels=14, cmap=cmap, **kw)
        ax.tricontour(tri, zs, levels=14, colors="k", linewidths=0.3, alpha=0.4)
        fig.colorbar(cf, ax=ax)
        ax.plot(xs, ys, "k.", ms=4)
        # diverged cases: shown, never interpolated through
        for r in bad:
            bx, by = fnum(r[a.x]), fnum(r[a.y])
            if bx is not None and by is not None:
                ax.plot(bx, by, "x", color="k", ms=9, mew=2)
        ax.set_title(title)
        ax.set_xlabel(a.x)
        ax.set_ylabel(a.y)

    fig.suptitle("PPPR regime map -- {}  (x = diverged)".format(a.group))
    fig.tight_layout()
    if a.save:
        fig.savefig(a.save, dpi=130)
        print("wrote " + a.save)
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
