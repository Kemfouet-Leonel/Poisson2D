"""Campagne reproductible : CSV, ordres, figures et performances.
Les temps sont mesures de nouveau a chaque repetition, pas extrapoles.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import platform
import statistics
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", nargs="+", type=int, default=[8, 12, 16, 24, 32, 48, 64])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--case", dest="case_id", type=int, choices=[1, 2], default=1)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.repeat < 1 or min(args.meshes) < 2 or len(set(args.meshes)) < 3:
        parser.error("Au moins 3 maillages distincts >=2 et repeat>=1 requis")
    out = ROOT / "results" / f"case_{args.case_id}"
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for n in sorted(set(args.meshes)):
        runs = []
        for rep in range(args.repeat):
            dest = out / f"n{n}" / f"repeat{rep}"
            dest.mkdir(parents=True, exist_ok=True)
            subprocess.run([str(ROOT / "bin/poisson2d"), str(n), str(args.case_id), str(dest)], check=True)
            with (dest / "metrics.csv").open() as f:
                row = {k: float(v) for k, v in next(csv.DictReader(f)).items()}
            runs.append(row)
        row = runs[-1].copy()
        for key in ["mesh_s", "assembly_s", "solve_s", "error_s"]:
            row[key] = statistics.median(r[key] for r in runs)
            row[key + "_min"] = min(r[key] for r in runs)
            row[key + "_max"] = max(r[key] for r in runs)
        row["repetitions"] = args.repeat
        for key in ["l2", "h1_semi", "h1", "linf_sample"]:
            row["p_" + key] = float("nan") if not records else float(
                np.log(records[-1][key] / row[key]) / np.log(records[-1]["h"] / row["h"]))
        row["dense_matrix_bytes"] = 8 * int(row["ndof"]) ** 2
        records.append(row)
    table = out / "convergence.csv"
    with table.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    environment = {
        "platform": platform.platform(), "python": sys.version,
        "numpy": np.__version__, "command": sys.argv,
        "compiler": subprocess.check_output(["gfortran", "--version"], text=True).splitlines()[0],
        "revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True).stdout.strip(),
        "note": "csr_bytes=stockage matrice seul, pas memoire RSS du processus"
    }
    (out / "environment.json").write_text(json.dumps(environment, indent=2))
    h = np.array([r["h"] for r in records])
    fig, ax = plt.subplots(figsize=(9, 6))
    for key in ["l2", "h1_semi", "h1", "linf_sample"]:
        ax.plot(np.log(h), np.log([r[key] for r in records]), "o-", label=key)
    ax.set(xlabel="log(h)", ylabel="log(erreur)", title="Convergence P1 — mesures reelles")
    ax.invert_xaxis(); ax.legend(); ax.grid(); fig.tight_layout()
    fig.savefig(out / "convergence.png", dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    dof = [r["ndof"] for r in records]
    for key in ["assembly_s", "solve_s", "error_s"]:
        axes[0].loglog(dof, [r[key] for r in records], "o-", label=key)
    axes[0].set(xlabel="Degres de liberte", ylabel="Temps median (s)")
    for key in ["csr_bytes", "dense_matrix_bytes"]:
        axes[1].loglog(dof, [r[key] for r in records], "o-", label=key)
    axes[1].set(xlabel="Degres de liberte", ylabel="Stockage matrice (octets)")
    for ax in axes: ax.legend(); ax.grid()
    fig.tight_layout(); fig.savefig(out / "performance.png", dpi=180); plt.close(fig)
    final = out / f"n{int(records[-1]['n'])}" / f"repeat{args.repeat-1}" / "solution.dat"
    x, y, uh, ue = np.loadtxt(final).T
    n = int(records[-1]["n"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    vmin = min(uh.min(), ue.min()); vmax = max(uh.max(), ue.max())
    for ax, values, title in zip(axes, [uh, ue], ["Approchee", "Exacte aux noeuds"]):
        image = ax.pcolormesh(x.reshape(n+1,n+1), y.reshape(n+1,n+1), values.reshape(n+1,n+1),
                              shading="nearest", cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set(title=title, xlabel="x", ylabel="y", aspect="equal")
        fig.colorbar(image, ax=ax)
    fig.tight_layout(); fig.savefig(out / "solution.png", dpi=180); plt.close(fig)
    if args.check:
        for r in records:
            if not all(np.isfinite(r[k]) and r[k] > 0 for k in ["l2", "h1", "h1_semi", "linf_sample"]):
                raise AssertionError("Erreur non finie ou non positive")
            if r["relres"] > 1.01e-11: raise AssertionError("Residu excessif")
            if abs(r["h1"]**2-r["l2"]**2-r["h1_semi"]**2) > 1e-12:
                raise AssertionError("Identite H1 invalide")
        # Test de non-regression pour CE probleme, pas theoreme universel.
        for key, low, high in [("l2", 1.7, 2.3), ("h1_semi", .8, 1.2)]:
            p = np.polyfit(np.log(h[-3:]), np.log([r[key] for r in records][-3:]), 1)[0]
            if not low < p < high:
                raise AssertionError(f"Ordre final {key} hors intervalle: {p}")
        print("Validation residu/normes/ordres : OK")
    print(f"Resultats : {out}")

if __name__ == "__main__":
    main()
