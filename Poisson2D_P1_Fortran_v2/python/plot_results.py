#!/usr/bin/env python3
# =============================================================================
# Script : python/plot_results.py
#
# Figure 1 :
#   Comparaison entre la solution approchée et la solution exacte aux nœuds.
#
# Figure 2 :
#   Logarithmes naturels des erreurs en fonction de log(h).
#
# Données utilisées :
#   results/convergence.csv
#   results/solution.dat
# =============================================================================

from pathlib import Path
import argparse
import csv

import numpy as np
import matplotlib

# Lire les arguments avant de charger pyplot pour pouvoir choisir
# un backend sans interface graphique avec --no-show.
parser = argparse.ArgumentParser(
    description="Tracer les solutions et la convergence de Poisson 2D."
)

parser.add_argument(
    "--case",
    dest="case_id",
    type=int,
    choices=[1, 2],
    default=1,
    help="Cas manufacturé à représenter : 1 ou 2."
)

parser.add_argument(
    "--solution",
    type=Path,
    default=None,
    help=(
        "Fichier solution.dat à utiliser pour la figure 1. "
        "Par défaut : dernière répétition du maillage le plus fin "
        "présent dans convergence.csv."
    )
)

parser.add_argument(
    "--no-show",
    action="store_true",
    help="Enregistrer les figures sans ouvrir de fenêtres."
)

args = parser.parse_args()

if args.no_show:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation


# Racine du projet, indépendante du dossier de lancement.
ROOT = Path(__file__).resolve().parents[1]


def read_convergence(path):
    """
    Lire le CSV produit par study.py.

    Les noms de colonnes évitent de dépendre de leur position.
    Les lignes sont triées du maillage grossier au maillage fin.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"Tableau de convergence introuvable : {path}\n"
            "Exécuter d'abord python/study.py."
        )

    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    if len(rows) < 2:
        raise ValueError(
            "Au moins deux maillages sont nécessaires "
            "pour représenter une convergence."
        )

    required = {
        "n",
        "h",
        "l2",
        "h1_semi",
        "h1",
        "linf_sample",
        "repetitions"
    }

    if not required.issubset(rows[0]):
        missing = required.difference(rows[0])

        raise ValueError(
            f"Colonnes manquantes dans {path.name} : "
            f"{sorted(missing)}"
        )

    # Ne convertir que les colonnes nécessaires.
    # Les colonnes p_* peuvent contenir un NaN intentionnel
    # pour le premier ordre, qui est indéfini.
    data = {
        key: np.array([float(row[key]) for row in rows])
        for key in required
    }

    for key, values in data.items():
        if not np.isfinite(values).all():
            raise ValueError(
                f"Valeur non finie dans la colonne {key}."
            )

    # Ordre décroissant en h : grossier -> fin.
    order = np.argsort(data["h"])[::-1]

    data = {
        key: values[order]
        for key, values in data.items()
    }

    if np.any(data["h"] <= 0):
        raise ValueError("Les pas h doivent être positifs.")

    if np.any(np.diff(data["h"]) >= 0):
        raise ValueError(
            "Les pas doivent être distincts après le tri."
        )

    for key in ("l2", "h1_semi", "h1", "linf_sample"):
        if np.any(data[key] <= 0):
            raise ValueError(
                f"Les erreurs {key} doivent être positives "
                "pour calculer leurs logarithmes."
            )

    return data


def read_solution(path):
    """
    Lire les colonnes :
        x, y, solution approchée, solution exacte.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"Fichier de solution introuvable : {path}"
        )

    data = np.loadtxt(path, ndmin=2)

    if data.shape[1] != 4:
        raise ValueError(
            "solution.dat doit contenir exactement quatre colonnes."
        )

    if not np.isfinite(data).all():
        raise ValueError(
            "La solution contient des valeurs NaN ou infinies."
        )

    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def build_triangulation(x, y):
    """
    Reconstruire le maillage triangulaire structuré du code Fortran.

    Les dimensions sont déduites des coordonnées.
    La numérotation attendue est ligne par ligne, x variant le plus vite.
    """
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    nx = x_unique.size - 1
    ny = y_unique.size - 1

    if nx < 1 or ny < 1:
        raise ValueError("Dimensions du maillage invalides.")

    if x.size != (nx + 1) * (ny + 1):
        raise ValueError(
            "Les coordonnées ne forment pas une grille complète."
        )

    expected_x, expected_y = np.meshgrid(
        x_unique,
        y_unique,
        indexing="xy"
    )

    if not (
        np.allclose(x, expected_x.ravel())
        and np.allclose(y, expected_y.ravel())
    ):
        raise ValueError(
            "L'ordre des nœuds ne correspond pas "
            "à la numérotation Fortran attendue."
        )

    triangles = []

    for j in range(ny):
        for i in range(nx):

            # Sommets inférieur gauche et supérieur gauche.
            p = j * (nx + 1) + i
            q = p + nx + 1

            # Même diagonale que dans domain.f90.
            triangles.append([p, p + 1, q])
            triangles.append([p + 1, q + 1, q])

    triangulation = Triangulation(
        x,
        y,
        np.asarray(triangles, dtype=int)
    )

    return triangulation, nx, ny


def plot_solutions(solution_file, output_file):
    """
    Figure 1 : les deux solutions en cartes de couleurs.

    Les niveaux sont communs pour éviter qu'une différence
    d'échelle de couleurs masque une différence de valeurs.
    """
    x, y, uh, u_exact = read_solution(solution_file)

    triangulation, nx, ny = build_triangulation(x, y)

    vmin = min(float(uh.min()), float(u_exact.min()))
    vmax = max(float(uh.max()), float(u_exact.max()))

    # Éviter des niveaux identiques si les champs sont constants.
    if vmax == vmin:
        padding = max(1.0, abs(vmin)) * 1.0e-6
        vmin -= padding
        vmax += padding

    levels = np.linspace(vmin, vmax, 61)

    fig, axes = plt.subplots(
        1,
        2,
        num=1,
        figsize=(13, 5),
        constrained_layout=True
    )

    fields = (uh, u_exact)

    titles = (
        r"Solution approchée $u_h$",
        r"Solution exacte aux nœuds $u$"
    )

    for ax, values, title in zip(axes, fields, titles):

        contour = ax.tricontourf(
            triangulation,
            values,
            levels=levels,
            cmap="viridis"
        )

        # Sur un maillage très fin, les arêtes masqueraient les couleurs.
        if max(nx, ny) <= 32:
            ax.triplot(
                triangulation,
                color="black",
                linewidth=0.2,
                alpha=0.25
            )

        ax.set_title(title)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_aspect("equal")

    colorbar = fig.colorbar(
        contour,
        ax=axes,
        shrink=0.9,
        pad=0.03
    )

    colorbar.set_label("Valeur de la solution")

    fig.suptitle(
        f"Poisson 2D — comparaison sur le maillage {nx} × {ny}"
    )

    fig.savefig(
        output_file,
        dpi=200,
        bbox_inches="tight"
    )

    print(f"Figure 1 enregistrée : {output_file}")


def plot_log_errors(data, output_file):
    """
    Figure 2 : log(erreur) en fonction de log(h).

    Les logarithmes sont calculés explicitement.
    Les axes de la figure restent donc linéaires.
    """
    h = data["h"]
    log_h = np.log(h)

    fig, ax = plt.subplots(
        num=2,
        figsize=(10, 6),
        constrained_layout=True
    )

    curves = (
        ("l2", "o", r"$\log(E_{L^2})$"),
        ("h1_semi", "s", r"$\log(E_{H^1,\mathrm{semi}})$"),
        ("h1", "D", r"$\log(E_{H^1})$"),
        (
            "linf_sample",
            "^",
            r"$\log(E_{\infty,\mathrm{échant.}})$"
        )
    )

    for key, marker, label in curves:
        ax.plot(
            log_h,
            np.log(data[key]),
            marker=marker,
            linewidth=1.8,
            markersize=6,
            label=label
        )

    # Repères graphiques ancrés sur le maillage le plus fin.
    # Ils ne modifient pas les mesures et ne prouvent pas les ordres.
    reference_l2 = (
        np.log(data["l2"][-1])
        + 2.0 * (log_h - log_h[-1])
    )

    reference_h1 = (
        np.log(data["h1_semi"][-1])
        + (log_h - log_h[-1])
    )

    ax.plot(
        log_h,
        reference_l2,
        "--",
        color="black",
        label="Repère : pente 2"
    )

    ax.plot(
        log_h,
        reference_h1,
        ":",
        color="gray",
        linewidth=2,
        label="Repère : pente 1"
    )

    # Le raffinement se lit de gauche à droite :
    # grand h à gauche, petit h à droite.
    ax.invert_xaxis()

    ax.set_xlabel(r"$\log(h)$")
    ax.set_ylabel(r"$\log(\mathrm{erreur})$")

    ax.set_title(
        "Évolution des erreurs avec le raffinement du maillage"
    )

    ax.grid(
        True,
        linestyle="--",
        alpha=0.45
    )

    ax.legend()

    fig.savefig(
        output_file,
        dpi=200,
        bbox_inches="tight"
    )

    print(f"Figure 2 enregistrée : {output_file}")


def main():
    """Lire les résultats du cas choisi et construire les deux figures."""
    results_dir = ROOT / "results" / f"case_{args.case_id}"

    convergence_file = results_dir / "convergence.csv"

    data = read_convergence(convergence_file)

    if args.solution is None:

        # Après tri, la dernière ligne correspond au plus petit h.
        n = int(data["n"][-1])
        repetitions = int(data["repetitions"][-1])

        if repetitions < 1:
            raise ValueError(
                "Le nombre de répétitions doit être positif."
            )

        solution_file = (
            results_dir
            / f"n{n}"
            / f"repeat{repetitions - 1}"
            / "solution.dat"
        )

    else:
        # Un chemin relatif explicite est interprété depuis
        # le répertoire courant du terminal.
        solution_file = args.solution.expanduser().resolve()

    print(f"Cas représenté : {args.case_id}")
    print(f"Solution utilisée : {solution_file}")
    print(f"Convergence utilisée : {convergence_file}")

    plot_solutions(
        solution_file,
        results_dir / "solution.png"
    )

    plot_log_errors(
        data,
        results_dir / "convergence.png"
    )

    # Sauvegarder les deux figures avant un éventuel affichage.
    if not args.no_show:
        plt.show()

    plt.close("all")


if __name__ == "__main__":
    main()