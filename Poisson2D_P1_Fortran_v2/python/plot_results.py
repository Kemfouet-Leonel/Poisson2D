#!/usr/bin/env python3
"""
Visualisation et étude des erreurs du solveur Poisson 2D.

Le script :
    1. lance le solveur Fortran sur plusieurs maillages ;
    2. lit les erreurs dans les fichiers metrics.csv ;
    3. affiche les erreurs dans le terminal ;
    4. enregistre un tableau de convergence CSV ;
    5. affiche deux fenêtres :
       - solutions approchée et exacte en cartes de couleurs 2D ;
       - logarithmes des erreurs en fonction de log(h).

Aucun fichier JSON n'est utilisé.

Exemple :
    python3 python/plot_results.py --case 1

Maillages personnalisés :
    python3 python/plot_results.py --case 1 --meshes 8 16 32 64 128
"""

# =============================================================================
# 1. Importer les bibliothèques
# =============================================================================

from pathlib import Path
import argparse
import csv
import subprocess

import numpy as np
import matplotlib


# =============================================================================
# 2. Lire les paramètres du terminal
# =============================================================================

parser = argparse.ArgumentParser(
    description="Visualisation et convergence de Poisson 2D."
)

parser.add_argument(
    "--case",
    type=int,
    choices=[1, 2],
    default=1,
    help="1 : solution trigonométrique ; 2 : solution polynomiale."
)

parser.add_argument(
    "--meshes",
    type=int,
    nargs="+",
    default=[8, 16, 32, 64],
    help="Nombres de subdivisions dans chaque direction."
)

parser.add_argument(
    "--no-show",
    action="store_true",
    help="Enregistrer les images sans ouvrir les fenêtres."
)

args = parser.parse_args()

# Il faut plusieurs maillages distincts pour étudier la convergence.
meshes = sorted(set(args.meshes))

if len(meshes) < 2:
    parser.error("Indiquer au moins deux maillages distincts.")

if min(meshes) < 2 or max(meshes) > 4096:
    parser.error("Les subdivisions doivent être comprises entre 2 et 4096.")

# Choisir un backend sans fenêtre uniquement si demandé.
if args.no_show:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation


# =============================================================================
# 3. Définir les chemins du projet
# =============================================================================

ROOT = Path(__file__).resolve().parents[1]

EXECUTABLE = ROOT / "bin" / "poisson2d"

OUTPUT = (
    ROOT
    / "results"
    / "visualisation"
    / f"case_{args.case}"
)

OUTPUT.mkdir(parents=True, exist_ok=True)

if not EXECUTABLE.is_file():
    raise FileNotFoundError(
        f"Exécutable absent : {EXECUTABLE}\n"
        "Compiler d'abord le programme avec la commande make."
    )


# =============================================================================
# 4. Définir la solution exacte
# =============================================================================

def exact_solution(x, y):
    """Évaluer la solution exacte du cas choisi."""
    if args.case == 1:
        return np.sin(np.pi * x) * np.sin(np.pi * y)

    return x * (1.0 - x) * y * (1.0 - y)


# =============================================================================
# 5. Lancer Fortran sur un maillage et lire ses résultats
# =============================================================================

def run_solver(n):
    """Calculer une solution et retourner les mesures d'erreur."""
    folder = OUTPUT / f"n{n}"
    folder.mkdir(parents=True, exist_ok=True)

    command = [
        str(EXECUTABLE),
        str(n),
        str(args.case),
        str(folder)
    ]

    # Capturer la sortie pour garder le tableau du terminal lisible.
    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True
    )

    # Conserver le journal de chaque calcul.
    log_text = result.stdout + result.stderr

    (folder / "run.log").write_text(
        log_text,
        encoding="utf-8"
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Échec du solveur pour n={n}.\n"
            f"Consulter : {folder / 'run.log'}\n\n"
            f"{log_text}"
        )

    metrics_file = folder / "metrics.csv"

    if not metrics_file.is_file():
        raise FileNotFoundError(
            f"Le solveur n'a pas produit : {metrics_file}"
        )

    with metrics_file.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)

    if len(rows) != 1:
        raise ValueError(
            f"Une ligne de mesures est attendue dans {metrics_file}."
        )

    row = rows[0]

    required = ("h", "l2", "h1_semi", "h1", "linf_sample")

    for name in required:
        if name not in row:
            raise ValueError(
                f"Colonne absente de {metrics_file.name} : {name}"
            )

    h = float(row["h"])

    errors = np.array([
        float(row["l2"]),
        float(row["h1_semi"]),
        float(row["h1"]),
        float(row["linf_sample"])
    ])

    if not np.isfinite(h) or h <= 0:
        raise ValueError(f"Pas invalide pour n={n}.")

    if not np.isclose(h, 1.0 / n):
        raise ValueError(f"Pas incompatible avec n={n}.")

    if not np.isfinite(errors).all() or np.any(errors <= 0):
        raise ValueError(
            f"Erreurs non finies ou non positives pour n={n}."
        )

    # Le nombre de nœuds vaut (n+1)^2 pour le maillage du projet.
    node_count = (n + 1) ** 2

    return h, node_count, errors


# =============================================================================
# 6. Effectuer la campagne et afficher les normes dans le terminal
# =============================================================================

print()
print(f"ÉTUDE DE CONVERGENCE — CAS {args.case}")
print()
print(
    f"{'n':>6} "
    f"{'Nœuds':>9} "
    f"{'h':>12} "
    f"{'Erreur L2':>15} "
    f"{'Erreur H1 semi':>17} "
    f"{'Erreur H1':>15} "
    f"{'Linf échant.':>15}"
)
print("-" * 100)

steps = []
node_counts = []
error_values = []

for n in meshes:
    h, node_count, errors = run_solver(n)

    steps.append(h)
    node_counts.append(node_count)
    error_values.append(errors)

    print(
        f"{n:6d} "
        f"{node_count:9d} "
        f"{h:12.5e} "
        f"{errors[0]:15.6e} "
        f"{errors[1]:17.6e} "
        f"{errors[2]:15.6e} "
        f"{errors[3]:15.6e}",
        flush=True
    )

print("-" * 100)
print("Linf échant. = maximum échantillonné, non supremum certifié.")
print()

steps = np.asarray(steps)
node_counts = np.asarray(node_counts)
error_values = np.asarray(error_values)


# =============================================================================
# 7. Calculer et afficher les ordres expérimentaux
# =============================================================================

# Le premier maillage n'a pas de prédécesseur.
orders = np.full_like(error_values, np.nan)

orders[1:, :] = (
    np.log(error_values[:-1, :] / error_values[1:, :])
    / np.log(steps[:-1] / steps[1:])[:, None]
)

print("ORDRES EXPÉRIMENTAUX ENTRE MAILLAGES SUCCESSIFS")
print()
print(
    f"{'Raffinement':>15} "
    f"{'p L2':>12} "
    f"{'p H1 semi':>12} "
    f"{'p H1':>12} "
    f"{'p Linf éch.':>14}"
)
print("-" * 70)

for i in range(1, len(meshes)):
    refinement = f"{meshes[i - 1]} -> {meshes[i]}"

    print(
        f"{refinement:>15} "
        f"{orders[i, 0]:12.5f} "
        f"{orders[i, 1]:12.5f} "
        f"{orders[i, 2]:12.5f} "
        f"{orders[i, 3]:14.5f}"
    )

print()


# =============================================================================
# 8. Enregistrer les mesures et les ordres dans un tableau CSV
# =============================================================================

table_file = OUTPUT / "convergence.csv"

with table_file.open(
    "w",
    newline="",
    encoding="utf-8"
) as stream:
    writer = csv.writer(stream)

    writer.writerow([
        "n",
        "ndof",
        "h",
        "l2",
        "h1_semi",
        "h1",
        "linf_sample",
        "p_l2",
        "p_h1_semi",
        "p_h1",
        "p_linf_sample"
    ])

    for i, n in enumerate(meshes):
        # Laisser les premiers ordres vides : ils sont indéfinis.
        order_row = (
            ["", "", "", ""]
            if i == 0
            else orders[i, :].tolist()
        )

        writer.writerow(
            [n, int(node_counts[i]), steps[i]]
            + error_values[i, :].tolist()
            + order_row
        )

print(f"Tableau enregistré : {table_file}")


# =============================================================================
# 9. Lire la solution du maillage le plus fin
# =============================================================================

finest_n = meshes[-1]
solution_file = OUTPUT / f"n{finest_n}" / "solution.dat"

data = np.loadtxt(solution_file, ndmin=2)

if data.shape != ((finest_n + 1) ** 2, 4):
    raise ValueError(
        "Dimensions inattendues dans solution.dat."
    )

if not np.isfinite(data).all():
    raise ValueError("La solution contient des valeurs non finies.")

x = data[:, 0]
y = data[:, 1]
u_approx = data[:, 2]
u_exact_nodes = data[:, 3]

if not np.allclose(
    u_exact_nodes,
    exact_solution(x, y),
    rtol=1.0e-10,
    atol=1.0e-12
):
    raise ValueError(
        "Les valeurs exactes du fichier ne correspondent pas au cas choisi."
    )


# =============================================================================
# 10. Reconstruire les triangles du maillage Fortran
# =============================================================================

# Vérifier la numérotation avant de reconstruire la connectivité.
grid = np.linspace(0.0, 1.0, finest_n + 1)
X_nodes, Y_nodes = np.meshgrid(grid, grid)

if not (
    np.allclose(x, X_nodes.ravel())
    and np.allclose(y, Y_nodes.ravel())
):
    raise ValueError(
        "La numérotation des nœuds ne correspond pas à domain.f90."
    )

triangles = []

for j in range(finest_n):
    for i in range(finest_n):
        p = j * (finest_n + 1) + i
        q = p + finest_n + 1

        triangles.append([p, p + 1, q])
        triangles.append([p + 1, q + 1, q])

triangulation = Triangulation(
    x,
    y,
    triangles=np.asarray(triangles, dtype=int)
)


# =============================================================================
# 11. Préparer une grille fine pour la solution analytique
# =============================================================================

fine_grid = np.linspace(0.0, 1.0, 301)
X, Y = np.meshgrid(fine_grid, fine_grid)

U_exact = exact_solution(X, Y)

# Échelle de couleurs commune aux deux sous-figures.
vmin = min(float(u_approx.min()), float(U_exact.min()))
vmax = max(float(u_approx.max()), float(U_exact.max()))

if vmax == vmin:
    vmax = vmin + 1.0e-12

levels = np.linspace(vmin, vmax, 61)


# =============================================================================
# 12. Fenêtre 1 : solutions approchée et exacte en 2D
# =============================================================================

fig_solution = plt.figure(
    num=1,
    figsize=(13, 6),
    constrained_layout=True
)

subfigures = fig_solution.subfigures(1, 2)

# Sous-figure gauche : approximation P1.
ax_approx = subfigures[0].subplots()

image_approx = ax_approx.tricontourf(
    triangulation,
    u_approx,
    levels=levels,
    cmap="viridis"
)

subfigures[0].suptitle("Solution approchée")

subfigures[0].colorbar(
    image_approx,
    ax=ax_approx,
    shrink=0.8,
    label="Valeur de u"
)

# Sous-figure droite : solution analytique.
ax_exact = subfigures[1].subplots()

image_exact = ax_exact.contourf(
    X,
    Y,
    U_exact,
    levels=levels,
    cmap="viridis"
)

subfigures[1].suptitle("Solution exacte")

subfigures[1].colorbar(
    image_exact,
    ax=ax_exact,
    shrink=0.8,
    label="Valeur de u"
)

for ax in (ax_approx, ax_exact):
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")

fig_solution.suptitle(
    f"Poisson 2D — cas {args.case} — "
    f"maillage {finest_n} × {finest_n}"
)

solution_image = OUTPUT / "solutions_2d.png"
fig_solution.savefig(solution_image, dpi=200)


# =============================================================================
# 13. Fenêtre 2 : logarithmes des erreurs en fonction de log(h)
# =============================================================================

fig_error, ax_error = plt.subplots(
    num=2,
    figsize=(10, 7),
    constrained_layout=True
)

labels = (
    r"Erreur $L^2$",
    r"Erreur $H^1$ semi-norme",
    r"Erreur $H^1$ complète",
    r"Erreur $L^\infty$ échantillonnée"
)

markers = ("o", "s", "D", "^")

log_h = np.log(steps)

for index, label in enumerate(labels):
    ax_error.plot(
        log_h,
        np.log(error_values[:, index]),
        marker=markers[index],
        linewidth=1.8,
        markersize=6,
        label=label
    )

# Repères graphiques de pentes 2 et 1.
# Ces repères ne remplacent pas les ordres calculés dans le tableau.
reference_l2 = (
    np.log(error_values[-1, 0])
    + 2.0 * (log_h - log_h[-1])
)

reference_h1 = (
    np.log(error_values[-1, 1])
    + (log_h - log_h[-1])
)

ax_error.plot(
    log_h,
    reference_l2,
    "--",
    color="black",
    label="Repère : pente 2"
)

ax_error.plot(
    log_h,
    reference_h1,
    ":",
    color="gray",
    linewidth=2,
    label="Repère : pente 1"
)

ax_error.set_xlabel(r"$\log(h)$")
ax_error.set_ylabel(r"$\log(E)$")

ax_error.set_title(
    "Convergence des erreurs en fonction du pas du maillage"
)

ax_error.grid(True, linestyle="--", alpha=0.4)
ax_error.legend()

# Le raffinement se lit de gauche à droite.
ax_error.invert_xaxis()

error_image = OUTPUT / "log_erreurs.png"
fig_error.savefig(error_image, dpi=200)


# =============================================================================
# 14. Afficher les deux fenêtres
# =============================================================================

print(f"Figure des solutions : {solution_image}")
print(f"Figure des erreurs   : {error_image}")

if not args.no_show:
    # Un seul appel affiche les deux figures déjà construites.
    plt.show()

plt.close("all")