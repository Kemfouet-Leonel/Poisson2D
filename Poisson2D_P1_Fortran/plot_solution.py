#!/usr/bin/env python3
# =============================================================================
# Projet : résolution de l'équation de Poisson 2D par éléments finis P1.
#
# Objectif de ce script
# --------------------
# Lire les fichiers produits par le programme Fortran, puis construire :
#
#   Figure 1 :
#       - la solution approchée u_h ;
#       - les valeurs de la solution exacte u aux nœuds du maillage.
#
#   Figure 2 :
#       - log(erreur_L2) en fonction de log(h) ;
#       - log(erreur_H1) en fonction de log(h) ;
#       - log(erreur_Linf) en fonction de log(h).
#
# Fichiers d'entrée
# -----------------
# results/solution.dat :
#       x   y   u_numerique   u_exacte
#
# results/error_convergence.dat :
#       n   h   erreur_L2   erreur_H1   erreur_Linf
#
# Fichiers de sortie
# ------------------
# results/solution.png
# results/log_convergence.png
#
# Hypothèse sur le maillage
# ------------------------
# Le maillage est rectangulaire structuré.
# Les nœuds sont numérotés ligne par ligne :
# tous les x d'une ligne de y, puis ceux de la ligne suivante.
# =============================================================================

# Path facilite la construction des chemins vers les fichiers.
from pathlib import Path

# NumPy sert à lire les tableaux et à effectuer les opérations numériques.
import numpy as np

# Matplotlib sert à construire et enregistrer les figures.
import matplotlib.pyplot as plt

# Triangulation représente les coordonnées et la connectivité du maillage.
from matplotlib.tri import Triangulation


# =============================================================================
# Lecture générique d'un tableau numérique.
# =============================================================================
def read_table(path, expected_columns):
    """
    Lire un fichier numérique et vérifier sa structure.

    Paramètres
    ----------
    path : str ou pathlib.Path
        Chemin du fichier à lire.

    expected_columns : int
        Nombre de colonnes attendu.

    Retour
    ------
    data : numpy.ndarray
        Tableau à deux dimensions contenant les données.
    """

    # Accepter aussi bien un chemin texte qu'un objet Path.
    path = Path(path)

    # Vérifier que le chemin désigne un fichier existant.
    if not path.is_file():
        raise FileNotFoundError(
            f"Fichier introuvable : {path}"
        )

    # Un fichier présent mais vide ne peut pas être exploité.
    if path.stat().st_size == 0:
        raise ValueError(
            f"Le fichier est vide : {path}"
        )

    # ndmin=2 conserve une structure à deux dimensions.
    # Sans cela, un fichier d'une seule ligne pourrait produire
    # un tableau 1D, incompatible avec des expressions data[:, 0].
    data = np.loadtxt(path, ndmin=2)

    # Vérifier que le fichier contient effectivement des données.
    if data.shape[0] == 0:
        raise ValueError(
            f"Aucune donnée numérique dans : {path}"
        )

    # Contrôler le nombre de colonnes avant leur extraction.
    if data.shape[1] != expected_columns:
        raise ValueError(
            f"{path.name} : {expected_columns} colonnes attendues, "
            f"mais {data.shape[1]} colonnes obtenues."
        )

    # NaN et infini peuvent provenir d'un problème de calcul.
    # On arrête le script plutôt que de produire une figure trompeuse.
    if not np.isfinite(data).all():
        raise ValueError(
            f"Présence de valeurs NaN ou infinies dans : {path}"
        )

    return data


# =============================================================================
# Lecture des coordonnées et des solutions.
# =============================================================================
def read_solution(path):
    """
    Extraire les quatre colonnes du fichier solution.dat.

    Retour
    ------
    x, y : tableaux des coordonnées des nœuds.
    uh : solution approchée aux nœuds.
    u_exact : solution exacte évaluée aux mêmes nœuds.
    """

    data = read_table(path, expected_columns=4)

    # Les indices de colonnes Python commencent à zéro.
    x = data[:, 0]
    y = data[:, 1]
    uh = data[:, 2]
    u_exact = data[:, 3]

    return x, y, uh, u_exact


# =============================================================================
# Détection des dimensions du maillage.
# =============================================================================
def detect_mesh(x, y):
    """
    Déterminer nx et ny à partir des coordonnées.

    Vérifier également que le tableau respecte la numérotation
    utilisée par le module Fortran mesh_module.
    """

    # Pour une grille structurée :
    # nombre de positions distinctes selon x = nx + 1.
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    nx = x_unique.size - 1
    ny = y_unique.size - 1

    # Il faut au moins un rectangle dans chaque direction.
    if nx < 1 or ny < 1:
        raise ValueError(
            "Le maillage doit contenir au moins un intervalle "
            "dans chaque direction."
        )

    # Une grille complète contient le produit du nombre de positions
    # selon x et selon y.
    expected_nodes = (nx + 1) * (ny + 1)

    if x.size != expected_nodes:
        raise ValueError(
            "Les données ne correspondent pas à une grille "
            "rectangulaire structurée complète. "
            f"Nœuds attendus : {expected_nodes}, "
            f"nœuds obtenus : {x.size}."
        )

    # Construire les coordonnées dans l'ordre attendu :
    # x varie d'abord, puis y.
    expected_x, expected_y = np.meshgrid(
        x_unique,
        y_unique,
        indexing="xy"
    )

    # ravel transforme les grilles 2D en vecteurs ligne par ligne.
    # Cette vérification évite de connecter les mauvais sommets.
    if not (
        np.allclose(x, expected_x.ravel())
        and np.allclose(y, expected_y.ravel())
    ):
        raise ValueError(
            "Les nœuds ne sont pas numérotés ligne par ligne "
            "comme dans le module Fortran mesh_module."
        )

    return nx, ny


# =============================================================================
# Reconstruction de la connectivité triangulaire.
# =============================================================================
def build_triangulation(nx, ny):
    """
    Construire les deux triangles de chaque rectangle élémentaire.

    Le découpage est identique à celui du code Fortran.

    Attention :
    les indices Python commencent à 0, tandis que les indices
    Fortran du projet commencent à 1.
    """

    triangles = []

    # Parcourir les rectangles du maillage.
    for j in range(ny):
        for i in range(nx):

            # Sommet inférieur gauche.
            n00 = j * (nx + 1) + i

            # Sommet inférieur droit.
            n10 = n00 + 1

            # Sommet supérieur gauche.
            n01 = n00 + nx + 1

            # Sommet supérieur droit.
            n11 = n01 + 1

            # Premier triangle : bas-gauche, bas-droit, haut-gauche.
            triangles.append([n00, n10, n01])

            # Second triangle : bas-droit, haut-droit, haut-gauche.
            triangles.append([n10, n11, n01])

    # Chaque ligne contient les trois indices d'un triangle.
    return np.asarray(triangles, dtype=int)


# =============================================================================
# Sauvegarde d'une figure et contrôle du fichier produit.
# =============================================================================
def save_figure(fig, output_file):
    """Enregistrer une figure et vérifier que l'image est non vide."""

    output_file = Path(output_file)

    # Créer le dossier de sortie si nécessaire.
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Enregistrer explicitement la figure reçue.
    # bbox_inches="tight" limite les marges inutiles.
    fig.savefig(
        output_file,
        dpi=200,
        bbox_inches="tight"
    )

    if not output_file.is_file():
        raise RuntimeError(
            f"L'image n'a pas été créée : {output_file}"
        )

    if output_file.stat().st_size == 0:
        raise RuntimeError(
            f"L'image créée est vide : {output_file}"
        )

    # flush=True rend le message immédiatement visible dans les logs CI.
    print(
        f"Image créée : {output_file} "
        f"({output_file.stat().st_size} octets)",
        flush=True
    )


# =============================================================================
# Figure 1 : comparaison des solutions.
# =============================================================================
def plot_solution(x, y, uh, u_exact, nx, ny, output_file):
    """
    Représenter les deux solutions sur la même triangulation.

    Une échelle de couleurs commune permet une comparaison cohérente.
    """

    # Construire la liste des triangles.
    triangles = build_triangulation(nx, ny)

    # Associer coordonnées et connectivité.
    triangulation = Triangulation(x, y, triangles)

    # Utiliser les extrema des deux solutions pour la même échelle.
    vmin = min(float(uh.min()), float(u_exact.min()))
    vmax = max(float(uh.max()), float(u_exact.max()))

    # Les niveaux doivent être distincts même si les champs sont constants.
    if vmax == vmin:
        padding = max(1.0, abs(vmin)) * 1.0e-6
        vmin -= padding
        vmax += padding

    # Définir les mêmes niveaux pour les deux cartes.
    levels = np.linspace(vmin, vmax, 61)

    # num=1 désigne explicitement la première figure.
    fig, axes = plt.subplots(
        1,
        2,
        num=1,
        figsize=(13, 5),
        constrained_layout=True
    )

    fields = [uh, u_exact]

    titles = [
        r"Solution approchée $u_h$",
        r"Solution exacte $u$"
    ]

    contour = None

    # Tracer successivement les deux champs.
    for ax, values, title in zip(axes, fields, titles):

        # Carte de couleurs sur les triangles.
        contour = ax.tricontourf(
            triangulation,
            values,
            levels=levels,
            cmap="viridis"
        )

        # Afficher les arêtes seulement pour les maillages modestes.
        # Sur un maillage très fin, elles masqueraient les couleurs.
        if max(nx, ny) <= 40:
            ax.triplot(
                triangulation,
                color="black",
                linewidth=0.25,
                alpha=0.25
            )

        ax.set_title(title)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")

        # Conserver la même échelle physique sur les deux axes.
        ax.set_aspect("equal")

    # Une seule barre de couleurs suffit puisque les niveaux sont communs.
    colorbar = fig.colorbar(
        contour,
        ax=axes,
        shrink=0.9,
        pad=0.03
    )

    colorbar.set_label("Valeur de la solution")

    fig.suptitle(
        rf"Poisson 2D : solutions sur le maillage {nx} × {ny}"
    )

    save_figure(fig, output_file)

    return fig


# =============================================================================
# Lecture des données de convergence.
# =============================================================================
def read_convergence_data(path):
    """
    Lire les erreurs et trier les mesures du maillage grossier au plus fin.

    Format :
        n  h  erreur_L2  erreur_H1  erreur_Linf
    """

    data = read_table(path, expected_columns=5)

    n = data[:, 0]
    h = data[:, 1]

    # Extraire simultanément les trois colonnes d'erreurs.
    errors = data[:, 2:5]

    if np.any(n <= 0):
        raise ValueError(
            "Les nombres de subdivisions doivent être positifs."
        )

    # Les logarithmes ne sont définis ici que pour des valeurs positives.
    if np.any(h <= 0) or np.any(errors <= 0):
        raise ValueError(
            "Les pas h et les erreurs doivent être strictement "
            "positifs pour calculer leurs logarithmes."
        )

    # Un seul pas permet d'afficher un point, mais pas une convergence.
    if np.unique(h).size < 2:
        raise ValueError(
            "Le tableau de convergence doit contenir "
            "au moins deux pas de maillage différents."
        )

    # argsort trie par valeurs croissantes.
    # [::-1] inverse cet ordre pour partir du plus grand h.
    order = np.argsort(h)[::-1]

    return (
        n[order],
        h[order],
        errors[order, 0],
        errors[order, 1],
        errors[order, 2]
    )


# =============================================================================
# Figure 2 : logarithmes des erreurs.
# =============================================================================
def plot_log_convergence(
    h,
    error_l2,
    error_h1,
    error_linf,
    output_file
):
    """
    Tracer log(erreur) en fonction de log(h).

    Les erreurs sont celles calculées par le programme Fortran.
    Les droites de référence ne modifient pas les mesures.
    """

    # np.log calcule le logarithme naturel.
    log_h = np.log(h)

    fig, ax = plt.subplots(
        num=2,
        figsize=(10, 6),
        constrained_layout=True
    )

    # Associer chaque erreur à un symbole et à une légende.
    curves = [
        (error_l2, "o", r"$\log(E_{L^2})$"),
        (error_h1, "s", r"$\log(E_{H^1})$"),
        (error_linf, "^", r"$\log(E_{\infty})$")
    ]

    for error, marker, label in curves:
        ax.plot(
            log_h,
            np.log(error),
            marker=marker,
            linewidth=1.8,
            markersize=6,
            label=label
        )

    # Si E(h) ≈ C*h^p, alors :
    #
    # log(E(h)) ≈ log(C) + p*log(h).
    #
    # Une droite de pente p sert donc de repère visuel.
    # On aligne les références sur le dernier point, le plus fin.
    reference_l2 = (
        np.log(error_l2[-1])
        + 2.0 * (log_h - log_h[-1])
    )

    reference_h1 = (
        np.log(error_h1[-1])
        + (log_h - log_h[-1])
    )

    ax.plot(
        log_h,
        reference_l2,
        "--",
        color="black",
        linewidth=1.3,
        label="Repère : pente 2"
    )

    ax.plot(
        log_h,
        reference_h1,
        ":",
        color="gray",
        linewidth=1.8,
        label="Repère : pente 1"
    )

    # Les grands h sont à gauche, les petits h à droite.
    # Le raffinement du maillage se lit donc de gauche à droite.
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

    save_figure(fig, output_file)

    return fig


# =============================================================================
# Fonction principale.
# =============================================================================
def main():
    """Lire les données, construire les figures et les sauvegarder."""

    # __file__ est le chemin de ce script.
    # Les données sont donc cherchées relativement au script,
    # et non au dossier depuis lequel Python est lancé.
    project_dir = Path(__file__).resolve().parent
    results_dir = project_dir / "results"

    solution_file = results_dir / "solution.dat"
    convergence_file = results_dir / "error_convergence.dat"

    solution_image = results_dir / "solution.png"
    convergence_image = results_dir / "log_convergence.png"

    print(
        f"Dossier du projet : {project_dir}",
        flush=True
    )

    # Lire et valider toutes les données avant de construire les figures.
    x, y, uh, u_exact = read_solution(solution_file)

    # Le maillage est déduit des coordonnées, sans nx=20 codé en dur.
    nx, ny = detect_mesh(x, y)

    n, h, error_l2, error_h1, error_linf = (
        read_convergence_data(convergence_file)
    )

    print(
        f"Maillage de la figure 1 : {nx} × {ny}",
        flush=True
    )

    print(
        f"Nombre de mesures de convergence : {n.size}",
        flush=True
    )

    # Construire et sauvegarder la première figure.
    plot_solution(
        x,
        y,
        uh,
        u_exact,
        nx,
        ny,
        solution_image
    )

    # Construire et sauvegarder la deuxième figure.
    plot_log_convergence(
        h,
        error_l2,
        error_h1,
        error_linf,
        convergence_image
    )

    print(
        "Les deux figures ont été sauvegardées.",
        flush=True
    )

    # Avec le backend Agg, les images sont générées sans fenêtre.
    # Avec un backend interactif, afficher les deux figures.
    if str(plt.get_backend()).lower() != "agg":
        plt.show()

    # Libérer les ressources associées aux figures.
    plt.close("all")


# Exécuter main seulement lorsque ce fichier est lancé directement.
# Importer ce module depuis un autre script ne lance pas les tracés.
if __name__ == "__main__":
    main()