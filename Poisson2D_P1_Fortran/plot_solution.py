# =============================================================================
# Script : plot_log_convergence.py
#
# Rôle :
#   Tracer le logarithme des erreurs en fonction du logarithme du pas h :
#
#       log(erreur) = f(log(h))
#
#   Ce type de graphique permet de visualiser l’ordre de convergence.
#
#   Pour une méthode éléments finis P1, on attend :
#
#       log(erreur_L2) ~ 2 log(h) + constante
#       log(erreur_H1) ~ 1 log(h) + constante
#
#   La pente de la droite donne donc l’ordre de convergence.
# =============================================================================

import os
import numpy as np
import matplotlib.pyplot as plt


def read_convergence_data(path):
    """
    Lit le fichier de convergence généré par le code Fortran.

    Format attendu :
        n   h   erreur_L2   erreur_H1   erreur_Linf

    Paramètre
    ---------
    path : str
        Chemin vers results/error_convergence.dat.

    Retour
    ------
    n, h, error_l2, error_h1, error_linf : numpy.ndarray
        Données de convergence.
    """
    data = np.loadtxt(path)

    n = data[:, 0]
    h = data[:, 1]
    error_l2 = data[:, 2]
    error_h1 = data[:, 3]
    error_linf = data[:, 4]

    return n, h, error_l2, error_h1, error_linf


def plot_log_convergence(h, error_l2, error_h1, error_linf, output_file):
    """
    Trace log(erreur) en fonction de log(h).

    Paramètres
    ----------
    h : numpy.ndarray
        Pas du maillage.

    error_l2 : numpy.ndarray
        Erreurs en norme L2.

    error_h1 : numpy.ndarray
        Erreurs en norme H1.

    error_linf : numpy.ndarray
        Erreurs en norme L_infini.

    output_file : str
        Chemin de l’image générée.
    """
    # Calcul des logarithmes naturels.
    log_h = np.log(h)
    log_error_l2 = np.log(error_l2)
    log_error_h1 = np.log(error_h1)
    log_error_linf = np.log(error_linf)

    # Création de la figure.
    plt.figure(figsize=(9, 6))

    # Courbe log(erreur_L2) en fonction de log(h).
    plt.plot(
        log_h,
        log_error_l2,
        marker="o",
        linewidth=2,
        label=r"$\log(\|u-u_h\|_{L^2})$"
    )

    # Courbe log(erreur_H1) en fonction de log(h).
    plt.plot(
        log_h,
        log_error_h1,
        marker="s",
        linewidth=2,
        label=r"$\log(\|u-u_h\|_{H^1})$"
    )

    # Courbe log(erreur_Linf) en fonction de log(h).
    plt.plot(
        log_h,
        log_error_linf,
        marker="^",
        linewidth=2,
        label=r"$\log(\|u-u_h\|_{L^\infty})$"
    )

    # Droite théorique de pente 2 pour la norme L2.
    # On l’aligne sur le premier point de la courbe L2.
    slope_2 = 2.0
    offset_l2 = log_error_l2[0] - slope_2 * log_h[0]

    plt.plot(
        log_h,
        slope_2 * log_h + offset_l2,
        linestyle="--",
        color="black",
        label=r"Référence : pente 2"
    )

    # Droite théorique de pente 1 pour la norme H1.
    # On l’aligne sur le premier point de la courbe H1.
    slope_1 = 1.0
    offset_h1 = log_error_h1[0] - slope_1 * log_h[0]

    plt.plot(
        log_h,
        slope_1 * log_h + offset_h1,
        linestyle=":",
        color="gray",
        label=r"Référence : pente 1"
    )

    # Inversion de l’axe x :
    # à gauche, maillage grossier ;
    # à droite, maillage fin.
    plt.gca().invert_xaxis()

    # Habillage du graphique.
    plt.xlabel(r"$\log(h)$")
    plt.ylabel(r"$\log(\mathrm{erreur})$")
    plt.title("Convergence logarithmique de la méthode éléments finis $P_1$")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    # Création du dossier results s’il n’existe pas.
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Sauvegarde de la figure.
    plt.savefig(output_file, dpi=300)

    # Affichage à l’écran.
    plt.show()


def main():
    """
    Fonction principale :
    lit les données de convergence et génère la figure logarithmique.
    """
    # Fichier généré par le programme Fortran.
    input_file = "results/error_convergence.dat"

    # Fichier image de convergence logarithmique.
    output_file = "results/log_convergence.png"

    # Lecture des données.
    n, h, error_l2, error_h1, error_linf = read_convergence_data(input_file)

    # Tracé de la convergence logarithmique.
    plot_log_convergence(h, error_l2, error_h1, error_linf, output_file)


# Point d’entrée du script.
if __name__ == "__main__":
    main()