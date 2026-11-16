# Poisson2D_P1_Fortran

Ce projet résout numériquement l’équation de Poisson en deux dimensions par la méthode des éléments finis \(P_1\), en langage Fortran 90. Il permet également d’étudier l’influence de la finesse du maillage sur l’erreur numérique et de visualiser la solution ainsi que la convergence à l’aide de Python.

Le problème résolu est le suivant :

\[
-\Delta u = f \quad \text{dans } \Omega = (0,1) \times (0,1),
\]

avec des conditions de Dirichlet homogènes sur toute la frontière :

\[
u = 0 \quad \text{sur } \partial\Omega.
\]

Le second membre est choisi non nul :

\[
f(x,y) = 2\pi^2 \sin(\pi x)\sin(\pi y).
\]

La solution exacte associée est :

\[
u(x,y) = \sin(\pi x)\sin(\pi y).
\]

## Objectifs du projet

Ce projet a pour objectifs de :

- résoudre l’équation de Poisson 2D par la méthode des éléments finis \(P_1\) ;
- construire un maillage triangulaire structuré du carré unité ;
- assembler la matrice de rigidité globale et le second membre global ;
- imposer les conditions de Dirichlet homogènes ;
- résoudre le système linéaire \(A U = b\) avec LAPACK ;
- calculer les erreurs en normes \(L^2\), \(H^1\) et \(L^\infty\) ;
- étudier la convergence de la méthode pour plusieurs maillages ;
- visualiser la solution approchée et la solution exacte en Python ;
- tracer \(\log(\text{erreur})\) en fonction de \(\log(h)\) afin d’observer l’ordre de convergence.

## Méthode numérique

Le domaine \(\Omega = (0,1)^2\) est discrétisé par un maillage triangulaire structuré. Chaque rectangle élémentaire est divisé en deux triangles.

Pour un maillage de taille \(n \times n\), on a :

\[
h = \frac{1}{n}.
\]

Le nombre de nœuds est :

\[
N_{\text{nœuds}} = (n+1)^2.
\]

Le nombre de triangles est :

\[
N_{\text{triangles}} = 2n^2.
\]

La formulation variationnelle du problème consiste à chercher \(u_h \in V_h\) telle que :

\[
\int_\Omega \nabla u_h \cdot \nabla v_h \, dx
=
\int_\Omega f v_h \, dx,
\]

pour toute fonction test \(v_h \in V_h\).

Avec des éléments finis \(P_1\), la solution approchée s’écrit :

\[
u_h(x,y) = \sum_{i=1}^{N} U_i \phi_i(x,y),
\]

où \(\phi_i\) sont les fonctions de base associées aux sommets du maillage.

Le système linéaire obtenu est :

\[
A U = b.
\]

La matrice \(A\) est la matrice de rigidité globale. Elle est symétrique, creuse et définie positive après application des conditions de Dirichlet.

## Arborescence du projet

```text
Poisson2D_Fortran/
├── README.md
├── Makefile
├── .gitignore
├── main.f90
├── params_module.f90
├── mesh_module.f90
├── exact_solution_module.f90
├── fem_module.f90
├── solver_module.f90
├── error_module.f90
├── plot_solution_convergence.py
├── bin/
│   └── poisson2d
└── results/
    ├── solution.dat
    ├── error_convergence.dat
    ├── solution.png
    └── log_convergence.png
```

## Rôle de chaque fichier

| Fichier | Rôle |
|---|---|
| `README.md` | Documentation complète du projet |
| `Makefile` | Automatise la compilation, l’exécution, l’étude de convergence et le nettoyage |
| `.gitignore` | Ignore les fichiers générés automatiquement, tels que les fichiers binaires, `.mod`, `.o` et les résultats |
| `main.f90` | Programme principal : lit les paramètres, construit le maillage, assemble et résout le système, calcule les erreurs et écrit les résultats |
| `params_module.f90` | Définit les paramètres du problème : dimensions du maillage, dimensions du domaine, nombre de nœuds, triangles et degrés de liberté |
| `mesh_module.f90` | Construit les nœuds, les triangles et la liste des nœuds situés sur la frontière |
| `exact_solution_module.f90` | Fournit la solution exacte, le second membre \(f\) et le gradient exact de la solution |
| `fem_module.f90` | Assemble la matrice de rigidité globale, le second membre global et applique les conditions de Dirichlet |
| `solver_module.f90` | Résout le système linéaire \(A U = b\) à l’aide de LAPACK |
| `error_module.f90` | Calcule les erreurs en normes \(L^2\), \(H^1\) et \(L^\infty\) |
| `plot_solution_convergence.py` | Visualise la solution et trace la convergence logarithmique des erreurs |

## Description détaillée des modules Fortran

### `params_module.f90`

Ce module contient les paramètres globaux du problème.

Il définit notamment :

- la précision des réels : `dp = kind(1.0d0)` ;
- le nombre de subdivisions selon \(x\) et \(y\) : `nx` et `ny` ;
- les dimensions du domaine : `lx` et `ly` ;
- le nombre de nœuds : `n_nodes` ;
- le nombre de triangles : `n_elements` ;
- le nombre de degrés de liberté : `n_dof`.

La subroutine principale est :

```fortran
call set_parameters(nx, ny)
```

Elle initialise tous les paramètres nécessaires à la construction du maillage et à la résolution du problème.

### `mesh_module.f90`

Ce module construit le maillage triangulaire structuré du carré unité.

Il fournit :

- les coordonnées des nœuds dans le tableau `nodes` ;
- la connectivité des triangles dans le tableau `elements` ;
- la liste des nœuds de bord dans le tableau `boundary_nodes`.

La subroutine principale est :

```fortran
call build_mesh()
```

Elle construit :

1. les coordonnées des nœuds ;
2. les triangles ;
3. les nœuds situés sur la frontière.

Les nœuds sont numérotés ligne par ligne. Pour un nœud d’indices \((i,j)\), son numéro global est :

\[
\text{index} = j(n+1) + i + 1.
\]

### `exact_solution_module.f90`

Ce module contient les données analytiques du problème.

Il fournit :

- la solution exacte :

\[
u(x,y) = \sin(\pi x)\sin(\pi y);
\]

- le second membre :

\[
f(x,y) = 2\pi^2\sin(\pi x)\sin(\pi y);
\]

- le gradient exact :

\[
\nabla u =
\begin{pmatrix}
\pi\cos(\pi x)\sin(\pi y) \\
\pi\sin(\pi x)\cos(\pi y)
\end{pmatrix}.
\]

Le gradient exact est utilisé pour calculer l’erreur en norme \(H^1\).

### `fem_module.f90`

Ce module construit le système linéaire éléments finis.

Pour chaque triangle \(K\), il calcule :

- la matrice élémentaire de rigidité :

\[
A^K_{ij}
=
\int_K \nabla\phi_i \cdot \nabla\phi_j \, dx;
\]

- le second membre élémentaire :

\[
b^K_i
=
\int_K f\phi_i \, dx.
\]

Ces contributions sont ensuite assemblées dans la matrice globale \(A\) et le vecteur global \(b\).

La subroutine principale est :

```fortran
call assemble_system()
```

Elle réalise les étapes suivantes :

1. initialisation de la matrice et du second membre ;
2. parcours de tous les triangles ;
3. calcul des matrices élémentaires ;
4. assemblage global ;
5. application des conditions de Dirichlet.

Les conditions de Dirichlet homogènes sont appliquées en annulant les lignes correspondant aux nœuds de bord, en plaçant 1 sur la diagonale et en mettant le second membre associé à zéro.

### `solver_module.f90`

Ce module résout le système linéaire :

\[
A U = b.
\]

Il utilise la routine LAPACK `dgesv`, qui résout un système linéaire dense en double précision par factorisation \(LU\) avec pivotage partiel.

La subroutine principale est :

```fortran
call solve_linear_system(global_matrix, rhs, solution)
```

Après l’appel, le tableau `solution` contient les valeurs de la solution approchée aux nœuds du maillage.

### `error_module.f90`

Ce module calcule l’erreur entre la solution numérique \(u_h\) et la solution exacte \(u\).

Trois normes sont calculées :

- la norme \(L^2\) :

\[
\|u-u_h\|_{L^2};
\]

- la norme \(H^1\) :

\[
\|u-u_h\|_{H^1};
\]

- la norme \(L^\infty\) :

\[
\|u-u_h\|_{L^\infty}.
\]

La subroutine principale est :

```fortran
call compute_errors(solution, error_l2, error_h1, error_linf)
```

Les erreurs sont calculées triangle par triangle, puis assemblées globalement.

### `main.f90`

Le programme principal orchestre toutes les étapes de la résolution.

Il réalise successivement :

1. la lecture du paramètre \(n\) depuis la ligne de commande ;
2. l’initialisation des paramètres ;
3. la construction du maillage ;
4. l’assemblage du système linéaire ;
5. la résolution du système ;
6. le calcul des erreurs ;
7. l’écriture de la solution dans `results/solution.dat` ;
8. l’ajout des erreurs dans `results/error_convergence.dat` ;
9. l’affichage des résultats dans le terminal.

Le maillage est carré :

\[
nx = ny = n.
\]

Le pas du maillage est :

\[
h = \frac{1}{n}.
\]

## Installation des dépendances pour la compilation du projet

```bash
sudo apt update
sudo apt install gfortran liblapack-dev libblas-dev python3 python3-pip
pip install numpy matplotlib
```

## Compilation du projet

Place-toi dans le dossier du projet, Puis compile :

```bash
make
```

Cette commande crée l’exécutable :

```text
bin/poisson2d
```

## Exécution du programme

Pour exécuter le programme avec un maillage \(20 \times 20\) :

```bash
make run
```

Tu peux aussi exécuter directement l’exécutable :

```bash
mkdir -p results
./bin/poisson2d 20
```

Le paramètre `20` correspond au nombre de subdivisions du maillage selon \(x\) et \(y\).

## Étude de convergence

Pour étudier l’effet de la finesse du maillage sur l’erreur, exécute :

```bash
make convergence
```

Cette commande teste automatiquement plusieurs maillages :

```text
n = 2
n = 5
n = 10
n = 20
n = 40
n = 80
n = 160
n = 320
```

Pour chaque maillage, le programme :

- construit un maillage plus fin ;
- résout le problème ;
- calcule les erreurs ;
- ajoute une ligne dans `results/error_convergence.dat`.

Le fichier `error_convergence.dat` contient les colonnes suivantes :

```text
n   h   erreur_L2   erreur_H1   erreur_Linf
```

## Visualisation des résultats

Après avoir exécuté le programme Fortran, lance le script Python :

```bash
python3 plot_solution_convergence.py
```

Ce script génère deux figures :

```text
results/solution.png
results/log_convergence.png
```

### Figure 1 : solution approchée et solution exacte

La première figure contient deux sous-graphiques :

- à gauche, la solution approchée \(u_h\) ;
- à droite, la solution exacte \(u\).

Les deux solutions sont représentées en couleur à l’aide de la colormap `viridis`. Les arêtes du maillage triangulaire sont également affichées.

### Figure 2 : convergence logarithmique

La deuxième figure trace :

\[
\log(\text{erreur})
\quad \text{en fonction de} \quad
\log(h).
\]

Elle contient les courbes suivantes :

- \(\log(\|u-u_h\|_{L^2})\) ;
- \(\log(\|u-u_h\|_{H^1})\) ;
- \(\log(\|u-u_h\|_{L^\infty})\) ;
- une droite de référence de pente 2 ;
- une droite de référence de pente 1.

## Résultats attendus

Pour une méthode éléments finis \(P_1\), on attend théoriquement :

\[
\|u-u_h\|_{L^2} \approx C h^2,
\]

\[
\|u-u_h\|_{H^1} \approx C h.
\]

Ainsi, lorsque le maillage est raffiné, c’est-à-dire lorsque \(h\) diminue, l’erreur doit diminuer.

| Maillage | Pas \(h\) | Comportement attendu |
|---:|---:|---|
| \(2 \times 2\) | 0.5 | Erreur élevée |
| \(5 \times 5\) | 0.2 | Erreur plus faible |
| \(10 \times 10\) | 0.1 | Erreur encore plus faible |
| \(20 \times 20\) | 0.05 | Erreur nettement réduite |
| \(40 \times 40\) | 0.025 | Erreur très faible |
| \(80 \times 80\) | 0.0125 | Erreur très petite |
| \(160 \times 160\) | 0.00625 | Convergence visible |
| \(320 \times 320\) | 0.003125 | Erreur minimale |

## Interprétation des résultats

Les résultats obtenus doivent montrer que :

- plus le maillage est fin, plus l’erreur est faible ;
- l’erreur en norme \(L^2\) décroît approximativement comme \(h^2\) ;
- l’erreur en norme \(H^1\) décroît approximativement comme \(h\) ;
- la solution approchée se rapproche visuellement de la solution exacte lorsque le maillage est raffiné ;
- la méthode éléments finis \(P_1\) converge donc vers la solution exacte.

Sur le graphique logarithmique, la pente de la courbe donne l’ordre de convergence :

\[
p = \frac{\log(E_2) - \log(E_1)}
{\log(h_2) - \log(h_1)}.
\]

Pour la norme \(L^2\), la pente attendue est proche de 2.  
Pour la norme \(H^1\), la pente attendue est proche de 1.

## Nettoyage du projet

Pour supprimer les fichiers générés par la compilation et l’exécution :

```bash
make clean
```

Cette commande supprime :

```text
bin/
*.mod
*.o
results/*.dat
results/*.png
```


## Conclusion

Ce projet met en œuvre une méthode éléments finis \(P_1\) pour résoudre l’équation de Poisson 2D avec second membre non nul et conditions de Dirichlet homogènes.

Le code est organisé en modules afin de séparer clairement :

- les paramètres du problème ;
- la construction du maillage ;
- la solution exacte ;
- l’assemblage éléments finis ;
- la résolution du système linéaire ;
- le calcul des erreurs ;
- la visualisation des résultats.

L’étude de la convergence permet de vérifier numériquement que plus le maillage est fin, plus l’erreur diminue. Les résultats obtenus confirment le comportement théorique attendu de la méthode éléments finis \(P_1\) : convergence d’ordre 2 en norme \(L^2\) et d’ordre 1 en norme \(H^1\).