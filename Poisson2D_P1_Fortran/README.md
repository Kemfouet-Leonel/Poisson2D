# Poisson2D_P1_Fortran

## Présentation du projet

Ce projet résout numériquement l’équation de Poisson en deux dimensions par la méthode des éléments finis P1 sur un maillage triangulaire structuré.

Le solveur est développé en Fortran et organisé en modules afin de séparer les différentes étapes du calcul : paramètres, maillage, assemblage, résolution et calcul des erreurs.

Un script Python permet de représenter la solution approchée, les valeurs de la solution exacte et les courbes de convergence.

Un pipeline GitHub Actions automatise la compilation, l’exécution sur plusieurs maillages, les contrôles des données, la génération des figures et la publication des résultats.

- Dossier du solveur : `Poisson2D_P1_Fortran/`
- Script de visualisation : `plot_solution.py`
- Workflow : `.github/workflows/ci-cd.yml`
- Artefact de résultats : `poisson2d-results`

## Sommaire

1. Objectifs
2. Problème mathématique
3. Méthode numérique
4. Arborescence
5. Description des fichiers
6. Installation des dépendances
7. Compilation et exécution
8. Étude de convergence
9. Visualisation
10. Résultats et interprétation
11. Pipeline CI/CD
12. Récupération des artefacts
13. Nettoyage et versionnement

## 1. Objectifs

Le projet vise à :

- résoudre un problème de Poisson 2D avec second membre non nul en utilisant les éléments finis $P_1$;
- construire un maillage triangulaire du carré unité ;
- assembler la matrice de rigidité et le second membre ;
- imposer des conditions de Dirichlet homogènes ;
- résoudre le système linéaire avec LAPACK ;
- comparer la solution approchée à une solution analytique connue ;
- étudier l’effet du raffinement sur les erreurs en normes $L^2, H^1, L^{\infty}$ ;
- tracer log(erreur) en fonction de log(h) afin d'observer l'ordre de convergence ;
- automatiser la chaîne de calcul avec GitHub Actions ;
- conserver les résultats sous forme d’artefacts.

## 2. Problème mathématique

### 2.1 Équation de Poisson

On considère le domaine :

\[
\Omega=(0,1)\times(0,1).
\]

Le problème consiste à chercher \(u\) telle que :

\[
-\Delta u=f \quad \text{dans }\Omega,
\]

avec :

\[
u=0 \quad \text{sur }\partial\Omega.
\]

### 2.2 Solution exacte

La solution analytique choisie est :

\[
u(x,y)=\sin(\pi x)\sin(\pi y).
\]

Elle satisfait les conditions de Dirichlet homogènes sur les quatre côtés du carré.

### 2.3 Second membre

Le second membre associé est :

\[
f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y).
\]

Cette expression vérifie :

\[
-\Delta u=f.
\]

### 2.4 Gradient exact

Le gradient de la solution est :

\[
\nabla u(x,y)=
\begin{pmatrix}
\pi\cos(\pi x)\sin(\pi y)\\
\pi\sin(\pi x)\cos(\pi y)
\end{pmatrix}.
\]

La connaissance de ce gradient permet de comparer les dérivées de la solution approchée à celles de la solution exacte.

## 3. Méthode numérique

### 3.1 Maillage

Le carré unité est subdivisé en rectangles, chacun étant découpé en deux triangles.

Pour un maillage carré :

\[
n_x=n_y=n.
\]

Le pas directionnel est :

\[
h=\frac{1}{n}.
\]

Le nombre de nœuds est :

\[
N=(n+1)^2.
\]

Le nombre de triangles est :

\[
N_T=2n^2.
\]

Le diamètre maximal des triangles vaut ici \(\sqrt{2}/n\). Le choix du paramètre \(h=1/n\) conserve les mêmes ordres de convergence.

### 3.2 Formulation variationnelle

On cherche \(u\in H_0^1(\Omega)\) telle que :

\[
\int_\Omega \nabla u\cdot\nabla v\,d\Omega
=
\int_\Omega fv\,d\Omega
\qquad
\forall v\in H_0^1(\Omega).
\]

### 3.3 Approximation P1

La solution approchée s’écrit :

\[
u_h(x,y)=\sum_{j=1}^{N}U_j\phi_j(x,y),
\]

où :

- \(U_j\) est la valeur approchée au nœud \(j\) ;
- \(\phi_j\) est une fonction de base affine par triangle.

L’assemblage conduit au système :

\[
AU=b.
\]

Les coefficients de rigidité sont :

\[
A_{ij}
=
\int_\Omega \nabla\phi_j\cdot\nabla\phi_i\,d\Omega.
\]

Le second membre est :

\[
b_i=\int_\Omega f\phi_i\,d\Omega.
\]

### 3.4 Assemblage élémentaire

Pour un triangle \(K\), les gradients des fonctions de base P1 sont constants :

\[
A_{ij}^{K}
=
|K|\,\nabla\phi_j\cdot\nabla\phi_i.
\]

Dans la version initiale du code, le second membre est approché au barycentre :

\[
b_i^K
\approx
\frac{|K|}{3}f(x_K,y_K).
\]

Cette quadrature est une approximation pour le second membre trigonométrique utilisé.

### 3.5 Conditions de Dirichlet

La version initiale impose les valeurs nulles au bord en :

1. annulant les lignes associées aux nœuds de bord ;
2. plaçant 1 sur leur diagonale ;
3. plaçant 0 dans le second membre correspondant.

Attention : cette modification des lignes seules ne conserve généralement pas la symétrie de la matrice globale.

La matrice de rigidité est structurellement creuse, mais l’implémentation actuelle la stocke dans un tableau dense.

## 4. Arborescence du dépôt

```text
Poisson2D/
├── .github/
│   └── workflows/
│       └── ci-cd.yml
│
├── .gitignore
│
└── Poisson2D_P1_Fortran/
    ├── README.md
    ├── Makefile
    ├── main.f90
    ├── params_module.f90
    ├── mesh_module.f90
    ├── exact_solution_module.f90
    ├── fem_module.f90
    ├── solver_module.f90
    ├── error_module.f90
    ├── plot_solution.py
    ├── bin/
    │   └── poisson2d
    └── results/
        ├── solution.dat
        ├── error_convergence.dat
        ├── solution.png
        └── log_convergence.png
```

Les dossiers `bin/` et `results/` contiennent des fichiers générés.

## 5. Description des fichiers

| Fichier | Rôle |
|---|---|
| `README.md` | Documentation du projet et des procédures |
| `Makefile` | Compilation, exécution selon les cibles définies et nettoyage |
| `.gitignore` | Exclusion des fichiers générés du suivi Git |
| `params_module.f90` | Paramètres numériques et dimensions |
| `mesh_module.f90` | Coordonnées, connectivité et frontière |
| `exact_solution_module.f90` | Solution exacte, source et gradient exact |
| `fem_module.f90` | Assemblage et conditions de Dirichlet |
| `solver_module.f90` | Résolution du système avec LAPACK |
| `error_module.f90` | Mesures d’erreur |
| `main.f90` | Orchestration du calcul |
| `plot_solution.py` | Lecture des résultats et génération des figures |
| `.github/workflows/ci-cd.yml` | Automatisation GitHub Actions |

### 5.1 `params_module.f90`

Ce module contient notamment :

```fortran
dp
nx
ny
lx
ly
n_nodes
n_elements
n_dof
```

La procédure :

```fortran
call set_parameters(n, n)
```

initialise les paramètres du maillage carré.

### 5.2 `mesh_module.f90`

Les principaux tableaux sont :

```fortran
nodes(2, n_nodes)
elements(3, n_elements)
boundary_nodes(n_boundary_nodes)
```

Ils contiennent respectivement :

- les coordonnées ;
- les trois sommets de chaque triangle ;
- les numéros des nœuds de bord.

Les principales procédures sont :

```fortran
build_mesh()
node_index(i, j)
build_boundary()
```

### 5.3 `exact_solution_module.f90`

Ce module fournit :

```fortran
exact_solution(x, y)
source_term(x, y)
exact_gradient(x, y, grad_u)
```

Ces fonctions définissent les données analytiques nécessaires à l’assemblage et à la comparaison des résultats.

### 5.4 `fem_module.f90`

Les tableaux principaux sont :

```fortran
global_matrix
rhs
```

La procédure :

```fortran
call assemble_system()
```

réalise :

1. l’allocation et l’initialisation ;
2. le parcours des triangles ;
3. le calcul des contributions locales ;
4. l’assemblage ;
5. l’application des conditions de Dirichlet.

### 5.5 `solver_module.f90`

La procédure :

```fortran
call solve_linear_system(global_matrix, rhs, solution)
```

Résoudre du système dense par la factorisation LU.


### 5.6 `error_module.f90`

L’appel est :

```fortran
call compute_errors(solution, error_l2, error_h1, error_linf)
```

Les sorties sont utilisées dans le tableau de convergence.

Leur interprétation doit correspondre à l’implémentation :

- `error_l2` : estimation de l’erreur L2 ;
- `error_h1` : erreur de gradient ou norme H1 selon le calcul retenu ;
- `error_linf` : maximum échantillonné dans la version initiale.

### 5.7 `main.f90`

Le programme principal :

1. lit le paramètre `n` ;
2. initialise le problème ;
3. construit le maillage ;
4. assemble le système ;
5. résout ;
6. calcule les erreurs ;
7. écrit la solution ;
8. ajoute une ligne au tableau de convergence ;
9. affiche les informations dans le terminal.

### 5.8 `plot_solution.py`

| Fonction | Rôle |
|---|---|
| `read_table` | Lire et vérifier un tableau |
| `read_solution` | Extraire les coordonnées et les solutions |
| `detect_mesh` | Déterminer les dimensions et vérifier la numérotation |
| `build_triangulation` | Reconstruire les triangles |
| `save_figure` | Sauvegarder et contrôler l’image |
| `plot_solution` | Construire la figure 1 |
| `read_convergence_data` | Lire et trier les mesures d’erreur |
| `plot_log_convergence` | Construire la figure 2 |
| `main` | Coordonner les opérations |

Les chemins sont construits depuis le dossier du script, et non depuis le répertoire courant du terminal.

## 6. Installation - Dépendances système

```bash
sudo apt update

sudo apt install -y \
    git \
    gfortran \
    liblapack-dev \
    python3 
```

## 7. Compilation et exécution

### 7.1 Compilation

Depuis `Poisson2D_P1_Fortran/` :

```bash
make
```

L’exécutable attendu est :

```text
bin/poisson2d
```

### 7.2 Exécution pour un maillage

```bash
mkdir -p results
./bin/poisson2d 20
```

Pour fixer les paramètres :

\[
n_x=n_y=20.
\]

La commande directe ci-dessus est la référence pour lancer un maillage précis. Le comportement de `make run` dépend des commandes présentes dans la cible `run` du Makefile.

### 7.3 Compilation depuis la racine du dépôt

Depuis `Poisson2D/` :

```bash
make -C Poisson2D_P1_Fortran
```

### 7.4 Ordre des sources

Les modules doivent être compilés avant les sources qui les utilisent.

Ordre utilisé :

```text
params_module.f90
mesh_module.f90
exact_solution_module.f90
fem_module.f90
solver_module.f90
error_module.f90
main.f90
```

## 8. Étude de la convergence

### 8.1 Maillages retenus pour la CI

L’automatisation utilise sept maillages :

```text
n = 2, 4, 6, 8, 10, 15, 20
```

| n | h | Nœuds | Triangles |
|---:|---:|---:|---:|
| 2 | 0.5 | 9 | 8 |
| 4 | 0.25 | 25 | 32 |
| 6 | 0.166667 | 49 | 72 |
| 8 | 0.125 | 81 | 128 |
| 10 | 0.1 | 121 | 200 |
| 15 | 0.066667 | 256 | 450 |
| 20 | 0.05 | 441 | 800 |

### 8.3 Enregistrement des erreurs

Le fichier :

```text
results/error_convergence.dat
```

contient :

```text
n   h   erreur_L2   erreur_H1   erreur_Linf
```

Le programme doit ajouter une ligne à chaque lancement :

```fortran
open(newunit=unit_out, file="results/error_convergence.dat", &
     status="unknown", position="append")

write(unit_out, *) n, 1.0_dp / real(n, dp), &
     error_l2, error_h1, error_linf

close(unit_out)
```

### 8.4 Enregistrement de la solution

Le fichier :

```text
results/solution.dat
```

contient :

```text
x   y   u_numerique   u_exacte
```

Il est remplacé à chaque lancement et correspond donc au dernier maillage calculé.

## 9. Visualisation

### 9.1 Lancement

```bash
python plot_solution.py
```

### 9.2 Mode sans fenêtre

```bash
MPLBACKEND=Agg python plot_solution.py
```

Ce mode est utilisé dans GitHub Actions.

### 9.3 Figure 1 : solutions

Fichier :

```text
results/solution.png
```

La figure contient :

- à gauche : la solution approchée ;
- à droite : les valeurs de la solution exacte.

### 9.4 Figure 2 : logarithmes des erreurs

Fichier :

```text
results/log_convergence.png
```

Les courbes représentent :

\[
\log(erreur)
 \;\ \text{en fonction de} \;\
\log(h).
\]

## 10. Résultats et interprétation

### 10.1 Résultats disponibles

Les résultats numériques effectivement calculés se trouvent dans :

```text
results/error_convergence.dat
```

Les images sont :

```text
results/solution.png
results/log_convergence.png
```

Après une exécution GitHub Actions ayant atteint l’archivage, ces fichiers sont disponibles dans l’artefact `poisson2d-results`.

### 10.2 Résultats théoriques de référence

Pour le problème régulier considéré, les estimations usuelles d’une approximation P1 correctement mise en œuvre donnent :

\[
\|u-u_h\|_{L^2(\Omega)}=O(h^2),
\]

\[
|u-u_h|_{H^1(\Omega)}=O(h).
\]

### 10.3 Ordre expérimental

Pour deux mesures :

\[
p=
\frac{\log(E_2/E_1)}
{\log(h_2/h_1)}.
\]

Si :

\[
E(h)\approx Ch^p,
\]

alors :

\[
\log(E(h))\approx\log(C)+p\log(h).
\]

La pente du graphique fournit donc une estimation de l’ordre.

## 11. Pipeline CI/CD

### 11.1 Objectif

Le pipeline automatise la chaîne :

```text
Récupération du dépôt
        ↓
Vérification des sources
        ↓
Installation des dépendances
        ↓
Compilation Fortran
        ↓
Exécution de sept maillages
        ↓
Contrôle des données
        ↓
Génération des figures
        ↓
Contrôle des images
        ↓
Archivage des résultats
        ↓
Release si push d'un tag v*
```

La partie publication distribue les résultats. Elle ne déploie pas une application web ou un service.

### 11.2 Emplacement

Le fichier est :

```text
.github/workflows/ci-cd.yml
```

à la racine du dépôt.

### 11.3 Déclencheurs

| Événement | Calcul et contrôles | Release |
|---|---|---|
| Push sur `main` | Oui | Non |
| Pull request vers `main` | Oui | Non |
| Lancement manuel | Oui | Non |
| Push d’un tag `v*` | Oui | Oui après réussite |

### 11.4 Jobs

Le workflow contient deux jobs :

- `build-test-visualize` : compile, exécute, vérifie et archive ;
- `release` : télécharge les résultats puis les publie.

Le second dépend de la réussite du premier.

### 11.5 Portée des contrôles

Le pipeline vérifie :

- les sources attendues ;
- la compilation ;
- l’exécutable ;
- les sept exécutions ;
- les dimensions des tableaux ;
- l’absence de NaN et d’infinis ;
- la positivité nécessaire aux logarithmes ;
- la correspondance entre les subdivisions et les pas ;
- la création et la lisibilité des figures.

Il ne valide pas encore automatiquement les pentes théoriques.

## 12. Récupération des résultats dans l'artefact

### 12.1 Où sont les résultats ?

Après une exécution locale :

```text
Poisson2D_P1_Fortran/results/
```

Après une exécution CI, les fichiers sont produits sur le runner GitHub, puis archivés sous le nom :

```text
poisson2d-results
```

### 12.2 Contenu

L’archive contient :

```text
solution.dat
error_convergence.dat
solution.png
log_convergence.png
```

L’exécutable Fortran n’est pas inclus dans l’artefact actuel.

## 16. Publication d’une release

### 16.1 Créer un tag

Depuis la racine du dépôt :

```bash
git switch main
git pull --ff-only origin main
```

Créer un nouveau tag :

```bash
git tag -a v1.0.0 \
    -m "Solveur Poisson 2D P1 avec chaîne CI/CD"

git push origin v1.0.0
```

Si ce tag existe déjà, choisir une nouvelle version.

## 13. Nettoyage et versionnement

### 13.1 Nettoyage

Depuis le dossier du solveur :

```bash
make clean
```

La cible utilisée supprime notamment :

```text
bin/
*.mod
*.o
results/*.dat
results/*.png
```

Attention : les résultats locaux sont supprimés. Les sauvegarder avant nettoyage.

### 13.2 Fichiers à ignorer

Exemple de `.gitignore` à la racine :

```gitignore
# Compilation
*.mod
*.o
Poisson2D_P1_Fortran/bin/

# Résultats générés
Poisson2D_P1_Fortran/results/*.dat
Poisson2D_P1_Fortran/results/*.png

# Python
Poisson2D_P1_Fortran/.venv/
__pycache__/
*.pyc
```

Les fichiers ignorés par Git peuvent tout de même être archivés par GitHub Actions.

Ajouter une règle ne retire pas automatiquement les fichiers déjà suivis.

### 13.3 Travail par branche

Le dépôt existe déjà : ne pas exécuter un nouveau `git init` dans le sous-dossier du solveur.

Exemple :

```bash
git switch -c docs/project-readme

git add Poisson2D_P1_Fortran/README.md
git commit -m "docs: documenter le projet et le pipeline"

git push -u origin docs/project-readme
```

Ouvrir ensuite une pull request vers `main`.

