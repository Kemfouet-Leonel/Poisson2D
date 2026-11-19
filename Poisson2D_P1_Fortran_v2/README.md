# Poisson2D_P1_Fortran_v2

Résolution numérique de l’équation de Poisson en deux dimensions par la méthode des éléments finis P1, avec stockage creux CSR, solveur itératif, validation des erreurs, tests de convergence, mesures de performance et automatisation GitHub Actions.

## Présentation

Ce projet est une évolution du premier solveur Poisson 2D écrit en Fortran.

La version initiale utilisait une matrice globale dense et une résolution directe avec LAPACK. La présente version introduit :

- un stockage creux CSR ;
- un gradient conjugué préconditionné par Jacobi ;
- une quadrature triangulaire à 25 points ;
- un calcul distinguant norme L2, semi-norme H1 et norme H1 complète ;
- deux problèmes à solution exacte connue ;
- des tests unitaires ;
- des tests de convergence directement dans le programme de tests Fortran ;
- une analyse Python des erreurs et des performances ;
- une visualisation des solutions ;
- un pipeline GitHub Actions ;
- un archivage des résultats et une publication par tags.

Le dépôt contient le projet dans le sous-dossier :

```text
Poisson2D_P1_Fortran_v2/
```

Le workflow GitHub Actions est situé à la racine du dépôt :

```text
.github/workflows/ci-v2.yml
```

## Sommaire

1. Objectifs du projet
2. Problèmes mathématiques
3. Méthode numérique
4. Stockage creux CSR
5. Résolution du système
6. Calcul des erreurs
7. Arborescence
8. Description des fichiers
9. Installation
10. Compilation
11. Tests unitaires et convergence
12. Exécution du solveur
13. Visualisation
14. Résultats et interprétation
15. Mesures de performance
16. Pipeline CI/CD


## 1. Objectifs du projet

Les objectifs sont les suivants :

- résoudre l’équation de Poisson 2D avec second membre non nul ;
- construire un maillage triangulaire du carré unité ;
- assembler la matrice et le second membre par éléments finis P1 ;
- imposer des conditions de Dirichlet homogènes ;
- éviter le stockage d’une matrice globale dense ;
- résoudre le système par une méthode itérative ;
- comparer la solution approchée à une solution analytique ;
- vérifier numériquement la convergence ;
- tester séparément les composants du programme ;
- mesurer les temps des différentes étapes ;
- générer des résultats traçables ;
- automatiser la compilation et les validations ;
- conserver les sorties dans les artefacts GitHub Actions.

## 2. Problèmes mathématiques

### 2.1 Domaine et équation

On considère :

\[
\Omega=(0,1)\times(0,1).
\]

On cherche une fonction \(u\) telle que :

\[
-\Delta u=f
\quad\text{dans }\Omega,
\]

avec :

\[
u=0
\quad\text{sur }\partial\Omega.
\]

Deux solutions manufacturées sont proposées.

Le principe d’une solution manufacturée consiste à choisir une fonction exacte, puis à calculer le second membre qui la rend solution de l’équation.

### 2.2 Cas 1 : solution trigonométrique

La solution exacte est :

\[
u(x,y)=\sin(\pi x)\sin(\pi y).
\]

Le second membre est :

\[
f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y).
\]

Le gradient est :

\[
\nabla u=
\begin{pmatrix}
\pi\cos(\pi x)\sin(\pi y)\\
\pi\sin(\pi x)\cos(\pi y)
\end{pmatrix}.
\]

### 2.3 Cas 2 : solution polynomiale

La solution exacte est :

\[
u(x,y)=x(1-x)y(1-y).
\]

Le second membre est :

\[
f(x,y)=2\left[x(1-x)+y(1-y)\right].
\]

Le gradient est :

\[
\nabla u=
\begin{pmatrix}
(1-2x)y(1-y)\\
x(1-x)(1-2y)
\end{pmatrix}.
\]

Les deux solutions satisfont les conditions de Dirichlet homogènes.

## 3. Méthode numérique

### 3.1 Maillage triangulaire

Le carré est subdivisé en \(n\times n\) rectangles.

Chaque rectangle est découpé en deux triangles.

On utilise :

\[
n_x=n_y=n,
\]

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

Le paramètre \(h\) désigne ici le pas directionnel. Le diamètre maximal des triangles vaut \(\sqrt{2}/n\), ce qui ne change pas les ordres de convergence.

### 3.2 Numérotation

Les nœuds sont numérotés ligne par ligne, en faisant varier x le plus rapidement.

Pour des indices commençant à zéro :

\[
\operatorname{index}(i,j)=j(n+1)+i+1.
\]

Les indices globaux Fortran commencent donc à 1.

La visualisation Python utilise la même connectivité, avec des indices commençant à 0.

### 3.3 Formulation variationnelle

On cherche \(u\in H_0^1(\Omega)\) telle que :

\[
\int_\Omega \nabla u\cdot\nabla v\,d\Omega
=
\int_\Omega fv\,d\Omega
\quad
\forall v\in H_0^1(\Omega).
\]

### 3.4 Approximation P1

La solution approchée est :

\[
u_h(x,y)=\sum_j U_j\phi_j(x,y),
\]

où les fonctions \(\phi_j\) sont continues et affines sur chaque triangle.

L’assemblage conduit à :

\[
AU=b.
\]

Les coefficients sont :

\[
A_{ij}
=
\int_\Omega\nabla\phi_j\cdot\nabla\phi_i\,d\Omega,
\]

\[
b_i=
\int_\Omega f\phi_i\,d\Omega.
\]

### 3.5 Matrice élémentaire

Les gradients P1 sont constants sur un triangle \(K\).

On obtient donc :

\[
A_{ij}^{K}
=
|K|\nabla\phi_j\cdot\nabla\phi_i.
\]

La matrice élémentaire est calculée directement à partir de l’aire et des gradients.

### 3.6 Quadrature du second membre

Le projet utilise une règle de Gauss-Legendre 5×5, transformée par Duffy vers le triangle.

Elle contient 25 points.

La convention est :

\[
\sum_q w_q=1.
\]

L’intégrale physique est approchée par :

\[
\int_K g
\approx
|K|\sum_q w_q g(x_q).
\]

Cette règle est exacte, à l’arrondi près, pour les polynômes triangulaires de degré total inférieur ou égal à 8.

Elle n’intègre pas analytiquement les fonctions trigonométriques : une erreur de quadrature subsiste.

### 3.7 Conditions de Dirichlet

Le code assemble uniquement les contributions du bloc intérieur.

Les couplages avec les nœuds de bord sont annulés symétriquement.

Une identité est ensuite placée sur les équations de bord :

\[
U_i=0.
\]

Cette méthode conserve la symétrie pour les conditions homogènes utilisées.

Une extension aux conditions non nulles nécessiterait aussi une correction du second membre.

## 4. Stockage creux CSR

### 4.1 Pourquoi utiliser une matrice creuse ?

Une fonction de base P1 n’interagit qu’avec les fonctions associées aux sommets voisins.

La matrice possède donc beaucoup de coefficients nuls.

Le stockage dense réserverait malgré tout :

\[
N^2
\]

coefficients.

En double précision, le tableau dense seul demanderait environ :

\[
8N^2\text{ octets}.
\]

Le stockage CSR conserve les entrées d’un motif creux, ce qui évite cette allocation globale dense.

### 4.2 Les trois tableaux

Le type `csr_t` contient :

```fortran
integer :: n
integer, allocatable :: row(:)
integer, allocatable :: col(:)
real(dp), allocatable :: val(:)
```

- `row(i)` indique le début de la ligne i ;
- `row(i+1)-1` indique sa fin ;
- `col(k)` donne la colonne de l’entrée k ;
- `val(k)` donne son coefficient.

Les indices sont ici basés sur 1.

### 4.3 Exemple

Considérons :

\[
A=
\begin{pmatrix}
4&-1&0\\
-1&4&-1\\
0&-1&4
\end{pmatrix}.
\]

Une représentation CSR est :

```text
val = [ 4, -1, -1, 4, -1, -1, 4 ]
col = [ 1,  2,  1, 2,  3,  2, 3 ]
row = [ 1,  3,  6, 8 ]
```

La deuxième ligne est stockée aux positions :

```text
row(2) : row(3)-1
```

soit :

```text
3 : 5
```

Les colonnes correspondantes sont `[1,2,3]`.

### 4.4 Application dans le projet

La procédure `pattern` :

1. parcourt les triangles ;
2. recense les sommets voisins de chaque sommet ;
3. élimine les doublons ;
4. trie les colonnes de chaque ligne ;
5. construit les bornes `row` ;
6. alloue `col` et `val`.

La procédure `add` ajoute une contribution élémentaire dans une entrée existante du motif.

La procédure `matvec` calcule :

\[
y=Ax
\]

en parcourant seulement les entrées stockées.

### 4.5 Limitation du constructeur

Le constructeur utilise une borne de sept voisins par sommet, diagonale comprise.

Cette borne est adaptée au maillage structuré de ce projet.

Elle ne doit pas être supposée valable pour un maillage triangulaire non structuré quelconque.

Le programme s’arrête explicitement si cette hypothèse est dépassée.

### 4.6 Signification de nnz

Dans les fichiers de résultats :

```text
nnz
```

désigne le nombre d’entrées stockées.

Certaines entrées structurelles peuvent rester nulles, notamment après traitement du bord.

Ce nombre n’est donc pas nécessairement le nombre exact de coefficients mathématiquement non nuls.

## 5. Résolution du système

Le système est résolu par gradient conjugué préconditionné par Jacobi.

Le préconditionneur est construit à partir de la diagonale :

\[
M=\operatorname{diag}(A).
\]

La méthode nécessite une matrice symétrique définie positive.

La résolution est déclarée réussie après vérification du vrai résidu :

\[
\frac{\|b-AU\|_2}{\|b\|_2}
\leq \mathrm{rtol}.
\]

Le résidu est recalculé explicitement avant le retour du solveur.

Les tolérances sont :

- `1e-11` pour le programme principal ;
- `1e-12` pour les tests de convergence.

Jacobi constitue une première approche simple. Il ne garantit pas un nombre d’itérations indépendant de la finesse du maillage.

## 6. Calcul des erreurs

### 6.1 Comparaison au même point

Aux points de quadrature :

\[
u_h(x_q,y_q)
=
\sum_{i=1}^{3}U_i\lambda_i(x_q,y_q).
\]

Cette valeur est comparée à :

\[
u(x_q,y_q).
\]

Il ne faut pas comparer une valeur numérique prise à un sommet avec une valeur exacte évaluée au barycentre.

### 6.2 Norme L2

\[
E_{L^2}
=
\left(
\int_\Omega (u_h-u)^2\,d\Omega
\right)^{1/2}.
\]

### 6.3 Semi-norme H1

\[
E_{H^1,\mathrm{semi}}
=
\left(
\int_\Omega
|\nabla u_h-\nabla u|^2\,d\Omega
\right)^{1/2}.
\]

### 6.4 Norme H1 complète

\[
E_{H^1}
=
\sqrt{
E_{L^2}^2+
E_{H^1,\mathrm{semi}}^2
}.
\]

### 6.5 Maximum échantillonné

La sortie :

```text
linf_sample
```

est le maximum sur une sous-grille de 15 points barycentriques par triangle.

Cette quantité n’est pas le supremum continu certifié sur tout le domaine.

Les légendes et les conclusions doivent conserver cette distinction.

## 7. Arborescence

```text
Poisson2D
└── Poisson2D_P1_Fortran_v2/
    ├── README.md
    ├── Makefile
    ├── .gitignore
    ├── src/
    │   ├── precision.f90
    │   ├── quadrature.f90
    │   ├── domain.f90
    │   ├── exact.f90
    │   ├── sparse.f90
    │   ├── fem.f90
    │   ├── solver.f90
    │   ├── errors.f90
    │   └── main.f90
    ├── tests/
    │   ├── unit.f90
    │   └── test_analysis.py
    ├── python/
    │   └── plot_results.py
    ├── bin/
    │   ├── poisson2d
    │   └── unit_tests
    ├── build/
    └── results/
        ├── tests/
        │   └── convergence_tests.csv
        ├── manual/
        │   ├── metrics.csv
        │   └── solution.dat
        ├── case_1/
        │   ├── convergence.csv
        │   ├── convergence.png
        │   ├── solution.png
        │   ├── performance.png
        └── case_2/
            └── ...
```

Les dossiers `bin/`, `build/` et `results/` sont générés.

Le dossier `manual/` est créé uniquement lorsqu’un calcul manuel utilise ce chemin.

## 8. Description des fichiers

### 8.1 `precision.f90`

Définit :

- le type réel `dp` ;
- la constante `pi`.

Ce module fournit une convention de précision commune.

### 8.2 `quadrature.f90`

Contient :

- la constante publique `nq=25` ;
- les points de Gauss-Legendre ;
- les poids ;
- la transformation de Duffy ;
- les coordonnées barycentriques.

La procédure principale est :

```fortran
call rule(lam, w)
```

### 8.3 `domain.f90`

Définit le type `mesh_t` :

```text
n          subdivisions
nv         nombre de nœuds
nt         nombre de triangles
xy         coordonnées
tri        connectivité
boundary   indicateur de frontière
```

Procédures :

```fortran
build_mesh(n, m)
geometry(x, area, g)
```

`geometry` utilise le déterminant signé pour calculer les gradients, et sa valeur absolue pour calculer l’aire.

### 8.4 `exact.f90`

Fournit :

```fortran
value(x, y, case_id)
source(x, y, case_id)
gradient(x, y, case_id)
```

Un identifiant de cas invalide provoque une erreur explicite.

### 8.5 `sparse.f90`

Définit `csr_t` et les opérations :

```fortran
pattern(m, a)
add(a, i, j, contribution)
matvec(a, x, y)
```

Il ne construit pas de matrice globale dense.

### 8.6 `fem.f90`

La procédure :

```fortran
call assemble(m, case_id, a, b)
```

réalise :

1. la construction du motif CSR ;
2. l’allocation du second membre ;
3. le calcul des matrices élémentaires ;
4. l’intégration de la charge ;
5. l’assemblage ;
6. le traitement de Dirichlet.

### 8.7 `solver.f90`

Contient :

```fortran
call pcg(a, b, x, rtol, maxit, it, relres)
```

Il contrôle :

- les dimensions ;
- les paramètres ;
- la positivité de la diagonale ;
- les produits nécessaires au PCG ;
- le résidu ;
- le nombre maximal d’itérations.

### 8.8 `errors.f90`

La procédure :

```fortran
call compute_errors(m, u, case_id, l2, h1semi, h1, linf_sample)
```

calcule les quatre mesures d’erreur.

### 8.9 `main.f90`

Le programme principal :

1. lit les arguments ;
2. vérifie le dossier de sortie ;
3. construit le maillage ;
4. assemble le système ;
5. résout ;
6. calcule les erreurs ;
7. mesure les temps ;
8. écrit les résultats.

Arguments :

```text
n  case_id  dossier_sortie
```

### 8.10 `tests/unit.f90`

Regroupe :

- tests de quadrature ;
- tests géométriques ;
- tests du CSR ;
- tests du solveur ;
- tests analytiques des normes ;
- tests de convergence des deux problèmes.

### 8.11 `tests/test_analysis.py`

Vérifie la fonction Python qui calcule les ordres.

Les tests utilisent notamment des pas non dyadiques pour éviter de supposer à tort un rapport de raffinement égal à 2.

### 8.12 `python/study.py`

Automatise les études :

- lancement du solveur ;
- répétitions ;
- lecture des métriques ;
- validation des données ;
- ordres expérimentaux ;
- statistiques de temps ;
- figures ;
- métadonnées ;
- rapport de validation.

### 8.13 `python/plot_results.py`

Relit les résultats sans relancer Fortran.

Il produit :

- figure 1 : solution approchée et solution exacte aux nœuds ;
- figure 2 : logarithmes des erreurs en fonction de log(h).

Par défaut, la solution choisie correspond à la dernière répétition du maillage le plus fin présent dans le tableau.

### 8.14 `Makefile`

Définit :

- l’ordre des modules ;
- la compilation optimisée ;
- la compilation diagnostique ;
- la cible de tests ;
- la cible d’exécution ;
- la cible d’étude ;
- le nettoyage.

## 9. Installation

Les commandes suivantes sont prévues pour Ubuntu ou Debian.

### 9.1 Dépendances système

```bash
sudo apt update

sudo apt install -y \
    git \
    make \
    gfortran \
    python3 \
    python3-pip \
    python3-venv
```

### 9.2 Environnement Python

```bash
python3 -m venv .venv

source .venv/bin/activate

python -m pip install \
    numpy==2.2.6 \
    matplotlib==3.10.3
```

Quitter l’environnement :

```bash
deactivate
```

## 10. Compilation

Depuis le dossier `Poisson2D_P1_Fortran_v2/` :

```bash
make clean
make test

python -m unittest discover -s tests -p 'test_*.py' -v

make
```

### Compilation diagnostique

`make test` active notamment :

```text
-O0
-g
-fcheck=all
-fbacktrace
-ffpe-trap=invalid,zero,overflow
```

### Compilation optimisée

`make` utilise notamment :

```text
-O3
-std=f2008
-Wall
-Wextra
```

Les fichiers `.mod` sont placés dans `build/`.

Le Makefile compile les modules dans l’ordre des dépendances.

La directive `.NOTPARALLEL` évite la compilation simultanée de recettes qui écriraient dans le même dossier de modules.

## 11. Tests unitaires et convergence

### 11.1 Lancer les tests Fortran

```bash
make test
```

Cette commande crée le dossier :

```text
results/tests/
```

puis exécute les tests.

### 11.2 Tests unitaires

| Groupe | Vérifications |
|---|---|
| Quadrature | Somme et positivité des poids, coordonnées barycentriques, monômes jusqu’au degré 8 |
| Maillage | Nombre de nœuds, triangles et nœuds de bord |
| Géométrie | Aire totale, somme des gradients, reproduction affine, orientation inverse |
| CSR | Bornes des lignes et symétrie |
| Solveur | Solution algébrique connue, vrai résidu, second membre nul |
| Normes | Références analytiques et identité H1 |

### 11.3 Références des normes

Pour :

\[
u=x(1-x)y(1-y),
\qquad u_h=0,
\]

on utilise :

\[
\|u\|_{L^2}=\frac{1}{30},
\]

\[
|u|_{H^1}=\sqrt{\frac{1}{45}},
\]

\[
\|u\|_{H^1}=\sqrt{\frac{7}{300}},
\]

\[
\max_\Omega u=\frac{1}{16}.
\]

### 11.4 Tests de convergence Fortran

Les deux problèmes sont résolus sur :

```text
n = 8, 16, 32, 64
```

Pour chaque calcul :

- contrôle du vrai résidu ;
- contrôle du bord ;
- contrôle des valeurs finies et positives ;
- contrôle de l’identité H1.

Pour chaque cas :

- vérification de la décroissance de L2, H1_semi et H1 ;
- calcul des ordres locaux ;
- contrôle d’un ordre global entre n=16 et n=64.

Seuils :

| Mesure | Intervalle |
|---|---|
| L2 | 1.7 < p < 2.3 |
| H1 semi-norme | 0.8 < p < 1.2 |
| H1 complète | 0.8 < p < 1.2 |

Ces seuils concernent les deux problèmes réguliers fournis.

Aucun seuil d’ordre n’est imposé au maximum échantillonné.

### 11.5 Rapport de tests

Le fichier :

```text
results/tests/convergence_tests.csv
```

contient les mesures et les ordres.

Le premier ordre de chaque cas est stocké à zéro avec :

```text
order_defined = 0
```

Il ne faut pas interpréter cette valeur comme un ordre mesuré.

Pour les mesures suivantes :

```text
order_defined = 1
```

### 11.6 Tests Python

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
```

Ils vérifient :

- l’ordre 1 pour une loi E=C*h ;
- l’ordre 2 pour E=C*h² ;
- les pas non dyadiques ;
- le rejet des entrées invalides ;
- le premier ordre indéfini.

## 12. Exécution du solveur

### 12.1 Cible pratique

```bash
make run
```

Valeurs par défaut :

```text
n = 32
case_id = 1
sortie = results/manual
```

Personnalisation :

```bash
make run N=64 CASE=2 OUT=results/polynomial
```

### 12.2 Exécution directe

```bash
mkdir -p results/manual

./bin/poisson2d 32 1 results/manual
```

Le dossier doit exister avant le lancement direct.

Le programme vérifie son accessibilité avant le calcul.

### 12.3 Sorties

```text
results/manual/
├── metrics.csv
└── solution.dat
```

Chaque lancement remplace les fichiers de son dossier.

Utiliser des dossiers différents pour conserver plusieurs calculs.

## 14. Visualisation

### 14.1 Relire les résultats

```bash
python python/plot_results.py
```

Sans fenêtre :

```bash
python python/plot_results.py --no-show
```

### 14.2 Figure 1

Fichier :

```text
results/case_1/solution.png
```

La figure représente :

- à gauche, la solution approchée ;
- à droite, les valeurs exactes aux nœuds.

Les deux panneaux utilisent :

- la même triangulation ;
- la colormap `viridis` ;
- les mêmes niveaux ;
- une barre de couleurs commune.

Le panneau exact visualise les valeurs analytiques aux nœuds, pas une évaluation indépendante sur une grille plus fine.

### 14.3 Figure 2

Fichier :

```text
results/case_1/convergence.png
```

Le graphique représente :

\[
X=\log(h),
\qquad
Y=\log(E).
\]

Il contient :

- L2 ;
- H1 semi-norme ;
- H1 complète ;
- maximum échantillonné ;
- repères de pentes 1 et 2.

Les logarithmes sont naturels.

L’axe horizontal est inversé pour lire le raffinement de gauche à droite.

Les droites de référence sont des repères, pas une validation des ordres.

## 15. Résultats et interprétation

### 15.1 Résultats réellement disponibles

Les résultats chiffrés doivent provenir des fichiers générés :

```text
results/tests/convergence_tests.csv
results/case_1/convergence.csv
results/case_2/convergence.csv
```

Cette documentation ne présente pas de valeurs Fortran inventées.

La compilation et la convergence de cette version doivent être confirmées par les exécutions locales et la CI.

### 15.2 Comportements de référence

Pour les solutions régulières considérées, on recherche :

\[
E_{L^2}=O(h^2),
\]

\[
E_{H^1,\mathrm{semi}}=O(h).
\]

Une décroissance observée indique un rapprochement selon la mesure considérée.

Cependant :

- une ressemblance visuelle ne prouve pas un ordre ;
- une décroissance stricte à chaque raffinement n’est pas universelle ;
- les maillages grossiers peuvent être hors du régime asymptotique ;
- les erreurs de quadrature et de résolution peuvent influencer les mesures.


## 16. Mesures de performance

### 16.1 Temps enregistrés

| Colonne | Étape |
|---|---|
| `mesh_s` | Construction du maillage |
| `assembly_s` | Motif CSR et assemblage |
| `solve_s` | Résolution PCG et allocations internes |
| `error_s` | Calcul des erreurs |

Ces temps muraux sont obtenus avec `system_clock`.

Ils excluent les écritures des fichiers et les tracés Python.

### 16.2 Répétitions

Le script conserve pour chaque étape :

- médiane ;
- minimum ;
- maximum.

Il n’y a pas de phase d’échauffement dédiée.

### 16.3 Stockage

Les colonnes :

```text
csr_bytes
dense_matrix_bytes
```

comparent :

- la taille des tableaux CSR ;
- une estimation théorique dense égale à 8*N².

`csr_bytes` ne mesure pas le pic mémoire du processus.

### 16.4 Figure de performance

```text
results/case_1/performance.png
```

Elle représente les temps et les tailles matricielles en fonction du nombre de degrés de liberté.

Aucun seuil de temps ne bloque la CI, car les runners partagés peuvent varier.

## 17. Pipeline CI/CD

### 17.1 Objectif

Le pipeline automatise :

```text
Récupération du dépôt
        ↓
Installation des dépendances
        ↓
Tests unitaires et convergence Fortran
        ↓
Tests Python
        ↓
Compilation optimisée
        ↓
Campagnes des deux problèmes
        ↓
Génération des figures
        ↓
Vérification des sorties
        ↓
Archivage
        ↓
Release si push d’un tag v* et réussite
```

### 17.2 Déclencheurs

| Événement | Validation | Release |
|---|---|---|
| Push main | Oui | Non |
| Pull request vers main | Oui | Non |
| Lancement manuel | Oui | Non |
| Push tag v* | Oui | Après réussite |

### 17.3 Deux jobs

- `scientific-tests` : tests, études, figures et artefact ;
- `release` : publication des résultats validés.


