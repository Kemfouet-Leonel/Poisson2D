# Poisson2D P1 — version creuse et validation

Version autonome à déposer dans `Poisson2D_P1_Fortran_v2/` à côté de l'ancienne version.
Le workflow fourni doit être copié à la RACINE du dépôt, dans `.github/workflows/ci-v2.yml`.
Ne pas laisser seulement une copie dans le sous-dossier : GitHub ne la détecterait pas.

## Installation et exécution

```bash
sudo apt install gfortran make python3-venv
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy==2.2.6 matplotlib==3.10.3
make test
make
mkdir -p results/manual
./bin/poisson2d 32 1 results/manual
python python/study.py --case 1 --check
python python/study.py --case 2 --check
```

Arguments : `n case_id dossier_sortie`. Le dossier doit déjà exister pour le lancement direct.
Cas 1 : u=sin(pi*x)sin(pi*y). Cas 2 : u=x(1-x)y(1-y).
Les deux ont Dirichlet homogène sur le carré unité et f=-Delta u non nul.

## Modules

- precision : double précision et pi.
- quadrature : 7 points Dunavant de degré 5, poids de somme 1 ; multiplier par aire physique.
- domain : maillage structuré triangulaire, géométrie et gradients signés.
- exact : deux solutions manufacturées et leurs gradients/sources.
- sparse : stockage CSR, motif et produit matrice-vecteur.
- fem : rigidité exacte pour P1, charge par quadrature ; suppression symétrique des couplages au bord et identité sur les nœuds de bord.
- solver : PCG-Jacobi, tolérance relative 1e-11 et contrôle du vrai résidu avant succès.
- errors : L2, semi-norme H1, norme H1 complète, maximum échantillonné explicitement nommé.
- main : arguments, chronométrages muraux, métriques et fichiers de solution.
- tests/unit.f90 : monômes jusqu'au degré 5, géométrie, symétrie CSR, solution algébrique connue, références analytiques des normes.
- python/study.py : répétitions, ordres locaux et régression finale, CSV, figures, métadonnées.

## Validation et limites

L2 et H1 sont des estimations par quadrature, pas des intégrales symboliques ni des bornes certifiées.
Dunavant degré 5 intègre exactement les polynômes de degré total <=5, PAS les fonctions trigonométriques.
Pour le test uh=0 et u=x(1-x)y(1-y), la semi-norme H1 au carré vaut 1/45 ; la norme L2 vaut 1/30. Son intégrande au carré est de degré 8 : le test L2 utilise une tolérance adaptée, sans prétendre à l'exactitude de la règle.
Le maximum `linf_sample` est mesuré sur une sous-grille de 15 points par triangle, PAS le supremum continu certifié.
Le constructeur CSR exploite au plus 7 voisins par sommet de CE maillage structuré. Il s'arrête explicitement si cette hypothèse est violée.
Les entrées structurelles nulles sont conservées ; nnz désigne donc le nombre d'entrées stockées.
Le traitement du bord ne prend en charge que Dirichlet homogène.
Jacobi n'est pas un préconditionneur optimal indépendant du maillage. Cette version n'implémente ni multigrille ni solveur creux direct.
Les intervalles d'ordres du script sont des tests de non-régression pour les deux problèmes fournis, pas un critère universel.
La norme H1 vérifie H1²=L2²+H1_semi², ce que les tests contrôlent.

## Sorties

Chaque lancement direct écrit metrics.csv et solution.dat dans son dossier, sans ajout ambigu.
La campagne crée results/case_1 ou results/case_2 avec :
- convergence.csv : mesures et ordres locaux p=log(Eancien/Enouveau)/log(hancien/hnouveau) ; premier ordre non défini (NaN).
- convergence.png : logarithmes naturels des erreurs.
- solution.png : cartes nodales communes aux deux solutions (ce n'est pas un rendu exact triangle par triangle).
- performance.png : temps médians et stockage.
- environment.json : environnement, compilateur, commande, commit.
- n*/repeat*/metrics.csv et solution.dat : sorties brutes de chaque répétition.

## Performances

mesh_s, assembly_s (motif inclus), solve_s et error_s sont des temps muraux via system_clock, hors écriture des fichiers.
Le script conserve médiane, minimum et maximum des répétitions. Les temps ne sont pas des tests CI bloquants : les runners fluctuent.
csr_bytes compte uniquement val/col/row, PAS le pic mémoire RSS, le maillage ou les vecteurs temporaires.
dense_matrix_bytes=8*N² est un comparatif de stockage théorique ; aucune matrice dense n'est allouée.
Le maillage peut être augmenté avec --meshes 16 24 32 48 64 96 128 ; mesurer avant de choisir les grands cas.
La limite d'argument n<=4096 est une garde technique, pas une garantie de faisabilité.

## CI

Copier `.github/workflows/ci-v2.yml` à la racine Git du dépôt ; son working-directory suppose le nom `Poisson2D_P1_Fortran_v2`.
La CI exécute les tests diagnostiques puis les deux études. Actions > exécution > Artifacts > poisson2d-v2-validation permet de récupérer les sorties disponibles, conservées 14 jours.
Ce workflow complète la CI scientifique, sans remplacer automatiquement l'ancien workflow ni sa release. Désactiver ou adapter l'ancien workflow si nécessaire pour ne pas confondre les validations.

## Statut de vérification de cette livraison

Les fichiers Python ont été vérifiés syntaxiquement. Les constantes de quadrature et références analytiques ont été contrôlées indépendamment en Python lors de la préparation.
Aucun compilateur Fortran n'était disponible dans l'environnement de préparation : la compilation et l'exécution Fortran restent à confirmer avec make test et la CI. Ne pas annoncer une convergence mesurée avant cette exécution.
