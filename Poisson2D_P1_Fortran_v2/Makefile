# =============================================================================
# Makefile : Poisson2D_P1_Fortran_v2
#
# Commandes :
#   make          : compiler le solveur optimisé
#   make test     : compiler et exécuter les tests en mode diagnostic
#   make run      : exécuter le cas 1 avec n=32 dans results/manual
#   make study    : étude de convergence du cas 1
#   make clean    : supprimer bin/ et build/, sans effacer les résultats
#
# Paramètres personnalisables :
#   make run N=64 CASE=2 OUT=results/polynomial
# =============================================================================

FC = gfortran
PYTHON ?= python3

FFLAGS = -O3 -std=f2008 -Wall -Wextra -ffree-line-length-none

TESTFLAGS = -O0 -g -std=f2008 -Wall -Wextra \
            -ffree-line-length-none \
            -fcheck=all \
            -fbacktrace \
            -ffpe-trap=invalid,zero,overflow

# Valeurs par défaut pour la cible run.
N ?= 32
CASE ?= 1
OUT ?= results/manual

# Sources communes, dans l'ordre des dépendances.
CORE = src/precision.f90 \
       src/quadrature.f90 \
       src/domain.f90 \
       src/exact.f90 \
       src/sparse.f90 \
       src/fem.f90 \
       src/solver.f90 \
       src/errors.f90

# Les recettes compilent plusieurs modules dans le même dossier build.
# Empêcher leur exécution simultanée par make -j.
.NOTPARALLEL:

#.PHONY: all test run study clean

all: bin/poisson2d

# Compilation du solveur.
bin/poisson2d: $(CORE) src/main.f90 Makefile
	mkdir -p bin build
	$(FC) $(FFLAGS) -Jbuild -Ibuild $(CORE) src/main.f90 -o $@

# Compilation des tests avec contrôles d'exécution.
bin/unit_tests: $(CORE) tests/unit.f90 Makefile
	mkdir -p bin build
	$(FC) $(TESTFLAGS) -Jbuild -Ibuild $(CORE) tests/unit.f90 -o $@

# Exécution des tests.
test: bin/unit_tests
	./bin/unit_tests

# Calcul individuel avec création préalable du dossier.
run: bin/poisson2d
	mkdir -p "$(OUT)"
	./bin/poisson2d "$(N)" "$(CASE)" "$(OUT)"

# Nettoyage des fichiers de compilation uniquement.
clean:
	rm -rf bin build