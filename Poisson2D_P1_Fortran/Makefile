# =============================================================================
# Makefile du projet Poisson2D_Fortran
#
# Commandes :
#   make             : compile le programme
#   make run         : compile puis exécute avec n = 20
#   make convergence : étudie plusieurs maillages
#   make clean       : supprime les fichiers générés
# =============================================================================

FC = gfortran
FFLAGS = -O2 -Wall -Wextra
LIBS = -llapack -lblas

MODULES = params_module.f90 \
          mesh_module.f90 \
          exact_solution_module.f90 \
          fem_module.f90 \
          solver_module.f90 \
          error_module.f90

MAIN = main.f90

EXEC = bin/poisson2d

all: $(EXEC)

$(EXEC): $(MODULES) $(MAIN)
	mkdir -p bin
	$(FC) $(FFLAGS) $(MODULES) $(MAIN) -o $(EXEC) $(LIBS)

run: $(EXEC)
	mkdir -p results
	./$(EXEC) 2
	./$(EXEC) 5
	./$(EXEC) 10
	./$(EXEC) 20
	./$(EXEC) 40
	./$(EXEC) 80
	./$(EXEC) 160
	./$(EXEC) 320

# Étude de convergence sur plusieurs maillages.
# Chaque exécution ajoute une ligne dans results/error_convergence.dat.
convergence: $(EXEC)
	mkdir -p results
	rm -f results/error_convergence.dat
	./$(EXEC) 2
	./$(EXEC) 5
	./$(EXEC) 10
	./$(EXEC) 20
	./$(EXEC) 40
	./$(EXEC) 80
	./$(EXEC) 160
	./$(EXEC) 320

clean:
	rm -rf bin *.mod *.o results/*.dat results/*.png

.PHONY: all run convergence clean