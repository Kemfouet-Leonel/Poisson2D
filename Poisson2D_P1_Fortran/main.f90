!==============================================================================
! Programme principal : poisson2d
!
! Rôle :
!   Résoudre l’équation de Poisson 2D par éléments finis P1,
!   puis étudier l’influence de la finesse du maillage sur l’erreur.
!
!   Le maillage est uniforme : nx = ny = n.
!
!   Le pas de maillage est :
!
!       h = 1 / n.
!
!   Pour chaque valeur de n, le programme :
!     - construit le maillage,
!     - assemble le système,
!     - résout le système,
!     - calcule les erreurs,
!     - ajoute les résultats dans results/error_convergence.dat.

!==============================================================================
program poisson2d

  use params_module
  use mesh_module
  use fem_module
  use solver_module
  use error_module

  implicit none

  ! Nombre de subdivisions du maillage selon x et y.
  ! Il est importé depuis params_module.
  integer :: n

  ! Indices de boucle et identifiant de fichier.
  integer :: i, unit_out

  ! Solution numérique.
  real(dp), allocatable :: solution(:)

  ! Erreurs calculées.
  real(dp) :: error_l2, error_h1, error_linf

  ! Variable pour lire l’argument de la ligne de commande.
  character(len=32) :: arg

  ! Valeur par défaut du maillage.
  n = 20

  ! Lecture du nombre de subdivisions depuis la ligne de commande.
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg)
     read(arg, *) n
  end if

  !---------------------------------------------------------------------------
  ! Étape 1 : initialisation des paramètres.
  !
  ! Le maillage est carré : nx = ny = n.
  !---------------------------------------------------------------------------
  call set_parameters(n, n)

  !---------------------------------------------------------------------------
  ! Étape 2 : construction du maillage.
  !---------------------------------------------------------------------------
  call build_mesh()

  !---------------------------------------------------------------------------
  ! Étape 3 : assemblage de la matrice et du second membre.
  !---------------------------------------------------------------------------
  call assemble_system()

  !---------------------------------------------------------------------------
  ! Étape 4 : résolution du système linéaire A U = b.
  !---------------------------------------------------------------------------
  allocate(solution(n_dof))

  call solve_linear_system(global_matrix, rhs, solution)

  !---------------------------------------------------------------------------
  ! Étape 5 : calcul des erreurs.
  !---------------------------------------------------------------------------
  call compute_errors(solution, error_l2, error_h1, error_linf)

  !---------------------------------------------------------------------------
  ! Étape 6 : écriture de la solution numérique et exacte.
  !
  ! Format du fichier solution.dat :
  !
  !   x   y   u_numerique   u_exacte
  !---------------------------------------------------------------------------
  open(newunit=unit_out, file="results/solution.dat", status="replace")

  do i = 1, n_nodes
     write(unit_out, *) nodes(1, i), nodes(2, i), solution(i), &
          exact_solution(nodes(1, i), nodes(2, i))
  end do

  close(unit_out)

  !---------------------------------------------------------------------------
  ! Étape 7 : ajout des erreurs dans le fichier de convergence.
  !
  ! Format du fichier error_convergence.dat :
  !
  !   n   h   erreur_L2   erreur_H1   erreur_Linf
  !---------------------------------------------------------------------------
  open(newunit=unit_out, file="results/error_convergence.dat", &
       status="unknown", position="append")

  write(unit_out, *) n, 1.0_dp / real(n, dp), &
       error_l2, error_h1, error_linf

  close(unit_out)

  !---------------------------------------------------------------------------
  ! Étape 8 : affichage des résultats.
  !---------------------------------------------------------------------------
  write(*, *) "----------------------------------------"
  write(*, *) "Résolution de -Delta u = f en 2D"
  write(*, *) "Maillage :", n, "x", n
  write(*, *) "Pas h =", 1.0_dp / real(n, dp)
  write(*, *) "Nombre de noeuds :", n_nodes
  write(*, *) "Nombre de triangles :", n_elements
  write(*, *) "Erreur L2 :", error_l2
  write(*, *) "Erreur H1 :", error_h1
  write(*, *) "Erreur Linf :", error_linf
  write(*, *) "----------------------------------------"

end program poisson2d