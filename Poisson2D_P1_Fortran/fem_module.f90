!==============================================================================
! Module : fem_module
! Rôle :
!   Construire le système linéaire issu de la méthode des éléments finis P1 :
!
!       A U = b.
!
!   La matrice A est la matrice de rigidité globale.
!   Le vecteur b est le second membre global.
!
!   Pour chaque triangle K, on calcule :
!
!       A^K(i,j) = intégrale de grad(phi_i) . grad(phi_j) sur K
!       b^K(i)   = intégrale de f * phi_i sur K
!
!   Puis on assemble ces contributions dans la matrice et le vecteur globaux.
!==============================================================================
module fem_module

  use params_module
  use mesh_module
  use exact_solution_module

  implicit none

  ! Matrice globale du système linéaire.
  real(dp), allocatable :: global_matrix(:, :)

  ! Second membre global.
  real(dp), allocatable :: rhs(:)

contains

  !--------------------------------------------------------------------------
  ! Subroutine : assemble_system
  ! Rôle :
  !   Assemble la matrice globale A et le second membre global b.
  !
  !   Étapes :
  !     1. Initialiser A et b à zéro.
  !     2. Parcourir tous les triangles.
  !     3. Calculer la matrice élémentaire et le second membre élémentaire.
  !     4. Ajouter les contributions aux positions globales.
  !     5. Appliquer les conditions de Dirichlet.
  !--------------------------------------------------------------------------
  subroutine assemble_system()
    integer :: e, i, j, a, b

    ! Coordonnées des trois sommets du triangle courant.
    real(dp) :: x1, y1, x2, y2, x3, y3

    ! Aire du triangle courant.
    real(dp) :: area

    ! Matrice contenant les dérivées des fonctions de base P1.
    real(dp) :: b_matrix(2, 3)

    ! Matrice élémentaire de rigidité.
    real(dp) :: local_matrix(3, 3)

    ! Second membre élémentaire.
    real(dp) :: local_rhs(3)

    ! Coordonnées du centre de gravité du triangle.
    real(dp) :: xc, yc

    ! Allocation de la matrice et du second membre globaux.
    allocate(global_matrix(n_dof, n_dof))
    allocate(rhs(n_dof))

    ! Initialisation à zéro.
    global_matrix = 0.0_dp
    rhs = 0.0_dp

    ! Boucle sur tous les triangles du maillage.
    do e = 1, n_elements

       ! Récupération des coordonnées des trois sommets.
       x1 = nodes(1, elements(1, e))
       y1 = nodes(2, elements(1, e))

       x2 = nodes(1, elements(2, e))
       y2 = nodes(2, elements(2, e))

       x3 = nodes(1, elements(3, e))
       y3 = nodes(2, elements(3, e))

       ! Aire du triangle par la formule du déterminant.
       area = 0.5_dp * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

       !---------------------------------------------------------------------
       ! Matrice B :
       !
       ! B contient les dérivées des trois fonctions de base P1.
       !
       ! B(1,i) = dérivée de phi_i selon x.
       ! B(2,i) = dérivée de phi_i selon y.
       !---------------------------------------------------------------------
       b_matrix(1, 1) = y2 - y3
       b_matrix(1, 2) = y3 - y1
       b_matrix(1, 3) = y1 - y2

       b_matrix(2, 1) = x3 - x2
       b_matrix(2, 2) = x1 - x3
       b_matrix(2, 3) = x2 - x1

       ! Normalisation par 2 * aire du triangle.
       b_matrix = b_matrix / (2.0_dp * area)

       !---------------------------------------------------------------------
       ! Matrice élémentaire de rigidité :
       !
       ! A^K(i,j) = aire * (grad phi_i . grad phi_j)
       !---------------------------------------------------------------------
       local_matrix = 0.0_dp

       do i = 1, 3
          do j = 1, 3
             local_matrix(i, j) = area * (b_matrix(1, i) * b_matrix(1, j) &
                  + b_matrix(2, i) * b_matrix(2, j))
          end do
       end do

       !---------------------------------------------------------------------
       ! Second membre élémentaire.
       !
       ! Avec une quadrature au centre de gravité :
       !
       ! b^K(i) = (aire / 3) * f(centre de gravité)
       !---------------------------------------------------------------------
       xc = (x1 + x2 + x3) / 3.0_dp
       yc = (y1 + y2 + y3) / 3.0_dp

       local_rhs = area / 3.0_dp * source_term(xc, yc)

       !---------------------------------------------------------------------
       ! Assemblage global.
       !
       ! Pour chaque contribution locale, on ajoute sa valeur
       ! dans la matrice et le vecteur globaux.
       !---------------------------------------------------------------------
       do i = 1, 3
          a = elements(i, e)

          rhs(a) = rhs(a) + local_rhs(i)

          do j = 1, 3
             b = elements(j, e)
             global_matrix(a, b) = global_matrix(a, b) + local_matrix(i, j)
          end do
       end do

    end do

    ! Application de la condition u = 0 sur la frontière.
    call apply_dirichlet()

  end subroutine assemble_system

  !--------------------------------------------------------------------------
  ! Subroutine : apply_dirichlet
  ! Rôle :
  !   Imposer la condition de Dirichlet homogène :
  !
  !       u = 0 sur la frontière.
  !
  !   Méthode utilisée :
  !     - la ligne correspondant à un nœud de bord est annulée,
  !     - un 1 est placé sur la diagonale,
  !     - le second membre associé est mis à zéro.
  !--------------------------------------------------------------------------
  subroutine apply_dirichlet()
    integer :: i, k

    do k = 1, n_boundary_nodes
       i = boundary_nodes(k)

       ! Annulation de la ligne i de la matrice.
       global_matrix(i, :) = 0.0_dp

       ! Mise à 1 sur la diagonale.
       global_matrix(i, i) = 1.0_dp

       ! Valeur imposée : u = 0 sur la frontière.
       rhs(i) = 0.0_dp
    end do

  end subroutine apply_dirichlet

end module fem_module