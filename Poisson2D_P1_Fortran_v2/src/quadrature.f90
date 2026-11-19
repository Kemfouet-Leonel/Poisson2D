!==============================================================================
! Fichier : src/quadrature.f90
!
! Rôle :
!   Fournir une quadrature sur un triangle à partir :
!     - d'une règle de Gauss-Legendre à 5 points sur [0,1] ;
!     - d'un produit tensoriel de cette règle ;
!     - d'une transformation de Duffy du carré vers le triangle.
!
! Nombre de points :
!   nq = 5 * 5 = 25.
!
! Convention des poids :
!   sum(w) = 1.
!
! Intégration sur un triangle physique K :
!
!   intégrale_K f ≈ aire(K) * sum_q w(q) * f(x_q,y_q).
!
! Les coordonnées lam(:,q) sont barycentriques :
!
!   x_q = sum_i lam(i,q) * x_i.
!
! La règle est exacte, à l'arrondi près, pour les polynômes
! triangulaires de degré total inférieur ou égal à 8.
!==============================================================================
module quadrature

  use precision, only: dp

  implicit none

  private

  integer, parameter, public :: nq = 25

  public :: rule

contains

  subroutine rule(lam, w)

    real(dp), intent(out) :: lam(3, nq)
    real(dp), intent(out) :: w(nq)

    ! Abscisses de Gauss-Legendre sur [-1,1].
    real(dp), parameter :: xi(5) = [ &
      -0.9061798459386639928_dp, &
      -0.5384693101056830910_dp, &
       0.0000000000000000000_dp, &
       0.5384693101056830910_dp, &
       0.9061798459386639928_dp ]

    ! Poids correspondants sur [-1,1].
    real(dp), parameter :: wi(5) = [ &
      0.2369268850561890875_dp, &
      0.4786286704993664680_dp, &
      0.5688888888888888889_dp, &
      0.4786286704993664680_dp, &
      0.2369268850561890875_dp ]

    real(dp) :: r, s
    real(dp) :: wr, ws

    integer :: i, j, q

    q = 0

    do i = 1, 5

      ! Passage de [-1,1] à [0,1].
      r = 0.5_dp * (xi(i) + 1.0_dp)
      wr = 0.5_dp * wi(i)

      do j = 1, 5

        s = 0.5_dp * (xi(j) + 1.0_dp)
        ws = 0.5_dp * wi(j)

        q = q + 1

        ! Transformation de Duffy :
        !
        !   x_ref = r
        !   y_ref = (1-r)*s.
        !
        ! Sur le triangle de référence :
        !   lambda_1 = 1 - x_ref - y_ref,
        !   lambda_2 = x_ref,
        !   lambda_3 = y_ref.
        lam(1, q) = (1.0_dp - r) * (1.0_dp - s)
        lam(2, q) = r
        lam(3, q) = (1.0_dp - r) * s

        ! Le jacobien de Duffy vaut (1-r).
        !
        ! Le facteur 2 permet de normaliser les poids à 1,
        ! puisque le triangle de référence a une aire de 1/2.
        w(q) = 2.0_dp * wr * ws * (1.0_dp - r)

      end do

    end do

  end subroutine rule

end module quadrature