!==============================================================================
! Fichier : src/errors.f90
!
! Rôle :
!   Calculer les mesures d'erreur entre la solution P1 et la solution exacte.
!
! Sorties :
!   l2          : estimation de ||u_h-u||_L2 ;
!   h1semi      : estimation de |u_h-u|_H1 ;
!   h1          : sqrt(l2**2 + h1semi**2) ;
!   linf_sample : maximum sur une sous-grille de points de chaque triangle.
!
! Important :
!   linf_sample n'est pas un supremum continu certifié.
!   Les intégrales sont calculées par quadrature.
!==============================================================================
module errors

  use precision, only: dp
  use domain, only: mesh_t, geometry
  use quadrature, only: nq, rule
  use exact, only: value, gradient

  implicit none

  private

  public :: compute_errors

contains

  subroutine compute_errors(m, u, case_id, l2, h1semi, h1, linf_sample)

    type(mesh_t), intent(in) :: m

    real(dp), intent(in) :: u(:)

    integer, intent(in) :: case_id

    real(dp), intent(out) :: l2
    real(dp), intent(out) :: h1semi
    real(dp), intent(out) :: h1
    real(dp), intent(out) :: linf_sample

    real(dp) :: lam(3, nq)
    real(dp) :: w(nq)

    real(dp) :: x(2, 3)
    real(dp) :: g(2, 3)
    real(dp) :: grad_num(2)
    real(dp) :: grad_diff(2)

    real(dp) :: area
    real(dp) :: pt(2)
    real(dp) :: bary(3)
    real(dp) :: diff

    integer :: ids(3)
    integer :: e, q, i, j

    if (size(u) /= m%nv) then
      error stop 'compute_errors : dimension de u incorrecte'
    end if

    call rule(lam, w)

    l2 = 0.0_dp
    h1semi = 0.0_dp
    linf_sample = 0.0_dp

    do e = 1, m%nt

      ids = m%tri(:, e)
      x = m%xy(:, ids)

      call geometry(x, area, g)

      ! Le gradient d'une fonction P1 est constant sur le triangle.
      grad_num = matmul(g, u(ids))

      do q = 1, nq

        pt = matmul(x, lam(:, q))

        ! Évaluer les deux solutions AU MÊME POINT.
        diff = dot_product(lam(:, q), u(ids)) &
             - value(pt(1), pt(2), case_id)

        grad_diff = grad_num &
                  - gradient(pt(1), pt(2), case_id)

        l2 = l2 + area * w(q) * diff**2

        h1semi = h1semi + area * w(q) &
               * dot_product(grad_diff, grad_diff)

      end do

      ! Échantillonnage sur une sous-grille barycentrique.
      !
      ! i+j+k=4 donne 15 points par triangle,
      ! incluant les sommets et des points sur les arêtes.
      do i = 0, 4

        do j = 0, 4-i

          bary = [ &
            real(i, dp), &
            real(j, dp), &
            real(4-i-j, dp) ] / 4.0_dp

          pt = matmul(x, bary)

          diff = dot_product(bary, u(ids)) &
               - value(pt(1), pt(2), case_id)

          linf_sample = max(linf_sample, abs(diff))

        end do

      end do

    end do

    ! Protection contre un éventuel arrondi négatif infinitésimal.
    l2 = sqrt(max(l2, 0.0_dp))
    h1semi = sqrt(max(h1semi, 0.0_dp))

    ! Norme H1 complète.
    h1 = sqrt(l2**2 + h1semi**2)

  end subroutine compute_errors

end module errors