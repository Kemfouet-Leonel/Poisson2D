!==============================================================================
! Module : error_module
! Rôle :
!   Calculer l’erreur entre la solution numérique uh et la solution exacte u.
!
!   Trois normes sont calculées :
!
!     - erreur L2 :
!       || uh - u ||_L2
!
!     - erreur H1 :
!       || uh - u ||_H1
!
!     - erreur L_infini :
!       max | uh - u |
!==============================================================================
module error_module

  use params_module
  use mesh_module
  use exact_solution_module

  implicit none

contains

  !--------------------------------------------------------------------------
  ! Subroutine : compute_errors
  ! Rôle :
  !   Calcule les erreurs L2, H1 et L_infini.
  !
  ! Entrée :
  !   uh : solution numérique obtenue par éléments finis.
  !
  ! Sorties :
  !   error_l2   : erreur en norme L2.
  !   error_h1   : erreur en norme H1.
  !   error_linf : erreur en norme L_infini.
  !--------------------------------------------------------------------------
  subroutine compute_errors(uh, error_l2, error_h1, error_linf)
    real(dp), intent(in) :: uh(:)

    real(dp), intent(out) :: error_l2, error_h1, error_linf

    integer :: e, i

    ! Coordonnées des sommets du triangle courant.
    real(dp) :: x1, y1, x2, y2, x3, y3

    ! Aire du triangle courant.
    real(dp) :: area

    ! Gradient numérique et gradient exact.
    real(dp) :: grad_num(2), grad_exact(2)

    ! Centre de gravité du triangle.
    real(dp) :: xc, yc

    ! Solution exacte au centre de gravité.
    real(dp) :: ue

    ! Différence entre solution numérique et solution exacte.
    real(dp) :: diff

    ! Matrice des dérivées des fonctions de base.
    real(dp) :: b_matrix(2, 3)

    ! Initialisation des erreurs.
    error_l2 = 0.0_dp
    error_h1 = 0.0_dp
    error_linf = 0.0_dp

    ! Boucle sur tous les triangles du maillage.
    do e = 1, n_elements

       ! Récupération des sommets du triangle.
       x1 = nodes(1, elements(1, e))
       y1 = nodes(2, elements(1, e))

       x2 = nodes(1, elements(2, e))
       y2 = nodes(2, elements(2, e))

       x3 = nodes(1, elements(3, e))
       y3 = nodes(2, elements(3, e))

       ! Aire du triangle.
       area = 0.5_dp * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

       ! Centre de gravité du triangle.
       xc = (x1 + x2 + x3) / 3.0_dp
       yc = (y1 + y2 + y3) / 3.0_dp

       ! Solution exacte au centre de gravité.
       ue = exact_solution(xc, yc)

       ! Différence entre la valeur numérique au premier sommet
       ! et la solution exacte au centre de gravité.
       diff = uh(elements(1, e)) - ue

       ! Contribution à l’erreur L2.
       error_l2 = error_l2 + area * diff**2

       ! Mise à jour de l’erreur L_infini.
       error_linf = max(error_linf, abs(diff))

       !---------------------------------------------------------------------
       ! Calcul du gradient numérique.
       !---------------------------------------------------------------------
       b_matrix(1, 1) = y2 - y3
       b_matrix(1, 2) = y3 - y1
       b_matrix(1, 3) = y1 - y2

       b_matrix(2, 1) = x3 - x2
       b_matrix(2, 2) = x1 - x3
       b_matrix(2, 3) = x2 - x1

       b_matrix = b_matrix / (2.0_dp * area)

       grad_num = 0.0_dp

       do i = 1, 3
          grad_num(1) = grad_num(1) + uh(elements(i, e)) * b_matrix(1, i)
          grad_num(2) = grad_num(2) + uh(elements(i, e)) * b_matrix(2, i)
       end do

       ! Gradient exact au centre de gravité.
       call exact_gradient(xc, yc, grad_exact)

       ! Contribution à l’erreur H1.
       error_h1 = error_h1 + area * ((grad_num(1) - grad_exact(1))**2 &
            + (grad_num(2) - grad_exact(2))**2)

    end do

    ! Les normes L2 et H1 sont définies comme des racines carrées.
    error_l2 = sqrt(error_l2)
    error_h1 = sqrt(error_h1)

  end subroutine compute_errors

end module error_module