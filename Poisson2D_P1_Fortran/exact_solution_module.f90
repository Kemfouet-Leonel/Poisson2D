!==============================================================================
! Module : exact_solution_module
! Rôle :
!   Fournir les données analytiques du problème :
!     - la solution exacte,
!     - le second membre f,
!     - le gradient exact de la solution.
!
!   Ces fonctions permettent de valider la solution numérique.
!==============================================================================
module exact_solution_module

  use params_module

  implicit none

contains

  !--------------------------------------------------------------------------
  ! Fonction : exact_solution
  ! Rôle :
  !   Retourne la solution exacte du problème :
  !
  !       u(x,y) = sin(pi x) * sin(pi y)
  !--------------------------------------------------------------------------
  function exact_solution(x, y) result(u)
    real(dp), intent(in) :: x, y
    real(dp) :: u

    u = sin(pi * x) * sin(pi * y)

  end function exact_solution

  !--------------------------------------------------------------------------
  ! Fonction : source_term
  ! Rôle :
  !   Retourne le second membre f(x,y).
  !
  !   Comme :
  !       u(x,y) = sin(pi x) sin(pi y),
  !
  !   on a :
  !       -Delta(u) = 2 pi^2 sin(pi x) sin(pi y).
  !
  !   Donc :
  !       f(x,y) = 2 pi^2 sin(pi x) sin(pi y).
  !--------------------------------------------------------------------------
  function source_term(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f

    f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)

  end function source_term

  !--------------------------------------------------------------------------
  ! Subroutine : exact_gradient
  ! Rôle :
  !   Calcule le gradient exact de la solution :
  !
  !       grad(u) = (du/dx, du/dy)
  !
  !   Il est utilisé pour calculer l’erreur en norme H1.
  !--------------------------------------------------------------------------
  subroutine exact_gradient(x, y, grad_u)
    real(dp), intent(in) :: x, y

    ! grad_u(1) = dérivée selon x.
    ! grad_u(2) = dérivée selon y.
    real(dp), intent(out) :: grad_u(2)

    grad_u(1) = pi * cos(pi * x) * sin(pi * y)
    grad_u(2) = pi * sin(pi * x) * cos(pi * y)

  end subroutine exact_gradient

end module exact_solution_module