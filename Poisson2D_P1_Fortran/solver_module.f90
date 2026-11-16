!==============================================================================
! Module : solver_module
! Rôle :
!   Résoudre le système linéaire :
!
!       A U = b.
!
!   Ici, on utilise la routine LAPACK dgesv.
!
!   dgesv résout un système linéaire dense réel en double précision
!   par factorisation LU avec pivotage partiel.
!==============================================================================
module solver_module

  use params_module

  implicit none

contains

  !--------------------------------------------------------------------------
  ! Subroutine : solve_linear_system
  ! Rôle :
  !   Résoudre A x = b.
  !
  ! Entrées :
  !   a : matrice du système.
  !   b : second membre.
  !
  ! Sortie :
  !   x : solution du système.
  !--------------------------------------------------------------------------
  subroutine solve_linear_system(a, b, x)
    real(dp), intent(inout) :: a(:, :)
    real(dp), intent(inout) :: b(:)
    real(dp), intent(out) :: x(:)

    integer :: n, info

    ! Tableau des pivots utilisé par LAPACK.
    integer, allocatable :: ipiv(:)

    ! Taille du système.
    n = size(b)

    allocate(ipiv(n))

    ! Appel à LAPACK :
    !
    ! dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
    !
    ! Ici :
    !   n    : taille du système,
    !   1    : un seul second membre,
    !   a    : matrice A,
    !   n    : dimension principale de A,
    !   ipiv : tableau des pivots,
    !   b    : second membre, remplacé par la solution,
    !   n    : dimension principale de b,
    !   info : code d’erreur.
    call dgesv(n, 1, a, n, ipiv, b, n, info)

    ! Si info /= 0, la factorisation ou la résolution a échoué.
    if (info /= 0) then
       write(*, *) "Erreur LAPACK : dgesv a échoué, info =", info
       stop 1
    end if

    ! Après l’appel, b contient la solution x.
    x = b

  end subroutine solve_linear_system

end module solver_module