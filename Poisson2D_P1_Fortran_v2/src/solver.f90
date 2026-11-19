!==============================================================================
! Fichier : src/solver.f90
!
! Rôle :
!   Résoudre A*x=b par gradient conjugué préconditionné par Jacobi.
!
! Hypothèses :
!   A doit être symétrique définie positive.
!
! Critère de convergence:
!   ||b-A*x||_2 / ||b||_2 <= rtol.
!
! Le second membre et la matrice ne sont pas modifiés.
!==============================================================================
module solver

  use precision, only: dp
  use sparse, only: csr_t, matvec

  implicit none

  private

  public :: pcg

contains

  subroutine pcg(a, b, x, rtol, maxit, it, relres)

    type(csr_t), intent(in) :: a

    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: rtol

    integer, intent(in) :: maxit

    real(dp), intent(out) :: x(:)
    real(dp), intent(out) :: relres

    integer, intent(out) :: it

    real(dp), allocatable :: r(:)
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: p(:)
    real(dp), allocatable :: ap(:)
    real(dp), allocatable :: diag(:)

    real(dp) :: rho
    real(dp) :: rho_new
    real(dp) :: alpha
    real(dp) :: beta
    real(dp) :: denominator
    real(dp) :: bnorm

    integer :: i, k, iteration

    if (size(b) /= a%n .or. size(x) /= a%n) then
      error stop 'PCG : dimensions incompatibles'
    end if

    if (rtol <= 0.0_dp .or. rtol >= 1.0_dp) then
      error stop 'PCG : rtol doit appartenir a ]0,1['
    end if

    if (maxit < 1) then
      error stop 'PCG : maxit doit etre positif'
    end if

    allocate(r(a%n))
    allocate(z(a%n))
    allocate(p(a%n))
    allocate(ap(a%n))
    allocate(diag(a%n))

    ! Extraire la diagonale pour le préconditionneur de Jacobi.
    diag = 0.0_dp

    do i = 1, a%n

      do k = a%row(i), a%row(i+1)-1

        if (a%col(k) == i) then
          diag(i) = a%val(k)
        end if

      end do

    end do

    if (any(diag <= 0.0_dp)) then
      error stop 'PCG : diagonale absente ou non positive'
    end if

    ! Initialisation avec x_0=0.
    x = 0.0_dp
    r = b

    it = 0

    ! norm2 évite les limites d'un calcul naïf sqrt(sum(b*b))
    ! pour certaines amplitudes extrêmes.
    bnorm = norm2(b)

    ! Cas nul ou extrêmement petit.
    if (bnorm <= tiny(1.0_dp)) then
      relres = 0.0_dp
      return
    end if

    relres = norm2(r) / bnorm

    ! Application de M^{-1}, M étant la diagonale de A.
    z = r / diag
    p = z

    rho = dot_product(r, z)

    if (rho <= 0.0_dp) then
      error stop 'PCG : produit r.M^-1.r non positif'
    end if

    do iteration = 1, maxit

      it = iteration

      call matvec(a, p, ap)

      denominator = dot_product(p, ap)

      if (denominator <= 0.0_dp) then
        error stop 'PCG : courbure non positive'
      end if

      alpha = rho / denominator

      x = x + alpha * p
      r = r - alpha * ap

      relres = norm2(r) / bnorm

      if (relres <= rtol) then

        ! Le résidu récursif peut dériver par arrondi.
        ! Recalculer explicitement b-A*x avant de conclure.
        call matvec(a, x, ap)

        r = b - ap
        relres = norm2(r) / bnorm

        if (relres <= rtol) return

        ! Si nécessaire, redémarrer depuis le vrai résidu.
        z = r / diag
        p = z

        rho = dot_product(r, z)

        if (rho <= 0.0_dp) then
          error stop 'PCG : redemarrage impossible'
        end if

        cycle

      end if

      z = r / diag

      rho_new = dot_product(r, z)

      if (rho_new <= 0.0_dp) then
        error stop 'PCG : produit preconditionne non positif'
      end if

      beta = rho_new / rho

      p = z + beta * p
      rho = rho_new

    end do

    ! Calculer une dernière fois le vrai résidu pour le diagnostic.
    call matvec(a, x, ap)

    relres = norm2(b-ap) / bnorm

    if (relres <= rtol) return

    write(*,*) 'PCG : iterations maximales = ', maxit
    write(*,*) 'PCG : residu relatif final = ', relres

    error stop 'PCG : convergence non atteinte'

  end subroutine pcg

end module solver