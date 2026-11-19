!==============================================================================
! Fichier : src/fem.f90
!
! Rôle :
!   Assembler le système éléments finis P1 :
!
!       A U = b.
!
! Méthode :
!   - rigidité calculée exactement sur chaque triangle P1 ;
!   - second membre intégré par quadrature ;
!   - matrice globale stockée en CSR ;
!   - Dirichlet homogène imposé en conservant la symétrie.
!==============================================================================
module fem

  use precision, only: dp
  use domain, only: mesh_t, geometry
  use quadrature, only: nq, rule
  use exact, only: source
  use sparse, only: csr_t, pattern, add

  implicit none

  private

  public :: assemble

contains

  subroutine assemble(m, case_id, a, b)

    type(mesh_t), intent(in) :: m

    integer, intent(in) :: case_id

    type(csr_t), intent(out) :: a

    real(dp), allocatable, intent(out) :: b(:)

    ! Coordonnées barycentriques et poids de quadrature.
    real(dp) :: lam(3, nq)
    real(dp) :: w(nq)

    ! Coordonnées des trois sommets du triangle.
    real(dp) :: x(2, 3)

    ! Gradients des fonctions de base P1.
    real(dp) :: g(2, 3)

    ! Matrice et second membre élémentaires.
    real(dp) :: ke(3, 3)
    real(dp) :: fe(3)

    real(dp) :: area
    real(dp) :: pt(2)
    real(dp) :: fq

    integer :: ids(3)
    integer :: e, q, i, j
    integer :: r, c

    ! Préparer la règle d'intégration et le motif CSR.
    call rule(lam, w)
    call pattern(m, a)

    allocate(b(m%nv))
    b = 0.0_dp

    do e = 1, m%nt

      ! Numéros globaux des trois sommets.
      ids = m%tri(:, e)

      ! Coordonnées physiques de ces sommets.
      x = m%xy(:, ids)

      ! Aire et gradients des fonctions de base.
      call geometry(x, area, g)

      ! Les gradients P1 sont constants sur le triangle :
      !
      ! ke(i,j) = aire * grad(phi_i) . grad(phi_j).
      ke = area * matmul(transpose(g), g)

      fe = 0.0_dp

      do q = 1, nq

        ! Passage des coordonnées barycentriques au point physique.
        pt = matmul(x, lam(:, q))

        fq = source(pt(1), pt(2), case_id)

        ! Au point q, phi_i = lambda_i.
        fe = fe + area * w(q) * fq * lam(:, q)

      end do

      ! Assemblage des degrés de liberté intérieurs.
      do i = 1, 3

        r = ids(i)

        ! Valeur imposée nulle sur le bord.
        if (m%boundary(r)) cycle

        b(r) = b(r) + fe(i)

        do j = 1, 3

          c = ids(j)

          ! Ne pas assembler les couplages vers les nœuds de bord.
          ! Cela préserve la symétrie pour Dirichlet homogène.
          if (.not. m%boundary(c)) then
            call add(a, r, c, ke(i, j))
          end if

        end do

      end do

    end do

    ! Les équations de bord deviennent U_i = 0.
    do r = 1, m%nv

      if (m%boundary(r)) then
        call add(a, r, r, 1.0_dp)
      end if

    end do

  end subroutine assemble

end module fem