!==============================================================================
! Fichier : tests/unit.f90
!
! Tests :
!   - poids et moments de quadrature ;
!   - maillage et frontière ;
!   - géométrie des triangles ;
!   - symétrie du stockage CSR ;
!   - solveur sur un système à solution connue ;
!   - second membre nul ;
!   - références analytiques des normes.
!==============================================================================
program unit_tests

  use precision, only: dp

  use quadrature, only: nq, rule
  use domain, only: mesh_t, build_mesh, geometry
  use sparse, only: csr_t, matvec
  use fem, only: assemble
  use solver, only: pcg
  use errors, only: compute_errors

  implicit none

  type(mesh_t) :: m
  type(csr_t) :: a

  real(dp) :: lam(3, nq)
  real(dp) :: w(nq)

  real(dp) :: q
  real(dp) :: ref
  real(dp) :: area
  real(dp) :: g(2, 3)

  real(dp) :: identity(2, 2)

  real(dp) :: l2
  real(dp) :: hs
  real(dp) :: h1
  real(dp) :: linf
  real(dp) :: rr

  real(dp), allocatable :: b(:)
  real(dp), allocatable :: u(:)
  real(dp), allocatable :: v(:)
  real(dp), allocatable :: av(:)
  real(dp), allocatable :: exact_u(:)

  integer :: i, j, k
  integer :: it
  integer :: idx
  integer :: e

  !---------------------------------------------------------------------------
  ! 1. Quadrature : poids.
  !---------------------------------------------------------------------------
  call rule(lam, w)

  call check( &
    abs(sum(w)-1.0_dp) < 1.0e-13_dp, &
    'somme des poids' &
  )

  call check( &
    all(w > 0.0_dp), &
    'positivite des poids' &
  )

  call check( &
    maxval(abs(sum(lam, dim=1)-1.0_dp)) < 1.0e-13_dp, &
    'somme des coordonnees barycentriques' &
  )

  call check( &
    all(lam >= 0.0_dp), &
    'points situes dans le triangle' &
  )

  !---------------------------------------------------------------------------
  ! 2. Quadrature : intégration des monômes x^i y^j.
  !
  ! Sur le triangle de référence :
  !
  ! intégrale x^i y^j = i! j! / (i+j+2)!.
  !---------------------------------------------------------------------------
  do i = 0, 8

    do j = 0, 8-i

      q = 0.5_dp * sum(w * lam(2,:)**i * lam(3,:)**j)

      ref = fact(i) * fact(j) / fact(i+j+2)

      call check( &
        abs(q-ref) < 1.0e-13_dp, &
        'exactitude des monomes de degre <=8' &
      )

    end do

  end do

  !---------------------------------------------------------------------------
  ! 3. Construction et comptage du maillage.
  !---------------------------------------------------------------------------
  call build_mesh(8, m)

  call check( &
    m%nv == 81 .and. m%nt == 128, &
    'nombre de noeuds et triangles' &
  )

  call check( &
    count(m%boundary) == 32, &
    'nombre de noeuds de bord' &
  )

  !---------------------------------------------------------------------------
  ! 4. Géométrie et fonctions P1.
  !---------------------------------------------------------------------------
  identity = reshape( &
    [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
    [2, 2] &
  )

  q = 0.0_dp

  do e = 1, m%nt

    call geometry(m%xy(:,m%tri(:,e)), area, g)

    q = q + area

    ! Puisque sum(phi_i)=1, la somme des gradients doit être nulle.
    call check( &
      maxval(abs(sum(g, dim=2))) < 1.0e-13_dp, &
      'partition de l''unite : gradients' &
    )

    ! Les fonctions P1 reproduisent exactement les coordonnées x et y.
    call check( &
      maxval(abs( &
        matmul(m%xy(:,m%tri(:,e)), transpose(g)) - identity &
      )) < 1.0e-13_dp, &
      'reproduction des coordonnees affines' &
    )

  end do

  call check( &
    abs(q-1.0_dp) < 1.0e-13_dp, &
    'aire totale du domaine' &
  )

  !---------------------------------------------------------------------------
  ! 5. Assemblage et symétrie.
  !---------------------------------------------------------------------------
  call assemble(m, 2, a, b)

  allocate(u(m%nv))
  allocate(v(m%nv))
  allocate(av(m%nv))
  allocate(exact_u(m%nv))

  do i = 1, a%n

    do k = a%row(i), a%row(i+1)-1

      j = a%col(k)

      ! Rechercher A(j,i) pour le comparer à A(i,j).
      ref = 0.0_dp

      do idx = a%row(j), a%row(j+1)-1

        if (a%col(idx) == i) then
          ref = a%val(idx)
        end if

      end do

      call check( &
        abs(a%val(k)-ref) < 1.0e-13_dp, &
        'symetrie CSR apres Dirichlet' &
      )

    end do

  end do

  !---------------------------------------------------------------------------
  ! 6. Test du solveur indépendant de la solution de l'EDP.
  !
  ! Choisir un vecteur x_exact, fabriquer b=A*x_exact,
  ! puis vérifier que le solveur retrouve ce vecteur.
  !---------------------------------------------------------------------------
  do i = 1, m%nv
    exact_u(i) = sin(real(i, dp))
  end do

  call matvec(a, exact_u, v)

  call check( &
    dot_product(exact_u, v) > 0.0_dp, &
    'energie positive pour le vecteur test' &
  )

  call pcg(a, v, u, 1.0e-12_dp, 1000, it, rr)

  call check( &
    maxval(abs(u-exact_u)) < 1.0e-9_dp, &
    'PCG : solution algebrique connue' &
  )

  ! Vérifier le résidu indépendamment du scalaire renvoyé par PCG.
  call matvec(a, u, av)

  call check( &
    norm2(v-av)/norm2(v) < 1.0e-12_dp, &
    'PCG : vrai residu' &
  )

  !---------------------------------------------------------------------------
  ! 7. Cas particulier du second membre nul.
  !---------------------------------------------------------------------------
  v = 0.0_dp

  call pcg(a, v, u, 1.0e-12_dp, 1000, it, rr)

  call check( &
    maxval(abs(u)) <= tiny(1.0_dp) .and. it == 0, &
    'PCG : second membre nul' &
  )

  !---------------------------------------------------------------------------
  ! 8. Références analytiques des normes.
  !
  ! On impose u_h=0 et on prend :
  !
  ! u = x(1-x)y(1-y).
  !
  ! Alors :
  !   ||u||_L2       = 1/30 ;
  !   |u|_H1         = sqrt(1/45) ;
  !   ||u||_H1       = sqrt(7/300) ;
  !   max(u)         = 1/16.
  !
  ! Le maillage n=8 contient le point central où le maximum est atteint.
  !---------------------------------------------------------------------------
  u = 0.0_dp

  call compute_errors(m, u, 2, l2, hs, h1, linf)

  call check( &
    abs(l2-1.0_dp/30.0_dp) < 1.0e-12_dp, &
    'norme L2 : reference polynomiale' &
  )

  call check( &
    abs(hs-sqrt(1.0_dp/45.0_dp)) < 1.0e-12_dp, &
    'semi norme H1 : reference polynomiale' &
  )

  call check( &
    abs(h1-sqrt(7.0_dp/300.0_dp)) < 1.0e-12_dp, &
    'norme H1 complete : reference polynomiale' &
  )

  call check( &
    abs(h1*h1-l2*l2-hs*hs) < 1.0e-14_dp, &
    'identite entre norme et semi-norme H1' &
  )

  call check( &
    abs(linf-1.0_dp/16.0_dp) < 1.0e-13_dp, &
    'maximum echantillonne polynomial' &
  )

  write(*,*) '--------------------------------------------------'
  write(*,*) 'Tous les tests unitaires ont reussi.'
  write(*,*) '--------------------------------------------------'

contains

  subroutine check(ok, message)

    logical, intent(in) :: ok
    character(len=*), intent(in) :: message

    if (.not. ok) then
      write(*,*) 'ECHEC : ', message
      error stop 1
    end if

  end subroutine check

  real(dp) function fact(n) result(f)

    integer, intent(in) :: n

    integer :: k

    f = 1.0_dp

    do k = 2, n
      f = f * real(k, dp)
    end do

  end function fact

end program unit_tests
