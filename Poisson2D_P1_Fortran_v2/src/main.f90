!==============================================================================
! Fichier : src/main.f90
!
! Valeurs par défaut :
!   n              = 20
!   case_id        = 1
!   dossier_sortie = results
!
! Exemples :
!   mkdir -p results
!   ./bin/poisson2d 32
!
!   mkdir -p results/manual
!   ./bin/poisson2d 32 1 results/manual
!
! Cas disponibles :
!   1 : u = sin(pi*x)*sin(pi*y)
!   2 : u = x*(1-x)*y*(1-y)
!==============================================================================
program poisson_v2

  use precision, only: dp, int64

  use domain, only: mesh_t, build_mesh
  use exact, only: value
  use sparse, only: csr_t
  use fem, only: assemble
  use solver, only: pcg
  use errors, only: compute_errors

  implicit none

  type(mesh_t) :: m
  type(csr_t) :: a

  real(dp), allocatable :: b(:)
  real(dp), allocatable :: u(:)

  integer :: n
  integer :: case_id
  integer :: ios
  integer :: unit
  integer :: i
  integer :: it

  integer(int64) :: c0, c1
  integer(int64) :: rate
  integer(int64) :: bytes

  real(dp) :: tm, ta, ts, te
  real(dp) :: l2, hs, h1, linf
  real(dp) :: rr
  real(dp) :: tol

  character(len=1024) :: arg
  character(len=1024) :: output
  character(len=1024) :: iomsg

  !---------------------------------------------------------------------------
  ! 1. Valeurs par défaut.
  !---------------------------------------------------------------------------
  n = 20
  case_id = 1

  output = 'results'
  tol = 1.0e-11_dp

  if (command_argument_count() > 3) then
    error stop 'Usage : poisson2d n case_id dossier_sortie'
  end if

  !---------------------------------------------------------------------------
  ! 2. Lecture des arguments.
  !---------------------------------------------------------------------------
  if (command_argument_count() >= 1) then

    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) n

    if (ios /= 0) then
      error stop 'Premier argument invalide : n doit etre entier'
    end if

  end if

  if (command_argument_count() >= 2) then

    call get_command_argument(2, arg)
    read(arg, *, iostat=ios) case_id

    if (ios /= 0) then
      error stop 'Deuxieme argument invalide : cas entier attendu'
    end if

  end if

  if (command_argument_count() >= 3) then
    call get_command_argument(3, output)
  end if

  if (n < 2 .or. n > 4096) then
    error stop 'n doit appartenir a [2,4096]'
  end if

  if (case_id /= 1 .and. case_id /= 2) then
    error stop 'Le cas manufacture doit etre 1 ou 2'
  end if

  if (len_trim(output) == 0) then
    error stop 'Le dossier de sortie ne peut pas etre vide'
  end if

  !---------------------------------------------------------------------------
  ! 3. Vérification anticipée du dossier et des permissions.
  !
  ! STATUS='UNKNOWN' permet de vérifier l'ouverture sans effacer
  ! immédiatement un ancien fichier de métriques.
  !---------------------------------------------------------------------------
  open( &
    newunit=unit, &
    file=trim(output)//'/metrics.csv', &
    status='unknown', &
    action='write', &
    iostat=ios, &
    iomsg=iomsg &
  )

  if (ios /= 0) then

    write(*,*) 'Impossible d''ecrire dans : ', trim(output)
    write(*,*) 'Diagnostic : ', trim(iomsg)
    write(*,*) 'Creer le dossier avant execution, par exemple :'
    write(*,*) 'mkdir -p results'
    write(*,*) 'Puis : ./bin/poisson2d 32 1 results'

    error stop 1

  end if

  close(unit)

  !---------------------------------------------------------------------------
  ! 4. Construction du maillage et chronométrage.
  !---------------------------------------------------------------------------
  call system_clock(count=c0, count_rate=rate)

  if (rate <= 0_int64) then
    error stop 'Horloge systeme indisponible'
  end if

  call build_mesh(n, m)

  call system_clock(count=c1)

  tm = real(c1-c0, dp) / real(rate, dp)
  c0 = c1

  !---------------------------------------------------------------------------
  ! 5. Assemblage CSR.
  ! Le temps comprend la construction du motif.
  !---------------------------------------------------------------------------
  call assemble(m, case_id, a, b)

  call system_clock(count=c1)

  ta = real(c1-c0, dp) / real(rate, dp)
  c0 = c1

  !---------------------------------------------------------------------------
  ! 6. Résolution.
  !---------------------------------------------------------------------------
  allocate(u(m%nv))

  call pcg( &
    a, b, u, &
    tol, &
    20*n, &
    it, rr &
  )

  call system_clock(count=c1)

  ts = real(c1-c0, dp) / real(rate, dp)
  c0 = c1

  !---------------------------------------------------------------------------
  ! 7. Calcul des erreurs.
  !---------------------------------------------------------------------------
  call compute_errors(m, u, case_id, l2, hs, h1, linf)

  call system_clock(count=c1)

  te = real(c1-c0, dp) / real(rate, dp)

  !---------------------------------------------------------------------------
  ! 8. Estimation du stockage des seuls tableaux CSR.
  !
  ! Ce n'est pas la mémoire totale ni le pic mémoire du processus.
  !---------------------------------------------------------------------------
  bytes = &
      int(size(a%val), int64) * int(storage_size(a%val(1))/8, int64) &
    + int(size(a%col), int64) * int(storage_size(a%col(1))/8, int64) &
    + int(size(a%row), int64) * int(storage_size(a%row(1))/8, int64)

  !---------------------------------------------------------------------------
  ! 9. Écriture des métriques.
  !---------------------------------------------------------------------------
  open( &
    newunit=unit, &
    file=trim(output)//'/metrics.csv', &
    status='replace', &
    action='write', &
    iostat=ios, &
    iomsg=iomsg &
  )

  if (ios /= 0) then
    write(*,*) trim(iomsg)
    error stop 'Impossible d''ouvrir metrics.csv'
  end if

  write(unit, '(a)') &
    'n,h,case_id,ndof,nnz,l2,h1_semi,h1,linf_sample,'// &
    'iterations,relres,mesh_s,assembly_s,solve_s,error_s,csr_bytes'

  write(unit, '(*(g0,:,","))') &
    n, 1.0_dp/real(n,dp), case_id, m%nv, size(a%val), &
    l2, hs, h1, linf, it, rr, tm, ta, ts, te, bytes

  close(unit)

  !---------------------------------------------------------------------------
  ! 10. Écriture de la solution.
  ! Colonnes : x, y, solution approchée, solution exacte.
  !---------------------------------------------------------------------------
  open( &
    newunit=unit, &
    file=trim(output)//'/solution.dat', &
    status='replace', &
    action='write', &
    iostat=ios, &
    iomsg=iomsg &
  )

  if (ios /= 0) then
    write(*,*) trim(iomsg)
    error stop 'Impossible d''ouvrir solution.dat'
  end if

  do i = 1, m%nv

    write(unit, '(4es25.16)') &
      m%xy(:,i), &
      u(i), &
      value(m%xy(1,i), m%xy(2,i), case_id)

  end do

  close(unit)

  !---------------------------------------------------------------------------
  ! 11. Résumé terminal.
  !---------------------------------------------------------------------------
  write(*,*) '--------------------------------------------------'
  write(*,*) 'Poisson 2D - elements finis P1 - stockage CSR'
  write(*,*) 'n                         : ', n
  write(*,*) 'Cas manufacture           : ', case_id
  write(*,*) 'Nombre de noeuds          : ', m%nv
  write(*,*) 'Nombre de triangles       : ', m%nt
  write(*,*) 'Entrees CSR stockees      : ', size(a%val)
  write(*,*) 'Iterations PCG            : ', it
  write(*,*) 'Vrai residu relatif       : ', rr
  write(*,*) 'Erreur L2                 : ', l2
  write(*,*) 'Erreur H1 semi-norme      : ', hs
  write(*,*) 'Erreur H1 complete        : ', h1
  write(*,*) 'Maximum echantillonne     : ', linf
  write(*,*) 'Temps maillage (s)        : ', tm
  write(*,*) 'Temps assemblage (s)      : ', ta
  write(*,*) 'Temps resolution (s)      : ', ts
  write(*,*) 'Temps erreurs (s)         : ', te
  write(*,*) 'Stockage CSR (octets)     : ', bytes
  write(*,*) 'Sorties                   : ', trim(output)
  write(*,*) '--------------------------------------------------'

end program poisson_v2