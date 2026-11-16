!==============================================================================
! Module : params_module
! Rôle :
!   Ce module contient les paramètres globaux du problème de Poisson 2D.
!   Il définit notamment :
!     - la précision des réels,
!     - les dimensions du maillage,
!     - les dimensions physiques du domaine,
!     - le nombre de nœuds, de triangles et de degrés de liberté.
!
!   Le problème résolu est :
!
!       -Delta(u) = f   dans Omega = (0,1) x (0,1)
!       u = 0           sur la frontière de Omega
!
!==============================================================================
module params_module

  implicit none

  integer, parameter :: dp = kind(1.0d0)  ! réel double précision
  real(dp), parameter :: pi = acos(-1.0_dp)  ! valeur de pi

  ! Nombre de subdivisions du maillage selon x et y.
  integer :: nx
  integer :: ny

  ! Nombre total de sommets du maillage.
  integer :: n_nodes

  ! Nombre total de triangles du maillage.
  integer :: n_elements

  ! Nombre total de degrés de liberté.
  integer :: n_dof

  ! Dimensions physiques du domaine.
  real(dp) :: lx
  real(dp) :: ly

contains

  !--------------------------------------------------------------------------
  ! Subroutine : set_parameters
  ! Rôle :
  !   Initialise les paramètres du problème à partir du nombre de
  !   subdivisions demandé dans les directions x et y.
  !
  ! Entrée :
  !   nx_in : nombre de subdivisions selon x
  !   ny_in : nombre de subdivisions selon y
  !--------------------------------------------------------------------------
  subroutine set_parameters(nx_in, ny_in)
    integer, intent(in) :: nx_in, ny_in

    ! Nombre de subdivisions dans chaque direction.
    nx = nx_in
    ny = ny_in

    ! Le domaine est le carré unité : [0,1] x [0,1].
    lx = 1.0_dp
    ly = 1.0_dp

    ! Si le maillage comporte nx subdivisions selon x,
    ! il possède nx + 1 nœuds selon x.
    ! De même pour y.
    n_nodes = (nx + 1) * (ny + 1)

    ! Chaque rectangle élémentaire est divisé en deux triangles.
    ! Il y a donc 2 * nx * ny triangles au total.
    n_elements = 2 * nx * ny

    ! En éléments finis P1, l’inconnue est associée à chaque sommet.
    n_dof = n_nodes

  end subroutine set_parameters

end module params_module