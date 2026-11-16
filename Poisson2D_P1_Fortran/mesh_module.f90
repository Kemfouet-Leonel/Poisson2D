!==============================================================================
! Module : mesh_module
! Rôle :
!   Construire le maillage triangulaire structuré du carré unité.
!
!   Le maillage est obtenu en découpant le domaine en rectangles,
!   puis chaque rectangle est divisé en deux triangles.
!
!   Ce module fournit :
!     - les coordonnées des nœuds,
!     - la connectivité des triangles,
!     - la liste des nœuds situés sur la frontière.
!==============================================================================
module mesh_module

  ! Le module dépend des paramètres globaux définis dans params_module.
  use params_module

  implicit none

  ! Tableau des coordonnées des nœuds.
  ! nodes(1, k) = abscisse x du nœud k.
  ! nodes(2, k) = ordonnée y du nœud k.
  real(dp), allocatable :: nodes(:, :)

  ! Tableau de connectivité des triangles.
  ! elements(1, e), elements(2, e), elements(3, e)
  ! sont les trois numéros globaux des sommets du triangle e.
  integer, allocatable :: elements(:, :)

  ! Liste des numéros globaux des nœuds situés sur la frontière.
  integer, allocatable :: boundary_nodes(:)

  ! Nombre de nœuds situés sur la frontière du domaine.
  integer :: n_boundary_nodes

contains

  !--------------------------------------------------------------------------
  ! Subroutine : build_mesh
  ! Rôle :
  !   Construit :
  !     1. les coordonnées des nœuds,
  !     2. la liste des triangles,
  !     3. la liste des nœuds de bord.
  !--------------------------------------------------------------------------
  subroutine build_mesh()
    integer :: i, j, k
    real(dp) :: hx, hy

    ! Allocation dynamique des tableaux du maillage.
    allocate(nodes(2, n_nodes))
    allocate(elements(3, n_elements))

    ! Pas du maillage selon x et y.
    ! Le maillage est uniforme.
    hx = lx / real(nx, dp)
    hy = ly / real(ny, dp)

    !------------------------------------------------------------------------
    ! Construction des coordonnées des nœuds.
    !
    ! Les nœuds sont numérotés ligne par ligne :
    !
    !   (0,1) --- (1,1)
    !     |          |
    !   (0,0) --- (1,0)
    !
    !------------------------------------------------------------------------
    k = 0

    do j = 0, ny
       do i = 0, nx
          k = k + 1

          nodes(1, k) = real(i, dp) * hx
          nodes(2, k) = real(j, dp) * hy
       end do
    end do

    !------------------------------------------------------------------------
    ! Construction des triangles.
    !
    ! Chaque rectangle élémentaire [i,i+1] x [j,j+1]
    ! est divisé en deux triangles.
    !------------------------------------------------------------------------
    k = 0

    do j = 0, ny - 1
       do i = 0, nx - 1

          ! Premier triangle du rectangle.
          k = k + 1
          elements(1, k) = node_index(i, j)
          elements(2, k) = node_index(i + 1, j)
          elements(3, k) = node_index(i, j + 1)

          ! Second triangle du rectangle.
          k = k + 1
          elements(1, k) = node_index(i + 1, j)
          elements(2, k) = node_index(i + 1, j + 1)
          elements(3, k) = node_index(i, j + 1)

       end do
    end do

    ! Identification des nœuds situés sur la frontière.
    call build_boundary()

  end subroutine build_mesh

  !--------------------------------------------------------------------------
  ! Fonction : node_index
  ! Rôle :
  !   Retourne le numéro global d’un nœud à partir de ses indices
  !   entiers (i, j) dans le maillage structuré.
  !
  !   i : indice selon x, avec 0 <= i <= nx
  !   j : indice selon y, avec 0 <= j <= ny
  !--------------------------------------------------------------------------
  function node_index(i, j) result(idx)
    integer, intent(in) :: i, j
    integer :: idx

    ! Numérotation ligne par ligne.
    idx = j * (nx + 1) + i + 1

  end function node_index

  !--------------------------------------------------------------------------
  ! Subroutine : build_boundary
  ! Rôle :
  !   Construit la liste des nœuds situés sur la frontière du carré.
  !
  !   Ces nœuds recevront la condition de Dirichlet homogène :
  !
  !       u = 0 sur la frontière.
  !--------------------------------------------------------------------------
  subroutine build_boundary()
    integer :: i, j, count

    ! Nombre de nœuds sur la frontière d’un rectangle maillé.
    n_boundary_nodes = 2 * (nx + ny)

    allocate(boundary_nodes(n_boundary_nodes))

    count = 0

    ! Bord inférieur : y = 0.
    do i = 0, nx
       count = count + 1
       boundary_nodes(count) = node_index(i, 0)
    end do

    ! Bord droit : x = 1.
    do j = 1, ny
       count = count + 1
       boundary_nodes(count) = node_index(nx, j)
    end do

    ! Bord supérieur : y = 1.
    do i = nx - 1, 0, -1
       count = count + 1
       boundary_nodes(count) = node_index(i, ny)
    end do

    ! Bord gauche : x = 0.
    do j = ny - 1, 1, -1
       count = count + 1
       boundary_nodes(count) = node_index(0, j)
    end do

  end subroutine build_boundary

end module mesh_module