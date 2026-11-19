module domain
  use precision
  implicit none
  type :: mesh_t
    integer :: n,nv,nt
    real(dp), allocatable :: xy(:,:)
    integer, allocatable :: tri(:,:)
    logical, allocatable :: boundary(:)
  end type
contains
  subroutine build_mesh(n,m)
    integer,intent(in):: n
    type(mesh_t),intent(out):: m
    integer:: i,j,k,p,q
    if(n<2) error stop 'n doit etre >=2'
    m%n=n; m%nv=(n+1)**2; m%nt=2*n*n
    allocate(m%xy(2,m%nv),m%tri(3,m%nt),m%boundary(m%nv))
    do j=0,n
      do i=0,n
        p=j*(n+1)+i+1
        m%xy(:,p)=[real(i,dp)/n,real(j,dp)/n]
        m%boundary(p)=i==0.or.i==n.or.j==0.or.j==n
      end do
    end do
    k=0
    do j=0,n-1
      do i=0,n-1
        p=j*(n+1)+i+1; q=p+n+1
        k=k+1; m%tri(:,k)=[p,p+1,q]
        k=k+1; m%tri(:,k)=[p+1,q+1,q]
      end do
    end do
  end subroutine

  subroutine geometry(x,area,g)
    real(dp),intent(in):: x(2,3)
    real(dp),intent(out):: area,g(2,3)
    real(dp):: det
    det=(x(1,2)-x(1,1))*(x(2,3)-x(2,1)) &
       -(x(1,3)-x(1,1))*(x(2,2)-x(2,1))
    if(abs(det)<=tiny(1.0_dp)) error stop 'Triangle degenere'
    area=abs(det)/2
    ! Le determinant SIGNE rend les gradients corrects dans les deux orientations.
    g(1,:)=[x(2,2)-x(2,3),x(2,3)-x(2,1),x(2,1)-x(2,2)]/det
    g(2,:)=[x(1,3)-x(1,2),x(1,1)-x(1,3),x(1,2)-x(1,1)]/det
  end subroutine
end module
