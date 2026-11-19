module sparse
  use precision
  use domain, only: mesh_t
  implicit none
  type :: csr_t
    integer :: n
    integer,allocatable :: row(:),col(:)
    real(dp),allocatable :: val(:)
  end type
contains
  subroutine pattern(m,a)
    type(mesh_t),intent(in):: m
    type(csr_t),intent(out):: a
    ! Ce maillage structure possede AU PLUS 7 voisins par sommet, diagonal inclus.
    ! Ce constructeur n'est pas destine a un maillage non structure quelconque.
    integer,allocatable :: neighbours(:,:),count(:)
    integer:: e,i,j,r,c,k,p,t
    allocate(neighbours(7,m%nv),count(m%nv))
    count=0; neighbours=0; a%n=m%nv
    do e=1,m%nt
      do i=1,3
        r=m%tri(i,e)
        do j=1,3
          c=m%tri(j,e)
          if(any(neighbours(1:count(r),r)==c)) cycle
          if(count(r)==7) error stop 'Valence >7: adapter constructeur CSR'
          count(r)=count(r)+1; neighbours(count(r),r)=c
        end do
      end do
    end do
    allocate(a%row(a%n+1)); a%row(1)=1
    do r=1,a%n
      a%row(r+1)=a%row(r)+count(r)
    end do
    allocate(a%col(a%row(a%n+1)-1),a%val(a%row(a%n+1)-1))
    a%val=0
    do r=1,a%n
      ! Trier les colonnes de chaque ligne (petit tri local).
      do i=2,count(r)
        t=neighbours(i,r); j=i-1
        do while(j>=1)
          if(neighbours(j,r)<=t) exit
          neighbours(j+1,r)=neighbours(j,r); j=j-1
        end do
        neighbours(j+1,r)=t
      end do
      p=a%row(r)
      do k=1,count(r)
        a%col(p)=neighbours(k,r); p=p+1
      end do
    end do
  end subroutine

  subroutine add(a,i,j,v)
    type(csr_t),intent(inout):: a
    integer,intent(in):: i,j
    real(dp),intent(in):: v
    integer:: k
    do k=a%row(i),a%row(i+1)-1
      if(a%col(k)==j) then
        a%val(k)=a%val(k)+v
        return
      end if
    end do
    error stop 'Entree absente du motif CSR'
  end subroutine

  subroutine matvec(a,x,y)
    type(csr_t),intent(in):: a
    real(dp),intent(in):: x(:)
    real(dp),intent(out):: y(:)
    integer:: i,k
    y=0
    do i=1,a%n
      do k=a%row(i),a%row(i+1)-1
        y(i)=y(i)+a%val(k)*x(a%col(k))
      end do
    end do
  end subroutine
end module
