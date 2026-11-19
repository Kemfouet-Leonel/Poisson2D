! Deux solutions manufacturees : trigonometrie et polynome biquadratique.
module exact
  use precision
  implicit none
contains
  pure real(dp) function value(x,y,case_id) result(u)
    real(dp),intent(in):: x,y
    integer,intent(in):: case_id
    if(case_id==1) then
      u=sin(pi*x)*sin(pi*y)
    else
      u=x*(1-x)*y*(1-y)
    end if
  end function
  pure real(dp) function source(x,y,case_id) result(f)
    real(dp),intent(in):: x,y
    integer,intent(in):: case_id
    if(case_id==1) then
      f=2*pi*pi*sin(pi*x)*sin(pi*y)
    else
      f=2*(x*(1-x)+y*(1-y))
    end if
  end function
  pure function gradient(x,y,case_id) result(g)
    real(dp),intent(in):: x,y
    integer,intent(in):: case_id
    real(dp):: g(2)
    if(case_id==1) then
      g=[pi*cos(pi*x)*sin(pi*y),pi*sin(pi*x)*cos(pi*y)]
    else
      g=[(1-2*x)*y*(1-y),x*(1-x)*(1-2*y)]
    end if
  end function
end module
