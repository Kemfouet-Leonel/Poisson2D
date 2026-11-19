module precision
  use iso_fortran_env, only: real64, int64
  implicit none
  integer, parameter :: dp=real64
  real(dp), parameter :: pi=acos(-1.0_dp)
end module
