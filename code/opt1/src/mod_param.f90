module mod_param
  implicit none

  private

  public :: wp
  public :: zero
  public :: one
  public :: maxchar
  public :: inpunit
  public :: outpunit


  integer, parameter :: wp=8
  integer, parameter :: maxchar=100
  integer, parameter :: inpunit=20
  integer, parameter :: outpunit=21
  real(wp), parameter :: zero=0.00e0_wp
  real(wp), parameter :: one=1.00e0_wp
  

end module mod_param
