program test_random
  use mod_param, only: wp
  use mod_random, only: init_rand,&
       generate_random_a_b
  implicit none

  integer :: i
  integer :: iseed
  integer :: nrand
  real(wp) :: a,b
  real(wp) :: xval

  write(*,*) "Enter seed: "
  read(*,*) iseed
  call init_rand(iseed)

  write(*,*) "Enter number of random numbers: "
  read(*,*) nrand

  write(*,*) "Limits a and b:"
  read(*,*) a, b

  do i = 1, nrand
     call generate_random_a_b(a,b,xval)
     write(*,*) i, xval
  end do
  


end program test_random
