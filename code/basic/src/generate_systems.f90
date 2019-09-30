program generate_systems
  use mod_param, only: wp,&
       maxchar
  use mod_random, only: init_rand
  use mod_matops, only: prettyprint
  use mod_system_generator, only: generateJacobisys
  use mod_io, only: write_system_to_file
  implicit none
  integer, parameter :: iseed=42
  real(wp), parameter :: too_small=1.0e-3_wp

  character(maxchar) :: filename
  integer :: n
  real(wp) :: xlow
  real(wp) :: xhigh
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: b(:,:)
  

  call init_rand(iseed)



  write(*,*) "Enter n: "
  read(*,*) n

  if(allocated(A)) deallocate(A)
  allocate(A(n,n))
  if(allocated(x)) deallocate(x)
  allocate(x(n,1))
  if(allocated(b)) deallocate(b)
  allocate(b(n,1))
  
  write(*,*) "Enter xlow, xhigh:"
  read(*,*) xlow, xhigh

  call generateJacobisys(n,xlow,xhigh,too_small,A,x,b)

  write(*,*) "A: "
  call prettyprint(n,n,A)
  write(*,*) "x: "
  call prettyprint(n,1,x)
  write(*,*) "b: "
  call prettyprint(n,1,b)

  write(*,*) "Enter filename to write to: "
  read(*,*) filename

  call write_system_to_file(filename,n,A,x,b)
  

  if(allocated(A)) deallocate(A)
  if(allocated(x)) deallocate(x)
  if(allocated(b)) deallocate(b)
  

end program generate_systems
