program benchmark_solver
  use mod_param, only: wp,&
       zero,&
       one
  use mod_random, only: init_rand
  use mod_system_generator, only: generateJacobisys
  use mod_solver, only: jacobi_method
  implicit none
  integer, parameter :: iseed=42
  real(wp), parameter :: too_small=1.0e-2_wp
  real(wp), parameter :: xlow=-1.0e2_wp
  real(wp), parameter :: xhigh=1.0e2_wp
  real(wp), parameter :: xtol=1.0e-6
  logical, parameter :: silent=.TRUE.

  ! Timing data
  real :: start_time, end_time
  real :: delta_time


  integer :: n
  integer :: niters
  real(wp) :: xnops
  real(wp) :: error
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: xguess(:,:)
  real(wp), allocatable :: ans(:,:)

  xnops = 0

  call init_rand(iseed)

  write(*,*) "Enter n: "
  read(*,*) n

  if(allocated(A)) deallocate(A)
  allocate(A(n,n))
  if(allocated(x)) deallocate(x)
  allocate(x(n,1))
  if(allocated(xguess)) deallocate(xguess)
  allocate(xguess(n,1))
  if(allocated(ans)) deallocate(ans)
  allocate(ans(n,1))    
  if(allocated(b)) deallocate(b)
  allocate(b(n,1))

  call generateJacobisys(n,xlow,xhigh,too_small,A,x,b)

  call cpu_time(start_time)
  call jacobi_method(silent,1000,n,A,b,xguess,xtol,ans,error,niters)
  call cpu_time(end_time)
  delta_time = end_time-start_time

  write(*,*) "Number_iterations: ", niters
  write(*,*) "Error: ", error
  xnops = (real(niters+1,wp))*(((2.0e0_wp)*real(n,wp)*real(n,wp)*real(n,wp))-(real(n,wp)*real(n,wp)))
  write(*,*) "Number_of_ops: ",xnops
  write(*,*) "Total_executation_time: ",delta_time
  write(*,*) "Operations per second: ", xnops/delta_time
  
  


  if(allocated(A)) deallocate(A)
  if(allocated(x)) deallocate(x)
  if(allocated(xguess)) deallocate(xguess)
  if(allocated(ans)) deallocate(ans)
  if(allocated(b)) deallocate(b)

end program benchmark_solver

  

  
