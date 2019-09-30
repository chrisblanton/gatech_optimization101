program test_solver
  use mod_param, only: wp,&
       maxchar,&
       zero
  use mod_solver, only: jacobi_method
  use mod_io, only: read_n_from_file,&
       read_system_from_file
  implicit none
  real(wp), parameter :: xtol=1.0e-6_wp
  logical, parameter :: silent=.FALSE.

  character(maxchar) :: filename
  integer :: n
  integer :: niters
  real(wp) :: error
  real(wp) :: xsum
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: xguess(:,:)
  real(wp), allocatable :: xtrue(:,:)  
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: ans(:,:)


  integer :: i, j

  write(*,*) "Enter filename: "
  read(*,*) filename

  call read_n_from_file(filename,n)

  if(allocated(A)) deallocate(A)
  allocate(A(n,n))
  if(allocated(xguess)) deallocate(xguess)
  allocate(xguess(n,1))
  if(allocated(xtrue)) deallocate(xtrue)
  allocate(xtrue(n,1))  
  if(allocated(ans)) deallocate(ans)
  allocate(ans(n,1))
  if(allocated(b)) deallocate(b)
  allocate(b(n,1))  

  call read_system_from_file(filename,n,A,xtrue,b)




  xguess = zero

  call jacobi_method(silent,1000,n,A,b,xguess,xtol,ans,error,niters)

  write(*,*) "Number_iterations: ", niters
  write(*,*) "Answer, xtrue, |Answer-xtrue| "
  xsum = zero
  do i = 1, n
     write(*,*) ans(i,1), xtrue(i,1), abs(ans(i,1)-xtrue(i,1))
     xsum = xsum + abs(ans(i,1)-xtrue(i,1))
  end do
  write(*,*) "Sum of Absolute Difference: ", xsum
  
  if(allocated(A)) deallocate(A)
  if(allocated(xguess)) deallocate(xguess)
  if(allocated(xtrue)) deallocate(xtrue)  
  if(allocated(b)) deallocate(b)
  if(allocated(ans)) deallocate(ans)


end program test_solver
