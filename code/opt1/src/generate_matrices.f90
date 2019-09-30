program generate_matrices
  use mod_param, only: wp, &
       zero
  use mod_random, only: init_rand
  use mod_matops, only: prettyprint,&
       generateDiagDomMat
  implicit none
  integer, parameter :: iseed=42

  integer :: i,j
  integer :: n
  real(wp) :: xlow, xhigh
  real(wp) :: xsum
  real(wp), allocatable :: A(:,:)


  write(*,*) "n: "
  read(*,*) n

  if(allocated(A)) deallocate(A)
  allocate(A(n,n))


  write(*,*) "Enter limits for non-diagonal entries (a,b): "
  read(*,*) xlow, xhigh

  call init_rand(iseed)
  call generateDiagDomMat(n, xlow, xhigh, A)
  call prettyprint(n,n,A)

  do i = 1, n
     xsum = zero
     do j = 1, n
        if (i  .ne. j) then
           xsum = xsum + abs(A(i,j))
        end if
     end do
     write(*,*) "|A(",i,",",i,")| - sum(| A(",i," j)|): ", abs(A(i,i))-xsum
  end do

     

  if(allocated(A)) deallocate(A)

  
end program generate_matrices
