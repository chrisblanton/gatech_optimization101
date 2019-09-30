program jacobi_allinone
  use mod_param, only: wp,&
       one,&
       zero
  use mod_matops, only: dgemm,&
       LDUdecompose,&
       prettyprint,&
       prettyentry
  implicit none

  integer :: i
  integer :: n
  real(wp) :: xnorm
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: L(:,:)
  real(wp), allocatable :: D(:,:)
  real(wp), allocatable :: U(:,:)
  real(wp), allocatable :: Dinverse(:,:)
  real(wp), allocatable :: Q(:,:)
  real(wp), allocatable :: C(:,:)
  real(wp), allocatable :: xk(:,:)
  real(wp), allocatable :: xkp1(:,:)

  write(*,*) "Enter n:"
  read(*,*) n
  
  if(allocated(A)) deallocate(A)
  if(allocated(b)) deallocate(b)
  if(allocated(L)) deallocate(L)
  if(allocated(D)) deallocate(D)
  if(allocated(U)) deallocate(U)
  if(allocated(Dinverse)) deallocate(Dinverse)
  if(allocated(Q)) deallocate(Q)
  if(allocated(C)) deallocate(C)
  if(allocated(xk)) deallocate(xk)
  if(allocated(xkp1)) deallocate(xkp1)

  allocate(A(n,n))
  allocate(b(n,1))
  allocate(L(n,n))
  allocate(D(n,n))
  allocate(U(n,n))
  allocate(Dinverse(n,n))
  allocate(Q(n,n))
  allocate(C(n,1))
  allocate(xk(n,1))
  allocate(xkp1(n,1))

  write(*,*) "A:"
  call prettyentry(n,n,A)
  write(*,*) "b: "
  call prettyentry(n,1,b)
  write(*,*) "Initial guess for x:"
  call prettyentry(n,1,xk)
  write(*,*) "The following were entered: "
  write(*,*) "A: "
  call prettyprint(n,n,A)
  write(*,*) "b: "
  call prettyprint(n,1,b)
  write(*,*) "xk: "
  call prettyprint(n,1,xk)

  call LDUdecompose(n,n,A,L,D,U)
  write(*,*) "L:"
  call prettyprint(n,n,L)
  write(*,*) "D:"
  call prettyprint(n,n,D)
  write(*,*) "U:"
  call prettyprint(n,n,U)
  write(*,*) "L+U:"
  call prettyprint(n,n,L+U)
  
 Dinverse = 0.0e0_wp
 do i = 1,n
    Dinverse(i,i) = (one)/D(i,i)
 end do

 call dgemm(n,n,n,one,one,-Dinverse,L+U,Q)
 call dgemm(n,1,n,one,one,Dinverse,b,C)
  
 i = 0 
 do 
    write(*,*) "k = ", i
    write(*,*) "x_{",i,"}:"
    call prettyprint(n,1,xk)
    xkp1 = C
    call dgemm(n,1,n,one,one,Q,xk,xkp1)
    write(*,*) "x_{",i+1,"}:"
    call prettyprint(n,1,xkp1)
    i = i + 1
    xk = xkp1
    read(*,*)
  end do




  if(allocated(A)) deallocate(A)
  if(allocated(b)) deallocate(b)
  if(allocated(L)) deallocate(L)
  if(allocated(D)) deallocate(D)
  if(allocated(U)) deallocate(U)
  if(allocated(Dinverse)) deallocate(Dinverse)
  if(allocated(Q)) deallocate(Q)
  if(allocated(C)) deallocate(C)
  if(allocated(xk)) deallocate(xk)
  if(allocated(xkp1)) deallocate(xkp1)

end program jacobi_allinone
