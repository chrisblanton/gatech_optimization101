program matmul
  use mod_param, only: &
       wp,&
       one
  use mod_matops, only: dgemm
  implicit none
  
  ! Variables
  integer :: m, n, k1, k2
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: B(:,:)
  real(wp), allocatable :: C(:,:)

  ! Dummy variables
  integer :: i, j


  do 
  write(*,*) "What are dimensions of A?"
  read(*,*) m, k1
  write(*,*) "What are dimensions of b?"
  read(*,*) k2, n
  if (k1 == k2) then
     exit
  else
     write(*,*) "A: ", m, k1
     write(*,*) "B: ", k2, n
     write(*,*) "Matrix mulplication condition not met. Try again."
  end if
  end do

  if(allocated(A)) deallocate(A)
  if(allocated(B)) deallocate(B)
  if(allocated(C)) deallocate(C)
  allocate(A(m,k1))
  allocate(B(k2,n))
  allocate(C(m,n))
  
  C(:,:) = 0.0e0_wp
  write(*,*) C(:,:)
  do i =  1, m
     write(*,*) "Enter row ",i, " of A: "
     read(*,*) A(i,:)
  end do

  do i =  1, k2
     write(*,*) "Enter row ",i, " of B: "
     read(*,*) B(i,:)
  end do

  call dgemm(m,n,k1,one,one,A,B,C)
  write(*,*) "C: "
  do i = 1, m
     write(*,*) C(i,:)
  end do


  if(allocated(A)) deallocate(A)
  if(allocated(B)) deallocate(B)
  if(allocated(C)) deallocate(C)




end program matmul
