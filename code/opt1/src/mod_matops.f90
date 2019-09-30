module mod_matops
  use mod_param, only: &
       wp, &
       zero,&
       one
  use mod_random, only: generate_random_a_b,&
       init_rand
  private

  public :: prettyprint
  public :: prettyentry
  public :: dgemm
  public :: dgemm_tiled
  public :: strictupperpart
  public :: strictlowerpart
  public :: diagonalpart
  public :: LDUdecompose
  public :: generateDiagDomMat
  public :: generateRandomMatrix


contains

  subroutine prettyprint(m,n,A)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: A(m,n)

    integer :: i 

    do i = 1, m
       write(*,*) A(i,:)
    end do

  end subroutine prettyprint

  subroutine prettyentry(m,n,A)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(inout) :: A(m,n)

    integer :: i
    do i = 1, m
       write(*,*) "Enter row ",i 
       read(*,*) A(i,:)
    end do

  end subroutine prettyentry

  ! Subroutine to calculate 
  ! C = alpha*A*B + beta*C
  ! where A is a m*k matrix
  !       B is a k*n matrix
  !       C is a m*n matrix
  subroutine dgemm(m,n,k,alpha,beta,A,B,C)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    real(wp), intent(in) :: A(m,k)
    real(wp), intent(in) :: B(k,n)
    real(wp), intent(inout) :: C(m,n)

    integer :: i, j, l

    
    do i = 1, m
       do j = 1, n
          do l = 1, k
             C(i,j) = (alpha*(A(i,l)*B(l,j)))+(beta*C(i,j))
          end do
       end do
    end do
  end subroutine dgemm


  subroutine dgemm_tiled(m,n,k,tile_size,alpha,beta,A,B,C)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: tile_size
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    real(wp), intent(in) :: A(m,k)
    real(wp), intent(in) :: B(k,n)
    real(wp), intent(inout) :: C(m,n)

    integer :: i0, j0, k0
    integer :: i1, j1, k1

    !We are implenenting a titled version to improve performance by
    !improving 
    i0lp: do i0 = 1, n, tile_size
      j0lp: do j0 = 1, m, tile_size
          k0lp:do k0 = 1, k, tile_size
             i1lp: do i1 = i0, min(i0+tile_size-1,m)
                j1lp: do j1 = j0, min(j0+tile_size-1,n)
                   k1lp: do k1 = k0, min(k0+tile_size-1,k)
                      C(i1,j1) = (alpha*(A(i1,k1)*B(k1,j1)))+(beta*C(i1,j1))
                   end do k1lp
                end do j1lp
             end do i1lp
          end do k0lp
       end do j0lp
    end do i0lp
  end subroutine dgemm_tiled

  ! Returns the strict upper of the Matrix A as U
  subroutine strictupperpart(m,n,A,U)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: A(m,n)
    real(wp), intent(out) :: U(m,n)

    integer :: i,j

    U = A
    do i = 1, m
       do j = 1, n
          if (i >= j) then
             U(i,j) = zero
          end if
       end do
    end do 
  end subroutine strictupperpart

  ! Returns the strict lower part of the Matrix A as L
  subroutine strictlowerpart(m,n,A,L)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: A(m,n)
    real(wp), intent(out) :: L(m,n)

    integer :: i,j

    L = A
    do i = 1, m
       do j = 1, n
          if (i <= j) then
             L(i,j) = zero
          end if
       end do
    end do 
  end subroutine strictlowerpart

  ! Returns the diagonal part of the Matrix A as D
    subroutine diagonalpart(m,n,A,D)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: A(m,n)
    real(wp), intent(out) :: D(m,n)

    integer :: i

    D = zero
    do i = 1, m
       D(i,i) = A(i,i)
    end do
  end subroutine diagonalpart

  ! Returns L,D,U of matrix A
  subroutine LDUdecompose(m,n,A,L,D,U)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: A(m,n)
    real(wp), intent(out) :: L(m,n)
    real(wp), intent(out) :: D(m,n)
    real(wp), intent(out) :: U(m,n)

    call strictlowerpart(m,n,A,L)
    call diagonalpart(m,n,A,D)
    call strictupperpart(m,n,A,U)
  end subroutine LDUdecompose

  subroutine generateDiagDomMat(n,xlow,xhigh,excludeSmallerThan,A)
    implicit none

    !integer, intent(in) :: iseed
    integer, intent(in) :: n
    real(wp), intent(in) :: xlow
    real(wp), intent(in) :: xhigh
    real(wp), intent(in) :: excludeSmallerThan
    real(wp), intent(out) :: A(n,n)

    integer :: i, j
    real(wp) :: xval
    real(wp) :: xsum

    !call init_rand(iseed)

    A = zero

    do i = 1, n
       do j = 1, n
          if (i .ne. j) then
             do 
                call generate_random_a_b(xlow,xhigh,xval)
                A(i,j) = xval
                if (abs(A(i,j)) > excludeSmallerThan) then
                   exit
                end if
             end do
          end if
       end do
    end do

    ! Strictly Diagonally Dominant Matrix  Definition
    ! A square matrix is diagonally domiant if
    ! |a_ii| > sum_{i=j} |a_ij| for all i
    do i = 1, n
       xsum = zero
       do j = 1, n
          if (i .ne. j) then
             xsum = xsum + abs(A(i,j))
          end if
       end do
       call generate_random_a_b(zero,max(abs(xlow),abs(xhigh)),xval)
       A(i,i) = xsum + xval
       ! Flip the sign to get negatives some of the time. 
       if (xval <= 0.5) then
          A(i,i) = -A(i,i)
       end if
    end do
  end subroutine generateDiagDomMat

  subroutine generateRandomMatrix(m,n,xlow,xhigh,excludeSmallerThan,A)
    implicit none

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: xlow
    real(wp), intent(in) :: xhigh
    real(wp), intent(in) :: excludeSmallerThan
    real(wp), intent(out) :: A(m,n)

    integer :: i, j
    real(wp) :: myrand

    A =  zero
    do i = 1, m
       do j = 1, n
          do 
             call generate_random_a_b(xlow,xhigh,myrand)
             A(i,j) = myrand
             if (abs(A(i,j)) > excludeSmallerThan) then
                exit
             end if
          end do
       end do
    end do
  end subroutine generateRandomMatrix


end module mod_matops
