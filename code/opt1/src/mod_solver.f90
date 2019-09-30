module mod_solver
  use mod_param, only: wp, &
       one, &
       zero
  use mod_matops, only: dgemm_tiled,&
       LDUdecompose, &
       prettyprint
  private


  public :: jacobi_method_tiled


contains

  ! In the Jacobi Method,
  ! x_k = Q_{k-1} + C
  ! where
  ! Q = -D^{-1}(L+U)
  ! C = D^{-1}b
  ! D is the diagonal matrix A,
  ! L is the stricly lower triangular part of A,
  ! U is the stictly upper triangular part of A,
  ! x_k is the current version of the answer vector,
  ! x_{k-1} is the previous iteration of the answer vector (or the initial guess),
  ! A is the matrix given,
  ! b is the vector given.
  ! All matrixes are n*n
  subroutine jacobi_method_tiled(silent,itermax,tile_size,n,A,b,xbegin,xtol,ans,error,niters)
    implicit none

    logical, intent(in) :: silent
    integer, intent(in) :: itermax
    integer, intent(in) :: tile_size
    integer, intent(in) :: n
    real(wp), intent(in) :: A(n,n)
    real(wp), intent(in) :: b(n,1)
    real(wp), intent(in) :: xbegin(n,1)
    real(wp), intent(in) :: xtol
    real(wp), intent(out) :: ans(n,1)
    real(wp), intent(out) :: error
    integer, intent(out) :: niters

    integer :: i
    integer :: j
    integer :: k
    real(wp) :: L(n,n)
    real(wp) :: D(n,n)
    real(wp) :: U(n,n)
    real(wp) :: Dinverse(n,n)
    real(wp) :: xk(n,1)
    real(wp) :: xkp1(n,1)
    real(wp) :: Q(n,n)
    real(wp) :: C(n,1)
    real(wp) :: xnorm
    
    call LDUdecompose(n,n,A,L,D,U)

    Dinverse = zero
    do i = 1, n
       Dinverse(i,i) = one/D(i,i)
    end do

    call dgemm_tiled(n,n,n,tile_size,one,one,-Dinverse,L+U,Q)
    call dgemm_tiled(n,1,n,tile_size,one,one,Dinverse,b,C)
    xk = xbegin

    k = 0 
    if (silent) then 
       do i = 1, itermax
          xkp1 = C
          call dgemm_tiled(n,1,n,tile_size,one,one,Q,xk,xkp1)
          xnorm = zero
          do j = 1, n
             xnorm = xnorm + abs(xkp1(j,1)-xk(j,1))
          end do
          !write(*,*) "k: ", i
          k = i
          !call prettyprint(n,1,xkp1)
          !write(*,*) "xnorm: ",xnorm
          if (xnorm <= xtol) then
             ans = xkp1
             error = xnorm
             niters = k
             exit
          end if
          xk = xkp1
       end do
    else
       write(*,*) "k: ", k
       call prettyprint(n,1,xk)
       do i = 1, itermax
          xkp1 = C
          write(*,*) "After copying C to xkp1"
          call prettyprint(n,1,xkp1)          
          call dgemm_tiled(n,1,n,tile_size,one,one,Q,xk,xkp1)
          write(*,*) "xkp1: "
          call prettyprint(n,1,xkp1)
          xnorm = zero
          do j = 1, n
             xnorm = xnorm + abs(xkp1(j,1)-xk(j,1))
          end do
          write(*,*) "k: ", i
          k = i
          call prettyprint(n,1,xkp1)
          write(*,*) "xnorm: ",xnorm
          write(*,*) "Press any key to continue"
          read(*,*)
          if (xnorm <= xtol) then
             write(*,*) "Convergence reached on iteration #",i
             ans = xkp1
             error = xnorm
             niters = k
             exit
          end if
          xk = xkp1
          write(*,*) "After copying"
          call prettyprint(n,1,xk)
       end do
    end if
    ans = xkp1
    error = xnorm
    !niters = itermax
  end subroutine jacobi_method_tiled


end module mod_solver
