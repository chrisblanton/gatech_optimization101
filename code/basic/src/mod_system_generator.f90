module mod_system_generator
  use mod_param, only: wp,&
       zero,&
       one
  use mod_matops, only: dgemm,&
       generateDiagDomMat,&
       generateRandomMatrix
  private

  public :: generateJacobisys


contains

  subroutine generateJacobisys(n,xlow,xhigh,excludeSmallerThan,A,x,b)
    implicit none

    integer, intent(in) :: n ! Size of A is n*n, x and b is n*1
    real(wp), intent(in) :: xlow ! The low for the non-diagonal parts
    real(wp), intent(in) :: xhigh ! The high for the non-diagonal
                                   ! The digaonal elments in A are
                                   !constructesuch that the magnitude
                                   !of the diagonal  than the sum
                                   !of the absolute value of the
                                   !non-diagonal row elments
    real(wp), intent(in) :: excludeSmallerThan
    real(wp), intent(out) :: A(n,n)
    real(wp), intent(out) :: x(n,1)
    real(wp), intent(out) :: b(n,1)


    call generateDiagDomMat(n,xlow,xhigh,excludeSmallerThan,A)

    call generateRandomMatrix(n,1,xlow,xhigh,excludeSmallerThan,x)
    b = zero
    call dgemm(n,1,n,one,one,A,x,b)

  end subroutine generateJacobisys

end module mod_system_generator
