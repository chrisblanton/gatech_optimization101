module mod_random
  use mod_param, only: wp
  private

  public :: init_rand
  public :: generate_random_a_b
  

contains

  subroutine init_rand(seed)
    implicit none

    integer, intent(in) :: seed

    integer :: n
    integer, allocatable :: iseed(:)

    call random_seed(size=n)
    allocate(iseed(n))
    iseed = seed
    call random_seed(put=iseed)

  end subroutine init_rand
  

  subroutine generate_random_a_b(a,b,ans)
    implicit none

    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp), intent(out) :: ans

    real(wp) :: myrand

    call random_number(myrand)
    ans =  ((b-a)*myrand)+a
    
  end subroutine generate_random_a_b

end module mod_random
