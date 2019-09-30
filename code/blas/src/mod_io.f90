module mod_io
  use mod_param, only: wp,&
       maxchar,&
       inpunit,&
       outpunit
  private


  public :: write_system_to_file
  public :: read_n_from_file
  public :: read_system_from_file


contains

  subroutine write_system_to_file(filename,n,A,x,b)
    implicit none

    character(maxchar), intent(in) :: filename
    integer, intent(in) :: n
    real(wp), intent(in) :: A(n,n)
    real(wp), intent(in) :: x(n,1)
    real(wp), intent(in) :: b(n,1)

    integer :: i

    open(unit=outpunit,file=filename)

    write(outpunit,*) "n: "
    write(outpunit,*) n
    write(outpunit,*) "A: "
    do i = 1, n
       write(outpunit,*) A(i,:)
    end do
    write(outpunit,*) "x: "
    do i = 1, n
       write(outpunit,*) x(i,1)
    end do
    write(outpunit,*) "b: "
    do i = 1, n
       write(outpunit,*) b(i,1)
    end do
    write(outpunit,*) ! Blank line just to be safe.

    close(outpunit)
    

  end subroutine write_system_to_file

  subroutine read_n_from_file(filename,n)
    implicit none

    character(maxchar), intent(in) :: filename
    integer, intent(out) :: n

    character(maxchar) :: junk

    open(unit=inpunit,file=filename)

    read(inpunit,*) junk ! "n:"
    read(inpunit,*) n

    close(inpunit)

  end subroutine read_n_from_file


  subroutine read_system_from_file(filename,n,A,x,b)
    implicit none

    character(maxchar), intent(in) :: filename
    integer, intent(in) :: n
    real(wp), intent(inout) :: A(n,n)
    real(wp), intent(inout) :: x(n,1)
    real(wp), intent(inout) :: b(n,1)

    character(maxchar) :: junk
    integer :: i
    integer :: m

    open(unit=inpunit,file=filename)

    read(inpunit,*) junk ! n:
    read(inpunit,*) m ! This is actually m, but we read as a check
    if (m .ne. n) then
       write(*,*) "*** Inconsient n in reading ",filename, " ***"
       stop 
    end if
    read(inpunit,*) junk ! A:
    do i = 1, n
       read(inpunit,*) A(i,:)
    end do
    read(inpunit,*) junk ! x:
    do i = 1, n
       read(inpunit,*) x(i,1)
    end do
    read(inpunit,*) junk ! b:
    do i = 1, n
       read(inpunit,*) b(i,1)
    end do

    close(inpunit)
    
    


  end subroutine read_system_from_file


end module mod_io
