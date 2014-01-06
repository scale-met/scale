module mod_debug

  use mod_precision
  use mod_const, only: &
       UNDEF => CONST_UNDEF
  implicit none
  private

  public :: CHECK


contains

  subroutine CHECK( line, v )
    integer,  intent(in) :: line
    real(RP), intent(in) :: v
    if ( .not. abs(v) .lt. abs(UNDEF) ) then
       write(*,*) "use uninitialized value at line ", line
       call abort
    end if
  end subroutine CHECK

end module mod_debug
