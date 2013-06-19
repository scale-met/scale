module mod_debug

  use mod_const, only: &
       UNDEF => CONST_UNDEF
  include 'inc_precision.h'


contains

  subroutine CHECK( line, v )
    integer,  intent(in) :: line
    real(RP), intent(in) :: v
    if ( abs(v) .ge. abs(UNDEF) ) then
       write(*,*) "use uninitialized value at line ", line
       call abort
    end if
  end subroutine CHECK

end module mod_debug
