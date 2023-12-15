module scale_da_read_mp_pawr_toshiba
  use iso_c_binding

  integer, parameter :: RDIM  = 800 ! max number of Range bins
  integer, parameter :: AZDIM = 320 ! max number of AZ angles
  integer, parameter :: ELDIM = 120 ! max number of AZ angles

  type, bind(c) :: c_mppawr_header
     character(c_char) :: data_name(32), site_name(32), sq_name(16)
     integer(c_int)    :: s_yr, s_mn, s_dy, s_hr, s_mi, s_sc
     integer(c_int)    :: e_yr, e_mn, e_dy, e_hr, e_mi, e_sc
     integer(c_int)    :: el_num, ray_num, range_num
     integer(c_int)    :: range_res
     real(c_double)    :: latitude, longitude, altitude
     real(c_float)     :: start_az, start_el, end_az, end_el
     real(c_float)     :: mesh_offset
  end type c_mppawr_header

  interface
     integer(kind=c_int) function read_toshiba_mpr_c(in_file, opt_verbose, hd, az, el, rtdat) bind(C, name="read_toshiba_mpr")
       use iso_c_binding
       import c_mppawr_header
       import RDIM, AZDIM, ELDIM

       character(c_char) :: in_file(*)
       integer(c_int), value :: opt_verbose
       type(c_mppawr_header) :: hd
       real(c_float) :: az(AZDIM, ELDIM)
       real(c_float) :: el(AZDIM, ELDIM)
       real(c_float) :: rtdat(RDIM, AZDIM, ELDIM)
     end function read_toshiba_mpr_c
  end interface

  public

contains

  function DA_read_mp_pawr_toshiba(fname, verbose, hd, az, el, rtdat)
    integer :: DA_read_mp_pawr_toshiba
    character(*), intent(in) :: fname
    integer(c_int), value :: verbose
    type(c_mppawr_header), intent(out) :: hd
    real(c_float), intent(out) :: az(AZDIM, ELDIM)
    real(c_float), intent(out) :: el(AZDIM, ELDIM)
    real(c_float), intent(out) :: rtdat(RDIM, AZDIM, ELDIM)
    character(c_char) :: c_fname*1025

    !write(*, *) "jitdt_read_toshiba_f#jitdt_read_toshiba"
    c_fname = trim(fname) // c_null_char

    DA_read_mp_pawr_toshiba = read_toshiba_mpr_c(c_fname, verbose, hd, az, el, rtdat)
  end function DA_read_mp_pawr_toshiba

end module scale_da_read_mp_pawr_toshiba

