module scale_da_read_pawr_toshiba
  use iso_c_binding

  integer, parameter :: RDIM = 600 ! max number of Range bins
  integer, parameter :: AZDIM = 320 ! max number of AZ angles
  integer, parameter :: ELDIM = 121 ! max number of AZ angles

  type, bind(c) :: c_pawr_header
     integer(kind=c_int) s_yr, s_mn, s_dy, s_hr, s_mi, s_sc
     integer(kind=c_int) e_yr, e_mn, e_dy, e_hr, e_mi, e_sc
     integer(kind=c_int) data_size
     integer(kind=c_int) total_step_num, el_num, total_el_num
     integer(kind=c_int) hit_num, sector_num, range_num, range_res, mesh_size
     real(kind=c_double) latitude, longitude, altitude
     real(kind=c_float)  start_angle, end_angle, mesh_lsb, mesh_offset
     real(kind=c_float)  tx_freq, tx_power, pulse_len_l, pulse_len_s
     real(kind=c_float)  ant_gain, beam_wid_h, beam_wid_v
     real(kind=c_float)  tx_loss, rx_loss, smin_h, smin_l
     real(kind=c_float)  prf_l, prf_h, zr_b, zr_beta
  end type c_pawr_header

  interface
     integer(kind=c_int) function read_toshiba_c(jitdt_place, hd, az, el, rtdat) bind(C, name="read_toshiba")
       use iso_c_binding
       import c_pawr_header
       import RDIM, AZDIM, ELDIM

       character(kind=c_char) :: jitdt_place(*)
       type(c_pawr_header) :: hd
       real(kind=c_float) :: az(AZDIM, ELDIM)
       real(kind=c_float) :: el(AZDIM, ELDIM)
       real(kind=c_float) :: rtdat(RDIM, AZDIM, ELDIM)
     end function read_toshiba_c
  end interface

  public

contains

  function DA_read_pawr_toshiba(fname, hd, az, el, rtdat)
    integer :: DA_read_pawr_toshiba
    character(*), intent(in) :: fname
    type(c_pawr_header), intent(out) :: hd
    real(kind=c_float), intent(out) :: az(AZDIM, ELDIM)
    real(kind=c_float), intent(out) :: el(AZDIM, ELDIM)
    real(kind=c_float), intent(out) :: rtdat(RDIM, AZDIM, ELDIM)
    character(kind=c_char) :: c_fname*1025

    !write(*, *) "jitdt_read_toshiba_f#jitdt_read_toshiba"
    c_fname = trim(fname) // c_null_char

    DA_read_pawr_toshiba = read_toshiba_c(c_fname, hd, az, el, rtdat)
  end function DA_read_pawr_toshiba

end module scale_da_read_pawr_toshiba
