!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Input and Output of the SDM variables
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-06-27 (S.Shima) [new] sdm_outasci is added
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_io

  implicit none
  private
  public :: sdm_outasci

contains
  subroutine sdm_outasci(otime,sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sd_dmpnskip)
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    integer, intent(in) :: sd_dmpnskip ! Base skip to store super droplets in text format

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp, ftmp2
    character(len=H_LONG) :: basename_sd_out
    character(len=H_LONG) :: basename_random
    character(len=H_LONG) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate

    !--- output restart file of Super Droplet
    write(basename_time,'(F15.3)') otime
    do n = 1, 15
       if( basename_time(n:n) == ' ' ) basename_time(n:n) = '0'
    enddo
    fid_sdm_o = IO_get_available_fid()

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,fmt3) 'SD_output', '_ASCII_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    open (fid_sdm_o, file = trim(basename_sd_out), & !action = "write", &
          access = "sequential", status = "replace", form = "formatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_ascii_out", "Write error"
      call PRC_MPIstop
    endif 

    write(fmt,'( "(", i2.2, "e16.8,i20,i10)" )')(sd_numasl+5)

    do m=1,sd_num,sd_dmpnskip

       write(fid_sdm_o,trim(fmt)) sd_x(m),              &
     &                           sd_y(m),              &
     &                           sd_z(m), sd_vz(m), sd_r(m),           &
     &                           (sd_asl(m,n),n=1,sd_numasl),sd_n(m),m

    end do

    close(fid_sdm_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (ASCII) of Super Droplet'

    return

  end subroutine sdm_outasci
end module m_sdm_io
