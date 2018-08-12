!-------------------------------------------------------------------------------
!> Module administration
!!
!! @par Description
!!         This module is for the management of process and region
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_adm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADM_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADM_setup
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_level, &
       PRC_RGN_local, &
       PRC_RGN_vlink
    implicit none

    integer :: glevel = -1 !> grid division level
    integer :: vlayer =  1 !> number of inner vertical layer

    namelist / ADMPARAM / &
        glevel, &
        vlayer

    integer :: nmax, dmd
    integer :: l, rgnid
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[adm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*,*) 'xxx Not found namelist ADMPARAM! STOP.'
       call PRC_abort
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=ADMPARAM)

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*,*) 'xxx [ADM_setup] glevel is not appropriate :', glevel
       call PRC_abort
    endif

    ADM_glevel  = glevel
    ADM_vlayer  = vlayer

    ADM_lall    = PRC_RGN_local

    nmax        = 2**(ADM_glevel-PRC_RGN_level)
    ADM_gall_1d = 1 + nmax + 1
    ADM_gmin    = 1 + 1
    ADM_gmax    = 1 + nmax

    ADM_gall    = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in = (   nmax+1 ) * (   nmax+1 )

    ADM_imin    = 1 + 1
    ADM_imax    = 1 + nmax
    ADM_iall    = 1 + nmax + 1
    ADM_jmin    = 1 + 1
    ADM_jmax    = 1 + nmax
    ADM_jall    = 1 + nmax + 1

    ADM_vlink   = PRC_RGN_vlink
    ADM_gall_pl = PRC_RGN_vlink + 1
    ADM_gmax_pl = PRC_RGN_vlink + 1

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif

    ! for physics grid
    KS = ADM_kmin
    KE = ADM_kmax
    KA = ADM_kall
    ! IS = 1
    IE = nmax + 1
    IA = IE
    ! JS = 1
    JE = nmax + 1
    JA = JE

    return
  end subroutine ADM_setup

end module mod_adm
