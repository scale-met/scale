!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  integer,                private, parameter :: QA = 1
  character(len=H_SHORT), private            :: QNAME(QA)
  character(len=H_MID),   private            :: QDESC(QA)
  character(len=H_SHORT), private            :: QUNIT(QA)

  data QNAME / 'QV' /

  data QDESC / 'Ratio of Water Vapor mass to total mass (Specific humidity)' /

  data QUNIT / 'kg/kg' /

  real(DP) :: t_npf = 21600.D0 ! duration time for new particle formation

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of other components
  subroutine USER_config
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist, &
       I_QV
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_HYDROMETEOR_regist( I_QV,               & ! (out)
                                   1, 0, 0,            & ! (in)
                                   QNAME, QDESC, QUNIT ) ! (in)

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_process, only: &
       PRC_abort
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    NAMELIST /PARAM_USER/ &
         t_npf

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    use scale_time, only: &
       TIME_NOWSEC, &
       TIME_STARTDAYSEC
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_KAJINO13_flag_npf
    implicit none
    !---------------------------------------------------------------------------
    real(DP) :: t_elaps

    t_elaps = TIME_NOWSEC - TIME_STARTDAYSEC

    if ( t_elaps > t_npf ) then ! no more new particle formation does not occur
       ATMOS_PHY_AE_KAJINO13_flag_npf = .false.
    end if

    return
  end subroutine USER_step

end module mod_user
