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
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  use scale_atmos_hydrometeor, only: &
     I_QV
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update
  public :: USER_resume

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: USER_resume_do = .false. !< do user resume step?
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: USER_do = .false. !< do user step?

  character(len=H_LONG), private :: LAND_RESTART_IN_BASENAME = ''
  logical, private :: READ_LAND_TEMP = .false.
  logical, private :: READ_LAND_WATER = .false.
  logical, private :: READ_LAND_SFC_TEMP = .false.
  logical, private :: READ_LAND_SFC_ALB_IR_dir = .false.
  logical, private :: READ_LAND_SFC_ALB_IR_dif = .false.
  logical, private :: READ_LAND_SFC_ALB_NIR_dir = .false.
  logical, private :: READ_LAND_SFC_ALB_NIR_dif = .false.
  logical, private :: READ_LAND_SFC_ALB_VIS_dir = .false.
  logical, private :: READ_LAND_SFC_ALB_VIS_dif = .false.
  logical, private :: READ_LAND_SFLX_MW = .false.
  logical, private :: READ_LAND_SFLX_MU = .false.
  logical, private :: READ_LAND_SFLX_MV = .false.
  logical, private :: READ_LAND_SFLX_SH = .false.
  logical, private :: READ_LAND_SFLX_LH = .false.
  logical, private :: READ_LAND_SFLX_GH = .false.
  logical, private :: READ_LAND_SFLX_evap = .false.
  character(len=H_LONG), private :: OCEAN_RESTART_IN_BASENAME = ''
  logical, private :: READ_OCEAN_TEMP = .false.
  logical, private :: READ_OCEAN_SFC_TEMP = .false.
  logical, private :: READ_OCEAN_SFC_ALB_IR_dir = .false.
  logical, private :: READ_OCEAN_SFC_ALB_IR_dif = .false.
  logical, private :: READ_OCEAN_SFC_ALB_NIR_dir = .false.
  logical, private :: READ_OCEAN_SFC_ALB_NIR_dif = .false.
  logical, private :: READ_OCEAN_SFC_ALB_VIS_dir = .false.
  logical, private :: READ_OCEAN_SFC_ALB_VIS_dif = .false.
  logical, private :: READ_OCEAN_SFC_Z0M = .false.
  logical, private :: READ_OCEAN_SFC_Z0H = .false.
  logical, private :: READ_OCEAN_SFC_Z0E = .false.
  logical, private :: READ_OCEAN_SFLX_MW = .false.
  logical, private :: READ_OCEAN_SFLX_MU = .false.
  logical, private :: READ_OCEAN_SFLX_MV = .false.
  logical, private :: READ_OCEAN_SFLX_SH = .false.
  logical, private :: READ_OCEAN_SFLX_LH = .false.
  logical, private :: READ_OCEAN_SFLX_WH = .false.
  logical, private :: READ_OCEAN_SFLX_evap = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of tracers
  subroutine USER_tracer_setup
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    ! if you want to add tracers, call the TRACER_regist subroutine.
    ! e.g.,
!    integer, parameter     :: NQ = 1
!    integer                :: QS
!    character(len=H_SHORT) :: NAME(NQ)
!    character(len=H_MID)   :: DESC(NQ)
!    character(len=H_SHORT) :: UNIT(NQ)
!
!    data NAME (/ 'name' /)
!    data DESC (/ 'tracer name' /)
!    data UNIT (/ 'kg/kg' /)
    !---------------------------------------------------------------------------

!    call TRACER_regist( QS,   & ! [OUT]
!                        NQ,   & ! [IN]
!                        NAME, & ! [IN]
!                        DESC, & ! [IN]
!                        UNIT  ) ! [IN]

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    use mod_land_admin, only: &
       LAND_do
    use mod_ocean_admin, only: &
       OCEAN_do
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       USER_resume_do, &
       LAND_RESTART_IN_BASENAME, &
       READ_LAND_TEMP, &
       READ_LAND_WATER, &
       READ_LAND_SFC_TEMP, &
       READ_LAND_SFC_ALB_IR_dir, &
       READ_LAND_SFC_ALB_IR_dif, &
       READ_LAND_SFC_ALB_NIR_dir, &
       READ_LAND_SFC_ALB_NIR_dif, &
       READ_LAND_SFC_ALB_VIS_dir, &
       READ_LAND_SFC_ALB_VIS_dif, &
       READ_LAND_SFLX_MW, &
       READ_LAND_SFLX_MU, &
       READ_LAND_SFLX_MV, &
       READ_LAND_SFLX_SH, &
       READ_LAND_SFLX_LH, &
       READ_LAND_SFLX_GH, &
       READ_LAND_SFLX_evap, &
       OCEAN_RESTART_IN_BASENAME, &
       READ_OCEAN_TEMP, &
       READ_OCEAN_SFC_TEMP, &
       READ_OCEAN_SFC_ALB_IR_dir, &
       READ_OCEAN_SFC_ALB_IR_dif, &
       READ_OCEAN_SFC_ALB_NIR_dir, &
       READ_OCEAN_SFC_ALB_NIR_dif, &
       READ_OCEAN_SFC_ALB_VIS_dir, &
       READ_OCEAN_SFC_ALB_VIS_dif, &
       READ_OCEAN_SFC_Z0M, &
       READ_OCEAN_SFC_Z0H, &
       READ_OCEAN_SFC_Z0E, &
       READ_OCEAN_SFLX_MW, &
       READ_OCEAN_SFLX_MU, &
       READ_OCEAN_SFLX_MV, &
       READ_OCEAN_SFLX_SH, &
       READ_OCEAN_SFLX_LH, &
       READ_OCEAN_SFLX_WH, &
       READ_OCEAN_SFLX_evap

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    call IO_filename_replace( LAND_RESTART_IN_BASENAME, 'LAND_RESTART_IN_BASENAME' )
    call IO_filename_replace( OCEAN_RESTART_IN_BASENAME, 'OCEAN_RESTART_IN_BASENAME' )

    if ( LAND_do .and. LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(LAND_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified in mod_user.'
       if( IO_L ) write(IO_FID_LOG,*) '*** no variables are overwritten.'
       READ_LAND_TEMP = .false.
       READ_LAND_WATER = .false.
       READ_LAND_SFC_TEMP = .false.
       READ_LAND_SFC_ALB_IR_dir = .false.
       READ_LAND_SFC_ALB_IR_dif = .false.
       READ_LAND_SFC_ALB_NIR_dir = .false.
       READ_LAND_SFC_ALB_NIR_dif = .false.
       READ_LAND_SFC_ALB_VIS_dir = .false.
       READ_LAND_SFC_ALB_VIS_dif = .false.
       READ_LAND_SFLX_MW = .false.
       READ_LAND_SFLX_MU = .false.
       READ_LAND_SFLX_MV = .false.
       READ_LAND_SFLX_SH = .false.
       READ_LAND_SFLX_LH = .false.
       READ_LAND_SFLX_GH = .false.
       READ_LAND_SFLX_evap = .false.
    end if

    if ( OCEAN_do .and. OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(OCEAN_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified in mod_user.'
       if( IO_L ) write(IO_FID_LOG,*) '*** no variables are overwritten.'
       READ_OCEAN_TEMP = .false.
       READ_OCEAN_SFC_TEMP = .false.
       READ_OCEAN_SFC_ALB_IR_dir = .false.
       READ_OCEAN_SFC_ALB_IR_dif = .false.
       READ_OCEAN_SFC_ALB_NIR_dir = .false.
       READ_OCEAN_SFC_ALB_NIR_dif = .false.
       READ_OCEAN_SFC_ALB_VIS_dir = .false.
       READ_OCEAN_SFC_ALB_VIS_dif = .false.
       READ_OCEAN_SFC_Z0M = .false.
       READ_OCEAN_SFC_Z0H = .false.
       READ_OCEAN_SFC_Z0E = .false.
       READ_OCEAN_SFLX_MW = .false.
       READ_OCEAN_SFLX_MU = .false.
       READ_OCEAN_SFLX_MV = .false.
       READ_OCEAN_SFLX_SH = .false.
       READ_OCEAN_SFLX_LH = .false.
       READ_OCEAN_SFLX_WH = .false.
       READ_OCEAN_SFLX_evap = .false.
    end if

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Finish reading namelist'

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    use scale_const, only: &
         I_SW  => CONST_I_SW, &
         I_LW  => CONST_I_LW
    use scale_file_cartesC, only: &
         FILE_CARTESC_read
    use mod_land_vars, only: &
         LAND_TEMP, &
         LAND_WATER, &
         LAND_SFC_TEMP, &
         LAND_SFC_albedo, &
         LAND_SFLX_MW, &
         LAND_SFLX_MU, &
         LAND_SFLX_MV, &
         LAND_SFLX_SH, &
         LAND_SFLX_LH, &
         LAND_SFLX_GH, &
         LAND_SFLX_QTRC
    use mod_ocean_vars, only: &
         OCEAN_TEMP, &
         OCEAN_SFC_TEMP, &
         OCEAN_SFC_albedo, &
         OCEAN_SFC_Z0M, &
         OCEAN_SFC_Z0H, &
         OCEAN_SFC_Z0E, &
         OCEAN_SFLX_MW, &
         OCEAN_SFLX_MU, &
         OCEAN_SFLX_MV, &
         OCEAN_SFLX_SH, &
         OCEAN_SFLX_LH, &
         OCEAN_SFLX_GH, &
         OCEAN_SFLX_QTRC
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_resume
    use scale_file_cartesC, only: &
         FILE_CARTESC_read
    use mod_land_vars, only: &
         LAND_TEMP, &
         LAND_WATER, &
         LAND_SFC_TEMP, &
         LAND_SFC_albedo, &
         LAND_SFLX_MW, &
         LAND_SFLX_MU, &
         LAND_SFLX_MV, &
         LAND_SFLX_SH, &
         LAND_SFLX_LH, &
         LAND_SFLX_GH, &
         LAND_SFLX_QTRC
    use mod_ocean_vars, only: &
         OCEAN_TEMP, &
         OCEAN_SFC_TEMP, &
         OCEAN_SFC_albedo, &
         OCEAN_SFC_Z0M, &
         OCEAN_SFC_Z0H, &
         OCEAN_SFC_Z0E, &
         OCEAN_SFLX_MW, &
         OCEAN_SFLX_MU, &
         OCEAN_SFLX_MV, &
         OCEAN_SFLX_SH, &
         OCEAN_SFLX_LH, &
         OCEAN_SFLX_GH, &
         OCEAN_SFLX_QTRC
    implicit none
    !---------------------------------------------------------------------------

!    if( IO_L ) write(IO_FID_LOG,*)
!    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (LAND) in mod_user ***'
    LOG_NEWLINE
    LOG_INFO("USER_resume",*) 'Start user_resume'

    if ( READ_LAND_TEMP ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_TEMP', 'LXY', & ! [IN]
                               LAND_TEMP(:,:,:), step=1                      ) ! [OUT]
    end if

    if ( READ_LAND_WATER ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_WATER', 'LXY', & ! [IN]
                               LAND_WATER(:,:,:), step=1                      ) ! [OUT]
    end if

    if ( READ_LAND_SFC_TEMP ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_TEMP', 'XY', & ! [IN]
                               LAND_SFC_TEMP(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_IR_dir ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_IR_dir', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct, I_R_IR), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_IR_dif ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_IR_dif', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse, I_R_IR), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_NIR_dir ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_NIR_dir', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct, I_R_NIR), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_NIR_dif ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_NIR_dif', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse, I_R_NIR), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_VIS_dir ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_VIS_dir', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_direct, I_R_VIS), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFC_ALB_VIS_dif ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFC_ALB_VIS_dif', 'XY', & ! [IN]
                               LAND_SFC_albedo(:,:,I_R_diffuse, I_R_VIS), step=1        ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_MW ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_MW', 'XY', & ! [IN]
                               LAND_SFLX_MW(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_MU ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_MU', 'XY', & ! [IN]
                               LAND_SFLX_MU(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_MV ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_MV', 'XY', & ! [IN]
                               LAND_SFLX_MV(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_SH ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_SH', 'XY', & ! [IN]
                               LAND_SFLX_SH(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_LH ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_LH', 'XY', & ! [IN]
                               LAND_SFLX_LH(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_GH ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_GH', 'XY', & ! [IN]
                               LAND_SFLX_GH(:,:), step=1                       ) ! [OUT]
    end if

    if ( READ_LAND_SFLX_evap ) then
       call FILE_CARTESC_read( LAND_RESTART_IN_BASENAME, 'LAND_SFLX_evap', 'XY', & ! [IN]
                               LAND_SFLX_QTRC(:,:,I_QV), step=1                  ) ! [OUT]
    endif

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (OCEAN) in mod_user ***'

    if ( READ_OCEAN_TEMP ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_TEMP', 'OXY', & ! [IN]
                               OCEAN_TEMP(:,:,:), step=1                       ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_TEMP ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_TEMP', 'XY', & ! [IN]
                               OCEAN_SFC_TEMP(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_IR_dir ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_IR_dir', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct, I_R_IR), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_IR_dif ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_IR_dif', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse, I_R_IR), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_NIR_dir ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_NIR_dir', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct, I_R_NIR), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_NIR_dif ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_NIR_dif', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse, I_R_NIR), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_VIS_dir ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_VIS_dir', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_direct, I_R_VIS), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_ALB_VIS_dif ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_ALB_VIS_dif', 'XY', & ! [IN]
                               OCEAN_SFC_albedo(:,:,I_R_diffuse, I_R_VIS), step=1         ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_Z0M ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_Z0M', 'XY', & ! [IN]
                               OCEAN_SFC_Z0M(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFC_Z0H ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_Z0H', 'XY', & ! [IN]
                               OCEAN_SFC_Z0H(:,:), step=1                        ) ! [OUT]
    endif


    if ( READ_OCEAN_SFC_Z0E ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFC_Z0E', 'XY', & ! [IN]
                               OCEAN_SFC_Z0E(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_MW ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_MW', 'XY', & ! [IN]
                               OCEAN_SFLX_MW(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_MU ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_MU', 'XY', & ! [IN]
                               OCEAN_SFLX_MU(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_MV ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_MV', 'XY', & ! [IN]
                               OCEAN_SFLX_MV(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_SH ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_SH', 'XY', & ! [IN]
                               OCEAN_SFLX_SH(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_LH ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_LH', 'XY', & ! [IN]
                               OCEAN_SFLX_LH(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_WH ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_WH', 'XY', & ! [IN]
                               OCEAN_SFLX_GH(:,:), step=1                        ) ! [OUT]
    endif

    if ( READ_OCEAN_SFLX_evap ) then
       call FILE_CARTESC_read( OCEAN_RESTART_IN_BASENAME, 'OCEAN_SFLX_evap', 'XY', & ! [IN]
                               OCEAN_SFLX_QTRC(:,:,I_QV), step=1                   ) ! [OUT]
    endif

    LOG_NEWLINE
    LOG_INFO("USER_resume",*) 'Finish user_resume'

    return
  end subroutine USER_resume

end module mod_user
