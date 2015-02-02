!-------------------------------------------------------------------------------
!> module administrator for restart
!!
!! @par Description
!!          Restart administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_admin_restart
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADMIN_restart_setup

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
  logical,               public :: RESTART_RUN              = .false.         !< is this run restart?
  logical,               public :: RESTART_OUTPUT           = .false.         !< output restart file?

  character(len=H_LONG), public :: RESTART_IN_BASENAME      = ''              !< basename of the restart file
  character(len=H_LONG), public :: RESTART_OUT_BASENAME     = ''              !< basename of the output file
  character(len=H_MID),  public :: RESTART_OUT_TITLE        = ''              !< title    of the output file
  character(len=H_MID),  public :: RESTART_OUT_DTYPE        = 'DEFAULT'       !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADMIN_restart_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_RESTART_OUTPUT,       &
       ATMOS_RESTART_IN_BASENAME,  &
       ATMOS_RESTART_OUT_BASENAME, &
       ATMOS_RESTART_OUT_TITLE,    &
       ATMOS_RESTART_OUT_DTYPE
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_RESTART_OUTPUT,       &
       ATMOS_DYN_RESTART_IN_BASENAME,  &
       ATMOS_DYN_RESTART_OUT_BASENAME, &
       ATMOS_DYN_RESTART_OUT_TITLE,    &
       ATMOS_DYN_RESTART_OUT_DTYPE
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_RESTART_OUTPUT,       &
       ATMOS_PHY_AE_RESTART_IN_BASENAME,  &
       ATMOS_PHY_AE_RESTART_OUT_BASENAME, &
       ATMOS_PHY_AE_RESTART_OUT_TITLE,    &
       ATMOS_PHY_AE_RESTART_OUT_DTYPE
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_RESTART_OUTPUT,       &
       ATMOS_PHY_CH_RESTART_IN_BASENAME,  &
       ATMOS_PHY_CH_RESTART_OUT_BASENAME, &
       ATMOS_PHY_CH_RESTART_OUT_TITLE,    &
       ATMOS_PHY_CH_RESTART_OUT_DTYPE
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_RESTART_OUTPUT,       &
       ATMOS_PHY_CP_RESTART_IN_BASENAME,  &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME, &
       ATMOS_PHY_CP_RESTART_OUT_TITLE,    &
       ATMOS_PHY_CP_RESTART_OUT_DTYPE
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_RESTART_OUTPUT,       &
       ATMOS_PHY_MP_RESTART_IN_BASENAME,  &
       ATMOS_PHY_MP_RESTART_OUT_BASENAME, &
       ATMOS_PHY_MP_RESTART_OUT_TITLE,    &
       ATMOS_PHY_MP_RESTART_OUT_DTYPE
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_RESTART_OUTPUT,       &
       ATMOS_PHY_RD_RESTART_IN_BASENAME,  &
       ATMOS_PHY_RD_RESTART_OUT_BASENAME, &
       ATMOS_PHY_RD_RESTART_OUT_TITLE,    &
       ATMOS_PHY_RD_RESTART_OUT_DTYPE
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_RESTART_OUTPUT,       &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,  &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME, &
       ATMOS_PHY_SF_RESTART_OUT_TITLE,    &
       ATMOS_PHY_SF_RESTART_OUT_DTYPE
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_RESTART_OUTPUT,       &
       ATMOS_PHY_TB_RESTART_IN_BASENAME,  &
       ATMOS_PHY_TB_RESTART_OUT_BASENAME, &
       ATMOS_PHY_TB_RESTART_OUT_TITLE,    &
       ATMOS_PHY_TB_RESTART_OUT_DTYPE
    use mod_ocean_vars, only: &
       OCEAN_RESTART_OUTPUT,       &
       OCEAN_RESTART_IN_BASENAME,  &
       OCEAN_RESTART_OUT_BASENAME, &
       OCEAN_RESTART_OUT_TITLE,    &
       OCEAN_RESTART_OUT_DTYPE
    use mod_land_vars, only: &
       LAND_RESTART_OUTPUT,       &
       LAND_RESTART_IN_BASENAME,  &
       LAND_RESTART_OUT_BASENAME, &
       LAND_RESTART_OUT_TITLE,    &
       LAND_RESTART_OUT_DTYPE
    use mod_urban_vars, only: &
       URBAN_RESTART_OUTPUT,       &
       URBAN_RESTART_IN_BASENAME,  &
       URBAN_RESTART_OUT_BASENAME, &
       URBAN_RESTART_OUT_TITLE,    &
       URBAN_RESTART_OUT_DTYPE

    implicit none

    NAMELIST / PARAM_RESTART / &
       RESTART_RUN,          &
       RESTART_OUTPUT,       &
       RESTART_IN_BASENAME,  &
       RESTART_OUT_BASENAME, &
       RESTART_OUT_TITLE,    &
       RESTART_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[RESTART] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_RESTART. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_RESTART)

    ! --- set default output switch
    ATMOS_RESTART_OUTPUT        = RESTART_OUTPUT
    ATMOS_DYN_RESTART_OUTPUT    = RESTART_OUTPUT
    ATMOS_PHY_AE_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_CH_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_CP_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_MP_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_RD_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_SF_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_TB_RESTART_OUTPUT = RESTART_OUTPUT
    OCEAN_RESTART_OUTPUT        = RESTART_OUTPUT
    LAND_RESTART_OUTPUT         = RESTART_OUTPUT
    URBAN_RESTART_OUTPUT        = RESTART_OUTPUT

    ! --- set default input filename
    if( RESTART_IN_BASENAME /= '' ) then
      ATMOS_RESTART_IN_BASENAME         = RESTART_IN_BASENAME
      ATMOS_DYN_RESTART_IN_BASENAME     = RESTART_IN_BASENAME
      ATMOS_PHY_AE_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_CH_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_CP_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_MP_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_RD_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_SF_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      ATMOS_PHY_TB_RESTART_IN_BASENAME  = RESTART_IN_BASENAME
      OCEAN_RESTART_IN_BASENAME         = RESTART_IN_BASENAME
      LAND_RESTART_IN_BASENAME          = RESTART_IN_BASENAME
      URBAN_RESTART_IN_BASENAME         = RESTART_IN_BASENAME
    end if

    ! --- set default output filename
    if( RESTART_OUT_BASENAME /= '' ) then
      ATMOS_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
      ATMOS_DYN_RESTART_OUT_BASENAME    = RESTART_OUT_BASENAME
      ATMOS_PHY_AE_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_CH_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_CP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_MP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_RD_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_SF_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      ATMOS_PHY_TB_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
      OCEAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
      LAND_RESTART_OUT_BASENAME         = RESTART_OUT_BASENAME
      URBAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
    end if

    ! --- set default output title
    if( RESTART_OUT_TITLE /= '' ) then
      ATMOS_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
      ATMOS_DYN_RESTART_OUT_TITLE       = RESTART_OUT_TITLE
      ATMOS_PHY_AE_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_CH_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_CP_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_MP_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_RD_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_SF_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      ATMOS_PHY_TB_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
      OCEAN_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
      LAND_RESTART_OUT_TITLE            = RESTART_OUT_TITLE
      URBAN_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
    end if

    ! --- set default output data type
    if( RESTART_OUT_DTYPE /= '' ) then
      ATMOS_RESTART_OUT_DTYPE           = RESTART_OUT_DTYPE
      ATMOS_DYN_RESTART_OUT_DTYPE       = RESTART_OUT_DTYPE
      ATMOS_PHY_AE_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_CH_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_CP_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_MP_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_RD_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_SF_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      ATMOS_PHY_TB_RESTART_OUT_DTYPE    = RESTART_OUT_DTYPE
      OCEAN_RESTART_OUT_DTYPE           = RESTART_OUT_DTYPE
      LAND_RESTART_OUT_DTYPE            = RESTART_OUT_DTYPE
      URBAN_RESTART_OUT_DTYPE           = RESTART_OUT_DTYPE
    end if

    return
  end subroutine ADMIN_restart_setup

end module mod_admin_restart
