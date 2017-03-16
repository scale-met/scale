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
  public :: ADMIN_restart

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
  logical,               public :: RESTART_RUN          = .false.   !< is this run restart?
  logical,               public :: RESTART_OUTPUT       = .false.   !< output restart file?

  character(len=H_LONG), public :: RESTART_IN_BASENAME  = ''        !< basename of the restart file
  character(len=H_LONG), public :: RESTART_OUT_BASENAME = ''        !< basename of the output file
  character(len=H_MID),  public :: RESTART_OUT_TITLE    = ''        !< title    of the output file
  character(len=H_MID),  public :: RESTART_OUT_DTYPE    = 'DEFAULT' !< REAL4 or REAL8

  integer, parameter,    public :: RESTART_OUT_ADDITIONAL_COPIES_MAX = 10                                  !< maximum number of additional copies of the outout file
  integer,               public :: RESTART_OUT_ADDITIONAL_COPIES = 0                                       !< number of additional copies of the outout file
  character(len=H_LONG), public :: RESTART_OUT_ADDITIONAL_BASENAME(RESTART_OUT_ADDITIONAL_COPIES_MAX) = '' !< basename of the additional output file

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
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_RESTART_OUTPUT,       &
       ATMOS_PHY_CH_RESTART_IN_BASENAME,  &
       ATMOS_PHY_CH_RESTART_OUT_BASENAME, &
       ATMOS_PHY_CH_RESTART_OUT_TITLE,    &
       ATMOS_PHY_CH_RESTART_OUT_DTYPE
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_RESTART_OUTPUT,       &
       ATMOS_PHY_AE_RESTART_IN_BASENAME,  &
       ATMOS_PHY_AE_RESTART_OUT_BASENAME, &
       ATMOS_PHY_AE_RESTART_OUT_TITLE,    &
       ATMOS_PHY_AE_RESTART_OUT_DTYPE
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
       RESTART_OUT_DTYPE,    &
       RESTART_OUT_ADDITIONAL_COPIES, &
       RESTART_OUT_ADDITIONAL_BASENAME

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[RESTART] / Origin[SCALE-RM]'

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

    !--- set default output switch
    ATMOS_RESTART_OUTPUT        = RESTART_OUTPUT
    ATMOS_DYN_RESTART_OUTPUT    = RESTART_OUTPUT
    ATMOS_PHY_CP_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_MP_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_RD_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_SF_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_TB_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_CH_RESTART_OUTPUT = RESTART_OUTPUT
    ATMOS_PHY_AE_RESTART_OUTPUT = RESTART_OUTPUT
    OCEAN_RESTART_OUTPUT        = RESTART_OUTPUT
    LAND_RESTART_OUTPUT         = RESTART_OUTPUT
    URBAN_RESTART_OUTPUT        = RESTART_OUTPUT

    !--- set default input filename
    if ( RESTART_IN_BASENAME /= '' ) then
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
    endif

    !--- set default output filename
    if ( RESTART_OUT_BASENAME /= '' ) then
       ATMOS_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
       ATMOS_DYN_RESTART_OUT_BASENAME    = RESTART_OUT_BASENAME
       ATMOS_PHY_CP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_MP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_RD_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_SF_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_TB_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_CH_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_AE_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       OCEAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
       LAND_RESTART_OUT_BASENAME         = RESTART_OUT_BASENAME
       URBAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
    endif

    !--- set default output title
    if ( RESTART_OUT_TITLE /= '' ) then
       ATMOS_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
       ATMOS_DYN_RESTART_OUT_TITLE       = RESTART_OUT_TITLE
       ATMOS_PHY_CP_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_MP_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_RD_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_SF_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_TB_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_CH_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       ATMOS_PHY_AE_RESTART_OUT_TITLE    = RESTART_OUT_TITLE
       OCEAN_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
       LAND_RESTART_OUT_TITLE            = RESTART_OUT_TITLE
       URBAN_RESTART_OUT_TITLE           = RESTART_OUT_TITLE
    endif

    !--- set default output data type
    if ( RESTART_OUT_DTYPE /= '' ) then
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
    endif

    return
  end subroutine ADMIN_restart_setup

  subroutine ADMIN_restart
    use mod_ocean_vars, only: &
       OCEAN_sw_restart => OCEAN_RESTART_OUTPUT, &
       OCEAN_vars_restart_create, &
       OCEAN_vars_restart_def_var, &
       OCEAN_vars_restart_enddef, &
       OCEAN_vars_restart_write_var, &
       OCEAN_vars_restart_close, &
       OCEAN_RESTART_OUT_BASENAME
    use mod_land_vars, only: &
       LAND_sw_restart => LAND_RESTART_OUTPUT, &
       LAND_vars_restart_create, &
       LAND_vars_restart_def_var, &
       LAND_vars_restart_enddef, &
       LAND_vars_restart_write_var, &
       LAND_vars_restart_close, &
       LAND_RESTART_OUT_BASENAME
    use mod_urban_vars, only: &
       URBAN_sw_restart => URBAN_RESTART_OUTPUT, &
       URBAN_vars_restart_create, &
       URBAN_vars_restart_def_var, &
       URBAN_vars_restart_enddef, &
       URBAN_vars_restart_write_var, &
       URBAN_vars_restart_close, &
       URBAN_RESTART_OUT_BASENAME
    use mod_atmos_vars, only: &
       ATMOS_sw_restart => ATMOS_RESTART_OUTPUT, &
       ATMOS_vars_restart_create, &
       ATMOS_vars_restart_def_var, &
       ATMOS_vars_restart_enddef, &
       ATMOS_vars_restart_write_var, &
       ATMOS_vars_restart_close, &
       ATMOS_RESTART_OUT_BASENAME
    use mod_admin_time, only: &
       TIME_DOATMOS_restart,  &
       TIME_DOLAND_restart,   &
       TIME_DOURBAN_restart,  &
       TIME_DOOCEAN_restart
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_RESTART_OUT_BASENAME
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_RESTART_OUT_BASENAME
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_RESTART_OUT_BASENAME
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_RESTART_OUT_BASENAME
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_RESTART_OUT_BASENAME
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_RESTART_OUT_BASENAME
    implicit none

    character(len=H_LONG) :: RESTART_OUT_BASENAME_ORIG
    integer :: irst

    ! restart files can be different for different models

    ! cread restart netCDF file
    if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_create
    if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_create
    if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_create
    if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_create

    ! define metadata in netCDF file
    if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_def_var
    if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_def_var
    if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_def_var
    if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_def_var

    ! exit define mode
    if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_enddef
    if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_enddef
    if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_enddef
    if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_enddef

    ! write variabes to netCDF file
    if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_write_var
    if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_write_var
    if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_write_var
    if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write_var

    ! clode the restart file
    if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_close
    if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_close
    if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_close
    if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_close

    ! additional restart files

    RESTART_OUT_BASENAME_ORIG = RESTART_OUT_BASENAME

    do irst = 1, min(RESTART_OUT_ADDITIONAL_COPIES, RESTART_OUT_ADDITIONAL_COPIES_MAX)
       RESTART_OUT_BASENAME = RESTART_OUT_ADDITIONAL_BASENAME(irst)
       ATMOS_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
       ATMOS_DYN_RESTART_OUT_BASENAME    = RESTART_OUT_BASENAME
       ATMOS_PHY_CP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_MP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_RD_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_SF_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_TB_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_CH_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       ATMOS_PHY_AE_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
       OCEAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
       LAND_RESTART_OUT_BASENAME         = RESTART_OUT_BASENAME
       URBAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME

       ! cread restart netCDF file
       if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_create
       if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_create
       if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_create
       if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_create

       ! define metadata in netCDF file
       if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_def_var
       if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_def_var
       if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_def_var
       if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_def_var

       ! exit define mode
       if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_enddef
       if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_enddef
       if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_enddef
       if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_enddef

       ! write variabes to netCDF file
       if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_write_var
       if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_write_var
       if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_write_var
       if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write_var

       ! close the restart file
       if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_close
       if(  LAND_sw_restart .AND. TIME_DOLAND_restart  ) call  LAND_vars_restart_close
       if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_close
       if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_close
    enddo

    RESTART_OUT_BASENAME = RESTART_OUT_BASENAME_ORIG
    ATMOS_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
    ATMOS_DYN_RESTART_OUT_BASENAME    = RESTART_OUT_BASENAME
    ATMOS_PHY_CP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_MP_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_RD_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_SF_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_TB_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_CH_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    ATMOS_PHY_AE_RESTART_OUT_BASENAME = RESTART_OUT_BASENAME
    OCEAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME
    LAND_RESTART_OUT_BASENAME         = RESTART_OUT_BASENAME
    URBAN_RESTART_OUT_BASENAME        = RESTART_OUT_BASENAME

  end subroutine ADMIN_restart

end module mod_admin_restart
