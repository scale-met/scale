!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for mod_atmos_phy_sf
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-27 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [mod]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_vars_setup
  public :: ATMOS_PHY_SF_vars_fillhalo
  public :: ATMOS_PHY_SF_vars_restart_read
  public :: ATMOS_PHY_SF_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public, save ::  ATMOS_PHY_SF_sw_restart = .false.

  real(RP), public, allocatable :: ATMOS_PHY_SF_DENS_t(:,:)   ! tendency DENS [    kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMZ_t(:,:)   ! tendency MOMZ [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMX_t(:,:)   ! tendency MOMX [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMY_t(:,:)   ! tendency MOMY [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOT_t(:,:)   ! tendency RHOT [K  *kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_QTRC_t(:,:,:) ! tendency QTRC [    kg/kg/s]

  ! surface flux (upward positive)
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MW   (:,:)   ! z-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MU   (:,:)   ! x-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MV   (:,:)   ! y-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SH   (:,:)   ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LH   (:,:)   ! latent   heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QTRC (:,:,:) ! tracer mass flux [kg/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LW_up(:,:)   ! surface upward longwave  flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SW_up(:,:)   ! surface upward shortwave flux [J/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_DENS  (:,:)   ! surface atmosphere density  [kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_PRES  (:,:)   ! surface atmosphere pressure [Pa]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_TEMP  (:,:)   ! surface skin temperature    [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo(:,:,:) ! surface albedo              [0-1]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_beta  (:,:)   ! evaporation efficiency      [0-1]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0    (:,:)   ! surface roughness length    [m]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo_land(:,:,:) ! surface albedo, land only [0-1]

!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QEMIS(:,:,:) ! tracer emission   flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QDEP (:,:,:) ! tracer deposition flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_VDEP (:,:,:) ! tracer deposition velocity [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: ATMOS_PHY_SF_RESTART_OUTPUT        = .false.
  character(len=H_LONG),  private, save      :: ATMOS_PHY_SF_RESTART_IN_BASENAME   = ''
  character(len=H_LONG),  private, save      :: ATMOS_PHY_SF_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),   private, save      :: ATMOS_PHY_SF_RESTART_OUT_TITLE     = 'ATMOS_PHY_SF restart'
  character(len=H_MID),   private, save      :: ATMOS_PHY_SF_RESTART_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 5       !< number of the variables
  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'ZMFLX', &
                  'XMFLX', &
                  'YMFLX', &
                  'SHFLX', &
                  'LHFLX'  /

  data VAR_DESC / 'z-momentum flux',    &
                  'x-momentum flux',    &
                  'y-momentum flux',    &
                  'sensible heat flux', &
                  'latent   heat flux'  /

  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'J/m2/s',  &
                  'J/m2/s'   /

  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_TEMP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,  &
       ATMOS_PHY_SF_RESTART_OUTPUT,       &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME, &
       ATMOS_PHY_SF_DEFAULT_SFC_TEMP,     &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo

    integer :: v, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS_PHY_SF VARS]/Categ[ATMOS]'

    ATMOS_PHY_SF_DEFAULT_SFC_TEMP   = UNDEF
    ATMOS_PHY_SF_DEFAULT_SFC_albedo = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_VARS)

    allocate( ATMOS_PHY_SF_DENS_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMZ_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMX_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMY_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOT_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_QTRC_t(IA,JA,QA) )

    allocate( ATMOS_PHY_SF_SFLX_MW   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MU   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MV   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_SH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_LH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_QTRC (IA,JA,QA) )

    allocate( ATMOS_PHY_SF_SFLX_LW_up(IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_SW_up(IA,JA)    )

    allocate( ATMOS_PHY_SF_SFC_DENS  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_PRES  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_TEMP  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_albedo(IA,JA,2)  )
    allocate( ATMOS_PHY_SF_SFC_beta  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_Z0    (IA,JA)    )

    allocate( ATMOS_PHY_SF_SFC_albedo_land(IA,JA,2) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_SF] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do v = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',v,'|',trim(VAR_NAME(v)),'|',VAR_DESC(v),'[',VAR_UNIT(v),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_PHY_SF_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_SF) : YES'
       ATMOS_PHY_SF_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_SF) : NO'
       ATMOS_PHY_SF_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_SF_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_SF_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: iq
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_SF_SFLX_MW(:,:), 1 )
    call COMM_vars8( ATMOS_PHY_SF_SFLX_MU(:,:), 2 )
    call COMM_vars8( ATMOS_PHY_SF_SFLX_MV(:,:), 3 )
    call COMM_vars8( ATMOS_PHY_SF_SFLX_SH  (:,:), 4 )
    call COMM_vars8( ATMOS_PHY_SF_SFLX_LH  (:,:), 5 )
    call COMM_wait ( ATMOS_PHY_SF_SFLX_MW(:,:), 1 )
    call COMM_wait ( ATMOS_PHY_SF_SFLX_MU(:,:), 2 )
    call COMM_wait ( ATMOS_PHY_SF_SFLX_MV(:,:), 3 )
    call COMM_wait ( ATMOS_PHY_SF_SFLX_SH  (:,:), 4 )
    call COMM_wait ( ATMOS_PHY_SF_SFLX_LH  (:,:), 5 )

    do iq = 1, QA
       call COMM_vars8( ATMOS_PHY_SF_SFLX_QTRC(:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( ATMOS_PHY_SF_SFLX_QTRC(:,:,iq), iq )
    enddo

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    character(len=H_SHORT) :: TRC_NAME

    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_SF) ***'

    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_PHY_SF_SFLX_MW(:,:),                                & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFLX_MU(:,:),                                & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(2), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFLX_MV(:,:),                                & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(3), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFLX_SH  (:,:),                                & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(4), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFLX_LH  (:,:),                                & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(5), 'XY', step=1 ) ! [IN]

       do iq = 1, QA
          TRC_NAME = 'SFLX_'//trim(AQ_NAME(iq))

          call FILEIO_read( ATMOS_PHY_SF_SFLX_QTRC(:,:,iq),                          & ! [OUT]
                            ATMOS_PHY_SF_RESTART_IN_BASENAME, TRC_NAME, 'XY', step=1 ) ! [IN]
       enddo

       call ATMOS_PHY_SF_vars_fillhalo

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_SF is not specified.'

       ATMOS_PHY_SF_SFLX_MW   (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_MU   (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_MV   (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_SH   (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_LH   (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_QTRC (:,:,:) = CONST_UNDEF

       ATMOS_PHY_SF_SFLX_LW_up(:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFLX_SW_up(:,:)   = CONST_UNDEF

       ATMOS_PHY_SF_SFC_DENS  (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFC_PRES  (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFC_TEMP  (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
       ATMOS_PHY_SF_SFC_albedo(:,:,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo
       ATMOS_PHY_SF_SFC_beta  (:,:)   = CONST_UNDEF
       ATMOS_PHY_SF_SFC_Z0    (:,:)   = CONST_UNDEF

       ATMOS_PHY_SF_SFC_albedo_land(:,:,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read


  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    character(len=H_SHORT) :: TRC_NAME
    character(len=H_MID)   :: TRC_DESC
    character(len=H_SHORT) :: TRC_UNIT

    integer :: iq
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_SF) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_SF_SFLX_MW(:,:), basename,       ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFLX_MU(:,:), basename,       ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFLX_MV(:,:), basename,       ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFLX_SH  (:,:), basename,       ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(4), VAR_DESC(4), VAR_UNIT(4), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFLX_LH  (:,:), basename,       ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(5), VAR_DESC(5), VAR_UNIT(5), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]

       do iq = 1, QA
          TRC_NAME = 'SFLX_'//trim(AQ_NAME(iq))
          TRC_DESC = 'surface flux of '//trim(AQ_NAME(iq))
          TRC_UNIT = 'kg/m2/s'

          call FILEIO_write( ATMOS_PHY_SF_SFLX_QTRC(:,:,iq), basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                             TRC_NAME, TRC_DESC, TRC_UNIT, 'XY',       ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

end module mod_atmos_phy_sf_vars
