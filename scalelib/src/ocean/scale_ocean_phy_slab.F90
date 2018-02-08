!-------------------------------------------------------------------------------
!> module OCEAN / Physics Slab model
!!
!! @par Description
!!          ocean physics module, slab model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_phy_slab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_ocean_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_SLAB_setup
  public :: OCEAN_PHY_SLAB

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
  real(RP),              private :: OCEAN_PHY_SLAB_DEPTH          = 10.0_RP !< water depth of slab ocean [m]
  real(RP),              private :: OCEAN_PHY_SLAB_HeatCapacity             !< heat capacity of slab ocean [J/K/m2]

  logical,               private :: OCEAN_PHY_SLAB_nudging        = .false. ! SST Nudging is used?
  real(DP),              private :: OCEAN_PHY_SLAB_nudging_tausec           ! Relaxation time [sec]
  logical,               private :: OCEAN_PHY_SLAB_fixedsst       = .false. ! SST is fixed?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_SLAB_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_file_limit, &
       FILE_EXTERNAL_INPUT_regist
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    real(DP)               :: OCEAN_PHY_SLAB_nudging_tau                                      = 0.0_DP  ! Relaxation time
    character(len=H_SHORT) :: OCEAN_PHY_SLAB_nudging_tau_unit                                 = "SEC"
    character(len=H_LONG)  :: OCEAN_PHY_SLAB_nudging_basename(FILE_EXTERNAL_INPUT_file_limit) = ''
    logical                :: OCEAN_PHY_SLAB_nudging_enable_periodic_year                     = .false.
    logical                :: OCEAN_PHY_SLAB_nudging_enable_periodic_month                    = .false.
    logical                :: OCEAN_PHY_SLAB_nudging_enable_periodic_day                      = .false.
    integer                :: OCEAN_PHY_SLAB_nudging_step_fixed                               = 0
    real(RP)               :: OCEAN_PHY_SLAB_nudging_offset                                   = 0.0_RP
    real(RP)               :: OCEAN_PHY_SLAB_nudging_defval                                  != UNDEF
    logical                :: OCEAN_PHY_SLAB_nudging_check_coordinates                        = .true.
    integer                :: OCEAN_PHY_SLAB_nudging_step_limit                               = 0

    NAMELIST / PARAM_OCEAN_PHY_SLAB / &
       OCEAN_PHY_SLAB_DEPTH,                         &
       OCEAN_PHY_SLAB_nudging,                       &
       OCEAN_PHY_SLAB_nudging_tau,                   &
       OCEAN_PHY_SLAB_nudging_tau_unit,              &
       OCEAN_PHY_SLAB_nudging_basename,              &
       OCEAN_PHY_SLAB_nudging_enable_periodic_year,  &
       OCEAN_PHY_SLAB_nudging_enable_periodic_month, &
       OCEAN_PHY_SLAB_nudging_enable_periodic_day,   &
       OCEAN_PHY_SLAB_nudging_step_fixed,            &
       OCEAN_PHY_SLAB_nudging_offset,                &
       OCEAN_PHY_SLAB_nudging_defval,                &
       OCEAN_PHY_SLAB_nudging_check_coordinates,     &
       OCEAN_PHY_SLAB_nudging_step_limit

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[OCEAN PHY] / Origin[SCALElib]'

    OCEAN_PHY_SLAB_nudging_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_PHY_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_PHY_SLAB)

    OCEAN_PHY_SLAB_HeatCapacity = DWATR * CL * OCEAN_PHY_SLAB_DEPTH

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Slab ocean depth [m]          : ', OCEAN_PHY_SLAB_DEPTH
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean heat capacity [J/K/m2]  : ', OCEAN_PHY_SLAB_HeatCapacity

    if ( OCEAN_PHY_SLAB_nudging ) then
       call CALENDAR_unit2sec( OCEAN_PHY_SLAB_nudging_tausec, OCEAN_PHY_SLAB_nudging_tau, OCEAN_PHY_SLAB_nudging_tau_unit )

       if( IO_L ) write(IO_FID_LOG,*) '*** Use nudging for OCEAN physics : ON'
       if( IO_L ) write(IO_FID_LOG,*) '*** Relaxation time Tau [sec]     : ', OCEAN_PHY_SLAB_nudging_tausec

       if ( OCEAN_PHY_SLAB_nudging_tausec == 0.0_RP ) then
          OCEAN_PHY_SLAB_fixedsst = .true.
          if( IO_L ) write(IO_FID_LOG,*) '*** Tau=0 means that SST is completely replaced by the external data.'
       endif

       if ( OCEAN_PHY_SLAB_nudging_basename(1) == '' ) then
          write(*,*) 'xxx OCEAN_PHY_SLAB_nudging_basename is necessary !!'
          call PRC_MPIstop
       endif
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Use nudging for OCEAN physics : OFF'
    endif

    if ( OCEAN_PHY_SLAB_nudging ) then
       call FILE_EXTERNAL_INPUT_regist( OCEAN_PHY_SLAB_nudging_basename(:),           & ! [IN]
                                        'OCEAN_TEMP',                                 & ! [IN]
                                        'OXY',                                        & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_enable_periodic_year,  & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_enable_periodic_month, & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_enable_periodic_day,   & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_step_fixed,            & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_offset,                & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_defval,                & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_check_coordinates,     & ! [IN]
                                        OCEAN_PHY_SLAB_nudging_step_limit             ) ! [IN]
    endif

    return
  end subroutine OCEAN_PHY_SLAB_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_SLAB( &
       OCEAN_TEMP_t,    &
       OCEAN_TEMP,      &
       OCEAN_SFLX_WH,   &
       OCEAN_SFLX_prec, &
       OCEAN_SFLX_evap, &
       dt               )
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    real(RP), intent(out) :: OCEAN_TEMP_t   (OKMAX,IA,JA)
    real(RP), intent(in)  :: OCEAN_TEMP     (OKMAX,IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
    real(DP), intent(in)  :: dt

    real(RP) :: OCEAN_TEMP_t_ndg(OKMAX,IA,JA)
    real(RP) :: OCEAN_TEMP_ref  (OKMAX,IA,JA)
    real(RP) :: rtau

    logical  :: error
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics step: Slab'

    if ( OCEAN_PHY_SLAB_nudging ) then

       call FILE_EXTERNAL_INPUT_update( 'OCEAN_TEMP', NOWDAYSEC, OCEAN_TEMP_ref(:,:,:), error )

       if ( error ) then
          write(*,*) 'xxx Requested data is not found!'
          call PRC_MPIstop
       endif

       ! if OCEAN_PHY_SLAB_nudging_tau < dt, Nudging acts as fixed boundary
       rtau = 1.0_RP / max(OCEAN_PHY_SLAB_nudging_tausec,dt)

       do j = JS, JE
       do i = IS, IE
       do k = OKS, OKE
          OCEAN_TEMP_t_ndg(k,i,j) = ( OCEAN_TEMP_ref(k,i,j) - OCEAN_TEMP(k,i,j) ) * rtau
       enddo
       enddo
       enddo

    else
       OCEAN_TEMP_t_ndg(:,:,:) = 0.0_RP
    endif

    do j = JS, JE
    do i = IS, IE
    do k = OKS, OKE
       if ( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
          OCEAN_TEMP_t(k,i,j) = OCEAN_TEMP_t_ndg(k,i,j)
       else
          OCEAN_TEMP_t(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo

    if ( .NOT. OCEAN_PHY_SLAB_fixedsst ) then ! heat flux from atm/ice at uppermost ocean layer
       do j = JS, JE
       do i = IS, IE
          if ( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
             OCEAN_TEMP_t(OKS,i,j) = OCEAN_TEMP_t(OKS,i,j) - OCEAN_SFLX_WH(i,j) / OCEAN_PHY_SLAB_HeatCapacity
          endif
       enddo
       enddo
    endif

    return
  end subroutine OCEAN_PHY_SLAB

end module scale_ocean_phy_slab
