!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Reference state
!!
!! @par Description
!!          Reference state of Atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2013-02-25 (H.Yashiro)  [mod] Separate ISA profile to scale_atmos_sub_profile
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_refstate
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
  public :: ATMOS_REFSTATE_setup
  public :: ATMOS_REFSTATE_read
  public :: ATMOS_REFSTATE_write
  public :: ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public, save :: ATMOS_REFSTATE_UPDATE_FLAG = .false.

  real(RP), public, allocatable :: ATMOS_REFSTATE_pres(:,:,:) !< refernce pressure [Pa]
  real(RP), public, allocatable :: ATMOS_REFSTATE_temp(:,:,:) !< refernce temperature [K]
  real(RP), public, allocatable :: ATMOS_REFSTATE_dens(:,:,:) !< refernce density [kg/m3]
  real(RP), public, allocatable :: ATMOS_REFSTATE_pott(:,:,:) !< refernce potential temperature [K]
  real(RP), public, allocatable :: ATMOS_REFSTATE_qv  (:,:,:) !< refernce vapor [kg/kg]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_REFSTATE_generate_isa
  private :: ATMOS_REFSTATE_generate_uniform
  private :: ATMOS_REFSTATE_generate_frominit

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_pres(:) !< 1D refernce pressure [Pa]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_temp(:) !< 1D refernce temperature [K]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_dens(:) !< 1D refernce density [kg/m3]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_pott(:) !< 1D refernce potential temperature [K]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_qv  (:) !< 1D refernce vapor [kg/kg]

  character(len=H_LONG),  private, save :: ATMOS_REFSTATE_IN_BASENAME  = ''                   !< basename of the input  file
  character(len=H_LONG),  private, save :: ATMOS_REFSTATE_OUT_BASENAME = ''                   !< basename of the output file
  character(len=H_MID) ,  private, save :: ATMOS_REFSTATE_OUT_TITLE    = 'SCALE-LES Refstate' !< title    of the output file
  character(len=H_MID) ,  private, save :: ATMOS_REFSTATE_OUT_DTYPE    = 'DEFAULT'            !< REAL4 or REAL8

  character(len=H_SHORT), private, save :: ATMOS_REFSTATE_TYPE         = 'UNIFORM'            !< profile type
  real(RP),               private, save :: ATMOS_REFSTATE_TEMP_SFC     = 300.0_RP             !< surface temperature           [K]
  real(RP),               private, save :: ATMOS_REFSTATE_RH           =   0.0_RP             !< surface & environment RH      [%]
  real(RP),               private, save :: ATMOS_REFSTATE_POTT_UNIFORM = 300.0_RP             !< uniform potential temperature [K]
  real(DP),               private, save :: ATMOS_REFSTATE_UPDATE_DT    = 0.0_DP

  real(DP),               private, save :: last_updated

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_REFSTATE_setup( &
       DENS, RHOT, QTRC )
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       CZ => GRID_CZ
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    NAMELIST / PARAM_ATMOS_REFSTATE / &
       ATMOS_REFSTATE_IN_BASENAME,  &
       ATMOS_REFSTATE_OUT_BASENAME, &
       ATMOS_REFSTATE_OUT_TITLE,    &
       ATMOS_REFSTATE_OUT_DTYPE,    &
       ATMOS_REFSTATE_TYPE,         &
       ATMOS_REFSTATE_POTT_UNIFORM, &
       ATMOS_REFSTATE_TEMP_SFC,     &
       ATMOS_REFSTATE_RH,           &
       ATMOS_REFSTATE_UPDATE_FLAG,  &
       ATMOS_REFSTATE_UPDATE_DT

    integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REFSTATE]/Categ[ATMOS]'

    allocate( ATMOS_REFSTATE_pres(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_temp(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_dens(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_pott(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_qv  (KA,IA,JA) )

    allocate( ATMOS_REFSTATE1D_pres(KA) )
    allocate( ATMOS_REFSTATE1D_temp(KA) )
    allocate( ATMOS_REFSTATE1D_dens(KA) )
    allocate( ATMOS_REFSTATE1D_pott(KA) )
    allocate( ATMOS_REFSTATE1D_qv  (KA) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_REFSTATE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_REFSTATE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_REFSTATE)


    ! input or generate reference profile
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       call ATMOS_REFSTATE_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found reference state file. Generate!'

       if ( ATMOS_REFSTATE_TYPE == 'ISA' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: ISA'
          call ATMOS_REFSTATE_generate_isa
          ATMOS_REFSTATE_UPDATE_FLAG = .false.
       elseif ( ATMOS_REFSTATE_TYPE == 'UNIFORM' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: UNIFORM POTT'
          call ATMOS_REFSTATE_generate_uniform
          ATMOS_REFSTATE_UPDATE_FLAG = .false.
       elseif ( ATMOS_REFSTATE_TYPE == 'INIT' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: make from initial data'
          call ATMOS_REFSTATE_generate_frominit( DENS, RHOT, QTRC ) ! (in)
       else
          write(*,*) 'xxx ATMOS_REFSTATE_TYPE must be "ISA" or "UNIFORM". Check!', trim(ATMOS_REFSTATE_TYPE)
          call PRC_MPIstop
       endif

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
       if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
       do k = KS, KE
          if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k),                    &
                                                    ATMOS_REFSTATE1D_pres(k), &
                                                    ATMOS_REFSTATE1D_temp(k), &
                                                    ATMOS_REFSTATE1D_dens(k), &
                                                    ATMOS_REFSTATE1D_pott(k), &
                                                    ATMOS_REFSTATE1D_qv  (k)
       enddo
       if( IO_L ) write(IO_FID_LOG,*) '####################################################'
    endif

    ! output reference profile
    call ATMOS_REFSTATE_write

    return
  end subroutine ATMOS_REFSTATE_setup

  !-----------------------------------------------------------------------------
  !> Read reference state profile
  subroutine ATMOS_REFSTATE_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input reference state profile ***'

    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_REFSTATE1D_pres(:),                           & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'PRES_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE1D_temp(:),                           & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'TEMP_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE1D_dens(:),                           & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'DENS_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE1D_pott(:),                           & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'POTT_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE1D_qv(:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'QV_ref',   'Z', step=1 ) ! [IN]

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** refstate file is not specified.'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_REFSTATE_read

  !-----------------------------------------------------------------------------
  !> Write reference state profile
  subroutine ATMOS_REFSTATE_write
    use scale_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output reference state profile ***'

       call FILEIO_write( ATMOS_REFSTATE1D_pres(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'PRES_ref', 'Reference profile of pres.', 'Pa', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE1D_temp(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'TEMP_ref', 'Reference profile of temp.', 'K', 'Z',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE1D_dens(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'DENS_ref', 'Reference profile of rho', 'kg/m3', 'Z',  ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE1D_pott(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'POTT_ref', 'Reference profile of theta', 'K', 'Z',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE1D_qv(:),   ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'QV_ref',   'Reference profile of qv', 'kg/kg', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (International Standard Atmosphere)
  subroutine ATMOS_REFSTATE_generate_isa
    use scale_const, only: &
       Pstd => CONST_Pstd
    use scale_grid_real, only: &
       REAL_CZ
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_atmos_profile, only: &
       PROFILE_isa => ATMOS_PROFILE_isa
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    implicit none

    real(RP) :: z(KA)

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc
    real(RP) :: pres_sfc
    real(RP) :: pott_sfc
    real(RP) :: qv_sfc
    real(RP) :: qc_sfc

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    pres_sfc = Pstd

    call COMM_horizontal_mean( z(:), REAL_CZ(:,:,:) )

    call PROFILE_isa( KA, KS, KE, & ! [IN]
                      pott_sfc,   & ! [IN]
                      pres_sfc,   & ! [IN]
                      z   (:),    & ! [IN]
                      pott(:)     ) ! [OUT]

    qv(:)  = 0.0_RP
    qc(:)  = 0.0_RP
    qv_sfc = 0.0_RP
    qc_sfc = 0.0_RP

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc, temp_sfc, pres_sfc )
    call SATURATION_pres2qsat_all( qsat(:),  temp(:),  pres(:)  )

    qv_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ATMOS_REFSTATE1D_pres(:) = pres(:)
    ATMOS_REFSTATE1D_temp(:) = temp(:)
    ATMOS_REFSTATE1D_dens(:) = dens(:)
    ATMOS_REFSTATE1D_pott(:) = pott(:)
    ATMOS_REFSTATE1D_qv  (:) = qv(:)

    call ATMOS_REFSTATE_calc3D

    return
  end subroutine ATMOS_REFSTATE_generate_isa

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (Uniform Potential Temperature)
  subroutine ATMOS_REFSTATE_generate_uniform
    use scale_const, only: &
       Pstd   => CONST_Pstd
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    implicit none

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc
    real(RP) :: pres_sfc
    real(RP) :: pott_sfc
    real(RP) :: qv_sfc
    real(RP) :: qc_sfc

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    pres_sfc = Pstd
    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    qv_sfc   = 0.0_RP
    qc_sfc   = 0.0_RP

    do k = KS, KE
       pott(k) = ATMOS_REFSTATE_POTT_UNIFORM
       qv  (k) = 0.0_RP
       qc  (k) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc, temp_sfc, pres_sfc )
    call SATURATION_pres2qsat_all( qsat(:),  temp(:),  pres(:)  )

    qv_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ATMOS_REFSTATE1D_pres(:) = pres(:)
    ATMOS_REFSTATE1D_temp(:) = temp(:)
    ATMOS_REFSTATE1D_dens(:) = dens(:)
    ATMOS_REFSTATE1D_pott(:) = pott(:)
    ATMOS_REFSTATE1D_qv  (:) = qv(:)

    call ATMOS_REFSTATE_calc3D

    return
  end subroutine ATMOS_REFSTATE_generate_uniform

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (Horizontal average from initial data)
  subroutine ATMOS_REFSTATE_generate_frominit( &
       DENS, RHOT, QTRC )
    use scale_time, only: &
       TIME_NOWSEC
    implicit none
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    !---------------------------------------------------------------------------

    last_updated = TIME_NOWSEC - ATMOS_REFSTATE_UPDATE_DT

    call ATMOS_REFSTATE_update( DENS, RHOT, QTRC ) ! (in)

    return
  end subroutine ATMOS_REFSTATE_generate_frominit

  !-----------------------------------------------------------------------------
  !> Update reference state profile (Horizontal average)
  subroutine ATMOS_REFSTATE_update( &
       DENS, RHOT, QTRC )
    use scale_time, only: &
       TIME_NOWSEC
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    implicit none
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: pott(KA,IA,JA)
    !---------------------------------------------------------------------------

    if ( TIME_NOWSEC - last_updated >= ATMOS_REFSTATE_UPDATE_DT ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** [REFSTATE] update reference state'

       call THERMODYN_temp_pres( temp(:,:,:),  & ! [OUT]
                                 pres(:,:,:),  & ! [OUT]
                                 DENS(:,:,:),  & ! [IN]
                                 RHOT(:,:,:),  & ! [IN]
                                 QTRC(:,:,:,:) ) ! [IN]

       pott(:,:,:) = RHOT(:,:,:) / DENS(:,:,:)

       call COMM_horizontal_mean( ATMOS_REFSTATE1D_temp(:), temp(:,:,:)      )
       call COMM_horizontal_mean( ATMOS_REFSTATE1D_pres(:), pres(:,:,:)      )
       call COMM_horizontal_mean( ATMOS_REFSTATE1D_dens(:), DENS(:,:,:)      )
       call COMM_horizontal_mean( ATMOS_REFSTATE1D_pott(:), pott(:,:,:)      )
       call COMM_horizontal_mean( ATMOS_REFSTATE1D_qv  (:), QTRC(:,:,:,I_QV) )

       call smoothing( ATMOS_REFSTATE1D_temp(:) )
       call smoothing( ATMOS_REFSTATE1D_pres(:) )
       call smoothing( ATMOS_REFSTATE1D_dens(:) )
       call smoothing( ATMOS_REFSTATE1D_pott(:) )
       call smoothing( ATMOS_REFSTATE1D_qv  (:) )

       call ATMOS_REFSTATE_calc3D

       last_updated = TIME_NOWSEC

    endif

    return
  end subroutine ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !> apply 1D reference to 3D (terrain-following)
  subroutine ATMOS_REFSTATE_calc3D
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       RovCP => CONST_RovCP, &
       P00   => CONST_PRE00
    use scale_grid_real, only: &
       PHI => REAL_PHI
    use scale_interpolation, only: &
       INTERP_vertical_z2xi
    implicit none

    real(RP) :: work(KA,IA,JA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    !--- temperature
    do j = 1, JA
    do i = 1, IA
       work(:,i,j) = ATMOS_REFSTATE1D_temp(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work               (:,:,:), & ! [IN]
                               ATMOS_REFSTATE_temp(:,:,:)  ) ! [OUT]

    !--- pressure
    do j = 1, JA
    do i = 1, IA
       work(:,i,j) = ATMOS_REFSTATE1D_pres(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work               (:,:,:), & ! [IN]
                               ATMOS_REFSTATE_pres(:,:,:)  ) ! [OUT]

    !--- density
    do j = 1, JA
    do i = 1, IA
       work(:,i,j) = ATMOS_REFSTATE1D_dens(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work               (:,:,:), & ! [IN]
                               ATMOS_REFSTATE_dens(:,:,:)  ) ! [OUT]

    !--- potential temperature
    do j = 1, JA
    do i = 1, IA
       work(:,i,j) = ATMOS_REFSTATE1D_pott(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work               (:,:,:), & ! [IN]
                               ATMOS_REFSTATE_pott(:,:,:)  ) ! [OUT]

    !--- water vapor
    do j = 1, JA
    do i = 1, IA
       work(:,i,j) = ATMOS_REFSTATE1D_qv(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work             (:,:,:), & ! [IN]
                               ATMOS_REFSTATE_qv(:,:,:)  ) ! [OUT]

    ! boundary condition
    do j = 1, JA
    do i = 1, IA
       ATMOS_REFSTATE_temp(KS-1,i,j) = ATMOS_REFSTATE_temp(KS,i,j)
       ATMOS_REFSTATE_temp(KE+1,i,j) = ATMOS_REFSTATE_temp(KE,i,j)

       ATMOS_REFSTATE_pres(KS-1,i,j) = ATMOS_REFSTATE_pres(KS+1,i,j) &
                                     - ATMOS_REFSTATE_dens(KS  ,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
       ATMOS_REFSTATE_pres(KE+1,i,j) = ATMOS_REFSTATE_pres(KE-1,i,j) &
                                     - ATMOS_REFSTATE_dens(KE  ,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )

       ATMOS_REFSTATE_dens(KS-1,i,j) = ATMOS_REFSTATE_pres(KS-1,i,j) / ( ATMOS_REFSTATE_temp(KS-1,i,j) * Rdry )
       ATMOS_REFSTATE_dens(KE+1,i,j) = ATMOS_REFSTATE_pres(KE+1,i,j) / ( ATMOS_REFSTATE_temp(KE+1,i,j) * Rdry )

       ATMOS_REFSTATE_pott(KS-1,i,j) = ATMOS_REFSTATE_temp(KS-1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KS-1,i,j) )**RovCP
       ATMOS_REFSTATE_pott(KE+1,i,j) = ATMOS_REFSTATE_temp(KE+1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KE+1,i,j) )**RovCP

       ATMOS_REFSTATE_qv  (KS-1,i,j) = ATMOS_REFSTATE_qv  (KS,i,j)
       ATMOS_REFSTATE_qv  (KE+1,i,j) = ATMOS_REFSTATE_qv  (KE,i,j)
    enddo
    enddo

    return
  end subroutine ATMOS_REFSTATE_calc3D

  !-----------------------------------------------------------------------------
  subroutine smoothing( &
       phi )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_grid, only: &
       FDZ  => GRID_FDZ,  &
       RCDZ => GRID_RCDZ
    implicit none

    real(RP), intent(inout) :: phi(KA)

    real(RP) :: dev (KA)
    real(RP) :: flux(KA)
    real(RP) :: fact(KA)

    integer, parameter :: iter_max = 100
    real(RP) :: sig0, sig1, zerosw
    logical  :: updated
    integer  :: k, iter
    !---------------------------------------------------------------------------

    dev(KS) = 0.0_RP
    dev(KE) = 0.0_RP

    flux(KS-1:KS+1) = 0.0_RP
    flux(KE-1:KE+1) = 0.0_RP

    fact(KS-1:KS+1) = 0.0_RP
    fact(KE-1:KE+1) = 0.0_RP

    do iter = 1, iter_max
       updated = .false.

       do k = KS+1, KE-1
          dev(k) = phi(k) - ( FDZ(k-1)*phi(k+1) + FDZ(k)*phi(k-1) ) / ( FDZ(k) + FDZ(k-1) )
       enddo

       do k = KS+2, KE-2
          sig0 = dev(k) * dev(k-1)
          sig1 = dev(k) * dev(k+1)
          ! if (sig0>0 .or. sig1>0) then flux(k) = 0.0
          flux(k) = dev(k) &
                  / ( 2.0_RP*RCDZ(k) + ( FDZ(k-1)*RCDZ(k+1) + FDZ(k)*RCDZ(k-1) ) / ( FDZ(k) + FDZ(k-1) ) ) &
                  * ( sign(0.5_RP ,sig0) + sign(0.5_RP ,sig1)          ) &
                  * ( sign(0.25_RP,sig0) + sign(0.25_RP,sig1) - 0.5_RP )
          updated = updated .or. ( sig0 < -EPS .and. sig1 < -EPS )
       enddo

       sig1 = dev(KS+1) * dev(KS+2)
       flux(KS+1) = dev(KS+1) &
                  / ( 2.0_RP*RCDZ(KS+1) + (FDZ(KS)*RCDZ(KS+2)+FDZ(KS+1)*RCDZ(KS))/(FDZ(KS+1)+FDZ(KS)) ) &
                  * ( 0.5_RP - sign(0.5_RP ,sig1) )
       updated = updated .or. ( sig1 < -EPS )

       sig0 = dev(KE-1) * dev(KE-2)
       flux(KE-1) = dev(KE-1) &
                  / ( 2.0_RP*RCDZ(KE-1) + (FDZ(KE-2)*RCDZ(KE)+FDZ(KE-1)*RCDZ(KE-2))/(FDZ(KE-1)+FDZ(KE-2)) ) &
                  * ( 0.5_RP - sign(0.5_RP ,sig0) )
       updated = updated .or. ( sig0 < -EPS )

       if ( .not. updated ) exit

       do k = KS+1, KE-1
          zerosw = 0.5_RP - sign( 0.5_RP, abs(flux(k))-EPS ) ! if flux(k) == 0 then fact(k) = 0.0
          fact(k) = flux(k) / ( flux(k) - flux(k+1) - flux(k-1) + zerosw )
       enddo

       do k = KS, KE
          phi(k) = phi(k) + ( flux(k+1) * fact(k+1)          &
                            - flux(k  ) * fact(k  ) * 2.0_RP &
                            + flux(k-1) * fact(k-1)          ) * RCDZ(k)
       enddo

       if ( iter == iter_max ) then
          if (IO_L) write(IO_FID_LOG,*) "*** [refstate smoothing] iteration not converged!", phi
!          call PRC_MPIstop
       endif
    enddo

    return
  end subroutine smoothing

end module scale_atmos_refstate
