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
  public :: ATMOS_REFSTATE_resume
  public :: ATMOS_REFSTATE_read
  public :: ATMOS_REFSTATE_write
  public :: ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public :: ATMOS_REFSTATE_UPDATE_FLAG = .false.

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
  private :: ATMOS_REFSTATE_generate_zero
  private :: ATMOS_REFSTATE_generate_frominit

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: ATMOS_REFSTATE_IN_BASENAME  = ''                   !< basename of the input  file
  character(len=H_LONG),  private :: ATMOS_REFSTATE_OUT_BASENAME = ''                   !< basename of the output file
  character(len=H_MID) ,  private :: ATMOS_REFSTATE_OUT_TITLE    = 'SCALE-RM RefState'  !< title    of the output file
  character(len=H_MID) ,  private :: ATMOS_REFSTATE_OUT_DTYPE    = 'DEFAULT'            !< REAL4 or REAL8

  character(len=H_SHORT), private :: ATMOS_REFSTATE_TYPE         = 'UNIFORM'            !< profile type
  real(RP),               private :: ATMOS_REFSTATE_TEMP_SFC     = 300.0_RP             !< surface temperature           [K]
  real(RP),               private :: ATMOS_REFSTATE_RH           =   0.0_RP             !< surface & environment RH      [%]
  real(RP),               private :: ATMOS_REFSTATE_POTT_UNIFORM = 300.0_RP             !< uniform potential temperature [K]
  real(DP),               private :: ATMOS_REFSTATE_UPDATE_DT    = 0.0_DP

  real(DP),               private :: last_updated

  real(RP), private, allocatable :: ATMOS_REFSTATE1D_pres(:) !< 1D refernce pressure [Pa]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_temp(:) !< 1D refernce temperature [K]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_dens(:) !< 1D refernce density [kg/m3]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_pott(:) !< 1D refernce potential temperature [K]
  real(RP), private, allocatable :: ATMOS_REFSTATE1D_qv  (:) !< 1D refernce vapor [kg/kg]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_REFSTATE_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_REFSTATE / &
       ATMOS_REFSTATE_IN_BASENAME,  &
       ATMOS_REFSTATE_OUT_BASENAME, &
       ATMOS_REFSTATE_OUT_TITLE,    &
       ATMOS_REFSTATE_OUT_DTYPE,    &
       ATMOS_REFSTATE_TYPE,         &
       ATMOS_REFSTATE_TEMP_SFC,     &
       ATMOS_REFSTATE_RH,           &
       ATMOS_REFSTATE_POTT_UNIFORM, &
       ATMOS_REFSTATE_UPDATE_FLAG,  &
       ATMOS_REFSTATE_UPDATE_DT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[REFSTATE] / Categ[ATMOS SHARE] / Origin[SCALElib]'

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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_REFSTATE)

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Input file of reference state   : ', trim(ATMOS_REFSTATE_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Input file of reference state   : Nothing, generate internally'
    endif

    ! input or generate reference profile
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       call ATMOS_REFSTATE_read
    else
       if ( ATMOS_REFSTATE_TYPE == 'ISA' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type               : ISA'
          if( IO_L ) write(IO_FID_LOG,*) '*** surface temperature      [K] : ', ATMOS_REFSTATE_TEMP_SFC
          if( IO_L ) write(IO_FID_LOG,*) '*** surface & environment RH [%] : ', ATMOS_REFSTATE_RH
          call ATMOS_REFSTATE_generate_isa
          ATMOS_REFSTATE_UPDATE_FLAG = .false.

       elseif ( ATMOS_REFSTATE_TYPE == 'UNIFORM' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type               : UNIFORM POTT'
          if( IO_L ) write(IO_FID_LOG,*) '*** potential temperature        : ', ATMOS_REFSTATE_POTT_UNIFORM
          call ATMOS_REFSTATE_generate_uniform
          ATMOS_REFSTATE_UPDATE_FLAG = .false.

       elseif ( ATMOS_REFSTATE_TYPE == 'ZERO' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type               : ZERO'
          call ATMOS_REFSTATE_generate_zero
          ATMOS_REFSTATE_UPDATE_FLAG = .false.

       elseif ( ATMOS_REFSTATE_TYPE == 'INIT' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type               : make from initial data'
          if( IO_L ) write(IO_FID_LOG,*) '*** Update state?                : ', ATMOS_REFSTATE_UPDATE_FLAG
          if( IO_L ) write(IO_FID_LOG,*) '*** Update interval [sec]        : ', ATMOS_REFSTATE_UPDATE_DT

       else
          write(*,*) 'xxx ATMOS_REFSTATE_TYPE must be "ISA" or "UNIFORM". Check! : ', trim(ATMOS_REFSTATE_TYPE)
          call PRC_MPIstop
       endif

    endif

    return
  end subroutine ATMOS_REFSTATE_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_REFSTATE_resume( &
       DENS, RHOT, QTRC )
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       CZ => GRID_CZ
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    integer :: k

    ! input or generate reference profile
    if ( ATMOS_REFSTATE_IN_BASENAME == '' ) then

       if ( ATMOS_REFSTATE_TYPE == 'INIT' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type               : make from initial data'
          if( IO_L ) write(IO_FID_LOG,*) '*** Update state?                : ', ATMOS_REFSTATE_UPDATE_FLAG
          if( IO_L ) write(IO_FID_LOG,*) '*** Update interval [sec]        : ', ATMOS_REFSTATE_UPDATE_DT
          call ATMOS_REFSTATE_generate_frominit( DENS, RHOT, QTRC ) ! (in)

       endif

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
       if( IO_L ) write(IO_FID_LOG,*) '   z*-coord.:    pressure: temperature:     density:   pot.temp.: water vapor'
       do k = KS, KE
          if( IO_L ) write(IO_FID_LOG,'(6F13.5)')   CZ(k),                    &
                                                    ATMOS_REFSTATE1D_pres(k), &
                                                    ATMOS_REFSTATE1D_temp(k), &
                                                    ATMOS_REFSTATE1D_dens(k), &
                                                    ATMOS_REFSTATE1D_pott(k), &
                                                    ATMOS_REFSTATE1D_qv  (k)
       enddo
       if( IO_L ) write(IO_FID_LOG,*) '####################################################'
    endif

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Reference state output? : ', trim(ATMOS_REFSTATE_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Reference state output? : NO'
    endif

    ! output reference profile
    call ATMOS_REFSTATE_write

    return
  end subroutine ATMOS_REFSTATE_resume

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

       ! 1D
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

       ! 3D
       call FILEIO_read( ATMOS_REFSTATE_pres(:,:,:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'PRES_ref3D', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_temp(:,:,:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'TEMP_ref3D', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_dens(:,:,:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'DENS_ref3D', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_pott(:,:,:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'POTT_ref3D', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_qv(:,:,:),                               & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'QV_ref3D',   'ZXY', step=1 ) ! [IN]

    else
       if( IO_L ) write(*,*) 'xxx refstate file is not specified.'
       call PRC_MPIstop
    endif

    call ATMOS_REFSTATE_calc3D


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

       ! 1D
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

       ! 3D
       call FILEIO_write( ATMOS_REFSTATE_pres(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'PRES_ref3D', 'Reference profile of pres.', 'Pa', 'ZXY',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE_temp(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'TEMP_ref3D', 'Reference profile of temp.', 'K', 'ZXY',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE_dens(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'DENS_ref3D', 'Reference profile of rho', 'kg/m3', 'ZXY',  ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE_pott(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'POTT_ref3D', 'Reference profile of theta', 'K', 'ZXY',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_REFSTATE_qv(:,:,:),   ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'QV_ref3D',   'Reference profile of qv', 'kg/kg', 'ZXY',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (International Standard Atmosphere)
  subroutine ATMOS_REFSTATE_generate_isa
    use scale_const, only: &
       EPSvap => CONST_EPSvap, &
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
       SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
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
    real(RP) :: psat_sfc

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
    call SATURATION_psat_all( psat_sfc, temp_sfc )
    call SATURATION_dens2qsat_all( qsat(:),  temp(:),  dens(:)  )

    psat_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * psat_sfc ! rh * e
    qv_sfc = EPSvap * psat_sfc / ( pres_sfc - (1.0_RP-EPSvap) * psat_sfc )
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
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
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
    real(RP) :: psat_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    pres_sfc = Pstd
    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    qv_sfc   = 0.0_RP
    qc_sfc   = 0.0_RP

    do k = 1, KA
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
    call SATURATION_psat_all( psat_sfc, temp_sfc )
    call SATURATION_dens2qsat_all( qsat(:),  temp(:),  pres(:)  )

    psat_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * psat_sfc ! rh * e
    qv_sfc = EPSvap * psat_sfc / ( pres_sfc - (1.0_RP - EPSvap) * psat_sfc )
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
  !> Generate reference state profile (None reference state)
  subroutine ATMOS_REFSTATE_generate_zero
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
       ATMOS_REFSTATE_dens(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_temp(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_pres(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_pott(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_qv  (k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_REFSTATE_generate_zero

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
    use scale_interpolation, only: &
       INTERP_vertical_xi2z
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_hydrometer, only: &
       I_QV
    implicit none
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: pott(KA,IA,JA)
    real(RP) :: work(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( TIME_NOWSEC - last_updated >= ATMOS_REFSTATE_UPDATE_DT ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** [REFSTATE] update reference state'

       call THERMODYN_temp_pres( temp(:,:,:),   & ! [OUT]
                                 pres(:,:,:),   & ! [OUT]
                                 DENS(:,:,:),   & ! [IN]
                                 RHOT(:,:,:),   & ! [IN]
                                 QTRC(:,:,:,:), & ! [IN]
                                 TRACER_CV(:),  & ! [IN]
                                 TRACER_CP(:),  & ! [IN]
                                 TRACER_MASS(:) ) ! [IN]


       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          pott(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       enddo
       enddo

       call INTERP_vertical_xi2z( temp(:,:,:), & ! [IN]
                                  work(:,:,:)  ) ! [OUT]

       call COMM_horizontal_mean( ATMOS_REFSTATE1D_temp(:), work(:,:,:) )

       call INTERP_vertical_xi2z( pres(:,:,:), & ! [IN]
                                  work(:,:,:)  ) ! [OUT]

       call COMM_horizontal_mean( ATMOS_REFSTATE1D_pres(:), work(:,:,:) )

       call INTERP_vertical_xi2z( DENS(:,:,:), & ! [IN]
                                  work(:,:,:)  ) ! [OUT]

       call COMM_horizontal_mean( ATMOS_REFSTATE1D_dens(:), work(:,:,:) )

       call INTERP_vertical_xi2z( pott(:,:,:), & ! [IN]
                                  work(:,:,:)  ) ! [OUT]

       call COMM_horizontal_mean( ATMOS_REFSTATE1D_pott(:), work(:,:,:) )

       if ( I_QV > 0 ) then
          call INTERP_vertical_xi2z( QTRC(:,:,:,I_QV), & ! [IN]
                                     work(:,:,:)       ) ! [OUT]

          call COMM_horizontal_mean( ATMOS_REFSTATE1D_qv(:), work(:,:,:) )
       else
          ATMOS_REFSTATE1D_qv(:) = 0.0_RP
       endif

       do k = KE-1, KS, -1 ! fill undefined value
          if( ATMOS_REFSTATE1D_dens(k) <= 0.0_RP ) ATMOS_REFSTATE1D_dens(k) = ATMOS_REFSTATE1D_dens(k+1)
          if( ATMOS_REFSTATE1D_temp(k) <= 0.0_RP ) ATMOS_REFSTATE1D_temp(k) = ATMOS_REFSTATE1D_temp(k+1)
          if( ATMOS_REFSTATE1D_pres(k) <= 0.0_RP ) ATMOS_REFSTATE1D_pres(k) = ATMOS_REFSTATE1D_pres(k+1)
          if( ATMOS_REFSTATE1D_pott(k) <= 0.0_RP ) ATMOS_REFSTATE1D_pott(k) = ATMOS_REFSTATE1D_pott(k+1)
          if( ATMOS_REFSTATE1D_qv  (k) <= 0.0_RP ) ATMOS_REFSTATE1D_qv  (k) = ATMOS_REFSTATE1D_qv  (k+1)
       enddo
       call smoothing( ATMOS_REFSTATE1D_pott(:) )
       if ( I_QV > 0 ) call smoothing( ATMOS_REFSTATE1D_qv(:) )

       call ATMOS_REFSTATE_calc3D

       last_updated = TIME_NOWSEC

    endif

    return
  end subroutine ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !> apply 1D reference to 3D (terrain-following) with re-calc hydrostatic balance
  subroutine ATMOS_REFSTATE_calc3D
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       P00   => CONST_PRE00
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid, only: &
       GRID_CZ, &
       GRID_FZ
    use scale_grid_real, only: &
       REAL_PHI, &
       REAL_CZ,  &
       REAL_FZ
    use scale_interpolation, only: &
       INTERP_vertical_z2xi
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_atmos_0D     => ATMOS_HYDROSTATIC_buildrho_atmos_0D,     &
       HYDROSTATIC_buildrho_atmos_rev_2D => ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D, &
       HYDROSTATIC_buildrho_atmos_rev_3D => ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D
    implicit none


    real(RP) :: dens(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: pott(KA,IA,JA)
    real(RP) :: qv  (KA,IA,JA)
    real(RP) :: qc  (KA,IA,JA)
    real(RP) :: dz  (KA,IA,JA)

    real(RP) :: dens_toa_1D
    real(RP) :: temp_toa_1D
    real(RP) :: pres_toa_1D
    real(RP) :: qc_1D
    real(RP) :: dz_1D

    real(RP) :: work(KA,IA,JA)
    real(RP) :: RovCP
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    !--- potential temperature
    do j = JSB, JEB
    do i = ISB, IEB
       work(:,i,j) = ATMOS_REFSTATE1D_pott(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work(:,:,:), & ! [IN]
                               pott(:,:,:)  ) ! [OUT]

    !--- water vapor
    do j = JSB, JEB
    do i = ISB, IEB
       work(:,i,j) = ATMOS_REFSTATE1D_qv(:)
    enddo
    enddo

    call INTERP_vertical_z2xi( work(:,:,:), & ! [IN]
                               qv  (:,:,:)  ) ! [OUT]



    !--- build up density to TOA (1D)
    qc_1D = 0.0_RP
    dz_1D = GRID_FZ(KE) - GRID_CZ(KE)

    call HYDROSTATIC_buildrho_atmos_0D( dens_toa_1D,               & ! [OUT]
                                        temp_toa_1D,               & ! [OUT]
                                        pres_toa_1D,               & ! [OUT]
                                        ATMOS_REFSTATE1D_pott(KE), & ! [IN]
                                        ATMOS_REFSTATE1D_qv  (KE), & ! [IN]
                                        qc_1D,                     & ! [IN]
                                        ATMOS_REFSTATE1D_dens(KE), & ! [IN]
                                        ATMOS_REFSTATE1D_pott(KE), & ! [IN]
                                        ATMOS_REFSTATE1D_qv  (KE), & ! [IN]
                                        qc_1D,                     & ! [IN]
                                        dz_1D,                     & ! [IN]
                                        KE+1                       ) ! [IN]

    ! build down density from TOA (3D)
    do j = JSB, JEB
    do i = ISB, IEB
       dz(KS,i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ! distance from surface to cell center
       do k = KS+1, KE
          dz(k,i,j) = REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) ! distance from cell center to cell center
       enddo
       dz(KE+1,i,j) = REAL_FZ(KE,i,j) - REAL_CZ(KE,i,j) ! distance from cell center to TOA
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       dens(KE+1,i,j) = dens_toa_1D
       temp(KE+1,i,j) = temp_toa_1D
       pres(KE+1,i,j) = pres_toa_1D
       pott(KE+1,i,j) = pott(KE,i,j)
       qv  (KE+1,i,j) = qv  (KE,i,j)
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       pott(KS-1,i,j) = pott(KS,i,j)
       qv  (KS-1,i,j) = qv  (KS,i,j)
    enddo
    enddo

    qc(:,:,:) = 0.0_RP

    call HYDROSTATIC_buildrho_atmos_rev_2D( dens(KE  ,:,:), & ! [OUT]
                                            temp(KE  ,:,:), & ! [OUT]
                                            pres(KE  ,:,:), & ! [OUT]
                                            pott(KE  ,:,:), & ! [IN]
                                            qv  (KE  ,:,:), & ! [IN]
                                            qc  (KE  ,:,:), & ! [IN]
                                            dens(KE+1,:,:), & ! [IN]
                                            pott(KE+1,:,:), & ! [IN]
                                            qv  (KE+1,:,:), & ! [IN]
                                            qc  (KE+1,:,:), & ! [IN]
                                            dz  (KE+1,:,:), & ! [IN]
                                            KE+1            ) ! [IN]

    call HYDROSTATIC_buildrho_atmos_rev_3D( dens(:,:,:), & ! [INOUT]
                                            temp(:,:,:), & ! [OUT]
                                            pres(:,:,:), & ! [OUT]
                                            pott(:,:,:), & ! [IN]
                                            qv  (:,:,:), & ! [IN]
                                            qc  (:,:,:), & ! [IN]
                                            dz  (:,:,:)  ) ! [IN]

!    call HYDROSTATIC_buildrho_atmos_rev_2D( dens(KS-1,:,:), & ! [OUT]
!                                            temp(KS-1,:,:), & ! [OUT]
!                                            pres(KS-1,:,:), & ! [OUT]
!                                            pott(KS-1,:,:), & ! [IN]
!                                            qv  (KS-1,:,:), & ! [IN]
!                                            qc  (KS-1,:,:), & ! [IN]
!                                            dens(KS  ,:,:), & ! [IN]
!                                            pott(KS  ,:,:), & ! [IN]
!                                            qv  (KS  ,:,:), & ! [IN]
!                                            qc  (KS  ,:,:), & ! [IN]
!                                            dz  (KS  ,:,:), & ! [IN]
!                                            KS              ) ! [IN]
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       ATMOS_REFSTATE_dens(k,i,j) = dens(k,i,j)
       ATMOS_REFSTATE_temp(k,i,j) = temp(k,i,j)
       ATMOS_REFSTATE_pres(k,i,j) = pres(k,i,j)
       ATMOS_REFSTATE_pott(k,i,j) = pott(k,i,j)
       ATMOS_REFSTATE_qv  (k,i,j) = qv  (k,i,j)
    enddo
    enddo
    enddo

    ! boundary condition
    do j = JSB, JEB
    do i = ISB, IEB

       ATMOS_REFSTATE_temp(1:KS-1,i,j) = temp(KS,i,j)
       ATMOS_REFSTATE_temp(KE+1:KA,i,j) = temp_toa_1D

       ATMOS_REFSTATE_qv  (1:KS-1,i,j) = qv  (KS,i,j)
       ATMOS_REFSTATE_qv  (KE+1:KA,i,j) = qv  (KE,i,j)

       ATMOS_REFSTATE_pres(1:KS-2,i,j) = UNDEF
       ATMOS_REFSTATE_pres(KS-1,i,j) = ATMOS_REFSTATE_pres(KS+1,i,j) &
                                     - ATMOS_REFSTATE_dens(KS  ,i,j) * ( REAL_PHI(KS-1,i,j) - REAL_PHI(KS+1,i,j) )
       ATMOS_REFSTATE_pres(KE+1,i,j) = ATMOS_REFSTATE_pres(KE-1,i,j) &
                                     - ATMOS_REFSTATE_dens(KE  ,i,j) * ( REAL_PHI(KE+1,i,j) - REAL_PHI(KE-1,i,j) )
       ATMOS_REFSTATE_pres(KE+2:KA,i,j) = UNDEF

       ATMOS_REFSTATE_dens(1:KS-2,i,j) = UNDEF
       ATMOS_REFSTATE_dens(KS-1,i,j) = ATMOS_REFSTATE_pres(KS-1,i,j) / ( ATMOS_REFSTATE_temp(KS-1,i,j) * Rdry )
       ATMOS_REFSTATE_dens(KE+1,i,j) = ATMOS_REFSTATE_pres(KE+1,i,j) / ( ATMOS_REFSTATE_temp(KE+1,i,j) * Rdry )
       ATMOS_REFSTATE_dens(KE+2:KA,i,j) = UNDEF

       ATMOS_REFSTATE_pott(1:KS-2,i,j) = UNDEF
       ATMOS_REFSTATE_pott(KS-1,i,j) = ATMOS_REFSTATE_temp(KS-1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KS-1,i,j) )**RovCP
       ATMOS_REFSTATE_pott(KE+1,i,j) = ATMOS_REFSTATE_temp(KE+1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KE+1,i,j) )**RovCP
       ATMOS_REFSTATE_pott(KE+2:KA,i,j) = UNDEF
    enddo
    enddo

    call COMM_vars8( ATMOS_REFSTATE_dens(:,:,:), 1 )
    call COMM_vars8( ATMOS_REFSTATE_temp(:,:,:), 2 )
    call COMM_vars8( ATMOS_REFSTATE_pres(:,:,:), 3 )
    call COMM_vars8( ATMOS_REFSTATE_pott(:,:,:), 4 )
    call COMM_vars8( ATMOS_REFSTATE_qv  (:,:,:), 5 )
    call COMM_wait ( ATMOS_REFSTATE_dens(:,:,:), 1, .false. )
    call COMM_wait ( ATMOS_REFSTATE_temp(:,:,:), 2, .false. )
    call COMM_wait ( ATMOS_REFSTATE_pres(:,:,:), 3, .false. )
    call COMM_wait ( ATMOS_REFSTATE_pott(:,:,:), 4, .false. )
    call COMM_wait ( ATMOS_REFSTATE_qv  (:,:,:), 5, .false. )

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
          ! if (sig0>0 .OR. sig1>0) then flux(k) = 0.0
          flux(k) = dev(k) &
                  / ( 2.0_RP*RCDZ(k) + ( FDZ(k-1)*RCDZ(k+1) + FDZ(k)*RCDZ(k-1) ) / ( FDZ(k) + FDZ(k-1) ) ) &
                  * ( sign(0.5_RP ,sig0) + sign(0.5_RP ,sig1)          ) &
                  * ( sign(0.25_RP,sig0) + sign(0.25_RP,sig1) - 0.5_RP )
          updated = updated .OR. ( sig0 < -EPS .AND. sig1 < -EPS )
       enddo

       sig1 = dev(KS+1) * dev(KS+2)
       flux(KS+1) = dev(KS+1) &
                  / ( 2.0_RP*RCDZ(KS+1) + (FDZ(KS)*RCDZ(KS+2)+FDZ(KS+1)*RCDZ(KS))/(FDZ(KS+1)+FDZ(KS)) ) &
                  * ( 0.5_RP - sign(0.5_RP ,sig1) )
       updated = updated .OR. ( sig1 < -EPS )

       sig0 = dev(KE-1) * dev(KE-2)
       flux(KE-1) = dev(KE-1) &
                  / ( 2.0_RP*RCDZ(KE-1) + (FDZ(KE-2)*RCDZ(KE)+FDZ(KE-1)*RCDZ(KE-2))/(FDZ(KE-1)+FDZ(KE-2)) ) &
                  * ( 0.5_RP - sign(0.5_RP ,sig0) )
       updated = updated .OR. ( sig0 < -EPS )

       if ( .NOT. updated ) exit

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
