!-------------------------------------------------------------------------------
!> module atmosphere / reference state
!!
!! @par Description
!!          Reference state of Atmosphere
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_refstate
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_REFSTATE_setup
  public :: ATMOS_REFSTATE_finalize
  public :: ATMOS_REFSTATE_read
  public :: ATMOS_REFSTATE_write
  public :: ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: ATMOS_REFSTATE_IN_BASENAME  = ''                   !< basename of the input  file
  logical,                private :: ATMOS_REFSTATE_IN_CHECK_COORDINATES = .true.
  character(len=H_LONG),  private :: ATMOS_REFSTATE_OUT_BASENAME = ''                   !< basename of the output file
  character(len=H_MID),   private :: ATMOS_REFSTATE_OUT_TITLE    = 'SCALE-RM RefState'  !< title    of the output file
  character(len=H_SHORT), private :: ATMOS_REFSTATE_OUT_DTYPE    = 'DEFAULT'            !< REAL4 or REAL8

  character(len=H_SHORT), private :: ATMOS_REFSTATE_TYPE         = 'UNIFORM'            !< profile type
  real(RP),               private :: ATMOS_REFSTATE_TEMP_SFC     = 300.0_RP             !< surface temperature           [K]
  real(RP),               private :: ATMOS_REFSTATE_RH           =   0.0_RP             !< surface & environment RH      [%]
  real(RP),               private :: ATMOS_REFSTATE_POTT_UNIFORM = 300.0_RP             !< uniform potential temperature [K]
  real(DP),               private :: ATMOS_REFSTATE_UPDATE_DT    =  -1.0_DP
  logical,                private :: ATMOS_REFSTATE_UPDATE_FLAG  = .false.

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
  subroutine ATMOS_REFSTATE_setup( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       CZ, FZ, REAL_CZ, REAL_FZ, REAL_PHI )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: CZ      (  KA)
    real(RP), intent(in) :: FZ      (0:KA)
    real(RP), intent(in) :: REAL_CZ (  KA,IA,JA)
    real(RP), intent(in) :: REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in) :: REAL_PHI(  KA,IA,JA)

    namelist / PARAM_ATMOS_REFSTATE / &
       ATMOS_REFSTATE_IN_BASENAME,  &
       ATMOS_REFSTATE_OUT_BASENAME, &
       ATMOS_REFSTATE_OUT_TITLE,    &
       ATMOS_REFSTATE_OUT_DTYPE,    &
       ATMOS_REFSTATE_TYPE,         &
       ATMOS_REFSTATE_TEMP_SFC,     &
       ATMOS_REFSTATE_RH,           &
       ATMOS_REFSTATE_POTT_UNIFORM, &
       ATMOS_REFSTATE_UPDATE_DT

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_REFSTATE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_REFSTATE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_REFSTATE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_REFSTATE_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_REFSTATE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_REFSTATE)

    allocate( ATMOS_REFSTATE_pres(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_temp(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_dens(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_pott(KA,IA,JA) )
    allocate( ATMOS_REFSTATE_qv  (KA,IA,JA) )
    ATMOS_REFSTATE_pres(:,:,:) = UNDEF
    ATMOS_REFSTATE_temp(:,:,:) = UNDEF
    ATMOS_REFSTATE_dens(:,:,:) = UNDEF
    ATMOS_REFSTATE_pott(:,:,:) = UNDEF
    ATMOS_REFSTATE_qv  (:,:,:) = UNDEF
    !$acc enter data create(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)

    allocate( ATMOS_REFSTATE1D_pres(KA) )
    allocate( ATMOS_REFSTATE1D_temp(KA) )
    allocate( ATMOS_REFSTATE1D_dens(KA) )
    allocate( ATMOS_REFSTATE1D_pott(KA) )
    allocate( ATMOS_REFSTATE1D_qv  (KA) )
    ATMOS_REFSTATE1D_pres(:) = UNDEF
    ATMOS_REFSTATE1D_temp(:) = UNDEF
    ATMOS_REFSTATE1D_dens(:) = UNDEF
    ATMOS_REFSTATE1D_pott(:) = UNDEF
    ATMOS_REFSTATE1D_qv  (:) = UNDEF
    !$acc enter data create(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    LOG_NEWLINE
    LOG_INFO("ATMOS_REFSTATE_setup",*) 'Reference state settings'
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       LOG_INFO_CONT(*) 'Input file of reference state  : ', trim(ATMOS_REFSTATE_IN_BASENAME)
    else
       LOG_INFO_CONT(*) 'Input file of reference state  : Nothing, generate internally'
    endif

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then
       LOG_INFO_CONT(*) 'Output file of reference state : ', trim(ATMOS_REFSTATE_OUT_BASENAME)
    else
       LOG_INFO_CONT(*) 'Output file of reference state : No output'
    endif

    ATMOS_REFSTATE_UPDATE_FLAG = .false.

    ! input or generate reference profile
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then

       call ATMOS_REFSTATE_read( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                                 REAL_PHI(:,:,:)                     ) ! [IN]

    else
       if ( ATMOS_REFSTATE_TYPE == 'ISA' ) then

          LOG_INFO_CONT(*) 'Reference type                 : ISA'
          LOG_INFO_CONT(*) 'Surface temperature      [K]   : ', ATMOS_REFSTATE_TEMP_SFC
          LOG_INFO_CONT(*) 'Surface & environment RH [%]   : ', ATMOS_REFSTATE_RH

          call ATMOS_REFSTATE_generate_isa( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                                            CZ(:), FZ(:),                       & ! [IN]
                                            REAL_CZ(:,:,:), REAL_FZ(:,:,:),     & ! [IN]
                                            REAL_PHI(:,:,:)                     ) ! [IN]

       elseif ( ATMOS_REFSTATE_TYPE == 'UNIFORM' ) then

          LOG_INFO_CONT(*) 'Reference type                 : UNIFORM POTT'
          LOG_INFO_CONT(*) 'Potential temperature          : ', ATMOS_REFSTATE_POTT_UNIFORM

          call ATMOS_REFSTATE_generate_uniform( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                                                CZ(:), FZ(:),                       & ! [IN]
                                                REAL_CZ(:,:,:), REAL_FZ(:,:,:),     & ! [IN]
                                                REAL_PHI(:,:,:)                     ) ! [IN]

       elseif ( ATMOS_REFSTATE_TYPE == 'ZERO' ) then

          LOG_INFO_CONT(*) 'Reference type                 : ZERO'

          call ATMOS_REFSTATE_generate_zero( KA, IA, JA )

       elseif ( ATMOS_REFSTATE_TYPE == 'INIT' ) then

          ATMOS_REFSTATE_UPDATE_FLAG = .true.

          LOG_INFO_CONT(*) 'Reference type                 : Generate from initial data'
          LOG_INFO_CONT(*) 'Update interval [sec]          : ', ATMOS_REFSTATE_UPDATE_DT

       else
          LOG_ERROR("ATMOS_REFSTATE_setup",*) 'ATMOS_REFSTATE_TYPE must be "ISA", "UNIFORM", "ZERO", or "INIT". Check! : ', trim(ATMOS_REFSTATE_TYPE)
          call PRC_abort
       endif

    endif

    last_updated = -1.0_RP

    return
  end subroutine ATMOS_REFSTATE_setup

  subroutine ATMOS_REFSTATE_finalize

    !$acc exit data delete(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)
    !$acc exit data delete(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    deallocate( ATMOS_REFSTATE_pres )
    deallocate( ATMOS_REFSTATE_temp )
    deallocate( ATMOS_REFSTATE_dens )
    deallocate( ATMOS_REFSTATE_pott )
    deallocate( ATMOS_REFSTATE_qv   )

    deallocate( ATMOS_REFSTATE1D_pres )
    deallocate( ATMOS_REFSTATE1D_temp )
    deallocate( ATMOS_REFSTATE1D_dens )
    deallocate( ATMOS_REFSTATE1D_pott )
    deallocate( ATMOS_REFSTATE1D_qv   )

    return
  end subroutine ATMOS_REFSTATE_finalize

  !-----------------------------------------------------------------------------
  !> Read reference state profile
  subroutine ATMOS_REFSTATE_read( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       REAL_PHI )
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_read, &
       FILE_CARTESC_close
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE
    real(RP), intent(in) :: REAL_PHI(KA,IA,JA)

    integer :: fid
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_REFSTATE_read",*) 'Input reference state profile '

    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then

       ! 1D
       call FILE_CARTESC_open( ATMOS_REFSTATE_IN_BASENAME, fid, single=.true. )
       call FILE_CARTESC_read( fid, 'PRES_ref', 'Z', ATMOS_REFSTATE1D_pres(:) )
       call FILE_CARTESC_read( fid, 'TEMP_ref', 'Z', ATMOS_REFSTATE1D_temp(:) )
       call FILE_CARTESC_read( fid, 'DENS_ref', 'Z', ATMOS_REFSTATE1D_dens(:) )
       call FILE_CARTESC_read( fid, 'POTT_ref', 'Z', ATMOS_REFSTATE1D_pott(:) )
       call FILE_CARTESC_read( fid, 'QV_ref',   'Z', ATMOS_REFSTATE1D_qv  (:) )
       call FILE_CARTESC_close( fid )

       ! 3D
       call FILE_CARTESC_open( ATMOS_REFSTATE_IN_BASENAME, fid )
       if ( ATMOS_REFSTATE_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( fid, atmos=.true. )
       end if
       call FILE_CARTESC_read( fid, 'PRES_ref3D', 'ZXY', ATMOS_REFSTATE_pres(:,:,:) )
       call FILE_CARTESC_read( fid, 'TEMP_ref3D', 'ZXY', ATMOS_REFSTATE_temp(:,:,:) )
       call FILE_CARTESC_read( fid, 'DENS_ref3D', 'ZXY', ATMOS_REFSTATE_dens(:,:,:) )
       call FILE_CARTESC_read( fid, 'POTT_ref3D', 'ZXY', ATMOS_REFSTATE_pott(:,:,:) )
       call FILE_CARTESC_read( fid, 'QV_ref3D',   'ZXY', ATMOS_REFSTATE_qv  (:,:,:) )
       call FILE_CARTESC_close( fid )

    else
       LOG_ERROR("ATMOS_REFSTATE_read",*) 'refstate file is not specified.'
       call PRC_abort
    endif

    !$acc update device(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)
    !$acc update device(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    call ATMOS_REFSTATE_fillhalo( KA, KS, KE, IA, IS, IE, JA, JS, JE, REAL_PHI(:,:,:) )

    return
  end subroutine ATMOS_REFSTATE_read

  !-----------------------------------------------------------------------------
  !> Write reference state profile
  subroutine ATMOS_REFSTATE_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write
    implicit none

    logical, save :: first = .true.
    !---------------------------------------------------------------------------

    if ( .not. first ) return
    first = .false.

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then

       !$acc update host(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)
       !$acc update host(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

       LOG_NEWLINE
       LOG_INFO("ATMOS_REFSTATE_write",*) 'Output reference state profile '

       ! 1D
       call FILE_CARTESC_write( ATMOS_REFSTATE1D_pres(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'PRES_ref', 'Reference profile of pres.', 'Pa', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE1D_temp(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'TEMP_ref', 'Reference profile of temp.', 'K', 'Z',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE1D_dens(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'DENS_ref', 'Reference profile of rho', 'kg/m3', 'Z',  ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE1D_pott(:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'POTT_ref', 'Reference profile of theta', 'K', 'Z',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE1D_qv(:),   ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'QV_ref',   'Reference profile of qv', 'kg/kg', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

       ! 3D
       call FILE_CARTESC_write( ATMOS_REFSTATE_pres(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'PRES_ref3D', 'Reference profile of pres.', 'Pa', 'ZXY',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE_temp(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'TEMP_ref3D', 'Reference profile of temp.', 'K', 'ZXY',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE_dens(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'DENS_ref3D', 'Reference profile of rho', 'kg/m3', 'ZXY',  ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE_pott(:,:,:), ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'POTT_ref3D', 'Reference profile of theta', 'K', 'ZXY',    ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_REFSTATE_qv(:,:,:),   ATMOS_REFSTATE_OUT_BASENAME, ATMOS_REFSTATE_OUT_TITLE,   & ! [IN]
                          'QV_ref3D',   'Reference profile of qv', 'kg/kg', 'ZXY',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (International Standard Atmosphere)
  subroutine ATMOS_REFSTATE_generate_isa( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       CZ, FZ, REAL_CZ, REAL_FZ, REAL_PHI )
    use scale_const, only: &
       EPSvap => CONST_EPSvap, &
       Pstd => CONST_Pstd
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_profile, only: &
       PROFILE_isa => ATMOS_PROFILE_isa
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: CZ      (  KA)
    real(RP), intent(in) :: FZ      (0:KA)
    real(RP), intent(in) :: REAL_CZ (  KA,IA,JA)
    real(RP), intent(in) :: REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in) :: REAL_PHI(  KA,IA,JA)

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

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    logical :: converged
    integer :: k
    !---------------------------------------------------------------------------

    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    pres_sfc = Pstd

    call PROFILE_isa( KA, KS, KE, & ! [IN]
                      pott_sfc, pres_sfc, & ! [IN]
                      CZ(:),              & ! [IN]
                      pott(:)             ) ! [OUT]

    qv(:)  = 0.0_RP
    qc(:)  = 0.0_RP
    qv_sfc = 0.0_RP
    qc_sfc = 0.0_RP

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),              & ! [IN]
                               pres_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                               CZ(:), FZ(:),                       & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),       & ! [WORK]
#endif
                               dens(:), temp(:), pres(:),          & ! [OUT]
                               temp_sfc,                           & ! [OUT]
                               converged                           ) ! [OUT]

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_REFSTATE_generate_isa",*) "not converged"
       call PRC_abort
    end if

    ! calc QV from RH
    call SATURATION_psat_all( temp_sfc, psat_sfc )
    call SATURATION_dens2qsat_all( KA, KS, KE, &
                                   temp(:),  dens(:), & ! [IN]
                                   qsat(:)            ) ! [OUT]

    psat_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * psat_sfc ! rh * e
    qv_sfc = EPSvap * psat_sfc / ( pres_sfc - (1.0_RP-EPSvap) * psat_sfc )
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),              & ! [IN]
                               pres_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                               CZ(:), FZ(:),                       & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),       & ! [WORK]
#endif
                               dens(:), temp(:), pres(:),          & ! [OUT]
                               temp_sfc,                           & ! [OUT]
                               converged                           )

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_REFSTATE_generate_isa",*) "not converged"
       call PRC_abort
    end if

    ATMOS_REFSTATE1D_pres(:) = pres(:)
    ATMOS_REFSTATE1D_temp(:) = temp(:)
    ATMOS_REFSTATE1D_dens(:) = dens(:)
    ATMOS_REFSTATE1D_pott(:) = pott(:)
    ATMOS_REFSTATE1D_qv  (:) = qv(:)

    !$acc update device(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    call ATMOS_REFSTATE_calc3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                CZ(:), FZ(:), REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:) )

    return
  end subroutine ATMOS_REFSTATE_generate_isa

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (Uniform Potential Temperature)
  subroutine ATMOS_REFSTATE_generate_uniform( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       CZ, FZ, REAL_CZ, REAL_FZ, REAL_PHI )
    use scale_const, only: &
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: CZ      (  KA)
    real(RP), intent(in) :: FZ      (0:KA)
    real(RP), intent(in) :: REAL_CZ (  KA,IA,JA)
    real(RP), intent(in) :: REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in) :: REAL_PHI(  KA,IA,JA)

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

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    logical :: converged
    integer :: k
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
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),              & ! [IN]
                               pres_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                               CZ(:), FZ(:),                       & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),       & ! [WORK]
#endif
                               dens(:), temp(:), pres(:),          & ! [OUT]
                               temp_sfc,                           & ! [OUT]
                               converged                           ) ! [OUT]
    if ( .not. converged ) then
       LOG_ERROR("ATMOS_REFSTATE_generate_uniform",*) "not converged"
       call PRC_abort
    end if

    ! calc QV from RH
    call SATURATION_psat_all( temp_sfc, psat_sfc )
    call SATURATION_dens2qsat_all( KA, KS, KE, &
                                   temp(:), dens(:), & ! [IN]
                                   qsat(:)           ) ! [OUT]

    psat_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * psat_sfc ! rh * e
    qv_sfc = EPSvap * psat_sfc / ( pres_sfc - (1.0_RP - EPSvap) * psat_sfc )
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),              & ! [IN]
                               pres_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                               CZ(:), FZ(:),                       & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),       & ! [WORK]
#endif
                               dens(:), temp(:), pres(:),          & ! [OUT]
                               temp_sfc,                           & ! [OUT]
                               converged                           ) ! [OUT]
    if ( .not. converged ) then
       LOG_ERROR("ATMOS_REFSTATE_generate_uniform",*) "not converged"
       call PRC_abort
    end if

    ATMOS_REFSTATE1D_pres(:) = pres(:)
    ATMOS_REFSTATE1D_temp(:) = temp(:)
    ATMOS_REFSTATE1D_dens(:) = dens(:)
    ATMOS_REFSTATE1D_pott(:) = pott(:)
    ATMOS_REFSTATE1D_qv  (:) = qv(:)

    !$acc update device(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    call ATMOS_REFSTATE_calc3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                CZ(:), FZ(:), REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:) )

    return
  end subroutine ATMOS_REFSTATE_generate_uniform

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (None reference state)
  subroutine ATMOS_REFSTATE_generate_zero( KA, IA, JA )
    implicit none
    integer, intent(in) :: KA, IA, JA

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do
    do k = 1, KA
       ATMOS_REFSTATE1D_pres(k) = 0.0_RP
       ATMOS_REFSTATE1D_temp(k) = 0.0_RP
       ATMOS_REFSTATE1D_dens(k) = 0.0_RP
       ATMOS_REFSTATE1D_pott(k) = 0.0_RP
       ATMOS_REFSTATE1D_qv  (k) = 0.0_RP
    end do

    !$acc update device(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

    !$omp parallel do collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_REFSTATE_dens(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_temp(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_pres(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_pott(k,i,j) = 0.0_RP
       ATMOS_REFSTATE_qv  (k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    !$acc update device(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)

    return
  end subroutine ATMOS_REFSTATE_generate_zero

  !-----------------------------------------------------------------------------
  !> Update reference state profile (Horizontal average)
  subroutine ATMOS_REFSTATE_update( &
       KA, KS, KE, IA, IS, IE, ISB, IEB, JA, JS, JE, JSB, JEB, &
       DENS, POTT, TEMP, PRES, QV,                          &
       CZ, FZ, FDZ, RCDZ, REAL_CZ, REAL_FZ, REAL_PHI, AREA, &
       nowsec                                               )
    use scale_statistics, only: &
       STATISTICS_horizontal_mean
    use scale_interp_vert, only: &
       INTERP_VERT_xi2z
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE, ISB, IEB
    integer,  intent(in) :: JA, JS, JE, JSB, JEB
    real(RP), intent(in) :: DENS    (KA,IA,JA)
    real(RP), intent(in) :: POTT    (KA,IA,JA)
    real(RP), intent(in) :: TEMP    (KA,IA,JA)
    real(RP), intent(in) :: PRES    (KA,IA,JA)
    real(RP), intent(in) :: QV      (KA,IA,JA)
    real(RP), intent(in) :: CZ      (  KA)
    real(RP), intent(in) :: FZ      (0:KA)
    real(RP), intent(in) :: FDZ     (  KA-1)
    real(RP), intent(in) :: RCDZ    (  KA)
    real(RP), intent(in) :: REAL_CZ (  KA,IA,JA)
    real(RP), intent(in) :: REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in) :: REAL_PHI(  KA,IA,JA)
    real(RP), intent(in) :: AREA    (     IA,JA)
    real(DP), intent(in) :: nowsec

    real(RP) :: work(KA,IA,JA)

    integer  :: k
    !---------------------------------------------------------------------------

    if ( .not. ATMOS_REFSTATE_UPDATE_FLAG ) return

    if ( last_updated < 0.0_RP .or. ( nowsec - last_updated >= ATMOS_REFSTATE_UPDATE_DT ) ) then

       !$acc data copyin(DENS, POTT, TEMP, PRES, QV, CZ, FZ, FDZ, RCDZ, REAL_CZ, REAL_FZ, REAL_PHI, AREA) create(work)

       LOG_INFO("ATMOS_REFSTATE_update",*) 'update reference state'

       call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              CZ(:), REAL_CZ(:,:,:), TEMP(:,:,:), & ! [IN]
                              work(:,:,:)                         ) ! [OUT]
       call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        work(:,:,:), area(:,:),  & ! [IN]
                                        ATMOS_REFSTATE1D_temp(:) ) ! [OUT]

       call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              CZ(:), REAL_CZ(:,:,:), PRES(:,:,:), & ! [IN]
                              work(:,:,:)                         ) ! [OUT]
       call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        work(:,:,:), area(:,:),  & ! [IN]
                                        ATMOS_REFSTATE1D_pres(:) ) ! [OUT]

       call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              CZ(:), REAL_CZ(:,:,:), DENS(:,:,:), & ! [IN]
                              work(:,:,:)                         ) ! [OUT]
       call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        work(:,:,:), area(:,:),  & ! [IN]
                                        ATMOS_REFSTATE1D_dens(:) ) ! [OUT]

       call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              CZ(:), REAL_CZ(:,:,:), POTT(:,:,:), & ! [IN]
                              work(:,:,:)                         ) ! [OUT]
       call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        work(:,:,:), area(:,:),  & ! [IN]
                                        ATMOS_REFSTATE1D_pott(:) ) ! [OUT]

       call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              CZ(:), REAL_CZ(:,:,:), QV(:,:,:), & ! [IN]
                              work(:,:,:)                       ) ! [OUT]
       call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                        work(:,:,:), area(:,:), & ! [IN]
                                        ATMOS_REFSTATE1D_qv(:)  ) ! [OUT]

       !$acc update host(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

       do k = KE-1, KS, -1 ! fill undefined value
          if( ATMOS_REFSTATE1D_dens(k) <= 0.0_RP ) ATMOS_REFSTATE1D_dens(k) = ATMOS_REFSTATE1D_dens(k+1)
          if( ATMOS_REFSTATE1D_temp(k) <= 0.0_RP ) ATMOS_REFSTATE1D_temp(k) = ATMOS_REFSTATE1D_temp(k+1)
          if( ATMOS_REFSTATE1D_pres(k) <= 0.0_RP ) ATMOS_REFSTATE1D_pres(k) = ATMOS_REFSTATE1D_pres(k+1)
          if( ATMOS_REFSTATE1D_pott(k) <= 0.0_RP ) ATMOS_REFSTATE1D_pott(k) = ATMOS_REFSTATE1D_pott(k+1)
          if( ATMOS_REFSTATE1D_qv  (k) <= 0.0_RP ) ATMOS_REFSTATE1D_qv  (k) = ATMOS_REFSTATE1D_qv  (k+1)
       enddo

       call ATMOS_REFSTATE_smoothing( KA, KS, KE, FDZ(:), RCDZ(:), ATMOS_REFSTATE1D_pott(:) )
       call ATMOS_REFSTATE_smoothing( KA, KS, KE, FDZ(:), RCDZ(:), ATMOS_REFSTATE1D_qv(:) )

       !$acc update device(ATMOS_REFSTATE1D_pres, ATMOS_REFSTATE1D_temp, ATMOS_REFSTATE1D_dens, ATMOS_REFSTATE1D_pott, ATMOS_REFSTATE1D_qv)

       call ATMOS_REFSTATE_calc3D( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                   CZ(:), FZ(:), REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:) )

       last_updated = nowsec


       LOG_NEWLINE
       LOG_INFO("ATMOS_REFSTATE_update",*) 'Generated reference state of atmosphere:'
       LOG_INFO_CONT(*) '============================================================================='
       LOG_INFO_CONT(*) '   z*-coord.:    pressure: temperature:     density:   pot.temp.: water vapor'
       do k = KS, KE
          LOG_INFO_CONT('(6F13.5)')   CZ(k),                    &
                                                    ATMOS_REFSTATE1D_pres(k), &
                                                    ATMOS_REFSTATE1D_temp(k), &
                                                    ATMOS_REFSTATE1D_dens(k), &
                                                    ATMOS_REFSTATE1D_pott(k), &
                                                    ATMOS_REFSTATE1D_qv  (k)
       enddo
       LOG_INFO_CONT(*) '============================================================================='

       ! output reference profile
       call ATMOS_REFSTATE_write

       !$acc end data

    endif

    if ( ATMOS_REFSTATE_UPDATE_DT < 0.0_RP ) ATMOS_REFSTATE_UPDATE_FLAG = .false.

    return
  end subroutine ATMOS_REFSTATE_update

  !-----------------------------------------------------------------------------
  !> apply 1D reference to 3D (terrain-following) with re-calc hydrostatic balance
  subroutine ATMOS_REFSTATE_calc3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       CZ, FZ, REAL_CZ, REAL_FZ, REAL_PHI )
    use scale_prc, only: &
       PRC_abort
    use scale_interp_vert, only: &
       INTERP_VERT_z2xi
    use scale_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho_atmos     => ATMOS_HYDROSTATIC_buildrho_atmos,     &
       HYDROSTATIC_buildrho_atmos_rev => ATMOS_HYDROSTATIC_buildrho_atmos_rev
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: CZ      (  KA)
    real(RP), intent(in) :: FZ      (0:KA)
    real(RP), intent(in) :: REAL_CZ (  KA,IA,JA)
    real(RP), intent(in) :: REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in) :: REAL_PHI(  KA,IA,JA)

    real(RP) :: dens(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: pott(KA,IA,JA)
    real(RP) :: qv  (KA,IA,JA)
    real(RP) :: qc  (KA,IA,JA)
    real(RP) :: dz  (KA,IA,JA)

    real(RP) :: pott1(IA,JA)
    real(RP) :: pott2(IA,JA)
    real(RP) :: dens1(IA,JA)
    real(RP) :: dens2(IA,JA)
    real(RP) :: temp1(IA,JA)
    real(RP) :: pres1(IA,JA)
    real(RP) :: qv1  (IA,JA)
    real(RP) :: qv2  (IA,JA)
    real(RP) :: qc1  (IA,JA)
    real(RP) :: qc2  (IA,JA)
    real(RP) :: dz1  (IA,JA)

    real(RP) :: dens_toa_1D
    real(RP) :: temp_toa_1D
    real(RP) :: pres_toa_1D
    real(RP) :: qc_1D
    real(RP) :: dz_1D

    real(RP) :: work(KA,IA,JA)

    logical :: converged
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(CZ, FZ, REAL_CZ, REAL_FZ, REAL_PHI) &
    !$acc      create(dens, temp, pres, pott, qv, qc, dz, work, pott1, pott2, dens1, dens2, temp1, pres1, qv1, qv2, qc1, qc2, dz1)

    !--- potential temperature
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       work(:,i,j) = ATMOS_REFSTATE1D_pott(:)
    enddo
    enddo
    !$acc end kernels

    call INTERP_VERT_z2xi( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                           REAL_CZ(:,:,:), CZ(:), & ! [IN]
                           work(:,:,:),           & ! [IN]
                           pott(:,:,:)            ) ! [OUT]

    !--- water vapor
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       work(:,i,j) = ATMOS_REFSTATE1D_qv(:)
    enddo
    enddo
    !$acc end kernels

    call INTERP_VERT_z2xi( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! [IN]
                           REAL_CZ(:,:,:), CZ(:), & ! [IN]
                           work(:,:,:),           & ! [IN]
                           qv(:,:,:)              ) ! [OUT]


    !--- build up density to TOA (1D)
    qc_1D = 0.0_RP
    dz_1D = FZ(KE) - CZ(KE)

    call HYDROSTATIC_buildrho_atmos( ATMOS_REFSTATE1D_pott(KE), & ! [IN]
                                     ATMOS_REFSTATE1D_qv  (KE), & ! [IN]
                                     qc_1D,                     & ! [IN]
                                     ATMOS_REFSTATE1D_dens(KE), & ! [IN]
                                     ATMOS_REFSTATE1D_pott(KE), & ! [IN]
                                     ATMOS_REFSTATE1D_qv  (KE), & ! [IN]
                                     qc_1D,                     & ! [IN]
                                     dz_1D,                     & ! [IN]
                                     KE+1,                      & ! [IN]
                                     dens_toa_1D,               & ! [OUT]
                                     temp_toa_1D,               & ! [OUT]
                                     pres_toa_1D,               & ! [OUT]
                                     converged                  ) ! [OUT]

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_REFSTATE_calc3D",*) "not converged"
       call PRC_abort
    end if

    ! build down density from TOA (3D)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       dz(KS-1,i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ! distance from surface to cell center
       do k = KS, KE-1
          dz(k,i,j) = REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) ! distance from cell center to cell center
       enddo
       dz(KE,i,j) = REAL_FZ(KE,i,j) - REAL_CZ(KE,i,j) ! distance from cell center to TOA
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       dens(KE+1,i,j) = dens_toa_1D
       temp(KE+1,i,j) = temp_toa_1D
       pres(KE+1,i,j) = pres_toa_1D
       pott(KE+1,i,j) = pott(KE,i,j)
       qv  (KE+1,i,j) = qv  (KE,i,j)
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       pott(KS-1,i,j) = pott(KS,i,j)
       qv  (KS-1,i,j) = qv  (KS,i,j)
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    qc(:,:,:) = 0.0_RP
    !$acc end kernels


    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       pott1(i,j) = pott(KE,i,j)
       qv1  (i,j) = qv  (KE,i,j)
       qc1  (i,j) = qc  (KE,i,j)
       dens2(i,j) = dens(KE+1,i,j)
       pott2(i,j) = pott(KE+1,i,j)
       qv2  (i,j) = qv  (KE+1,i,j)
       qc2  (i,j) = qc  (KE+1,i,j)
       dz1  (i,j) = dz  (KE,i,j)
    end do
    end do
    !$acc end kernels
    call HYDROSTATIC_buildrho_atmos_rev( IA, IS, IE, JA, JS, JE, &
                                         pott1(:,:), qv1(:,:), qc1(:,:),             & ! [IN]
                                         dens2(:,:), pott2(:,:), qv2(:,:), qc2(:,:), & ! [IN]
                                         dz1(:,:), KE+1,                             & ! [IN]
                                         dens1(:,:), temp1(:,:), pres1(:,:)          ) ! [OUT]
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       dens(KE,i,j) = dens1(i,j)
       temp(KE,i,j) = temp1(i,j)
       pres(KE,i,j) = pres1(i,j)
    end do
    end do
    !$acc end kernels

    call HYDROSTATIC_buildrho_atmos_rev( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                         pott(:,:,:), qv(:,:,:), qc(:,:,:), & ! [IN]
                                         dz  (:,:,:),                       & ! [IN]
                                         dens(:,:,:),                       & ! [INOUT]
                                         temp(:,:,:), pres(:,:,:)           ) ! [OUT]

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_REFSTATE_dens(k,i,j) = dens(k,i,j)
       ATMOS_REFSTATE_temp(k,i,j) = temp(k,i,j)
       ATMOS_REFSTATE_pres(k,i,j) = pres(k,i,j)
       ATMOS_REFSTATE_pott(k,i,j) = pott(k,i,j)
       ATMOS_REFSTATE_qv  (k,i,j) = qv  (k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    call ATMOS_REFSTATE_fillhalo( KA, KS, KE, IA, IS, IE, JA, JS, JE, REAL_PHI(:,:,:) )

    return
  end subroutine ATMOS_REFSTATE_calc3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_fillhalo( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       REAL_PHI )
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       P00   => CONST_PRE00, &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: REAL_PHI(KA,IA,JA)

    real(RP) :: RovCP
    integer :: k, i, j

    RovCP = Rdry / CPdry

    !$acc kernels
    !$acc loop collapse(2)
    do j = JS, JE
    do i = IS, IE
       !$acc loop seq
       do k = 1, KS-1
          ATMOS_REFSTATE_temp(k, i,j) = ATMOS_REFSTATE_temp(KS,i,j)
       end do
       !$acc loop seq
       do k = KE+1, KA
          ATMOS_REFSTATE_temp(k,i,j) = ATMOS_REFSTATE_temp(KE,i,j)
       end do

       !$acc loop seq
       do k = 1, KS-1
          ATMOS_REFSTATE_qv  (k, i,j) = ATMOS_REFSTATE_qv  (KS,i,j)
       end do
       !$acc loop seq
       do k = KE+1, KA
          ATMOS_REFSTATE_qv  (k,i,j) = ATMOS_REFSTATE_qv  (KE,i,j)
       end do

       !$acc loop seq
       do k = 1, KS-2
          ATMOS_REFSTATE_pres(k, i,j) = UNDEF
       end do
       ATMOS_REFSTATE_pres(KS-1,   i,j) = ATMOS_REFSTATE_pres(KS+1,i,j) &
                                        - ATMOS_REFSTATE_dens(KS  ,i,j) * ( REAL_PHI(KS-1,i,j) - REAL_PHI(KS+1,i,j) )
       ATMOS_REFSTATE_pres(KE+1,   i,j) = ATMOS_REFSTATE_pres(KE-1,i,j) &
                                        - ATMOS_REFSTATE_dens(KE  ,i,j) * ( REAL_PHI(KE+1,i,j) - REAL_PHI(KE-1,i,j) )
       !$acc loop seq
       do k = KE+2, KA
          ATMOS_REFSTATE_pres(k,i,j) = UNDEF
       end do

       !$acc loop seq
       do k = 1, KS-2
          ATMOS_REFSTATE_dens(k, i,j) = UNDEF
       end do
       ATMOS_REFSTATE_dens(KS-1,   i,j) = ATMOS_REFSTATE_pres(KS-1,i,j) / ( ATMOS_REFSTATE_temp(KS-1,i,j) * Rdry )
       ATMOS_REFSTATE_dens(KE+1,   i,j) = ATMOS_REFSTATE_pres(KE+1,i,j) / ( ATMOS_REFSTATE_temp(KE+1,i,j) * Rdry )
       !$acc loop seq
       do k = KE+2, KA
          ATMOS_REFSTATE_dens(k,i,j) = UNDEF
       end do

       !$acc loop seq
       do k = 1, KS-2
          ATMOS_REFSTATE_pott(k, i,j) = UNDEF
       end do
       ATMOS_REFSTATE_pott(KS-1,   i,j) = ATMOS_REFSTATE_temp(KS-1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KS-1,i,j) )**RovCP
       ATMOS_REFSTATE_pott(KE+1,   i,j) = ATMOS_REFSTATE_temp(KE+1,i,j) * ( P00 / ATMOS_REFSTATE_pres(KE+1,i,j) )**RovCP
       !$acc loop seq
       do k = KE+2, KA
          ATMOS_REFSTATE_pott(k,i,j) = UNDEF
       end do
    enddo
    enddo
    !$acc end kernels

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

    !$acc update host(ATMOS_REFSTATE_pres, ATMOS_REFSTATE_temp, ATMOS_REFSTATE_dens, ATMOS_REFSTATE_pott, ATMOS_REFSTATE_qv)

    return
  end subroutine ATMOS_REFSTATE_fillhalo

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_REFSTATE_smoothing( &
       KA, KS, KE, &
       FDZ, RCDZ, &
       phi )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: FDZ (KA-1)
    real(RP), intent(in) :: RCDZ(KA)

    real(RP), intent(inout) :: phi(KA)

    real(RP) :: dev (KA)
    real(RP) :: correction(KA)
    real(RP) :: fact(KA)

    integer, parameter :: iter_max = 100
    real(RP) :: sig0, sig1, zerosw
    logical  :: updated
    integer  :: k, iter
    !---------------------------------------------------------------------------

    dev(KS) = 0.0_RP
    dev(KE) = 0.0_RP

    correction(KS-1:KS+1) = 0.0_RP
    correction(KE-1:KE+1) = 0.0_RP

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
          if ( sig0 < -EPS .and. sig1 < -EPS ) then
             correction(k) = dev(k) &
                     / ( 2.0_RP*RCDZ(k) + ( FDZ(k-1)*RCDZ(k+1) + FDZ(k)*RCDZ(k-1) ) / ( FDZ(k) + FDZ(k-1) ) )
             updated = .true.
          else
             correction(k) = 0.0_RP
          end if
       enddo

       sig1 = dev(KS+1) * dev(KS+2)
       if ( sig1 < -EPS ) then
          correction(KS+1) = dev(KS+1) &
                     / ( 2.0_RP*RCDZ(KS+1) + (FDZ(KS)*RCDZ(KS+2)+FDZ(KS+1)*RCDZ(KS))/(FDZ(KS+1)+FDZ(KS)) )
          updated = .true.
       else
          correction(KS+1) = 0.0_RP
       end if

       sig0 = dev(KE-1) * dev(KE-2)
       if ( sig0 < -EPS ) then
          correction(KE-1) = dev(KE-1) &
                     / ( 2.0_RP*RCDZ(KE-1) + (FDZ(KE-2)*RCDZ(KE)+FDZ(KE-1)*RCDZ(KE-2))/(FDZ(KE-1)+FDZ(KE-2)) )
          updated = .true.
       else
          correction(KE-1) = 0.0_RP
       end if

       if ( .NOT. updated ) exit

       do k = KS+1, KE-1
          zerosw = 0.5_RP - sign( 0.5_RP, abs(correction(k))-EPS ) ! if correction(k) == 0 then fact(k) = 0.0
          fact(k) = correction(k) / ( correction(k) - correction(k+1) - correction(k-1) + zerosw )
       enddo

       do k = KS, KE
          phi(k) = phi(k) + ( correction(k+1) * fact(k+1)          &
                            - correction(k  ) * fact(k  ) * 2.0_RP &
                            + correction(k-1) * fact(k-1)          ) * RCDZ(k)
       enddo

       if ( iter == iter_max ) then
          LOG_INFO("ATMOS_REFSTATE_smoothing",*) "iteration not converged!", phi
       endif
    enddo

    return
  end subroutine ATMOS_REFSTATE_smoothing

end module scale_atmos_refstate
