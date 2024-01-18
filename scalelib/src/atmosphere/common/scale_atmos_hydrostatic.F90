!-------------------------------------------------------------------------------
!> module atmosphere / hydrostatic barance
!!
!! @par Description
!!          make hydrostatic profile in the model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_const, only: &
     GRAV    => CONST_GRAV,    &
     Rdry    => CONST_Rdry,    &
     Rvap    => CONST_Rvap,    &
     CVdry   => CONST_CVdry,   &
     CVvap   => CONST_CVvap,   &
     CPdry   => CONST_CPdry,   &
     CPvap   => CONST_CPvap,   &
     LAPS    => CONST_LAPS,    &
     LAPSdry => CONST_LAPSdry, &
     P00     => CONST_PRE00,   &
     EPSTvap => CONST_EPSTvap
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_HYDROSTATIC_setup
  public :: ATMOS_HYDROSTATIC_buildrho
  public :: ATMOS_HYDROSTATIC_buildrho_real
  public :: ATMOS_HYDROSTATIC_buildrho_atmos
  public :: ATMOS_HYDROSTATIC_buildrho_atmos_rev
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp_atmos

  public :: ATMOS_HYDROSTATIC_barometric_law_mslp
  public :: ATMOS_HYDROSTATIC_barometric_law_pres

  interface ATMOS_HYDROSTATIC_buildrho
#ifdef _OPENACC
     module procedure ATMOS_HYDROSTATIC_buildrho_1D_cpu
#endif
     module procedure ATMOS_HYDROSTATIC_buildrho_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_3D
  end interface ATMOS_HYDROSTATIC_buildrho

  interface ATMOS_HYDROSTATIC_buildrho_real
     module procedure ATMOS_HYDROSTATIC_buildrho_real_3D
  end interface ATMOS_HYDROSTATIC_buildrho_real

  interface ATMOS_HYDROSTATIC_buildrho_atmos
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_0D
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_3D
  end interface ATMOS_HYDROSTATIC_buildrho_atmos

  interface ATMOS_HYDROSTATIC_buildrho_atmos_rev
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D
  end interface ATMOS_HYDROSTATIC_buildrho_atmos_rev

  interface ATMOS_HYDROSTATIC_buildrho_bytemp
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_3D
  end interface ATMOS_HYDROSTATIC_buildrho_bytemp

  interface ATMOS_HYDROSTATIC_buildrho_bytemp_atmos
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D
  end interface ATMOS_HYDROSTATIC_buildrho_bytemp_atmos

  interface ATMOS_HYDROSTATIC_barometric_law_mslp
     module procedure ATMOS_HYDROSTATIC_barometric_law_mslp_0D
     module procedure ATMOS_HYDROSTATIC_barometric_law_mslp_2D
  end interface ATMOS_HYDROSTATIC_barometric_law_mslp

  interface ATMOS_HYDROSTATIC_barometric_law_pres
     module procedure ATMOS_HYDROSTATIC_barometric_law_pres_0D
     module procedure ATMOS_HYDROSTATIC_barometric_law_pres_2D
  end interface ATMOS_HYDROSTATIC_barometric_law_pres

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_HYDROSTATIC_buildrho_atmos_2D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: itelim = 100 !< itelation number limit
  real(RP), private            :: criteria                            !< convergence judgement criteria
  logical,  private            :: HYDROSTATIC_uselapserate  = .false. !< use lapse rate?
  integer,  private            :: HYDROSTATIC_buildrho_real_kref = 1
  integer,  private            :: HYDROSTATIC_barometric_law_mslp_kref = 1 !< reference layer for MSLP calculation
  !$acc declare create(criteria, HYDROSTATIC_uselapserate, HYDROSTATIC_barometric_law_mslp_kref)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_HYDROSTATIC_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_EPS
    implicit none

    namelist / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_uselapserate, &
       HYDROSTATIC_buildrho_real_kref, &
       HYDROSTATIC_barometric_law_mslp_kref

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_HYDROSTATIC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_HYDROSTATIC)

    criteria = sqrt( CONST_EPS )

    LOG_NEWLINE
    LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Use lapse rate for estimation of surface temperature? : ', HYDROSTATIC_uselapserate
    LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Buildrho conversion criteria                          : ', criteria

    !$acc update device(criteria, HYDROSTATIC_uselapserate, HYDROSTATIC_barometric_law_mslp_kref)

    return
  end subroutine ATMOS_HYDROSTATIC_setup

  !-----------------------------------------------------------------------------
  !> Build up density from surface (1D)
#ifdef _OPENACC
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_1D_cpu( &
       KA, KS, KE, &
       pott, qv, qc,                       &
       pres_sfc, pott_sfc, qv_sfc, qc_sfc, &
       cz, fz,                             &
       dens, temp, pres,                   &
       temp_sfc,                           &
       converged                           )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: pott(KA)        !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA)        !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA)        !< liquid water          [kg/kg]
    real(RP), intent(in)  :: pres_sfc        !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc        !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc          !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc          !< surface liquid water          [kg/kg]
    real(RP), intent(in)  :: cz  (KA)
    real(RP), intent(in)  :: fz  (0:KA)
    real(RP), intent(out) :: dens(KA)        !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA)        !< temperature           [K]
    real(RP), intent(out) :: pres(KA)        !< pressure              [Pa]
    real(RP), intent(out) :: temp_sfc        !< surface temperature           [K]
    logical,  intent(out) :: converged

    real(RP) :: dz(KA)
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)

    call ATMOS_HYDROSTATIC_buildrho_1D( &
         KA, KS, KE, &
         pott, qv, qc,                       &
         pres_sfc, pott_sfc, qv_sfc, qc_sfc, &
         cz, fz,                             &
         dz, work1, work2,                   &
         dens, temp, pres,                   &
         temp_sfc,                           &
         converged                           )

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_1D_cpu
#endif
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_1D( &
       KA, KS, KE, &
       pott, qv, qc,                       &
       pres_sfc, pott_sfc, qv_sfc, qc_sfc, &
       cz, fz,                             &
#ifdef _OPENACC
       dz, work1, work2,                   &
#endif
       dens, temp, pres,                   &
       temp_sfc,                           &
       converged                           )
    !$acc routine seq
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: pott(KA)        !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA)        !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA)        !< liquid water          [kg/kg]
    real(RP), intent(in)  :: pres_sfc        !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc        !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc          !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc          !< surface liquid water          [kg/kg]
    real(RP), intent(in)  :: cz  (KA)
    real(RP), intent(in)  :: fz  (0:KA)
    real(RP), intent(out) :: dens(KA)        !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA)        !< temperature           [K]
    real(RP), intent(out) :: pres(KA)        !< pressure              [Pa]
    real(RP), intent(out) :: temp_sfc        !< surface temperature           [K]
    logical,  intent(out) :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: dz(KA)
    real(RP), intent(out) :: work1(KA)
    real(RP), intent(out) :: work2(KA)
#else
    real(RP) :: dz(KA)
#endif

    real(RP) :: dens_sfc
    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: CPtot_sfc
    real(RP) :: CPovCV_sfc
    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPtot
    real(RP) :: CPovCV

    integer  :: k
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CV_VAPOR * qv_sfc                    &
               + CV_WATER * qc_sfc
    CPtot_sfc  = CPdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CP_VAPOR * qv_sfc                    &
               + CP_WATER * qc_sfc
    CPovCV_sfc = CPtot_sfc / CVtot_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CV_VAPOR * qv(KS)                    &
           + CV_WATER * qc(KS)
    CPtot  = CPdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CP_VAPOR * qv(KS)                    &
           + CP_WATER * qc(KS)
    CPovCV = CPtot / CVtot

    ! density at surface
    dens_sfc   = P00 / Rtot_sfc / pott_sfc * ( pres_sfc/P00 )**(CVtot_sfc/CPtot_sfc)
    temp_sfc   = pres_sfc / ( dens_sfc * Rtot_sfc )

    dz(KS-1) = CZ(KS) - FZ(KS-1)
    do k = KS, KE-1
       dz(k) = CZ(k+1) - CZ(k)
    end do

    ! make density at lowermost cell center
    if ( HYDROSTATIC_uselapserate ) then

       temp(KS) = pott_sfc - LAPSdry * dz(KS-1) ! use dry lapse rate
       pres(KS) = P00 * ( temp(KS)/pott(KS) )**(CPtot/Rtot)
       dens(KS) = P00 / Rtot / pott(KS) * ( pres(KS)/P00 )**(CVtot/CPtot)

       converged = .true.

    else ! use itelation

       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( pott(KS), qv(KS), qc(KS),           & ! [IN]
                                                 dens_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                                                 dz(KS-1), KS-1,                     & ! [IN]
                                                 dens(KS), temp(KS), pres(KS),       & ! [OUT]
                                                 converged                           ) ! [OUT]

    endif

    if ( converged ) then
       !--- from lowermost atmosphere to top of atmosphere
       call ATMOS_HYDROSTATIC_buildrho_atmos_1D( KA, KS, KE, &
                                                 pott(:), qv(:), qc(:), & ! [IN]
                                                 dz(:),                 & ! [IN]
#ifdef _OPENACC
                                                 work1(:), work2(:),    & ! [WORK]
#endif
                                                 dens(:),               & ! [INOUT]
                                                 temp(:), pres(:),      & ! [OUT]
                                                 converged              ) ! [OUT]
       if ( converged ) then
          ! fill dummy
          dens(   1:KS-1) = 0.0_RP
          dens(KE+1:KA  ) = 0.0_RP
       end if
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_1D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pott, qv, qc,       &
       pres_sfc, pott_sfc, &
       qv_sfc, qc_sfc,     &
       cz, fz, area,       &
       dens, temp, pres,   &
       temp_sfc            )
    use scale_prc, only: &
       PRC_abort
    use scale_statistics, only: &
       STATISTICS_horizontal_mean
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: pott(KA,IA,JA)  !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA)  !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA)  !< liquid water          [kg/kg]
    real(RP), intent(in)  :: pres_sfc(IA,JA) !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc(IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc  (IA,JA) !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (IA,JA) !< surface liquid water          [kg/kg]
    real(RP), intent(in)  :: cz(  KA,IA,JA)
    real(RP), intent(in)  :: fz(0:KA,IA,JA)
    real(RP), intent(in)  :: area(IA,JA)

    real(RP), intent(out) :: dens(KA,IA,JA)  !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA,IA,JA)  !< temperature           [K]
    real(RP), intent(out) :: pres(KA,IA,JA)  !< pressure              [Pa]
    real(RP), intent(out) :: temp_sfc(IA,JA) !< surface temperature           [K]

    real(RP) :: dz(KA,IA,JA), dz_top(IA,JA)

    ! TOA
    real(RP) :: pott_toa(IA,JA)
    real(RP) :: qv_toa  (IA,JA)
    real(RP) :: qc_toa  (IA,JA)
    real(RP) :: dens_toa(IA,JA)
    real(RP) :: temp_toa(IA,JA)
    real(RP) :: pres_toa(IA,JA)
    ! k = KE
    real(RP) :: pott_ke(IA,JA)
    real(RP) :: qv_ke  (IA,JA)
    real(RP) :: qc_ke  (IA,JA)
    real(RP) :: dens_ke(IA,JA)
    real(RP) :: temp_ke(IA,JA)
    real(RP) :: pres_ke(IA,JA)

    real(RP) :: dens_mean

#ifdef _OPENACC
    real(RP) :: work1(KA,IA,JA)
    real(RP) :: work2(KA,IA,JA)
    real(RP) :: work3(KA,IA,JA)
#endif

    logical  :: converged, flag
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(pott, qv, qc, pres_sfc, pott_sfc, qv_sfc, qc_sfc, cz, fz, area) &
    !$acc      copyout(dens, temp, pres, temp_sfc) &
    !$acc      create(dz, dz_top, pott_toa, qv_toa, qc_toa, dens_toa, temp_toa, pres_toa, pott_ke, qv_ke, qc_ke, dens_ke, temp_ke, pres_ke)


    converged = .true.

    !--- from surface to lowermost atmosphere
    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_1D( KA, KS, KE, &
                                           pott(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                           pres_sfc(i,j), pott_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                           cz(:,i,j), fz(:,i,j),                                   & ! [IN]
#ifdef _OPENACC
                                           work1(:,i,j), work2(:,i,j), work3(:,i,j),               & ! [WORK]
#endif
                                           dens(:,i,j), temp(:,i,j), pres(:,i,j),                  & ! [OUT]
                                           temp_sfc(i,j),                                          & ! [OUT]
                                           flag                                                    ) ! [OUT]
       converged = converged .and. flag
    end do
    end do
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_3D",*) "not converged"
       call PRC_abort
    end if

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE-1
          dz(k,i,j) = CZ(k+1,i,j) - CZ(k,i,j)
       end do
       dz_top(i,j) = FZ(KE,i,j) - CZ(KE,i,j) ! distance from cell center to TOA

       ! value at TOA
       pott_toa(i,j) = pott(KE,i,j)
       qv_toa  (i,j) = qv  (KE,i,j)
       qc_toa  (i,j) = qc  (KE,i,j)
       ! value at k=KE
       dens_ke(i,j) = DENS(KE,i,j)
       pott_ke(i,j) = pott(KE,i,j)
       qv_ke  (i,j) = qv  (KE,i,j)
       qc_ke  (i,j) = qc  (KE,i,j)
    end do
    end do
    !$acc end kernels

    call ATMOS_HYDROSTATIC_buildrho_atmos_2D( IA, IS, IE, JA, JS, JE, &
                                              pott_toa(:,:), qv_toa(:,:), qc_toa(:,:),            & ! [IN]
                                              dens_ke(:,:), pott_ke(:,:), qv_ke(:,:), qc_ke(:,:), & ! [IN]
                                              dz_top(:,:), KE+1,                                  & ! [IN]
                                              dens_toa(:,:), temp_toa(:,:), pres_toa(:,:)         ) ! [OUT]

    call STATISTICS_horizontal_mean( IA, IS, IE, JA, JS, JE, &
                                     dens_toa(:,:), area(:,:), & ! [IN]
                                     dens_mean                 ) ! [OUT]

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       dens_toa(i,j) = dens_mean
    enddo
    enddo
    !$acc end kernels

    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D( IA, IS, IE, JA, JS, JE, &
                                                  pott_ke(:,:), qv_ke(:,:), qc_ke(:,:),                   & ! [IN]
                                                  dens_toa(:,:), pott_toa(:,:), qv_toa(:,:), qc_toa(:,:), & ! [IN]
                                                  dz_top(:,:), KE+1,                                      & ! [IN]
                                                  dens_ke(:,:), temp_ke(:,:), pres_ke(:,:)                ) ! [OUT]

    !$omp parallel do OMP_SCHEDULE_
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       DENS(KE,i,j) = dens_ke(i,j)
    end do
    end do
    !$acc end kernels

    !--- from top of atmosphere to lowermost atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                                  pott(:,:,:), qv(:,:,:), qc(:,:,:), & ! [IN]
                                                  dz(:,:,:),                         & ! [IN]
                                                  dens(:,:,:),                       & ! [INOUT]
                                                  temp(:,:,:), pres(:,:,:)           ) ! [OUT]

    !$acc end data

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_3D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D), not to reverse from TOA
  subroutine ATMOS_HYDROSTATIC_buildrho_real_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pott, qv, qc, &
       cz,           &
       pres,         &
       dens, temp    )
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: cz  (KA,IA,JA)

    real(RP), intent(inout) :: pres(KA,IA,JA) !< pressure              [Pa]

    real(RP), intent(out)   :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]

    real(RP) :: dz(KA,IA,JA)

    real(RP) :: pott_toa(IA,JA) !< surface potential temperature [K]
    real(RP) :: qv_toa  (IA,JA) !< surface water vapor           [kg/kg]
    real(RP) :: qc_toa  (IA,JA) !< surface liquid water          [kg/kg]

    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPtot

    integer  :: kref(IA,JA)
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(pott, qv, qc, cz) copy(pres), copyout(dens, temp) create(dz, pott_toa, qv_toa, qc_toa, kref)

    !--- from surface to lowermost atmosphere

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE-1
          dz(k,i,j) = CZ(k+1,i,j) - CZ(k,i,j) ! distance from cell center to cell center
       enddo
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       pott_toa(i,j) = pott(KE,i,j)
       qv_toa  (i,j) = qv  (KE,i,j)
       qc_toa  (i,j) = qc  (KE,i,j)
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       kref(i,j) = HYDROSTATIC_buildrho_real_kref + KS - 1
       !$acc loop seq
       do k = kref(i,j), KE
          if ( pres(k,i,j) .ne. UNDEF ) then
             kref(i,j) = k
             exit
          end if
       end do
    end do
    end do
    !$acc end kernels

    ! calc density at reference level
    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(Rtot,CVtot,CPtot,k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       k = kref(i,j)
       Rtot = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
            + Rvap  * qv(k,i,j)
       CVtot = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
             + CV_VAPOR * qv(k,i,j)                          &
             + CV_WATER * qc(k,i,j)
       CPtot = CPdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
             + CP_VAPOR * qv(k,i,j)                          &
             + CP_WATER * qc(k,i,j)
       dens(k,i,j) = P00 / ( Rtot * pott(k,i,j) ) * ( pres(k,i,j)/P00 )**(CVtot/CPtot)
    enddo
    enddo
    !$acc end kernels

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                              pott(:,:,:), qv(:,:,:), qc(:,:,:), & ! [IN]
                                              dz  (:,:,:),                       & ! [IN]
                                              dens(:,:,:),                       & ! [INOUT]
                                              temp(:,:,:), pres(:,:,:),          & ! [OUT]
                                              kref = kref                        ) ! [IN]
    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                                  pott(:,:,:), qv(:,:,:), qc(:,:,:), & ! [IN]
                                                  dz  (:,:,:),                       & ! [IN]
                                                  dens(:,:,:),                       & ! [INOUT]
                                                  temp(:,:,:), pres(:,:,:),          & ! [OUT]
                                                  kref = kref                        ) ! [IN]

!!$    call ATMOS_HYDROSTATIC_buildrho_atmos_2D( IA, IS, IE, JA, JS, JE, &
!!$                                              pott_toa(:,:), qv_toa(:,:), qc_toa(:,:),    & ! [IN]
!!$                                              dens(KE,:,:),                               & ! [IN]
!!$                                              pott(KE,:,:), qv(KE,:,:), qc(KE,:,:),       & ! [IN]
!!$                                              dz(KE+1,:,:), KE+1,                         & ! [IN]
!!$                                              dens_toa(:,:), temp_toa(:,:), pres_toa(:,:) ) ! [OUT]

    ! density at TOA
    !$acc kernels
    dens(   1:KS-1,:,:) = 0.0_RP ! fill dummy
!!$    dens(KE+2:KA  ,:,:) = 0.0_RP ! fill dummy
    dens(KE+1:KA  ,:,:) = 0.0_RP ! fill dummy
    !$acc end kernels

    !$acc end data

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_real_3D

  !-----------------------------------------------------------------------------
  !> Build up density (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
       pott_L2, qv_L2, qc_L2,          &
       dens_L1, pott_L1, qv_L1, qc_L1, &
       dz, k,                          &
       dens_L2, temp_L2, pres_L2,      &
       converged                       )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    !$acc routine seq
    implicit none
    real(RP), intent(in)  :: pott_L2 !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2   !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2   !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dens_L1 !< density               at layer 1 [Pa]
    real(RP), intent(in)  :: pott_L1 !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1   !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1   !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dz      !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k       !< for monitor

    real(RP), intent(out) :: dens_L2 !< density               at layer 2 [kg/m3]
    real(RP), intent(out) :: temp_L2 !< temperature           at layer 2 [K]
    real(RP), intent(out) :: pres_L2 !< pressure              at layer 2 [Pa]
    logical,  intent(out) :: converged

    real(RP) :: Rtot_L1  , Rtot_L2
    real(RP) :: CVtot_L1 , CVtot_L2
    real(RP) :: CPtot_L1 , CPtot_L2
    real(RP) :: CPovCV_L1, CPovCV_L2

    real(RP) :: pres_L1
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    !---------------------------------------------------------------------------

    Rtot_L1   = Rdry  * ( 1.0_RP - qv_L1 - qc_L1 ) &
              + Rvap  * qv_L1
    CVtot_L1  = CVdry * ( 1.0_RP - qv_L1 - qc_L1 ) &
              + CV_VAPOR * qv_L1                   &
              + CV_WATER * qc_L1
    CPtot_L1  = CPdry * ( 1.0_RP - qv_L1 - qc_L1 ) &
              + CP_VAPOR * qv_L1                   &
              + CP_WATER * qc_L1
    CPovCV_L1 = CPtot_L1 / CVtot_L1

    Rtot_L2   = Rdry  * ( 1.0_RP - qv_L2 - qc_L2 ) &
              + Rvap  * qv_L2
    CVtot_L2  = CVdry * ( 1.0_RP - qv_L2 - qc_L2 ) &
              + CV_VAPOR * qv_L2                   &
              + CV_WATER * qc_L2
    CPtot_L2  = CPdry * ( 1.0_RP - qv_L2 - qc_L2 ) &
              + CP_VAPOR * qv_L2                   &
              + CP_WATER * qc_L2
    CPovCV_L2 = CPtot_L2 / CVtot_L2

    dens_s  = 0.0_RP
    dens_L2 = dens_L1 ! first guess

    pres_L1 = P00 * ( dens_L1 * Rtot_L1 * pott_L1 / P00 )**CPovCV_L1

    converged = .false.
    do ite = 1, itelim
       if ( abs(dens_L2-dens_s) <= criteria ) then
          converged = .true.
          exit
       endif

       dens_s = dens_L2
       pres_L2 = P00 * ( dens_s  * Rtot_L2 * pott_L2 / P00 )**CPovCV_L2

       dhyd = + ( pres_L1 - pres_L2 ) / dz           & ! dp/dz
              - GRAV * 0.5_RP * ( dens_L1 + dens_s )   ! rho*g

       dgrd = - pres_L2  * CPovCV_L2  / dens_s / dz &
              - GRAV * 0.5_RP

       dens_L2 = max( dens_s - dhyd/dgrd, dens_s * 0.1_RP )

       if ( dens_L2*0.0_RP /= 0.0_RP ) exit
    enddo

    if ( converged ) then
       pres_L2 = P00 * ( dens_L2 * Rtot_L2 * pott_L2 / P00 )**CPovCV_L2
       temp_L2 = pres_L2 / ( dens_L2 * Rtot_L2 )
#ifndef _OPENACC
    else
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_0D",*) 'iteration not converged!', &
                  k,dens_L2,ite,dens_s,dhyd,dgrd
       call PRC_abort
#endif
    endif

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_0D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D( &
       KA, KS, KE, &
       pott, qv, qc, &
       dz,           &
#ifdef _OPENACC
       Rtot, CPovCV, &
#endif
       dens,         &
       temp, pres,   &
       converged     )
    !$acc routine seq
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: temp(KA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    logical,  intent(out)   :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: Rtot  (KA)
    real(RP), intent(out) :: CPovCV(KA)
#else
    real(RP) :: Rtot  (KA)
    real(RP) :: CPovCV(KA)
#endif
    real(RP) :: CVtot
    real(RP) :: CPtot

    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot  (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                 + Rvap  * qv(k)
       CVtot     = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CV_VAPOR * qv(k)                   &
                 + CV_WATER * qc(k)
       CPtot     = CPdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CP_VAPOR * qv(k)                   &
                 + CP_WATER * qc(k)
       CPovCV(k) = CPtot / CVtot
    enddo

    pres(KS) = P00 * ( dens(KS) * Rtot(KS) * pott(KS) / P00 )**CPovCV(KS)

    do k = KS+1, KE

       dens_s  = 0.0_RP
       dens(k) = dens(k-1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)
          pres(k) = P00 * ( dens_s * Rtot(k) * pott(k) / P00 )**CPovCV(k)

          dhyd = + ( pres(k-1) - pres(k) ) / dz(k-1)    & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s ) ! rho*g

          dgrd = - CPovCV(k) * pres(k) / dens_s / dz(k-1) &
                 - GRAV * 0.5_RP

          dens(k) = dens_s - dhyd/dgrd

          if ( dens(k)*0.0_RP /= 0.0_RP ) exit
       enddo

       if ( .NOT. converged ) then
#ifndef _OPENACC
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
#endif
          exit
       endif

       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)

    enddo

    if ( converged ) then
       do k = KS, KE
          temp(k) = pres(k) / ( dens(k) * Rtot(k) )
       enddo

       dens(KE+1:KA  ) = dens(KE)
       pres(KE+1:KA  ) = pres(KE)
       temp(KE+1:KA  ) = temp(KE)
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from upermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D( &
       KA, KS, KE, &
       pott, qv, qc, &
       dz,           &
#ifdef _OPENACC
       Rtot, CPovCV, &
#endif
       dens,         &
       temp, pres,   &
       converged     )
    !$acc routine seq
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: temp(KA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    logical,  intent(out)   :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: Rtot  (KA)
    real(RP), intent(out) :: CPovCV(KA)
#else
    real(RP) :: Rtot  (KA)
    real(RP) :: CPovCV(KA)
#endif
    real(RP) :: CVtot
    real(RP) :: CPtot

    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot  (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                 + Rvap  * qv(k)
       CVtot     = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CV_VAPOR * qv(k)                   &
                 + CV_WATER * qc(k)
       CPtot     = CPdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CP_VAPOR * qv(k)                   &
                 + CP_WATER * qc(k)
       CPovCV(k) = CPtot / CVtot
    enddo

    pres(KE) = P00 * ( dens(KE) * Rtot(KE) * pott(KE) / P00 )**CPovCV(KE)

    do k = KE-1, KS, -1

       dens_s  = 0.0_RP
       dens(k) = dens(k+1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)
          pres(k) = P00 * ( dens_s * Rtot(k) * pott(k) / P00 )**CPovCV(k)

          dhyd = - ( pres(k+1) - pres(k) ) / dz(k)        & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k+1) + dens_s ) ! rho*g

          dgrd = + CPovCV(k) * pres(k) / dens_s / dz(k) &
                 - GRAV * 0.5_RP

          dens(k) = dens_s - dhyd/dgrd

          if ( dens(k)*0.0_RP /= 0.0_RP ) exit
       enddo

       if ( .NOT. converged ) then
#ifndef _OPENACC
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
#endif
          exit
       endif

       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)

    enddo

    if ( converged ) then
       do k = KS, KE
          temp(k) = pres(k) / ( dens(k) * Rtot(k) )
       enddo

!!$       dens(   1:KS-1) = dens(KS)
!!$       pres(   1:KS-1) = pres(KS)
!!$       temp(   1:KS-1) = temp(KS)
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D

  !-----------------------------------------------------------------------------
  !> Build up density (2D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_2D( &
       IA, IS, IE, JA, JS, JE, &
       pott_L2, qv_L2, qc_L2,          &
       dens_L1, pott_L1, qv_L1, qc_L1, &
       dz, k,                          &
       dens_L2, temp_L2, pres_L2       )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: pott_L2(IA,JA) !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2  (IA,JA) !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2  (IA,JA) !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dens_L1(IA,JA) !< density               at layer 1 [Pa]
    real(RP), intent(in)  :: pott_L1(IA,JA) !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1  (IA,JA) !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1  (IA,JA) !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dz     (IA,JA) !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k              !< for monitor

    real(RP), intent(out) :: dens_L2(IA,JA) !< density               at layer 2 [kg/m3]
    real(RP), intent(out) :: temp_L2(IA,JA) !< temperature           at layer 2 [K]
    real(RP), intent(out) :: pres_L2(IA,JA) !< pressure              at layer 2 [Pa]

    logical :: converged, flag
    integer :: i, j
    !---------------------------------------------------------------------------

    converged = .true.

    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels copyin(pott_L2, qv_l2, qc_l2, dens_l1, pott_L1, qv_L1, qc_L1, dz) copyout(dens_L2, temp_L2, pres_L2)
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
            pott_L2(i,j), qv_L2(i,j), qc_L2(i,j),               & ! [IN]
            dens_L1(i,j), pott_L1(i,j), qv_L1(i,j), qc_L1(i,j), & ! [IN]
            dz(i,j), k,                                         & ! [IN]
            dens_L2(i,j), temp_L2(i,j), pres_L2(i,j),           & ! [OUT]
            flag                                                ) ! [OUT]
       converged = converged .and. flag
    enddo
    enddo
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_2D",*) "not converged"
       call PRC_abort
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_2D

  !-----------------------------------------------------------------------------
  !> Build up density (2D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D( &
       IA, IS, IE, JA, JS, JE, &
       pott_L1, qv_L1, qc_L1,          &
       dens_L2, pott_L2, qv_L2, qc_L2, &
       dz, k,                          &
       dens_L1, temp_L1, pres_L1       )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: pott_L1(IA,JA) !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1  (IA,JA) !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1  (IA,JA) !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dens_L2(IA,JA) !< density               at layer 2 [Pa]
    real(RP), intent(in)  :: pott_L2(IA,JA) !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2  (IA,JA) !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2  (IA,JA) !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dz     (IA,JA) !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k              !< for monitor

    real(RP), intent(out) :: dens_L1(IA,JA) !< density               at layer 1 [kg/m3]
    real(RP), intent(out) :: temp_L1(IA,JA) !< temperature           at layer 1 [K]
    real(RP), intent(out) :: pres_L1(IA,JA) !< pressure              at layer 1 [Pa]

    logical :: converged, flag
    integer :: i, j
    !---------------------------------------------------------------------------

    converged = .true.

    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels copyin(pott_L1, qv_L1, qc_L1, dens_L2, pott_L2, qv_L2, qc_L2, dz) copyout(dens_L1, temp_L1, pres_L1)
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
            pott_L1(i,j), qv_L1(i,j), qc_L1(i,j),               & ! [IN]
            dens_L2(i,j), pott_L2(i,j), qv_L2(i,j), qc_L2(i,j), & ! [IN]
            -dz(i,j), k,                                        & ! [IN]
            dens_L1(i,j), temp_L1(i,j), pres_L1(i,j),           & ! [OUT]
            flag                                                ) ! [OUT]
       converged = converged .and. flag
    enddo
    enddo
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D",*) "not converged"
       call PRC_abort
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pott, qv, qc, &
       dz,           &
       dens,         &
       temp, pres,   &
       kref          )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA) !< distance between the layer (center) [m]

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]

    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]

    integer, intent(in), optional, target :: kref(IA,JA)

    integer, pointer :: kref_(:,:)

#ifdef _OPENACC
    real(RP) :: work1(KA,IA,JA)
    real(RP) :: work2(KA,IA,JA)
#endif

    logical :: converged, flag
    integer :: i, j
    !---------------------------------------------------------------------------

    if ( present(kref) ) then
       kref_ => kref
    else
       allocate( kref_(IA,JA) )
       !$acc enter data create(kref_)
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          kref_(i,j) = KS
       end do
       end do
       !$acc end kernels
    end if

    converged = .true.

    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels copyin(pott, qv, qc, dz, kref_) copy(dens) copyout(temp, pres)
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       if ( kref_(i,j) < KE ) then
          call ATMOS_HYDROSTATIC_buildrho_atmos_1D( KA, kref_(i,j), KE, &
                                                    pott(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                    dz(:,i,j),                         & ! [IN]
#ifdef _OPENACC
                                                    work1(:,i,j), work2(:,i,j),        & ! [WORK]
#endif
                                                    dens(:,i,j),                       & ! [INOUT]
                                                    temp(:,i,j), pres(:,i,j),          & ! [OUT]
                                                    flag                               ) ! [OUT]
          converged = converged .and. flag
       end if
    enddo
    enddo
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_3D",*) "not converged"
       call PRC_abort
    end if

    if ( .not. present(kref) ) then
       !$acc exit data delete(kref_)
       deallocate( kref_ )
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pott, qv, qc, &
       dz,           &
       dens,         &
       temp, pres,   &
       kref          )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA) !< distance between the layer (center) [m]
    integer,  intent(in), optional, target :: kref(IA,JA)

    integer, pointer :: kref_(:,:)
    logical :: converged, flag
    integer :: i, j
#ifdef _OPENACC
    real(RP) :: work1(KA,IA,JA)
    real(RP) :: work2(KA,IA,JA)
#endif
    !---------------------------------------------------------------------------

    converged = .true.

    if ( present(kref) ) then
       kref_ => kref

       !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
       !$acc kernels copyin(pott, qv, qc, dz, kref_) copy(dens) copyout(temp, pres)
       !$acc loop independent collapse(2) &
       !$acc private(flag) reduction(.and.: converged)
       do j = JS, JE
       do i = IS, IE
          if ( kref_(i,j) > KS ) then
             call ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D( KA, KS, kref_(i,j), &
                                                           pott(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                           dz(:,i,j),                         & ! [IN]
#ifdef _OPENACC
                                                           work1(:,i,j), work2(:,i,j),        & ! [WORK]
#endif
                                                           dens(:,i,j),                       & ! [INOUT]
                                                           temp(:,i,j), pres(:,i,j),          & ! [OUT]
                                                           flag                               ) ! [OUT]
             converged = converged .and. flag
          end if
       enddo
       enddo
       !$acc end kernels

    else
       !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
       !$acc kernels copyin(pott, qv, qc, dz) copy(dens) copyout(temp, pres)
       !$acc loop independent collapse(2) &
       !$acc private(flag) reduction(.and.: converged)
       do j = JS, JE
       do i = IS, IE
          call ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D( KA, KS, KE,                        &
                                                        pott(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                        dz(:,i,j),                         & ! [IN]
#ifdef _OPENACC
                                                        work1(:,i,j), work2(:,i,j),        & ! [WORK]
#endif
                                                        dens(:,i,j),                       & ! [INOUT]
                                                        temp(:,i,j), pres(:,i,j),          & ! [OUT]
                                                        flag                               ) ! [OUT]
          converged = converged .and. flag
       enddo
       enddo
       !$acc end kernels
    end if

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D",*) "not converged"
       call PRC_abort
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_1D( &
       KA, KS, KE, &
       temp, qv, qc,                       &
       pres_sfc, temp_sfc, qv_sfc, qc_sfc, &
       cz, fz,                             &
#ifdef _OPENACC
       dz, work1, work2, work3,            &
#endif
       dens, pott, pres, pott_sfc,         &
       converged                           )
    !$acc routine seq
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)  :: pres_sfc !< surface pressure              [Pa]
    real(RP), intent(in)  :: temp_sfc !< surface temperature           [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface liquid water          [kg/kg]
    real(RP), intent(in)  :: cz(  KA)
    real(RP), intent(in)  :: fz(0:KA)

    real(RP), intent(out) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out) :: pott(KA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure              [Pa]
    real(RP), intent(out) :: pott_sfc !< surface potential temperature [K]
    logical,  intent(out) :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: dz(KA)
    real(RP), intent(out) :: work1(KA)
    real(RP), intent(out) :: work2(KA)
    real(RP), intent(out) :: work3(KA)
#else
    real(RP) :: dz(KA)
#endif

    real(RP) :: dens_sfc

    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: CPtot_sfc
    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPtot

    real(RP) :: RovCP_sfc
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite

    integer :: k
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CV_VAPOR * qv_sfc                    &
               + CV_WATER * qc_sfc
    CPtot_sfc  = CPdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CP_VAPOR * qv_sfc                    &
               + CP_WATER * qc_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CV_VAPOR * qv(KS)                    &
           + CV_WATER * qc(KS)
    CPtot  = CPdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CP_VAPOR * qv(KS)                    &
           + CP_WATER * qc(KS)

    ! density at surface
    RovCP_sfc = Rtot_sfc / CPtot_sfc
    dens_sfc  = pres_sfc / ( Rtot_sfc * temp_sfc )
    pott_sfc  = temp_sfc * ( P00/pres_sfc )**RovCP_sfc

    ! make density at lowermost cell center
    dens_s   = 0.0_RP
    dens(KS) = dens_sfc ! first guess

    dz(KS-1) = CZ(KS) - FZ(KS-1)
    do k = KS, KE-1
       dz(k) = CZ(k+1) - CZ(k)
    end do

    converged = .false.
    do ite = 1, itelim
       if ( abs(dens(KS)-dens_s) <= criteria ) then
          converged = .true.
          exit
       endif

       dens_s = dens(KS)

       dhyd = + ( pres_sfc - dens_s * Rtot * temp(KS) ) / dz(KS-1) & ! dp/dz
              - GRAV * 0.5_RP * ( dens_sfc + dens_s )                ! rho*g

       dgrd = - Rtot * temp(KS) / dz(KS-1) &
              - GRAV * 0.5_RP

       dens(KS) = dens_s - dhyd/dgrd

       if ( dens(KS)*0.0_RP /= 0.0_RP ) exit
    enddo

    if ( converged ) then
       !--- from lowermost atmosphere to top of atmosphere
       call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( KA, KS, KE, &
                                                        temp(:), qv(:), qc(:), & ! [IN]
                                                        dz  (:),               & ! [IN]
#ifdef _OPENACC
                                                        work1(:), work2(:), work3(:), & ! [WORK]
#endif
                                                        dens(:),               & ! [INOUT]
                                                        pott(:), pres(:),      & ! [OUT]
                                                        converged              ) ! [OUT]
    end if

#ifndef _OPENACC
    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_1D",*) 'iteration not converged!', &
                  dens(KS),ite,dens_s,dhyd,dgrd
       call PRC_abort
    endif
#endif

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_1D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, qv, qc,                       &
       pres_sfc, temp_sfc, qv_sfc, qc_sfc, &
       cz, fz,                             &
       dens, pott, pres,                   &
       pott_sfc                            )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp(KA,IA,JA)  !< temperature          [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA)  !< water vapor          [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA)  !< liquid water         [kg/kg]
    real(RP), intent(in)  :: pres_sfc(IA,JA) !< surface pressure     [Pa]
    real(RP), intent(in)  :: temp_sfc(IA,JA) !< surface temperature  [K]
    real(RP), intent(in)  :: qv_sfc  (IA,JA) !< surface water vapor  [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (IA,JA) !< surface liquid water [kg/kg]
    real(RP), intent(in)  :: cz(KA,IA,JA)
    real(RP), intent(in)  :: fz(KA,IA,JA)

    real(RP), intent(out) :: dens(KA,IA,JA)  !< density               [kg/m3]
    real(RP), intent(out) :: pott(KA,IA,JA)  !< potential temperature [K]
    real(RP), intent(out) :: pres(KA,IA,JA)  !< pressure              [Pa]
    real(RP), intent(out) :: pott_sfc(IA,JA) !< surface potential temperature [K]

    logical :: converged, flag

#ifdef _OPENACC
    real(RP) :: work1(KA,IA,JA)
    real(RP) :: work2(KA,IA,JA)
    real(RP) :: work3(KA,IA,JA)
    real(RP) :: work4(KA,IA,JA)
#endif

    integer  :: i, j
    !---------------------------------------------------------------------------

    converged = .true.

    !--- from surface to lowermost atmosphere
    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels copyin(temp, qv, qc, pres_sfc, temp_sfc, qv_sfc, qc_sfc, cz, fz) copyout(dens, pott, pres, pott_sfc)
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_bytemp_1D( KA, KS, KE, &
                                                  temp(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                                  pres_sfc(i,j), temp_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                                  cz(:,i,j), fz(:,i,j),                                   & ! [IN]
#ifdef _OPENACC
                                                  work1(:,i,j), work2(:,i,j), work3(:,i,j), work4(:,i,j), & ! [WORK]
#endif
                                                  dens(:,i,j), pott(:,i,j), pres(:,i,j),                  & ! [OUT]
                                                  pott_sfc(i,j),                                          & ! [OUT]
                                                  flag                                                    ) ! [OUT]
       converged = converged .and. flag
    enddo
    enddo
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_3D",*) "not converged"
       call PRC_abort
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_3D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( &
       KA, KS, KE, &
       temp, qv, qc, &
       dz,           &
#ifdef _OPENACC
       Rtot,         &
       CVtot,        &
       CPtot,        &
#endif
       dens,         &
       pott, pres,   &
       converged     )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    !$acc routine seq
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: temp(KA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: pott(KA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    logical,  intent(out)   :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: Rtot  (KA)
    real(RP), intent(out) :: CVtot (KA)
    real(RP), intent(out) :: CPtot (KA)
#else
    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)
    real(RP) :: CPtot (KA)
#endif

    real(RP) :: RovCP
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                + Rvap  * qv(k)
       CVtot(k) = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CV_VAPOR * qv(k)                   &
                + CV_WATER * qc(k)
       CPtot(k) = CPdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CP_VAPOR * qv(k)                   &
                + CP_WATER * qc(k)
    enddo

    pres(KS) = dens(KS) * Rtot(KS) * temp(KS)

    do k = KS+1, KE

       dens_s  = 0.0_RP
       dens(k) = dens(k-1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)

          dhyd = + ( pres(k-1) - dens_s * Rtot(k) * temp(k) ) / dz(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )            ! rho*g

          dgrd = - Rtot(k) * temp(k) / dz(k-1) &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if ( dens(k)*0.0_RP /= 0.0_RP ) exit
       enddo

       if ( .NOT. converged ) then
#ifndef _OPENACC
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
#endif
          exit
       endif

       pres(k) = dens(k) * Rtot(k) * temp(k)

    enddo

    if ( converged ) then
       do k = KS, KE
          RovCP   = Rtot(k) / CPtot(k)
          pott(k) = temp(k) * ( P00 / pres(k) )**RovCP
       enddo
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from upermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_rev_1D( &
       KA, KS, KE, &
       temp, qv, qc, &
       dz,           &
#ifdef _OPENACC
       Rtot,         &
       CVtot,        &
       CPtot,        &
#endif
       dens,         &
       pott, pres,   &
       converged     )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    !$acc routine seq
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: temp(KA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: pott(KA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    logical,  intent(out)   :: converged

#ifdef _OPENACC
    real(RP), intent(out) :: Rtot  (KA)
    real(RP), intent(out) :: CVtot (KA)
    real(RP), intent(out) :: CPtot (KA)
#else
    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)
    real(RP) :: CPtot (KA)
#endif

    real(RP) :: RovCP
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                + Rvap  * qv(k)
       CVtot(k) = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CV_VAPOR * qv(k)                   &
                + CV_WATER * qc(k)
       CPtot(k) = CPdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CP_VAPOR * qv(k)                   &
                + CP_WATER * qc(k)
    enddo

    pres(KE) = dens(KE) * Rtot(KE) * temp(KE)

    do k = KE-1, KE, -1

       dens_s  = 0.0_RP
       dens(k) = dens(k+1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)

          dhyd = - ( pres(k+1) - dens_s * Rtot(k) * temp(k) ) / dz(k) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k+1) + dens_s )               ! rho*g

          dgrd = + Rtot(k) * temp(k) / dz(k) &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if ( dens(k)*0.0_RP /= 0.0_RP ) exit
       enddo

       if ( .NOT. converged ) then
#ifdef _OPENACC
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_rev_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
#endif
          exit
       endif

       pres(k) = dens(k) * Rtot(k) * temp(k)

    enddo

    if ( converged ) then
       do k = KS, KE
          RovCP   = Rtot(k) / CPtot(k)
          pott(k) = temp(k) * ( P00 / pres(k) )**RovCP
       enddo
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_rev_1D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, qv, qc, &
       dz,           &
       dens,         &
       pott, pres    )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)    :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA)

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]

    logical :: converged, flag

#ifdef _OPENACC
    real(RP) :: work1(KA,IA,JA)
    real(RP) :: work2(KA,IA,JA)
    real(RP) :: work3(KA,IA,JA)
#endif
    integer  :: i, j
    !---------------------------------------------------------------------------

    converged = .true.

    !$omp parallel do OMP_SCHEDULE_ collapse(2) private(flag) reduction(.and.:converged)
    !$acc kernels copyin(temp, qv, qc, dz) copy(dens) copyout(pott, pres)
    !$acc loop independent collapse(2) &
    !$acc private(flag) reduction(.and.: converged)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( KA, KS, KE, &
                                                        temp(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                        dz(:,i,j),                         & ! [IN]
#ifdef _OPENACC
                                                        work1(:,i,j), work2(:,i,j), work3(:,i,j), & ! [WORK]
#endif
                                                        dens(:,i,j),                       & ! [INOUT]
                                                        pott(:,i,j), pres(:,i,j),          & ! [OUT]
                                                        flag                               ) ! [OUT]
       converged = converged .and. flag
    enddo
    enddo
    !$acc end kernels

    if ( .not. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D",*) "iteration not converged"
       call PRC_abort
    end if

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D( &
       KA, KS, KE, &
       pres, temp, qv, &
       cz,             &
       mslp            )
    !$acc routine seq
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: pres(KA) !< pressure          [Pa]
    real(RP), intent(in)  :: temp(KA) !< air temperature   [K]
    real(RP), intent(in)  :: qv  (KA) !< specific humidity [kg/kg]
    real(RP), intent(in)  :: cz  (KA) !< height from MSL   [m]

    real(RP), intent(out) :: mslp !< mean sea-level pressure [Pa]

    ! work
    integer  :: kref
    real(RP) :: vtemp
    !---------------------------------------------------------------------------

    kref = HYDROSTATIC_barometric_law_mslp_kref + KS - 1

    ! virtual temperature
    vtemp = temp(kref) * (1.0_RP + EPSTvap * qv(kref) )

    ! barometric law assuming constant lapse rate
    mslp = pres(kref) * ( ( vtemp + LAPS * cz(kref) ) / vtemp ) ** ( GRAV / ( Rdry * LAPS ) )

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (2D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_2D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pres, temp, qv, &
       cz,             &
       mslp            )
    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: pres(KA,IA,JA) !< pressure          [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< air temperature   [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< specific humidity [kg/kg]
    real(RP), intent(in)  :: cz  (KA,IA,JA) !< height from MSL   [m]

    real(RP), intent(out) :: mslp(IA,JA) !< mean sea-level pressure [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels copyin(pres, temp, qv, cz) copyout(mslp)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_barometric_law_mslp_0D( KA, KS, KE, &
                                                      pres(:,i,j), temp(:,i,j), qv(:,i,j), & ! [IN]
                                                      cz(:,i,j),                           & ! [IN]
                                                      mslp(i,j)                            ) ! [OUT]
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_2D

  !-----------------------------------------------------------------------------
  !> Calculate surface pressure from barometric law (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_HYDROSTATIC_barometric_law_pres_0D( &
       mslp, temp, &
       dz,         &
       pres        )
    !$acc routine seq
    implicit none
    real(RP), intent(in)  :: mslp !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: temp !< surface air temperature [K]
    real(RP), intent(in)  :: dz   !< surface height from MSL [m]

    real(RP), intent(out) :: pres !< surface pressure        [Pa]

    ! work
    real(RP) :: TM
    !---------------------------------------------------------------------------

    TM = temp + LAPS * dz * 0.5_RP ! column-mean air temperature

    ! barometric law
    pres = mslp / exp( GRAV * dz / ( Rdry * TM ) )

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_pres_0D

  !-----------------------------------------------------------------------------
  !> Calculate surface pressure from barometric law (2D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_pres_2D( &
       IA, IS, IE, JA, JS, JE, &
       mslp, temp, &
       dz,         &
       pres        )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: mslp(IA,JA) !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: temp(IA,JA) !< surface air temperature [K]
    real(RP), intent(in)  :: dz  (IA,JA) !< surface height from MSL [m]

    real(RP), intent(out) :: pres(IA,JA) !< surface pressure        [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$acc kernels copyin(mslp, temp, dz) copyout(pres)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_barometric_law_pres_0D( mslp(i,j), temp(i,j), dz(i,j), & ! [IN]
                                                      pres(i,j)                      ) ! [OUT]
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_pres_2D

end module scale_atmos_hydrostatic
