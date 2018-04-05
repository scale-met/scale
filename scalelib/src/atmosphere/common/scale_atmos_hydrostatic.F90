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
  use scale_stdio
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
     P00     => CONST_PRE00
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

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_uselapserate, &
       HYDROSTATIC_buildrho_real_kref

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[HYDROSTATIC] / Categ[ATMOS] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_HYDROSTATIC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_HYDROSTATIC)

    criteria = sqrt( CONST_EPS )

    LOG_NEWLINE
    LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Use lapse rate for estimation of surface temperature? : ', HYDROSTATIC_uselapserate
    LOG_INFO("ATMOS_HYDROSTATIC_setup",*) 'Buildrho conversion criteria                          : ', criteria

    return
  end subroutine ATMOS_HYDROSTATIC_setup

  !-----------------------------------------------------------------------------
  !> Build up density from surface (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_1D( &
       KA, KS, KE, &
       pott, qv, qc,                       &
       pres_sfc, pott_sfc, qv_sfc, qc_sfc, &
       cz, fz,                             &
       dens, temp, pres,                   &
       temp_sfc                            )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)  :: pres_sfc !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface liquid water          [kg/kg]
    real(RP), intent(in)  :: cz  (KA)
    real(RP), intent(in)  :: fz  (0:KA)

    real(RP), intent(out) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA) !< temperature           [K]
    real(RP), intent(out) :: pres(KA) !< pressure              [Pa]
    real(RP), intent(out) :: temp_sfc !< surface temperature           [K]

    real(RP) :: dens_sfc

    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: CPtot_sfc
    real(RP) :: CPovCV_sfc
    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPtot
    real(RP) :: CPovCV

    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    real(RP) :: dz(KA)

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

    else ! use itelation

       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( pott(KS), qv(KS), qc(KS),           & ! [IN]
                                                 dens_sfc, pott_sfc, qv_sfc, qc_sfc, & ! [IN]
                                                 dz(KS-1), KS-1,                     & ! [IN]
                                                 dens(KS), temp(KS), pres(KS)        ) ! [OUT]

    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_1D( KA, KS, KE, &
                                              pott(:), qv(:), qc(:), & ! [IN]
                                              dz(:),                 & ! [IN]
                                              dens(:),               & ! [INOUT]
                                              temp(:), pres(:)       ) ! [OUT]
    ! fill dummy
    dens(   1:KS-1) = 0.0_RP
    dens(KE+1:KA  ) = 0.0_RP

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_1D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       pott, qv, qc,       &
       pres_sfc, pott_sfc, &
       qv_sfc, qc_sfc,     &
       cz, fz,             &
       dens, temp, pres,   &
       temp_sfc            )
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
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

    real(RP) :: dens_mean

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_1D( KA, KS, KE, &
                                           pott(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                           pres_sfc(i,j), pott_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                           cz(:,i,j), fz(:,i,j),                                   & ! [IN]
                                           dens(:,i,j), temp(:,i,j), pres(:,i,j),                  & ! [OUT]
                                           temp_sfc(i,j)                                           ) ! [OUT]
    end do
    end do

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
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
    end do
    end do

    call ATMOS_HYDROSTATIC_buildrho_atmos_2D( IA, IS, IE, JA, JS, JE, &
                                              pott_toa(:,:), qv_toa(:,:), qc_toa(:,:),            & ! [IN]
                                              dens(KE,:,:), pott(KE,:,:), qv(KE,:,:), qc(KE,:,:), & ! [IN]
                                              dz_top(:,:), KE+1,                                  & ! [IN]
                                              dens_toa(:,:), temp_toa(:,:), pres_toa(:,:)         ) ! [OUT]

    call COMM_horizontal_mean( dens_mean, dens_toa(:,:) )
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       dens_toa(i,j) = dens_mean
    enddo
    enddo

    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D( IA, IS, IE, JA, JS, JE, &
                                                  pott(KE,:,:), qv(KE,:,:), qc(KE,:,:),                   & ! [IN]
                                                  dens_toa(:,:), pott_toa(:,:), qv_toa(:,:), qc_toa(:,:), & ! [IN]
                                                  dz_top(:,:), KE+1,                                      & ! [IN]
                                                  dens(KE,:,:), temp(KE,:,:), pres(KE,:,:)                ) ! [OUT]

    !--- from top of atmosphere to lowermost atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                                  pott(:,:,:), qv(:,:,:), qc(:,:,:), & ! [IN]
                                                  dz(:,:,:),                         & ! [IN]
                                                  dens(:,:,:),                       & ! [INOUT]
                                                  temp(:,:,:), pres(:,:,:)           ) ! [OUT]

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

    integer  :: kref
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE-1
          dz(k,i,j) = CZ(k+1,i,j) - CZ(k,i,j) ! distance from cell center to cell center
       enddo
    enddo
    enddo

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       pott_toa(i,j) = pott(KE,i,j)
       qv_toa  (i,j) = qv  (KE,i,j)
       qc_toa  (i,j) = qc  (KE,i,j)
    enddo
    enddo

    kref = HYDROSTATIC_buildrho_real_kref + KS - 1

    ! calc density at reference level
    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(Rtot,CVtot,CPtot)
    do j = JS, JE
    do i = IS, IE
       Rtot = Rdry  * ( 1.0_RP - qv(kref,i,j) - qc(kref,i,j) ) &
            + Rvap  * qv(kref,i,j)
       CVtot = CVdry * ( 1.0_RP - qv(kref,i,j) - qc(kref,i,j) ) &
             + CV_VAPOR * qv(kref,i,j)                          &
             + CV_WATER * qc(kref,i,j)
       CPtot = CPdry * ( 1.0_RP - qv(kref,i,j) - qc(kref,i,j) ) &
             + CP_VAPOR * qv(kref,i,j)                          &
             + CP_WATER * qc(kref,i,j)
       dens(kref,i,j) = P00 / ( Rtot * pott(kref,i,j) ) * ( pres(kref,i,j)/P00 )**(CVtot/CPtot)
    enddo
    enddo

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
    dens(   1:KS-1,:,:) = 0.0_RP ! fill dummy
!!$    dens(KE+2:KA  ,:,:) = 0.0_RP ! fill dummy
    dens(KE+1:KA  ,:,:) = 0.0_RP ! fill dummy

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
       dens_L2, temp_L2, pres_L2       )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
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

    real(RP) :: Rtot_L1  , Rtot_L2
    real(RP) :: CVtot_L1 , CVtot_L2
    real(RP) :: CPtot_L1 , CPtot_L2
    real(RP) :: CPovCV_L1, CPovCV_L2

    real(RP) :: pres_L1
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged
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

       dens_L2 = dens_s - dhyd/dgrd

       if ( dens_L2*0.0_RP /= 0.0_RP ) exit
    enddo

    if ( .NOT. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_0D",*) 'iteration not converged!', &
                  k,dens_L2,ite,dens_s,dhyd,dgrd
       call PRC_abort
    endif

    pres_L2 = P00 * ( dens_L2 * Rtot_L2 * pott_L2 / P00 )**CPovCV_L2
    temp_L2 = pres_L2 / ( dens_L2 * Rtot_L2 )

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_0D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D( &
       KA, KS, KE, &
       pott, qv, qc, &
       dz,           &
       dens,         &
       temp, pres,   &
       kref          )
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

    integer, intent(in), optional :: kref

    real(RP) :: Rtot  (KA)
    real(RP) :: CPovCV(KA)
    real(RP) :: CVtot
    real(RP) :: CPtot

    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer :: kref_
    integer :: k
    !---------------------------------------------------------------------------

    if ( present(kref) ) then
       kref_ = kref
    else
       kref_ = KS
    end if

    do k = kref_, KE
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

    pres(kref_) = P00 * ( dens(kref_) * Rtot(kref_) * pott(kref_) / P00 )**CPovCV(kref_)

    do k = kref_+1, KE

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
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
       endif

       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)

    enddo

    do k = kref_, KE
       temp(k) = pres(k) / ( dens(k) * Rtot(k) )
    enddo

    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)
    pres(   1:KS-1) = pres(KS)
    pres(KE+1:KA  ) = pres(KE)
    temp(   1:KS-1) = temp(KS)
    temp(KE+1:KA  ) = temp(KE)

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from upermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D( &
       KA, KS, KE, &
       pott, qv, qc, &
       dz,           &
       dens,         &
       temp, pres,   &
       kref          )
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

    integer, intent(in), optional :: kref

    real(RP) :: Rtot  (KA)
    real(RP) :: CPovCV(KA)
    real(RP) :: CVtot
    real(RP) :: CPtot

    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer :: kref_
    integer :: k
    !---------------------------------------------------------------------------

    if ( present(kref) ) then
       kref_ = kref
    else
       kref_ = KE
    end if

    do k = KS, kref_
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

    pres(kref_) = P00 * ( dens(kref_) * Rtot(kref_) * pott(kref_) / P00 )**CPovCV(kref_)

    do k = kref_-1, KS, -1

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
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
       endif

       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)

    enddo

    do k = KS, kref_
       temp(k) = pres(k) / ( dens(k) * Rtot(k) )
    enddo

    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)
    pres(   1:KS-1) = pres(KS)
    pres(KE+1:KA  ) = pres(KE)
    temp(   1:KS-1) = temp(KS)
    temp(KE+1:KA  ) = temp(KE)

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

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
            pott_L2(i,j), qv_L2(i,j), qc_L2(i,j),               & ! [IN]
            dens_L1(i,j), pott_L1(i,j), qv_L1(i,j), qc_L1(i,j), & ! [IN]
            dz(i,j), k,                                         & ! [IN]
            dens_L2(i,j), temp_L2(i,j), pres_L2(i,j)            ) ! [OUT]
    enddo
    enddo

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

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
            pott_L1(i,j), qv_L1(i,j), qc_L1(i,j),               & ! [IN]
            dens_L2(i,j), pott_L2(i,j), qv_L2(i,j), qc_L2(i,j), & ! [IN]
            -dz(i,j), k,                                        & ! [IN]
            dens_L1(i,j), temp_L1(i,j), pres_L1(i,j)            ) ! [OUT]
    enddo
    enddo

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

    integer, intent(in), optional :: kref

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_1D( KA, KS, KE, &
                                                 pott(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                 dz(:,i,j),                         & ! [IN]
                                                 dens(:,i,j),                       & ! [INOUT]
                                                 temp(:,i,j), pres(:,i,j),          & ! [OUT]
                                                 kref=kref                          ) ! [IN]
    enddo
    enddo

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
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA) !< distance between the layer (center) [m]
    integer, intent(in), optional :: kref

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_atmos_rev_1D( KA, KS, KE, &
                                                     pott(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                     dz(:,i,j),                         & ! [IN]
                                                     dens(:,i,j),                       & ! [INOUT]
                                                     temp(:,i,j), pres(:,i,j),          & ! [OUT]
                                                     kref=kref                          ) ! [IN]
    enddo
    enddo

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
       dens, pott, pres, pott_sfc          )
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
    logical  :: converged

    real(RP) :: dz(KA)

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

    if ( .NOT. converged ) then
       LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_1D",*) 'iteration not converged!', &
                  dens(KS),ite,dens_s,dhyd,dgrd
       call PRC_abort
    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( KA, KS, KE, &
                                                     temp(:), qv(:), qc(:), & ! [IN]
                                                     dz  (:),               & ! [IN]
                                                     dens(:),               & ! [INOUT]
                                                     pott(:), pres(:)       ) ! [OUT]

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

    integer  :: i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_bytemp_1D( KA, KS, KE, &
                                                  temp(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                                  pres_sfc(i,j), temp_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                                  cz(:,i,j), fz(:,i,j),                                   & ! [IN]
                                                  dens(:,i,j), pott(:,i,j), pres(:,i,j),                  & ! [OUT]
                                                  pott_sfc(i,j)                                           ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_3D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( &
       KA, KS, KE, &
       temp, qv, qc, &
       dz,           &
       dens,         &
       pott, pres    )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: temp(KA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: pott(KA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]

    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)
    real(RP) :: CPtot (KA)

    real(RP) :: RovCP
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

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
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
       endif

       pres(k) = dens(k) * Rtot(k) * temp(k)

    enddo

    do k = KS, KE
       RovCP   = Rtot(k) / CPtot(k)
       pott(k) = temp(k) * ( P00 / pres(k) )**RovCP
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from upermost atmosphere (1D)
!OCL SERIAL
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_rev_1D( &
       KA, KS, KE, &
       temp, qv, qc, &
       dz,           &
       dens,         &
       pott, pres    )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_VAPOR, &
       CV_WATER, &
       CP_VAPOR, &
       CP_WATER
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)    :: temp(KA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA)

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]

    real(RP), intent(out)   :: pott(KA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]

    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)
    real(RP) :: CPtot (KA)

    real(RP) :: RovCP
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

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
          LOG_ERROR("ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_rev_1D",*) 'iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_abort
       endif

       pres(k) = dens(k) * Rtot(k) * temp(k)

    enddo

    do k = KS, KE
       RovCP   = Rtot(k) / CPtot(k)
       pott(k) = temp(k) * ( P00 / pres(k) )**RovCP
    enddo

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

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( KA, KS, KE, &
                                                        temp(:,i,j), qv(:,i,j), qc(:,i,j), & ! [IN]
                                                        dz(:,i,j),                         & ! [IN]
                                                        dens(:,i,j),                       & ! [INOUT]
                                                        pott(:,i,j), pres(:,i,j)           ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D( &
       pres, temp, &
       dz,         &
       mslp        )
    implicit none
    real(RP), intent(in)  :: pres !< surface pressure        [Pa]
    real(RP), intent(in)  :: temp !< surface air temperature [K]
    real(RP), intent(in)  :: dz   !< surface height from MSL [m]

    real(RP), intent(out) :: mslp !< mean sea-level pressure [Pa]

    ! work
    real(RP) :: TM
    !---------------------------------------------------------------------------

    TM = temp + LAPS * dz * 0.5_RP ! column-mean air temperature

    ! barometric law
    mslp = pres * exp( GRAV * dz / ( Rdry * TM ) )

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (2D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_2D( &
       IA, IS, IE, JA, JS, JE, &
       pres, temp, &
       dz,         &
       mslp        )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: pres(IA,JA) !< surface pressure        [Pa]
    real(RP), intent(in)  :: temp(IA,JA) !< surface air temperature [K]
    real(RP), intent(in)  :: dz  (IA,JA) !< surface height from MSL [m]

    real(RP), intent(out) :: mslp(IA,JA) !< mean sea-level pressure [Pa]

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_barometric_law_mslp_0D( pres(i,j), temp(i,j), dz(i,j), & ! [IN]
                                                      mslp(i,j)                      ) ! [OUT]
    enddo
    enddo

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
    do j = JS, JE
    do i = IS, IE
       call ATMOS_HYDROSTATIC_barometric_law_pres_0D( mslp(i,j), temp(i,j), dz(i,j), & ! [IN]
                                                      pres(i,j)                      ) ! [OUT]
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_pres_2D

end module scale_atmos_hydrostatic
