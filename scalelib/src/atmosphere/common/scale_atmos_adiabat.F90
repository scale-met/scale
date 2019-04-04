!-------------------------------------------------------------------------------
!> module atmosphere / adiabat
!!
!! @par Description
!!          Moist adiabatic process for calculation of CAPE, CIN
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_adiabat
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
  public :: ATMOS_ADIABAT_setup
  public :: ATMOS_ADIABAT_cape
  public :: ATMOS_ADIABAT_liftparcel

  interface ATMOS_ADIABAT_cape
     module procedure ATMOS_ADIABAT_cape_1D
     module procedure ATMOS_ADIABAT_cape_3D
  end interface ATMOS_ADIABAT_cape

  interface ATMOS_ADIABAT_liftparcel
     module procedure ATMOS_ADIABAT_liftparcel_1D
     module procedure ATMOS_ADIABAT_liftparcel_3D
  end interface ATMOS_ADIABAT_liftparcel

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine ATMOS_ADIABAT_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_setup
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_setup

    call ATMOS_HYDROMETEOR_setup
    call ATMOS_SATURATION_setup

    return
  end subroutine ATMOS_ADIABAT_setup

  !-----------------------------------------------------------------------------
  !> Calc CAPE and CIN
  !> Type of parcel method: Pseudo-adiabatic ascend from lowermost layer of the model
  !> Reference: Emanuel(1994)
!OCL SERIAL
  subroutine ATMOS_ADIABAT_cape_1D( &
       KA, KS, KE, &
       Kstr,                         &
       DENS, TEMP, PRES,             &
       QV, QC, Qdry, Rtot, CPtot,    &
       CZ, FZ,                       &
       CAPE, CIN, LCL, LFC, LNB,     &
       DENS_p, TEMP_p, BUOY_p, QV_p, &
       converged                    )
    use scale_const, only: &
       GRAV => CONST_GRAV
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    integer,  intent(in)  :: Kstr
    real(RP), intent(in)  :: DENS (KA)
    real(RP), intent(in)  :: TEMP (KA)
    real(RP), intent(in)  :: PRES (KA)
    real(RP), intent(in)  :: QV   (KA)
    real(RP), intent(in)  :: QC   (KA)
    real(RP), intent(in)  :: Qdry (KA)
    real(RP), intent(in)  :: Rtot (KA)
    real(RP), intent(in)  :: CPtot(KA)
    real(RP), intent(in)  :: CZ   (  KA)
    real(RP), intent(in)  :: FZ   (0:KA)

    real(RP), intent(out) :: CAPE
    real(RP), intent(out) :: CIN
    real(RP), intent(out) :: LCL
    real(RP), intent(out) :: LFC
    real(RP), intent(out) :: LNB
    real(RP), intent(out) :: DENS_p(KA)
    real(RP), intent(out) :: TEMP_p(KA)
    real(RP), intent(out) :: BUOY_p(KA)
    real(RP), intent(out) :: QV_p  (KA)
    logical,  intent(out) :: converged

    real(RP) :: BUOY_pf(KA)
    integer  :: kLCL, kLFC, kLNB

    integer  :: k
    !---------------------------------------------------------------------------

    ! lift parcel
    call ATMOS_ADIABAT_liftparcel_1D( KA, KS, KE, &
                                      Kstr,                          & ! [IN]
                                      dens(:), temp(:), pres(:),     & ! [IN]
                                      qv(:), qc(:),                  & ! [IN]
                                      qdry(:), Rtot(:), CPtot(:),    & ! [IN]
                                      dens_p(:), temp_p(:), qv_p(:), & ! [OUT]
                                      kLCL, converged                ) ! [OUT]

    if ( .not. converged ) return

    ! parcel buoyancy
    do k = KS, KE
       BUOY_p(k) = GRAV * ( DENS(k) - DENS_p(k) ) / DENS_p(k)
    end do

    do k = KS, KE-1
       BUOY_pf(k) = 0.5_RP * ( BUOY_p(k+1) + BUOY_p(k) )
    end do
    BUOY_p(KE) = BUOY_p(KE-1)

    LCL  = 0.0_RP
    LFC  = 0.0_RP
    LNB  = 0.0_RP
    CAPE = 0.0_RP
    CIN  = 0.0_RP

    kLFC = -1
    kLNB = -1

    do k = kLCL, KE
       if ( BUOY_p(k) >= 0.0_RP ) then
          kLFC = k
          exit
       endif
    enddo

    do k = KE, Kstr, -1
       if ( BUOY_p(k) >= 0.0_RP ) then
          kLNB = k
          exit
       endif
    enddo

    if ( kLCL >= Kstr ) then
       LCL = CZ(kLCL)
    endif

    if ( kLFC >= Kstr ) then
       LFC = CZ(kLFC)
    endif

    if ( kLNB >= Kstr ) then
       LNB = CZ(kLNB)
    endif

    if ( kLFC >= Kstr .AND. kLNB > Kstr ) then
       do k = kLFC, kLNB
          CAPE = CAPE + BUOY_pf(k-1) * ( FZ(k)-FZ(k-1) )
       enddo
    endif

    if ( kLFC >= Kstr ) then
       do k = Kstr+1, kLFC
          CIN = CIN + BUOY_pf(k-1) * ( FZ(k)-FZ(k-1) )
       enddo
    endif

    return
  end subroutine ATMOS_ADIABAT_cape_1D

  subroutine ATMOS_ADIABAT_cape_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Kstr, &
       DENS, TEMP, PRES,            &
       QV, QC, Qdry, Rtot, CPtot,   &
       CZ, FZ,                      &
       CAPE, CIN, LCL, LFC, LNB,    &
       DENS_p, TEMP_p, BUOY_p, QV_p )
    use scale_prc, only: &
       PRC_abort
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    integer,  intent(in)  :: Kstr
    real(RP), intent(in)  :: DENS (KA,IA,JA)
    real(RP), intent(in)  :: TEMP (KA,IA,JA)
    real(RP), intent(in)  :: PRES (KA,IA,JA)
    real(RP), intent(in)  :: QV   (KA,IA,JA)
    real(RP), intent(in)  :: QC   (KA,IA,JA)
    real(RP), intent(in)  :: Qdry (KA,IA,JA)
    real(RP), intent(in)  :: Rtot (KA,IA,JA)
    real(RP), intent(in)  :: CPtot(KA,IA,JA)
    real(RP), intent(in)  :: CZ   (  KA,IA,JA)
    real(RP), intent(in)  :: FZ   (0:KA,IA,JA)

    real(RP), intent(out) :: CAPE(IA,JA)
    real(RP), intent(out) :: CIN (IA,JA)
    real(RP), intent(out) :: LCL (IA,JA)
    real(RP), intent(out) :: LFC (IA,JA)
    real(RP), intent(out) :: LNB (IA,JA)

    real(RP), intent(out), optional, target :: DENS_p(KA,IA,JA)
    real(RP), intent(out), optional, target :: TEMP_p(KA,IA,JA)
    real(RP), intent(out), optional, target :: BUOY_p(KA,IA,JA)
    real(RP), intent(out), optional, target :: QV_p  (KA,IA,JA)

    logical :: converged, error

    real(RP), pointer :: P_DENS(:,:,:), P_TEMP(:,:,:), P_BUOY(:,:,:), P_QV(:,:,:)

    integer :: i, j
    !---------------------------------------------------------------------------

    error = .false.

    if ( present(DENS_p) ) then
       P_DENS => DENS_p
    else
       allocate( P_DENS(KA,IA,JA) )
    end if
    if ( present(TEMP_p) ) then
       P_TEMP => TEMP_p
    else
       allocate( P_TEMP(KA,IA,JA) )
    end if
    if ( present(BUOY_p) ) then
       P_BUOY => BUOY_p
    else
       allocate( P_BUOY(KA,IA,JA) )
    end if
    if ( present(QV_p) ) then
       P_QV => QV_p
    else
       allocate( P_QV(KA,IA,JA) )
    end if


    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(converged)
    do j = JS, JE
    do i = IS, IE

       call ATMOS_ADIABAT_cape_1D( KA, KS, KE, &
                                   Kstr, &
                                   DENS(:,i,j), TEMP(:,i,j), PRES(:,i,j), & ! [IN]
                                   QV(:,i,j), QC(:,i,j), Qdry(:,i,j),     & ! [IN]
                                   Rtot(:,i,j), CPtot(:,i,j),             & ! [IN]
                                   CZ(:,i,j), FZ(:,i,j),                  & ! [IN]
                                   CAPE(i,j), CIN(i,j),                   & ! [OUT]
                                   LCL(i,j), LFC(i,j), LNB(i,j),          & ! [OUT]
                                   P_DENS(:,i,j), P_TEMP(:,i,j),          & ! [OUT]
                                   P_BUOY(:,i,j), P_QV(:,i,j),            & ! [OUT]
                                   converged                              ) ! [OUT]

       if ( .not. converged ) then
          LOG_ERROR("ATMOS_ADIABAT_cape_3D",*) '[liftparcel] not converged! ', i, j
          error = .true.
       end if

    end do
    end do

    if ( error ) call PRC_abort

    call FILE_HISTORY_in( P_DENS(:,:,:), 'DENS_parcel', 'density profile in lifting parcel',     'kg/m3' )
    call FILE_HISTORY_in( P_TEMP(:,:,:), 'TEMP_parcel', 'temperature profile in lifting parcel', 'K'     )
    call FILE_HISTORY_in( P_BUOY(:,:,:), 'BUOY_parcel', 'buoyancy profile in lifting parcel',    'm/s2'  )
    call FILE_HISTORY_in( P_QV  (:,:,:), 'QV_parcel',   'humidity profile in lifting parcel',    'kg/kg' )

    if ( .not. present(DENS_p) ) deallocate( P_DENS )
    if ( .not. present(TEMP_p) ) deallocate( P_TEMP )
    if ( .not. present(BUOY_p) ) deallocate( P_BUOY )
    if ( .not. present(QV_p  ) ) deallocate( P_QV   )

    return
  end subroutine ATMOS_ADIABAT_cape_3D

  !-----------------------------------------------------------------------------
  !> Calc temperature profile with lifting parcel
  !> Method: Pseudo-adiabatic ascend from lowermost layer of the model
  !> Reference: Emanuel(1994)
!OCL SERIAL
  subroutine ATMOS_ADIABAT_liftparcel_1D( &
       KA, KS, KE, &
       Kstr,                       &
       DENS, TEMP, PRES, QV, QC,   &
       QDRY, Rtot, CPtot,          &
       DENS_p1D, TEMP_p1D, QV_p1D, &
       kLCL, converged             )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_entr => ATMOS_HYDROMETEOR_entr, &
       CP_WATER
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_moist_conversion_pres_liq
    implicit none
    integer, intent(in) :: KA, KS, KE

    integer,  intent(in) :: Kstr
    real(RP), intent(in) :: DENS (KA)
    real(RP), intent(in) :: TEMP (KA)
    real(RP), intent(in) :: PRES (KA)
    real(RP), intent(in) :: QV   (KA)
    real(RP), intent(in) :: QC   (KA)
    real(RP), intent(in) :: QDRY (KA)
    real(RP), intent(in) :: Rtot (KA)
    real(RP), intent(in) :: CPtot(KA)

    real(RP), intent(out) :: DENS_p1D(KA)
    real(RP), intent(out) :: TEMP_p1D(KA)
    real(RP), intent(out) :: QV_p1D  (KA)
    integer,  intent(out) :: kLCL
    logical,  intent(out) :: converged

    real(RP) :: ENTR_p
    real(RP) :: QV_p, QC_p, Qdry_p
    real(RP) :: Rtot_p, CPtot_p
    real(RP) :: TEMP_p

    real(RP) :: fact

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, Kstr
       DENS_p1D(k) = DENS(k)
       TEMP_p1D(k) = TEMP(k)
       QV_p1D  (k) = QV  (k)
    enddo

    ! vapor + liquid water at the start point
    QV_p    = QV   (Kstr)
    QC_p    = QC   (Kstr)
    Qdry_p  = Qdry (Kstr)
    Rtot_p  = Rtot (Kstr)
    CPtot_p = CPtot(Kstr)

    call HYDROMETEOR_entr( TEMP(Kstr), PRES(Kstr), & ! [IN]
                           QV_p, QC_p, Qdry_p,     & ! [IN]
                           Rtot_p, CPtot_p,        & ! [IN]
                           ENTR_p                  ) ! [OUT]

    kLCL = -1
    do k = Kstr, KE

       call ATMOS_SATURATION_moist_conversion_pres_liq( &
            PRES(k), Entr_p, Qdry_p,     & ! [IN]
            QV_p, QC_p, Rtot_p, CPtot_p, & ! [INOUT]
            TEMP_p,                      & ! [OUT]
            converged                    ) ! [OUT]

       if ( .NOT. converged ) then
          exit
       endif

       if ( QC_p > 0.0_RP .and. kLCL == -1 ) kLCL = k

       ! remove condensed water
       fact = 1.0_RP / ( 1.0_RP - QC_p )
       CPtot_p = ( CPtot_p - CP_WATER * QC_p ) * fact
       Rtot_p = Rtot_p * fact
       Qdry_p = Qdry_p * fact
       QV_p = QV_p * fact
       ! Entr_p = Entr_p - QC_p * CP_WATER * log( TEMP_p / TEM00 )
       QC_p = 0.0_RP

       DENS_p1D(k) = PRES(k) / ( Rtot_p * TEMP_p )
       TEMP_p1D(k) = TEMP_p
       QV_p1D  (k) = QV_p

    enddo

    return
  end subroutine ATMOS_ADIABAT_liftparcel_1D

  subroutine ATMOS_ADIABAT_liftparcel_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Kstr,                       &
       DENS, TEMP, PRES, QV, QC,   &
       QDRY, Rtot, CPtot,          &
       DENS_p3D, TEMP_p3D, QV_p3D, &
       kLCL, converged             )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_entr => ATMOS_HYDROMETEOR_entr
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_moist_conversion_pres_liq
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    integer,  intent(in) :: Kstr
    real(RP), intent(in) :: DENS (KA,IA,JA)
    real(RP), intent(in) :: TEMP (KA,IA,JA)
    real(RP), intent(in) :: PRES (KA,IA,JA)
    real(RP), intent(in) :: QV   (KA,IA,JA)
    real(RP), intent(in) :: QC   (KA,IA,JA)
    real(RP), intent(in) :: QDRY (KA,IA,JA)
    real(RP), intent(in) :: Rtot (KA,IA,JA)
    real(RP), intent(in) :: CPtot(KA,IA,JA)

    real(RP), intent(out) :: DENS_p3D(KA,IA,JA)
    real(RP), intent(out) :: TEMP_p3D(KA,IA,JA)
    real(RP), intent(out) :: QV_p3D  (KA,IA,JA)
    integer,  intent(out) :: kLCL(IA,JA)
    logical,  intent(out) :: converged

    logical :: error

    integer :: i, j
    !---------------------------------------------------------------------------

    error = .false.

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       call ATMOS_ADIABAT_liftparcel_1D( KA, KS, KE, &
                                         Kstr,                                            & ! [IN]
                                         DENS(:,i,j), TEMP(:,i,j), PRES(:,i,j),           & ! [IN]
                                         QV(:,i,j), QC(:,i,j),                            & ! [IN]
                                         QDRY(:,i,j), Rtot(:,i,j), CPtot(:,i,j),          & ! [IN]
                                         DENS_p3D(:,i,j), TEMP_p3D(:,i,j), QV_p3D(:,i,j), & ! [OUT]
                                         kLCL(i,j), converged                             ) ! [OUT]
       if ( .not. converged ) then
          LOG_ERROR("ATMOS_ADIABAT_liftparcel_3D",*) '[liftparcel] not converged! ', i, j
          error = .true.
       end if

    enddo
    enddo

    if ( error ) call PRC_abort

    return
  end subroutine ATMOS_ADIABAT_liftparcel_3D

end module scale_atmos_adiabat
