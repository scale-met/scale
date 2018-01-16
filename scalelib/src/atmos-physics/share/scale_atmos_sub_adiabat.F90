!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Adiabatic process
!!
!! @par Description
!!          Moist adiabatic process for calculation of CAPE, CIN
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_adiabat
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
  public :: ATMOS_ADIABAT_cape
  public :: ATMOS_ADIABAT_liftparcel

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
  !> Calc CAPE and CIN
  !> Type of parcel method: Pseudo-adiabatic ascend from lowermost layer of the model
  !> Reference: Emanuel(1994)
  subroutine ATMOS_ADIABAT_cape( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Kstr, &
       DENS, TEMP, PRES,          &
       QV, QC, Qdry, Rtot, CPtot, &
       CZ, FZ,                    &
       CAPE, CIN, LCL, LFC, LNB   )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    use scale_history, only: &
       HIST_in
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

    real(RP) :: DENS_p (KA,IA,JA)
    real(RP) :: TEMP_p (KA,IA,JA)
    real(RP) :: QV_p   (KA,IA,JA)
    real(RP) :: QL_p   (KA,IA,JA)
    real(RP) :: QI_p   (KA,IA,JA)
    real(RP) :: QC_p   (KA,IA,JA)
    real(RP) :: BUOY_p (KA,IA,JA)
    real(RP) :: BUOY_pf(KA,IA,JA)
    real(RP) :: QSAT_p (KA,IA,JA)

    integer :: kLCL
    integer :: kLFC
    integer :: kLNB

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! lift parcel
    call ATMOS_ADIABAT_liftparcel( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         Kstr,                                     & ! [IN]
         dens(:,:,:), temp(:,:,:), pres(:,:,:),    & ! [IN]
         qv(:,:,:), qc(:,:,:),                     & ! [IN]
         qdry(:,:,:), Rtot(:,:,:), CPtot(:,:,:),   & ! [IN]
         dens_p(:,:,:), temp_p(:,:,:), qv_p(:,:,:) ) ! [OUT]

    ! parcel buoyancy
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          BUOY_p(k,i,j) = GRAV * ( DENS(k,i,j) - DENS_p(k,i,j) ) / DENS_p(k,i,j)
       enddo
       BUOY_p(   1:KS-1,i,j) = 0.0_RP
       BUOY_p(KE+1:KA  ,i,j) = 0.0_RP

       do k = KS, KE-1
          BUOY_pf(k,i,j) = 0.5_RP * ( BUOY_p(k+1,i,j) + BUOY_p(k,i,j) )
       enddo
       BUOY_p(KE,i,j) = BUOY_p(KE-1,i,j)

       BUOY_p(   1:KS-1,i,j) = 0.0_RP
       BUOY_p(KE+1:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    ! saturation point profile (lifted parcel)
    call SATURATION_dens2qsat_liq( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                   TEMP_p(:,:,:), DENS_p(:,:,:), & ! [IN]
                                   QSAT_p(:,:,:)                 ) ! [OUT]

    call HIST_in( DENS_p(:,:,:), 'DENS_parcel', 'density profile in lifting parcel',     'kg/m3' )
    call HIST_in( TEMP_p(:,:,:), 'TEMP_parcel', 'temperature profile in lifting parcel', 'K'     )
    call HIST_in( BUOY_p(:,:,:), 'BUOY_parcel', 'buoyancy profile in lifting parcel',    'm/s2'  )
    call HIST_in( QSAT_p(:,:,:), 'QSAT_parcel', 'saturation profile in lifting parcel',  'kg/kg' )
    call HIST_in( QV_p  (:,:,:), 'QV_parcel',   'humidity profile in lifting parcel',    'kg/kg' )

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(kLFC,kLCL,kLNB,k)
    do j = JS, JE
    do i = IS, IE

       LCL (i,j) = 0.0_RP
       LFC (i,j) = 0.0_RP
       LNB (i,j) = 0.0_RP
       CAPE(i,j) = 0.0_RP
       CIN (i,j) = 0.0_RP

       kLFC = -1
       kLCL = -1
       kLNB = -1

       do k = Kstr, KE
          if ( QV_p(Kstr,i,j) >= QSAT_p(k,i,j) ) then
             kLCL = k
             exit
          endif
       enddo

       do k = kLCL, KE
          if ( BUOY_p(k,i,j) >= 0.0_RP ) then
             kLFC = k
             exit
          endif
       enddo

       do k = KE, Kstr, -1
          if ( BUOY_p(k,i,j) >= 0.0_RP ) then
             kLNB = k
             exit
          endif
       enddo

       if ( kLCL >= Kstr ) then
          LCL(i,j) = CZ(kLCL,i,j)
       endif

       if ( kLFC >= Kstr ) then
          LFC(i,j) = CZ(kLFC,i,j)
       endif

       if ( kLNB >= Kstr ) then
          LNB(i,j) = CZ(kLNB,i,j)
       endif

       if ( kLFC >= Kstr .AND. kLNB > Kstr ) then
          do k = kLFC, kLNB
             CAPE(i,j) = CAPE(i,j) + BUOY_pf(k-1,i,j) * ( FZ(k,i,j)-FZ(k-1,i,j) )
          enddo
       endif

       if ( kLFC >= Kstr ) then
          do k = Kstr+1, kLFC
             CIN (i,j) = CIN (i,j) + BUOY_pf(k-1,i,j) * ( FZ(k,i,j)-FZ(k-1,i,j) )
          enddo
       endif

    enddo
    enddo

    return
  end subroutine ATMOS_ADIABAT_cape

  !-----------------------------------------------------------------------------
  !> Calc temperature profile with lifting parcel
  !> Method: Pseudo-adiabatic ascend from lowermost layer of the model
  !> Reference: Emanuel(1994)
  subroutine ATMOS_ADIABAT_liftparcel( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Kstr,                      &
       DENS, TEMP, PRES, QV, QC,  &
       QDRY, Rtot, CPtot,         &
       DENS_p3D, TEMP_p3D, QV_p3D )
    use scale_process, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_entr => ATMOS_HYDROMETEOR_entr, &
       CP_WATER
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

    real(RP) :: ENTR_p
    real(RP) :: QV_p, QC_p, Qdry_p
    real(RP) :: Rtot_p, CPtot_p
    real(RP) :: TEMP_p

    logical  :: converged, error
    real(RP) :: fact

    integer :: k, i, j
    !---------------------------------------------------------------------------

    error = .false.

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(ENTR_p,QV_p,QC_p,Qdry_p,Rtot_p,CPtot_p,TEMP_p)
    do j = JS, JE
    do i = IS, IE

       do k = 1, Kstr
          DENS_p3D(k,i,j) = DENS(k,i,j)
          TEMP_p3D(k,i,j) = TEMP(k,i,j)
          QV_p3D  (k,i,j) = QV(k,i,j)
       enddo

       ! vapor + liquid water at the start point
       QV_p    = QV(Kstr,i,j)
       QC_p    = QC(Kstr,i,j)
       Qdry_p  = Qdry (Kstr,i,j)
       Rtot_p  = Rtot (Kstr,i,j)
       CPtot_p = CPtot(Kstr,i,j)

       call HYDROMETEOR_entr( &
            TEMP(Kstr,i,j), PRES(Kstr,i,j), & ! [IN]
            QV_p, QC_p, Qdry_p,             & ! [IN]
            Rtot_p, CPtot_p,                & ! [IN]
            ENTR_p                          ) ! [OUT]

       do k = Kstr, KE

          call ATMOS_SATURATION_moist_conversion_pres_liq( &
               PRES(k,i,j), Entr_p, Qdry_p, & ! [IN]
               QV_p, QC_p, Rtot_p, CPtot_p, & ! [INOUT]
               TEMP_p,                      & ! [OUT]
               converged                    ) ! [OUT]

          if ( .NOT. converged ) then
             error = .true.
             write(*,*) 'xxx [moist_conversion] not converged!', k,i,j
             exit
          endif

          ! remove condensed water
          fact = 1.0_RP / ( 1.0_RP - QC_p )
          CPtot_p = ( CPtot_p - CP_WATER * QC_p ) * fact
          Rtot_p = Rtot_p * fact
          Qdry_p = Qdry_p * fact
          QV_p = QV_p * fact
          ! Entr_p = Entr_p - QC_p * CP_WATER * log( TEMP_p / TEM00 )
          QC_p = 0.0_RP

          DENS_p3D(k,i,j) = PRES(k,i,j) / ( Rtot_p * TEMP_p )
          TEMP_p3D(k,i,j) = TEMP_p
          QV_p3D  (k,i,j) = QV_p

       enddo

    enddo
    enddo

    if ( error ) call PRC_abort

    return
  end subroutine ATMOS_ADIABAT_liftparcel

end module scale_atmos_adiabat
