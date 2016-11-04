!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Adiabatic process
!!
!! @par Description
!!          Moist adiabatic process for calculation of CAPE, CIN
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
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
       Kstr, &
       DENS, &
       TEMP, &
       PRES, &
       QTRC, &
       CZ,   &
       FZ,   &
       CAPE, &
       CIN,  &
       LCL,  &
       LFC,  &
       LNB   )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       HYDROMETEOR_entr => ATMOS_HYDROMETEOR_entr
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    use scale_history, only: &
       HIST_in
    implicit none

    integer,  intent(in)  :: Kstr
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: TEMP(KA,IA,JA)
    real(RP), intent(in)  :: PRES(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ  (  KA,IA,JA)
    real(RP), intent(in)  :: FZ  (0:KA,IA,JA)
    real(RP), intent(out) :: CAPE(IA,JA)
    real(RP), intent(out) :: CIN (IA,JA)
    real(RP), intent(out) :: LCL (IA,JA)
    real(RP), intent(out) :: LFC (IA,JA)
    real(RP), intent(out) :: LNB (IA,JA)

    real(RP) :: DENS_p (KA,IA,JA)
    real(RP) :: TEMP_p (KA,IA,JA)
    real(RP) :: QTRC_p (KA,IA,JA,QA)
    real(RP) :: ENTR_p (KA,IA,JA)
    real(RP) :: BUOY_p (KA,IA,JA)
    real(RP) :: BUOY_pf(KA,IA,JA)
    real(RP) :: QSAT_p (KA,IA,JA)

    integer :: kLCL(IA,JA)
    integer :: kLFC(IA,JA)
    integer :: kLNB(IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! entropy at start point
    call HYDROMETEOR_entr( ENTR_p(Kstr,:,:),   & ! [OUT]
                          TEMP  (Kstr,:,:),   & ! [IN]
                          PRES  (Kstr,:,:),   & ! [IN]
                          QTRC  (Kstr,:,:,:), & ! [IN]
                          TRACER_R(:)         ) ! [IN]

    ! lift parcel
    call ATMOS_ADIABAT_liftparcel( Kstr,             & ! [IN]
                                   TEMP  (:,:,:),    & ! [IN]
                                   PRES  (:,:,:),    & ! [IN]
                                   QTRC  (:,:,:,:),  & ! [IN]
                                   ENTR_p(Kstr,:,:), & ! [IN]
                                   DENS_p(:,:,:),    & ! [OUT]
                                   TEMP_p(:,:,:),    & ! [OUT]
                                   QTRC_p(:,:,:,:)   ) ! [OUT]

    ! entropy profile (lifted parcel)
    call HYDROMETEOR_entr( ENTR_p(:,:,:),   & ! [OUT]
                          TEMP_p(:,:,:),   & ! [IN]
                          PRES  (:,:,:),   & ! [IN]
                          QTRC_p(:,:,:,:), & ! [IN]
                          TRACER_R(:)      ) ! [IN]

    ! parcel buoyancy
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
    call SATURATION_dens2qsat_liq( QSAT_p(:,:,:), & ! [OUT]
                                   TEMP_p(:,:,:), & ! [IN]
                                   DENS_p(:,:,:)  ) ! [IN]

    call HIST_in( DENS_p(:,:,:), 'DENS_parcel', 'density profile in lifting parcel',     'kg/m3' )
    call HIST_in( TEMP_p(:,:,:), 'TEMP_parcel', 'temperature profile in lifting parcel', 'K'     )
    call HIST_in( ENTR_p(:,:,:), 'ENTR_parcel', 'entropy profile in lifting parcel',     'J/K'   )
    call HIST_in( BUOY_p(:,:,:), 'BUOY_parcel', 'buoyancy profile in lifting parcel',    'm/s2'  )
    call HIST_in( QSAT_p(:,:,:), 'QSAT_parcel', 'saturation profile in lifting parcel',  'kg/kg' )

    ! detect layer number of LNB, LFC, LCL
    kLCL(:,:) = -1
    do j = JS, JE
    do i = IS, IE
       do k = Kstr, KE
          if ( QTRC_p(Kstr,i,j,I_QV) >= QSAT_p(k,i,j) ) then
             kLCL(i,j) = k
             exit
          endif
       enddo
    enddo
    enddo

    kLFC(:,:) = -1
    do j = JS, JE
    do i = IS, IE
       do k = kLCL(i,j), KE
          if ( BUOY_p(k,i,j) >= 0.0_RP ) then
             kLFC(i,j) = k
             exit
          endif
       enddo
    enddo
    enddo

    kLNB(:,:) = -1
    do j = JS, JE
    do i = IS, IE
       do k = KE, Kstr, -1
          if ( BUOY_p(k,i,j) >= 0.0_RP ) then
             kLNB(i,j) = k
             exit
          endif
       enddo
    enddo
    enddo

    LCL(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       if ( kLCL(i,j) >= Kstr ) then
          LCL(i,j) = CZ(kLCL(i,j),i,j)
       endif
    enddo
    enddo

    LFC(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       if ( kLFC(i,j) >= Kstr ) then
          LFC(i,j) = CZ(kLFC(i,j),i,j)
       endif
    enddo
    enddo

    LNB(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       if ( kLNB(i,j) >= Kstr ) then
          LNB(i,j) = CZ(kLNB(i,j),i,j)
       endif
    enddo
    enddo

    CAPE(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       if ( kLFC(i,j) >= Kstr .AND. kLNB(i,j) > Kstr ) then
          do k = kLFC(i,j), kLNB(i,j)
             CAPE(i,j) = CAPE(i,j) + BUOY_pf(k-1,i,j) * ( FZ(k,i,j)-FZ(k-1,i,j) )
          enddo
       endif
    enddo
    enddo

    CIN(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       if ( kLFC(i,j) >= Kstr ) then
          do k = Kstr+1, kLFC(i,j)
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
       Kstr,   &
       TEMP,   &
       PRES,   &
       QTRC,   &
       ENTR_p, &
       DENS_p, &
       TEMP_p, &
       QTRC_p  )
    use scale_const, only: &
       EPS   => CONST_EPS,   &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       CPvap => CONST_CPvap, &
       LHV0  => CONST_LHV0,  &
       PSAT0 => CONST_PSAT0, &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       I_QC, &
       HYDROMETEOR_entr => ATMOS_HYDROMETEOR_entr
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    implicit none

    integer,  intent(in)  :: Kstr
    real(RP), intent(in)  :: TEMP  (KA,IA,JA)
    real(RP), intent(in)  :: PRES  (KA,IA,JA)
    real(RP), intent(in)  :: QTRC  (KA,IA,JA,QA)
    real(RP), intent(in)  :: ENTR_p(IA,JA)
    real(RP), intent(out) :: DENS_p(KA,IA,JA)
    real(RP), intent(out) :: TEMP_p(KA,IA,JA)
    real(RP), intent(out) :: QTRC_p(KA,IA,JA,QA)

    real(RP) :: qsat_p(KA,IA,JA)
    real(RP) :: qtot_p(IA,JA)

    real(RP) :: qdry_p, Rtot, CPtot
    real(RP) :: pres_dry, pres_vap
    real(RP) :: TEMP_unsat
    real(RP) :: qsat_unsat

    real(RP) :: TEMP_ite
    real(RP) :: qsat_ite
    real(RP) :: QTRC_ite(QA)
    real(RP) :: ENTR_ite
    real(RP) :: TEMP_prev
    real(RP) :: dENTR_dT

    real(RP), parameter :: criteria = 1.E-8_RP
    integer,  parameter :: itelim   = 100
    integer             :: ite

    real(RP), parameter :: TEMMIN = 0.1_RP

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = 1, Kstr
       TEMP_p(k,i,j) = TEMP(k,i,j)
       do iqw = 1, QA
          QTRC_p(k,i,j,iqw) = QTRC(k,i,j,iqw)
       enddo
    enddo
    enddo
    enddo

    ! vapor + cloud water at the start point
    do j = JS, JE
    do i = IS, IE
       qtot_p(i,j) = QTRC(Kstr,i,j,I_QV) + QTRC(Kstr,i,j,I_QC)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = Kstr, KE
       TEMP_p(k,i,j) = TEMP(k,i,j) ! first guess

       ! T1: unsaturated temperature, S = U1(PRES, TEMP_unsat, qtot_p)
       qdry_p = 1.0_RP - qtot_p(i,j)
       Rtot   = Rdry  * qdry_p + Rvap  * qtot_p(i,j)
       CPtot  = CPdry * qdry_p + CPvap * qtot_p(i,j)

       ! dry air + vapor
       pres_dry = max( PRES(k,i,j) * qdry_p      * Rdry / Rtot, EPS )
       pres_vap = max( PRES(k,i,j) * qtot_p(i,j) * Rvap / Rtot, EPS )

       TEMP_unsat = TEM00 * exp( ( ENTR_p(i,j) + qdry_p      * Rdry * log( pres_dry / PRE00 ) &
                                               + qtot_p(i,j) * Rvap * log( pres_vap / PSAT0 ) &
                                               - qtot_p(i,j) * LHV0 / TEM00                   ) / CPtot )

       call SATURATION_pres2qsat_liq( qsat_unsat, & ! [OUT]
                                      TEMP_unsat, & ! [IN]
                                      PRES(k,i,j) ) ! [IN]

       ! T2: saturated temperature, S = U2(PRES, TEMP_ite, QTRC_ite)
       if ( qtot_p(i,j) > qsat_unsat ) then

          TEMP_ite = TEM00 * exp( ( ENTR_p(i,j) + Rdry * log( PRES(k,i,j) / PRE00 ) ) / CPdry )

          do ite = 1, itelim

             call SATURATION_pres2qsat_liq( qsat_ite,   & ! [OUT]
                                            TEMP_ite,   & ! [IN]
                                            PRES(k,i,j) ) ! [IN]

             QTRC_ite(:)    = 0.0_RP ! Pseudo-adiabatic: no cloud water
             QTRC_ite(I_QV) = min( qtot_p(i,j), qsat_ite )

             call HYDROMETEOR_entr( ENTR_ite,    & ! [OUT]
                                   TEMP_ite,    & ! [IN]
                                   PRES(k,i,j), & ! [IN]
                                   QTRC_ite(:), & ! [IN]
                                   TRACER_R(:)  ) ! [IN]

             qdry_p   = 1.0_RP - QTRC_ite(I_QV)
             CPtot    = CPdry * qdry_p + CPvap * QTRC_ite(I_QV)

             dENTR_dT  = CPtot                    / TEMP_ite    &
                       - QTRC_ite(I_QV) * LHV0    / TEMP_ite**2 &
                       + QTRC_ite(I_QV) * LHV0**2 / TEMP_ite**3 / Rvap
             dENTR_dT  = max( dENTR_dT, EPS )

             TEMP_prev = TEMP_ite
             TEMP_ite  = TEMP_ite - ( ENTR_ite - ENTR_p(i,j) ) / dENTR_dT
             TEMP_ite  = max( TEMP_ite, TEMMIN )

             if( abs(TEMP_ite-TEMP_prev) < criteria ) exit

          enddo

          TEMP_p(k,i,j) = TEMP_ite

       endif

       ! parcel satulation point
       call SATURATION_pres2qsat_liq( qsat_p(k,i,j), & ! [OUT]
                                      TEMP_p(k,i,j), & ! [IN]
                                      PRES  (k,i,j)  ) ! [IN]

       ! update parcel vapor : remove condensed water
       qtot_p(i,j) = min( qtot_p(i,j), qsat_p(k,i,j) )

       QTRC_p(k,i,j,I_QV) = qtot_p(i,j)
       QTRC_p(k,i,j,I_QC) = 0.0_RP

       qdry_p = 1.0_RP - qtot_p(i,j)
       Rtot   = Rdry * qdry_p + Rvap * qtot_p(i,j)

       DENS_p(k,i,j) = PRES(k,i,j) / ( Rtot * TEMP_p(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_ADIABAT_liftparcel

end module scale_atmos_adiabat
