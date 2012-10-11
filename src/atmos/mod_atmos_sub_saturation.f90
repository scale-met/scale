!-------------------------------------------------------------------------------
!> module Saturation Adjustment
!!
!! @par Description
!!          Saturation adjustment module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)   [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro) [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_saturation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use mod_const, only : &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     LHV00  => CONST_LH00,   &
     LHS00  => CONST_LHS00,  &
     LHV0   => CONST_LH0,    &
     LHS0   => CONST_LHS0,   &
     EPSvap => CONST_EPSvap, &
     PSAT0  => CONST_PSAT0,  &
     TEM00  => CONST_TEM00
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_SATURATION_psat_water
  public :: ATMOS_SATURATION_psat_ice
  public :: ATMOS_SATURATION_qsat_sfc
  public :: ATMOS_SATURATION_qsat_water
  public :: ATMOS_SATURATION_qsat_ice

  public :: ATMOS_SATURATION_dqsw_dtem_rho  
  public :: ATMOS_SATURATION_dqsi_dtem_rho  
  public :: ATMOS_SATURATION_dqsw_dtem_dpre 
  public :: ATMOS_SATURATION_dqsi_dtem_dpre 

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  real(RP), private, parameter :: TEM_MIN   = 10.E0_RP

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  ! psat : Clasius-Clapeyron: based on CPV, CPL constant
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_psat_water( psat, temp )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LHV00 / Rvap

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat(k,i,j) = PSAT0                                 &
                   * ( TEM * RTEM00 )**CPovRvap            &
                   * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_psat_water

  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_psat_ice( psat, temp )
    implicit none

    real(RP), intent(out) :: psat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat(k,i,j) = PSAT0                                 &
                   * ( TEM * RTEM00 )**CPovRvap            &
                   * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_psat_ice

  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_qsat_sfc( qsat, temp, pres )
    implicit none

    real(RP), intent(out) :: qsat(1,IA,JA)
    real(RP), intent(in)  :: temp(1,IA,JA)
    real(RP), intent(in)  :: pres(1,IA,JA)
    
    real(RP) :: psat
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LHV00 / Rvap

    k = 1
    do j = JS, JE
    do i = IS, IE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0                                 &
            * ( TEM * RTEM00 )**CPovRvap            &
            * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_qsat_sfc

  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_qsat_water( qsat, temp, pres )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: pres(KA,IA,JA)
    
    real(RP) :: psat
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LHV00 / Rvap

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0                                 &
            * ( TEM * RTEM00 )**CPovRvap            &
            * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_qsat_water

  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_qsat_ice( qsat, temp, pres )
    implicit none

    real(RP), intent(out) :: qsat(KA,IA,JA)
    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(in)  :: pres(KA,IA,JA)
    
    real(RP) :: psat
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEM = max( temp(k,i,j), TEM_MIN )

       psat = PSAT0                                 &
            * ( TEM * RTEM00 )**CPovRvap            &
            * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )

       qsat(k,i,j) = EPSvap * psat / ( pres(k,i,j) - ( 1.0_RP-EPSvap ) * psat )
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_qsat_ice

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsw_dtem_rho( dqsdtem, temp, dens )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LHV00 / Rvap

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                  &
                   * ( TEM * RTEM00 )**CPovRvap            &
                   * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          lhv(k)  = LHV0 + ( CPvap-CL ) * ( temp(k,i,j)-TEM00 )

          dqsdtem(k,i,j) = psat(k) / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                         * ( lhv(k) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsi_dtem_rho( dqsdtem, temp, dens )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: dens   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovRvap            &
                  * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          lhv(k) = LHS0 + ( CPvap-CI ) * ( temp(k,i,j)-TEM00 )

          dqsdtem(k,i,j) = psat(k) / ( dens(k,i,j) * Rvap * temp(k,i,j) * temp(k,i,j) ) &
                         * ( lhv(k) / ( Rvap * temp(k,i,j) ) - 1.0_RP )
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsw_dtem_dpre( dqsdtem, dqsdpre, temp, pres )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: den1(KA), den2(KA) ! denominator
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LHV00 / Rvap

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovRvap            &
                  * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          den1(k) = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) ) &
                  * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) )
          den2(k) = den1(k) * Rvap * temp(k,i,j) * temp(k,i,j)
          lhv (k) = LHV0 + ( CPvap-CL ) * ( temp(k,i,j)-TEM00 )
       enddo

       do k = KS, KE
          dqsdpre(k,i,j) = - EPSvap * psat(k) / den1(k)
          dqsdtem(k,i,j) =   EPSvap * psat(k) / den2(k) * lhv(k) * pres(k,i,j)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SATURATION_dqsi_dtem_dpre( dqsdtem, dqsdpre, temp, pres )
    implicit none

    real(RP), intent(out) :: dqsdtem(KA,IA,JA)
    real(RP), intent(out) :: dqsdpre(KA,IA,JA)
    real(RP), intent(in)  :: temp   (KA,IA,JA)
    real(RP), intent(in)  :: pres   (KA,IA,JA)

    real(RP) :: psat(KA) ! saturation vapor pressure
    real(RP) :: lhv (KA) ! latent heat for condensation

    real(RP) :: den1(KA), den2(KA) ! denominator
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_satadjust')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_satadjust')
#endif

    RTEM00   = 1.0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          TEM = max( temp(k,i,j), TEM_MIN )

          psat(k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovRvap            &
                  * exp( LHovRvap * ( RTEM00 - 1.0_RP/TEM ) )
       enddo

       do k = KS, KE
          den1(k) = ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) ) &
                  * ( pres(k,i,j) - (1.0_RP-EPSvap) * psat(k) )
          den2(k) = den1(k) * Rvap * temp(k,i,j) * temp(k,i,j)
          lhv (k) = LHS0 + ( CPvap-CI ) * ( temp(k,i,j)-TEM00 )
       enddo

       do k = KS, KE
          dqsdpre(k,i,j) = - EPSvap * psat(k) / den1(k)
          dqsdtem(k,i,j) =   EPSvap * psat(k) / den2(k) * lhv(k) * pres(k,i,j)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_satadjust')
#endif
    call TIME_rapend  ('SUB_satadjust')

    return
  end subroutine ATMOS_SATURATION_dqsi_dtem_dpre

end module mod_atmos_saturation
!-------------------------------------------------------------------------------
