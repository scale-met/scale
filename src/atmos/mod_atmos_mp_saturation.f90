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
module mod_mp_saturation
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
     LH00   => CONST_LH00,   &
     LHS00  => CONST_LHS00,  &
     LH0    => CONST_LH0,    &
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
  public :: MP_SATURATION_psat_water
  public :: MP_SATURATION_psat_ice
  public :: MP_SATURATION_qsat_water
  public :: MP_SATURATION_qsat_ice

  public :: MP_SATURATION_dqsw_dtem_rho  
  public :: MP_SATURATION_dqsi_dtem_rho  
  public :: MP_SATURATION_dqsw_dtem_dpre 
  public :: MP_SATURATION_dqsi_dtem_dpre 

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
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
  subroutine MP_SATURATION_psat_water ( temp, psat )
    implicit none

    real(RP), intent(in)  :: temp(IJA,KA)
    real(RP), intent(out) :: psat(IJA,KA)

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = 1, KA
    do ij = 1, IJA
       TEM = max( temp(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovRvap            &
                  * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_psat_water

  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_psat_ice ( temp, psat )
    implicit none

    real(RP), intent(in)  :: temp(IJA,KA)
    real(RP), intent(out) :: psat(IJA,KA)

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = 1, KA
    do ij = 1, IJA
       TEM = max( temp(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_psat_ice

  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_qsat_water( temp, pres, qsat )
    implicit none

    real(RP), intent(in)  :: temp(IJA,KA)
    real(RP), intent(in)  :: pres(IJA,KA)
    real(RP), intent(out) :: qsat(IJA,KA)
    
    real(RP) :: psat
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = 1, KA
    do ij = 1, IJA
       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )

       qsat(ij,k) = EPSvap * psat / ( pres(ij,k) - ( 1.E0_RP-EPSvap ) * psat )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_qsat_water

  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_qsat_ice ( temp, pres, qsat )
    implicit none

    real(RP), intent(in)  :: temp(IJA,KA)
    real(RP), intent(in)  :: pres(IJA,KA)
    real(RP), intent(out) :: qsat(IJA,KA)
    
    real(RP) :: psat
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = 1, KA
    do ij = 1, IJA

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )

       qsat(ij,k) = EPSvap * psat / ( pres(ij,k) - ( 1.E0_RP-EPSvap ) * psat )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_qsat_ice

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_dqsw_dtem_rho( temp, dens, dqsdtem )
    implicit none

    real(RP), intent(in)  :: temp   (IJA,KA)
    real(RP), intent(in)  :: dens   (IJA,KA)
    real(RP), intent(out) :: dqsdtem(IJA,KA)

    real(RP) :: psat(IJA) ! saturation vapor pressure
    real(RP) :: lhv (IJA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = KS, KE
!    do k  = 1, KA

       do ij = 1, IJA
          TEM = max( temp(ij,k), TEM_MIN )

          psat(ij) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )
       enddo

       do ij = 1, IJA
          lhv(ij)  = LH0 + ( CPvap-CL ) * ( temp(ij,k)-TEM00 )

          dqsdtem(ij,k) = psat(ij) / ( dens(ij,k) * Rvap * temp(ij,k) * temp(ij,k) ) &
                      * ( lhv(ij) / ( Rvap * temp(ij,k) ) - 1.E0_RP )
       enddo

    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_dqsi_dtem_rho( temp, dens, dqsdtem )
    implicit none

    real(RP), intent(in)  :: temp   (IJA,KA)
    real(RP), intent(in)  :: dens   (IJA,KA)
    real(RP), intent(out) :: dqsdtem(IJA,KA)

    real(RP) :: psat(IJA) ! saturation vapor pressure
    real(RP) :: lhv (IJA) ! latent heat for condensation

    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

!    do k  = 1, KA
    do k  = KS, KE

       do ij = 1, IJA
          TEM = max( temp(ij,k), TEM_MIN )

          psat(ij) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )
       enddo

       do ij = 1, IJA
          lhv(ij)  = LHS0 + ( CPvap-CI ) * ( temp(ij,k)-TEM00 )

          dqsdtem(ij,k) = psat(ij) / ( dens(ij,k) * Rvap * temp(ij,k) * temp(ij,k) ) &
                      * ( lhv(ij) / ( Rvap * temp(ij,k) ) - 1.E0_RP )
       enddo

    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_dqsw_dtem_dpre( temp, pres, dqsdtem, dqsdpre )
    implicit none

    real(RP), intent(in)  :: temp   (IJA,KA)
    real(RP), intent(in)  :: pres   (IJA,KA)
    real(RP), intent(out) :: dqsdtem(IJA,KA)
    real(RP), intent(out) :: dqsdpre(IJA,KA)

    real(RP) :: psat(IJA) ! saturation vapor pressure
    real(RP) :: lhv (IJA) ! latent heat for condensation

    real(RP) :: den1(IJA), den2(IJA) ! denominator
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = 1, KA

       do ij = 1, IJA
          TEM = max( temp(ij,k), TEM_MIN )

          psat(ij) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )
       enddo

       do ij = 1, IJA
          den1(ij) = ( pres(ij,k) - (1.E0_RP-EPSvap) * psat(ij) ) &
                   * ( pres(ij,k) - (1.E0_RP-EPSvap) * psat(ij) )
          den2(ij) = den1(ij) * Rvap * temp(ij,k) * temp(ij,k)
          lhv(ij)  = LH0 + ( CPvap-CL ) * ( temp(ij,k)-TEM00 )
       enddo

       do ij = 1, IJA
          dqsdpre(ij,k) = - EPSvap * psat(ij) / den1(ij)
          dqsdtem(ij,k) =   EPSvap * psat(ij) / den2(ij) * lhv(ij) * pres(ij,k)
       enddo

    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine MP_SATURATION_dqsi_dtem_dpre( temp, pres, dqsdtem, dqsdpre )
    implicit none

    real(RP), intent(in)  :: temp   (IJA,KA)
    real(RP), intent(in)  :: pres   (IJA,KA)
    real(RP), intent(out) :: dqsdtem(IJA,KA)
    real(RP), intent(out) :: dqsdpre(IJA,KA)

    real(RP) :: psat(IJA) ! saturation vapor pressure
    real(RP) :: lhv (IJA) ! latent heat for condensation

    real(RP) :: den1(IJA), den2(IJA) ! denominator
    real(RP) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("satadjust")
#endif

    RTEM00   = 1.E0_RP / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = 1, KA

       do ij = 1, IJA
          TEM = max( temp(ij,k), TEM_MIN )

          psat(ij) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.E0_RP/TEM ) )
       enddo

       do ij = 1, IJA
          den1(ij) = ( pres(ij,k) - (1.E0_RP-EPSvap) * psat(ij) ) &
                   * ( pres(ij,k) - (1.E0_RP-EPSvap) * psat(ij) )
          den2(ij) = den1(ij) * Rvap * temp(ij,k) * temp(ij,k)
          lhv(ij)  = LHS0 + ( CPvap-CI ) * ( temp(ij,k)-TEM00 )
       enddo

       do ij = 1, IJA
          dqsdpre(ij,k) = - EPSvap * psat(ij) / den1(ij)
          dqsdtem(ij,k) =   EPSvap * psat(ij) / den2(ij) * lhv(ij) * pres(ij,k)
       enddo

    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("satadjust")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine MP_SATURATION_dqsi_dtem_dpre

end module mod_mp_saturation
!-------------------------------------------------------------------------------
