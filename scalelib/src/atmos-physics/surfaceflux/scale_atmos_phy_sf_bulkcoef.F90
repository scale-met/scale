!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface bulk coefficient
!!
!! @par Description
!!          Calculation of Bulk coefficient at the surface
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-18 (T.Yamaura)   [new]
!! @li      2014-05-20 (H.Yashiro)   [mod] move from coupler to atmos_phy_sf
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_sf_bulkcoef
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
  public :: ATMOS_PHY_SF_bulkcoef_setup

  abstract interface
     subroutine bc( &
          ATM_Uabs, ATM_POTT, Z1, &
          SFC_TEMP,               &
          Z0M, Z0H, Z0E,          &
          Cm, Ch, Ce ,            &
          R10M, R2H, R2E          )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(in)  :: ATM_Uabs(IA,JA) ! absolute velocity at Z1       [m/s]
       real(RP), intent(in)  :: ATM_POTT(IA,JA) ! potential temperature at z1, based on the local surface pressure [K]
       real(RP), intent(in)  :: Z1      (IA,JA) ! height of lowermost atmosphere grid (cell center) [m]
       real(RP), intent(in)  :: SFC_TEMP(IA,JA) ! potential temperature at surface skin [K]
       real(RP), intent(in)  :: Z0M     (IA,JA) ! roughness length for momentum [m]
       real(RP), intent(in)  :: Z0H     (IA,JA) ! roughness length for heat     [m]
       real(RP), intent(in)  :: Z0E     (IA,JA) ! roughness length for moisture [m]
       real(RP), intent(out) :: Cm      (IA,JA) ! bulk coefficient for momentum
       real(RP), intent(out) :: Ch      (IA,JA) ! bulk coefficient for heat
       real(RP), intent(out) :: Ce      (IA,JA) ! bulk coefficient for moisture
       real(RP), intent(out) :: R10M    (IA,JA)
       real(RP), intent(out) :: R2H     (IA,JA)
       real(RP), intent(out) :: R2E     (IA,JA)
     end subroutine bc
  end interface

  procedure(bc), pointer :: ATMOS_PHY_SF_bulkcoef => NULL()

  public :: ATMOS_PHY_SF_bulkcoef

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_PHY_SF_bulkcoef_uno
  private :: ATMOS_PHY_SF_bulkcoef_beljaars
  private :: fm_unstable
  private :: fh_unstable
  private :: fm_stable
  private :: fh_stable

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: ATMOS_PHY_SF_BULKCOEF_TYPE = 'BH91'

  real(RP), private, parameter :: ATMOS_PHY_SF_BULKCOEF_Cm_max = 2.5E-3_RP ! maximum limit of bulk coefficient for momentum [NIL]
  real(RP), private, parameter :: ATMOS_PHY_SF_BULKCOEF_Ch_max = 1.0E+0_RP ! maximum limit of bulk coefficient for heat     [NIL]
  real(RP), private, parameter :: ATMOS_PHY_SF_BULKCOEF_Ce_max = 1.0E+0_RP ! maximum limit of bulk coefficient for moisture [NIL]
  real(RP), private            :: ATMOS_PHY_SF_BULKCOEF_Cm_min = 1.0E-5_RP ! minimum limit of bulk coefficient for momentum [NIL]
  real(RP), private            :: ATMOS_PHY_SF_BULKCOEF_Ch_min = 1.0E-5_RP ! minimum limit of bulk coefficient for heat     [NIL]
  real(RP), private            :: ATMOS_PHY_SF_BULKCOEF_Ce_min = 1.0E-5_RP ! minimum limit of bulk coefficient for moisture [NIL]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_bulkcoef_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_BULKCOEF / &
       ATMOS_PHY_SF_BULKCOEF_TYPE,   &
       ATMOS_PHY_SF_BULKCOEF_Cm_min, &
       ATMOS_PHY_SF_BULKCOEF_Ch_min, &
       ATMOS_PHY_SF_BULKCOEF_Ce_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Bulk coefficient parameter'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_BULKCOEF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_BULKCOEF. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_BULKCOEF)

    select case( ATMOS_PHY_SF_BULKCOEF_TYPE )
    case ('U95')
       ATMOS_PHY_SF_bulkcoef => ATMOS_PHY_SF_bulkcoef_uno
    case ('BH91')
       ATMOS_PHY_SF_bulkcoef => ATMOS_PHY_SF_bulkcoef_beljaars
    case default
       write(*,*) 'xxx invalid bulk scheme (', trim(ATMOS_PHY_SF_BULKCOEF_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_SF_bulkcoef_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_bulkcoef_uno( &
       ATM_Uabs, ATM_POTT, Z1, &
       SFC_TEMP,               &
       Z0M, Z0H, Z0E,          &
       Cm, Ch, Ce ,            &
       R10M, R2H, R2E          )
    use scale_const, only: &
       GRAV   => CONST_GRAV,  &
       KARMAN => CONST_KARMAN
    implicit none

    real(RP), intent(in)  :: ATM_Uabs(IA,JA) ! absolute velocity at Z1       [m/s]
    real(RP), intent(in)  :: ATM_POTT(IA,JA) ! potential temperature at z1, based on the local surface pressure [K]
    real(RP), intent(in)  :: Z1      (IA,JA) ! height of lowermost atmosphere grid (cell center) [m]
    real(RP), intent(in)  :: SFC_TEMP(IA,JA) ! potential temperature at surface skin [K]
    real(RP), intent(in)  :: Z0M     (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in)  :: Z0H     (IA,JA) ! roughness length for heat     [m]
    real(RP), intent(in)  :: Z0E     (IA,JA) ! roughness length for moisture [m]
    real(RP), intent(out) :: Cm      (IA,JA) ! bulk coefficient for momentum
    real(RP), intent(out) :: Ch      (IA,JA) ! bulk coefficient for heat
    real(RP), intent(out) :: Ce      (IA,JA) ! bulk coefficient for moisture
    real(RP), intent(out) :: R10M    (IA,JA)
    real(RP), intent(out) :: R2H     (IA,JA)
    real(RP), intent(out) :: R2E     (IA,JA)

    ! specific parameters
    real(RP), parameter :: tPrn = 0.74_RP ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  =  9.4_RP ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp =  4.7_RP ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm =  7.4_RP ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh =  5.3_RP ! Louis factor d for heat (Louis 1979)

    real(RP) :: RiB (IA,JA) ! bulk Richardson number   [no unit]
    real(RP) :: RiBT(IA,JA)
    real(RP) :: C0  (IA,JA) ! initial drag coefficient [no unit]
    real(RP) :: fm  (IA,JA)
    real(RP) :: fh  (IA,JA)
    real(RP) :: Psi (IA,JA)
    real(RP) :: t0th(IA,JA)
    real(RP) :: q0qe(IA,JA)

    real(RP) :: C0_10m(IA,JA)
    real(RP) :: fm_10m(IA,JA)
    real(RP) :: C0_2m (IA,JA)
    real(RP) :: fm_2m (IA,JA)
    real(RP) :: fh_2m (IA,JA)
    real(RP) :: Psi_2m(IA,JA)

    real(RP) :: tmp1, tmp2, tmp3, Psi1
    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       C0    (i,j)  = ( KARMAN / log( Z1(i,j)/Z0M(i,j) ) )**2
       C0_10m(i,j)  = ( KARMAN / log( 10.0_RP/Z0M(i,j) ) )**2
       C0_2m (i,j)  = ( KARMAN / log(  2.0_RP/Z0M(i,j) ) )**2

       RiBT(i,j)  = GRAV * Z1(i,j) * ( ATM_POTT(i,j) - SFC_TEMP(i,j)    ) &
                                   / ( ATM_POTT(i,j) * ATM_Uabs(i,j)**2 )

       RiB (i,j)  = RiBT(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       if ( RiBT(i,j) < 0.0_RP ) then ! unstable condition
          tmp1 = C0(i,j) * LFb * sqrt( Z1(i,j)/Z0M(i,j) ) * sqrt( abs(RiB(i,j)) )

          fm(i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdm * tmp1 )
          fh(i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdh * tmp1 )
       else                        ! stable condition
          fm(i,j) = 1.0_RP / ( 1.0_RP + LFbp * RiB(i,j) )**2
          fh(i,j) = fm(i,j)
       endif
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Psi(i,j) = tPrn * sqrt(fm(i,j)) / fh(i,j) * log( Z1(i,j)/Z0M(i,j) )

       t0th(i,j) = 1.0_RP / ( 1.0_RP + tPrn * log( Z0M(i,j)/Z0H(i,j) ) / Psi(i,j) )
       q0qe(i,j) = 1.0_RP / ( 1.0_RP + tPrn * log( Z0M(i,j)/Z0E(i,j) ) / Psi(i,j) )

       RiB (i,j) = RiB(i,j) * t0th(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       if ( RiBT(i,j) < 0.0_RP ) then ! unstable condition
          tmp1 = C0    (i,j) * LFb * sqrt( Z1(i,j)/Z0M(i,j) ) * sqrt( abs(RiB(i,j)) )
          tmp2 = C0_10m(i,j) * LFb * sqrt( 10.0_RP/Z0M(i,j) ) * sqrt( abs(RiB(i,j)) )
          tmp3 = C0_2m (i,j) * LFb * sqrt(  2.0_RP/Z0M(i,j) ) * sqrt( abs(RiB(i,j)) )

          fm    (i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdm * tmp1 )
          fh    (i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdh * tmp1 )

          fm_10m(i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdm * tmp2 )
          fm_2m (i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdm * tmp3 )
          fh_2m (i,j) = 1.0_RP - LFb * RiB(i,j) / ( 1.0_RP + LFdh * tmp3 )
       else                            ! stable condition
          fm    (i,j) = 1.0_RP / ( 1.0_RP + LFbp * RiB(i,j) )**2
          fh    (i,j) = fm(i,j)

          fm_10m(i,j) = fm(i,j)
          fm_2m (i,j) = fm(i,j)
          fh_2m (i,j) = fm(i,j)
       endif
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Psi1        = tPrn * sqrt(fm   (i,j)) / fh   (i,j) * log( Z1(i,j)/Z0M(i,j) )
       Psi_2m(i,j) = tPrn * sqrt(fm_2m(i,j)) / fh_2m(i,j) * log(  2.0_RP/Z0M(i,j) )

       t0th(i,j) = 1.0_RP / ( 1.0_RP + tPrn * log( Z0M(i,j)/Z0H(i,j) ) / Psi1 )
       q0qe(i,j) = 1.0_RP / ( 1.0_RP + tPrn * log( Z0M(i,j)/Z0E(i,j) ) / Psi1 )
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Cm(i,j) = C0(i,j) * fm(i,j)
       Ch(i,j) = C0(i,j) * fh(i,j) * t0th(i,j) / tPrn
       Ce(i,j) = C0(i,j) * fh(i,j) * q0qe(i,j) / tPrn

       Cm(i,j) = min( max( Cm(i,j), ATMOS_PHY_SF_BULKCOEF_Cm_min ), ATMOS_PHY_SF_BULKCOEF_Cm_max )
       Ch(i,j) = min( max( Ch(i,j), ATMOS_PHY_SF_BULKCOEF_Ch_min ), ATMOS_PHY_SF_BULKCOEF_Ch_max )
       Ce(i,j) = min( max( Ce(i,j), ATMOS_PHY_SF_BULKCOEF_Ce_min ), ATMOS_PHY_SF_BULKCOEF_Ce_max )
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       R10M(i,j) = log( 10.0_RP/Z0M(i,j) ) / log( Z1(i,j)/Z0M(i,j) ) &
                 / sqrt( fm_10m(i,j) / fm(i,j) )

       tmp1 = log( Z0M(i,j)/Z0H(i,j) ) / log( Z1(i,j)/Z0M(i,j) )

       R2H(i,j) = ( tmp1 + Psi_2m(i,j) / tPrn ) &
                / ( tmp1 + Psi   (i,j) / tPrn )

       tmp2 = log( Z0M(i,j)/Z0E(i,j) ) / log( Z1(i,j)/Z0M(i,j) )

       R2E(i,j) = ( tmp2 + Psi_2m(i,j) / tPrn ) &
                / ( tmp2 + Psi   (i,j) / tPrn )

       R10M(i,j) = min( max( R10M(i,j), 0.0_RP ), 1.0_RP )
       R2H (i,j) = min( max( R2H (i,j), 0.0_RP ), 1.0_RP )
       R2E (i,j) = min( max( R2E (i,j), 0.0_RP ), 1.0_RP )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_bulkcoef_uno

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_bulkcoef_beljaars( &
       ATM_Uabs, ATM_POTT, Z1, &
       SFC_TEMP,               &
       Z0M, Z0H, Z0E,          &
       Cm, Ch, Ce ,            &
       R10M, R2H, R2E          )
    use scale_const, only: &
       GRAV   => CONST_GRAV,  &
       KARMAN => CONST_KARMAN
    implicit none

    real(RP), intent(in)  :: ATM_Uabs(IA,JA) ! absolute velocity at Z1       [m/s]
    real(RP), intent(in)  :: ATM_POTT(IA,JA) ! potential temperature at z1, based on the local surface pressure [K]
    real(RP), intent(in)  :: Z1      (IA,JA) ! height of lowermost atmosphere grid (cell center) [m]
    real(RP), intent(in)  :: SFC_TEMP(IA,JA) ! potential temperature at surface skin [K]
    real(RP), intent(in)  :: Z0M     (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in)  :: Z0H     (IA,JA) ! roughness length for heat     [m]
    real(RP), intent(in)  :: Z0E     (IA,JA) ! roughness length for moisture [m]
    real(RP), intent(out) :: Cm      (IA,JA) ! bulk coefficient for momentum
    real(RP), intent(out) :: Ch      (IA,JA) ! bulk coefficient for heat
    real(RP), intent(out) :: Ce      (IA,JA) ! bulk coefficient for moisture
    real(RP), intent(out) :: R10M    (IA,JA)
    real(RP), intent(out) :: R2H     (IA,JA)
    real(RP), intent(out) :: R2E     (IA,JA)

    ! specific parameters
    real(RP), parameter :: RiB_min      = 1.E-4_RP
    real(RP), parameter :: res_criteria = 1.E-6_RP
    real(RP), parameter :: dL           = 1.E-10_RP ! delta Obukhov length [m]

    real(RP) :: RiB0(IA,JA) ! bulk Richardson number [no unit]
    real(RP) :: L   (IA,JA) ! Obukhov length [m]
    real(RP) :: Lpls(IA,JA) ! Obukhov length [m]

    real(RP) :: LogZ1Z0M(IA,JA)
    real(RP) :: LogZ1Z0H(IA,JA)
    real(RP) :: LogZ1Z0E(IA,JA)

    real(RP) ::  CmUS(IA,JA),  CmS(IA,JA)
    real(RP) ::  ChUS(IA,JA),  ChS(IA,JA)
    real(RP) ::  CeUS(IA,JA),  CeS(IA,JA)
    real(RP) :: dCmUS(IA,JA), dCmS(IA,JA), dCm(IA,JA)
    real(RP) :: dChUS(IA,JA), dChS(IA,JA), dCh(IA,JA)
    real(RP) :: dCeUS(IA,JA), dCeS(IA,JA), dCe(IA,JA)
    real(RP) :: res  (IA,JA), dres

    real(RP) :: fmUS_Z1L, fmUS_Z0ML, fhUS_Z1L, fhUS_Z0HL, fhUS_Z0EL
    real(RP) :: fmS_Z1L,  fmS_Z0ML,  fhS_Z1L,  fhS_Z0HL,  fhS_Z0EL

    real(RP) :: fmUS_10mL, fmUS_2mL, fhUS_2mL
    real(RP) :: fmS_10mL,  fmS_2mL,  fhS_2mL
    real(RP) :: fm    (IA,JA)
    real(RP) :: fh    (IA,JA)
    real(RP) :: fm_10m(IA,JA)
    real(RP) :: fm_2m (IA,JA)
    real(RP) :: fh_2m (IA,JA)
    real(RP) :: Psi   (IA,JA)
    real(RP) :: Psi_2m(IA,JA)

    integer, parameter :: itelim = 100 ! maximum iteration number

    real(RP) :: KARMAN2
    real(RP) :: tmp1, tmp2, sw
    integer  :: i, j, ite
    !---------------------------------------------------------------------------

    KARMAN2  = KARMAN * KARMAN

    do j = JS, JE
    do i = IS, IE
       RiB0(i,j)  = GRAV * Z1(i,j) * ( ATM_POTT(i,j) - SFC_TEMP(i,j)    ) &
                                   / ( ATM_POTT(i,j) * ATM_Uabs(i,j)**2 )
       RiB0(i,j) = max( abs(RiB0(i,j)), RiB_min )
    enddo
    enddo

    ! The initial Obukhov length is assumed by bulk Richardson number.
    do j = JS, JE
    do i = IS, IE
       L(i,j) = Z1(i,j) / RiB0(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       LogZ1Z0M(i,j) = log( Z1(i,j) / Z0M(i,j) )
       LogZ1Z0H(i,j) = log( Z1(i,j) / Z0H(i,j) )
       LogZ1Z0E(i,j) = log( Z1(i,j) / Z0E(i,j) )
    enddo
    enddo

    do ite = 1, itelim

       !---< calc coefficient >---

       ! unstable condition
       do j = JS, JE
       do i = IS, IE
          fmUS_Z1L  = fm_unstable( Z1 (i,j),L(i,j) )
          fmUS_Z0ML = fm_unstable( Z0M(i,j),L(i,j) )

          fhUS_Z1L  = fh_unstable( Z1 (i,j),L(i,j) )
          fhUS_Z0HL = fh_unstable( Z0H(i,j),L(i,j) )
          fhUS_Z0EL = fh_unstable( Z0E(i,j),L(i,j) )

          CmUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML )**2
          ChUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML ) / ( LogZ1Z0H(i,j) - fhUS_Z1L + fhUS_Z0HL )
          CeUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML ) / ( LogZ1Z0E(i,j) - fhUS_Z1L + fhUS_Z0EL )
       enddo
       enddo

       ! stable condition
       do j = JS, JE
       do i = IS, IE
          fmS_Z1L  = fm_stable( Z1 (i,j),L(i,j) )
          fmS_Z0ML = fm_stable( Z0M(i,j),L(i,j) )

          fhS_Z1L  = fh_stable( Z1 (i,j),L(i,j) )
          fhS_Z0HL = fh_stable( Z0H(i,j),L(i,j) )
          fhS_Z0EL = fh_stable( Z0E(i,j),L(i,j) )

          CmS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML )**2
          ChS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML ) / ( LogZ1Z0H(i,j) - fhS_Z1L + fhS_Z0HL )
          CeS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML ) / ( LogZ1Z0E(i,j) - fhS_Z1L + fhS_Z0EL )
       enddo
       enddo

       ! select unstable / stable
       do j = JS, JE
       do i = IS, IE
          sw = 0.5_RP - sign(0.5_RP,L(i,j)) ! if unstable, sw = 1

          Cm(i,j) = (        sw ) * CmUS(i,j) &
                  + ( 1.0_RP-sw ) * CmS (i,j)
          Ch(i,j) = (        sw ) * ChUS(i,j) &
                  + ( 1.0_RP-sw ) * ChS (i,j)
          Ce(i,j) = (        sw ) * CeUS(i,j) &
                  + ( 1.0_RP-sw ) * CeS (i,j)
       enddo
       enddo

       !---< calc derivative >---

       do j = JS, JE
       do i = IS, IE
          Lpls(i,j) = L(i,j) + dL
       enddo
       enddo

       ! unstable condition
       do j = JS, JE
       do i = IS, IE
          fmUS_Z1L  = fm_unstable( Z1 (i,j),Lpls(i,j) )
          fmUS_Z0ML = fm_unstable( Z0M(i,j),Lpls(i,j) )

          fhUS_Z1L  = fh_unstable( Z1 (i,j),Lpls(i,j) )
          fhUS_Z0HL = fh_unstable( Z0H(i,j),Lpls(i,j) )
          fhUS_Z0EL = fh_unstable( Z0E(i,j),Lpls(i,j) )

          dCmUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML )**2
          dChUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML ) / ( LogZ1Z0H(i,j) - fhUS_Z1L + fhUS_Z0HL )
          dCeUS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmUS_Z1L + fmUS_Z0ML ) / ( LogZ1Z0E(i,j) - fhUS_Z1L + fhUS_Z0EL )
       enddo
       enddo

       ! stable condition
       do j = JS, JE
       do i = IS, IE
          fmS_Z1L  = fm_stable( Z1 (i,j),Lpls(i,j) )
          fmS_Z0ML = fm_stable( Z0M(i,j),Lpls(i,j) )

          fhS_Z1L  = fh_stable( Z1 (i,j),Lpls(i,j) )
          fhS_Z0HL = fh_stable( Z0H(i,j),Lpls(i,j) )
          fhS_Z0EL = fh_stable( Z0E(i,j),Lpls(i,j) )

          dCmS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML )**2
          dChS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML ) / ( LogZ1Z0H(i,j) - fhS_Z1L + fhS_Z0HL )
          dCeS(i,j) = KARMAN2 / ( LOGZ1Z0M(i,j) - fmS_Z1L + fmS_Z0ML ) / ( LogZ1Z0E(i,j) - fhS_Z1L + fhS_Z0EL )
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          sw = 0.5_RP - sign(0.5_RP,Lpls(i,j)) ! if unstable, sw = 1

          dCm(i,j) = (        sw ) * dCmUS(i,j) &
                   + ( 1.0_RP-sw ) * dCmS (i,j)
          dCh(i,j) = (        sw ) * dChUS(i,j) &
                   + ( 1.0_RP-sw ) * dChS (i,j)
          dCe(i,j) = (        sw ) * dCeUS(i,j) &
                   + ( 1.0_RP-sw ) * dCeS (i,j)
       enddo
       enddo

       ! update Obukhov length
       do j = JS, JE
       do i = IS, IE
          ! residual
          res(i,j) = L(i,j) - Z1(i,j) * Cm(i,j)**1.5_RP / ( KARMAN * Ch(i,j) * RiB0(i,j) )
          ! d(residual)
          dres = Lpls(i,j) - Z1(i,j) * dCm(i,j)**1.5_RP / ( KARMAN * dCh(i,j) * RiB0(i,j) )

          L(i,j) = L(i,j) - res(i,j) / ( dres - res(i,j) ) * dL
       enddo
       enddo

       if( maxval( abs( res(:,:) ) ) < res_criteria ) exit
    enddo

    ! result
!    RiB(:,:) = Z1(:,:) * Cm(:,:)**1.5_RP / ( L(:,:) * KARMAN * Ch(:,:) )

    do j = JS, JE
    do i = IS, IE
       ! unstable condition
       fmUS_Z1L  = fm_unstable( Z1(i,j),L(i,j) )
       fhUS_Z1L  = fh_unstable( Z1(i,j),L(i,j) )
       fmUS_10mL = fm_unstable( 10.0_RP,L(i,j) )
       fmUS_2mL  = fm_unstable(  2.0_RP,L(i,j) )
       fhUS_2mL  = fh_unstable(  2.0_RP,L(i,j) )
       ! stable condition
       fmS_Z1L   = fm_stable  ( Z1(i,j),L(i,j) )
       fhS_Z1L   = fh_stable  ( Z1(i,j),L(i,j) )
       fmS_10mL  = fm_stable  ( 10.0_RP,L(i,j) )
       fmS_2mL   = fm_stable  (  2.0_RP,L(i,j) )
       fhS_2mL   = fh_stable  (  2.0_RP,L(i,j) )

       ! select unstable / stable
       sw = 0.5_RP - sign(0.5_RP,L(i,j)) ! if unstable, sw = 1

       fm    (i,j) = (        sw ) * fmUS_Z1L  &
                   + ( 1.0_RP-sw ) * fmS_Z1L
       fh    (i,j) = (        sw ) * fhUS_Z1L  &
                   + ( 1.0_RP-sw ) * fhS_Z1L
       fm_10m(i,j) = (        sw ) * fmUS_10mL &
                   + ( 1.0_RP-sw ) * fmS_10mL
       fm_2m (i,j) = (        sw ) * fmUS_2mL  &
                   + ( 1.0_RP-sw ) * fmS_2mL
       fh_2m (i,j) = (        sw ) * fhUS_2mL  &
                   + ( 1.0_RP-sw ) * fhS_2mL
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       R10M(i,j) = log( 10.0_RP/Z0M(i,j) ) / LogZ1Z0M(i,j) &
                 / sqrt( fm_10m(i,j) / fm(i,j) )
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Psi   (i,j) = sqrt(fm   (i,j)) / fh   (i,j) * LogZ1Z0M(i,j)
       Psi_2m(i,j) = sqrt(fm_2m(i,j)) / fh_2m(i,j) * log( 2.0_RP/Z0M(i,j) )

       tmp1 = log( Z0M(i,j)/Z0H(i,j) ) / LogZ1Z0M(i,j)

       R2H(i,j) = ( tmp1 + Psi_2m(i,j) ) &
                / ( tmp1 + Psi   (i,j) )

       tmp2 = log( Z0M(i,j)/Z0E(i,j) ) / LogZ1Z0M(i,j)

       R2E(i,j) = ( tmp2 + Psi_2m(i,j) ) &
                / ( tmp2 + Psi   (i,j) )
    enddo
    enddo


    if ( ite > itelim ) then
       if( IO_L ) write(IO_FID_LOG,*) 'Error: reach maximum iteration in the function of ATMOS_PHY_SF_bulkcoef_beljaars.'
    endif

    do j = JS, JE
    do i = IS, IE
       Cm(i,j) = min( max( Cm(i,j), ATMOS_PHY_SF_BULKCOEF_Cm_min ), ATMOS_PHY_SF_BULKCOEF_Cm_max )
       Ch(i,j) = min( max( Ch(i,j), ATMOS_PHY_SF_BULKCOEF_Ch_min ), ATMOS_PHY_SF_BULKCOEF_Ch_max )
       Ce(i,j) = min( max( Ce(i,j), ATMOS_PHY_SF_BULKCOEF_Ce_min ), ATMOS_PHY_SF_BULKCOEF_Ce_max )
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       R10M(i,j) = min( max( R10M(i,j), 0.0_RP ), 1.0_RP )
       R2H (i,j) = min( max( R2H (i,j), 0.0_RP ), 1.0_RP )
       R2E (i,j) = min( max( R2E (i,j), 0.0_RP ), 1.0_RP )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_bulkcoef_beljaars

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in unstable condition
  function fm_unstable( Z, L ) result(fmUS)
    use scale_const, only: &
       PI  => CONST_PI,  &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L
    real(RP)             :: fmUS

    real(RP) :: R, sqsqR
    !---------------------------------------------------------------------------

    R     = Z / min(L,-EPS) ! should be negative
    sqsqR = ( 1.0_RP - 16.0_RP * R )**0.25_RP

    ! Paulson (1974); Dyer (1974)
    fmUS = log( ( 1.0_RP+sqsqR )**2 * ( 1.0_RP+sqsqR*sqsqR ) * 0.125_RP ) &
         - 2.0_RP * atan( sqsqR )                                         &
         + 0.5_RP * PI

    return
  end function fm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, L ) result(fhUS)
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L
    real(RP)             :: fhUS

    real(RP) :: R, sqR
    !---------------------------------------------------------------------------

    R   = Z / min(L,-EPS) ! should be negative
    sqR = sqrt( 1.0_RP - 16.0_RP * R )

    ! Paulson (1974); Dyer (1974)
    fhUS = 2.0_RP * log( ( 1.0_RP+sqR ) * 0.5_RP )

    return
  end function fh_unstable

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in stable condition
  function fm_stable( Z, L ) result(fmS)
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L
    real(RP)             :: fmS

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a =   1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c =   5.0_RP
    real(RP), parameter :: d =  0.35_RP

    real(RP) :: R
    !---------------------------------------------------------------------------

    R = Z / max(L,EPS) ! should be positive

    ! Holtslag and DeBruin (1988)
    fmS = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fm_stable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in stable condition
  function fh_stable( Z, L ) result(fhS)
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L
    real(RP)             :: fhS

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a =   1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c =   5.0_RP
    real(RP), parameter :: d =  0.35_RP

    real(RP) :: R
    !---------------------------------------------------------------------------

    R = Z / max(L,EPS) ! should be positive

    ! Beljaars and Holtslag (1991)
    fhS = 1.0_RP - ( 1.0_RP + 2.0_RP/3.0_RP * a*R )**1.5_RP - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fh_stable

end module scale_atmos_phy_sf_bulkcoef
