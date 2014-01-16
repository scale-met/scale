!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Typical vertical profile
!!
!! @par Description
!!          Generate typical vertical profile of atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-25 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_profile
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PROFILE_isa

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
  !> Generate profile (International Standard Atmosphere)
  subroutine ATMOS_PROFILE_isa( &
       pott,     &
       temp_sfc, &
       pres_sfc  )
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       RovCP => CONST_RovCP, &
       P00   => CONST_PRE00
    use mod_grid, only: &
       CZ => GRID_CZ
    implicit none

    real(RP), intent(out) :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: temp_sfc !< surface temperature   [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure      [Pa]

    ! ISA profile
    integer,  parameter :: nref = 8
    real(RP), parameter :: CZ_isa(nref) = (/     0.0_RP, &
                                             11000.0_RP, &
                                             20000.0_RP, &
                                             32000.0_RP, &
                                             47000.0_RP, &
                                             51000.0_RP, &
                                             71000.0_RP, &
                                             84852.0_RP  /)
    real(RP), parameter :: GAMMA(nref)  = (/ -6.5E-3_RP, &
                                                 0.0_RP, &
                                              1.0E-3_RP, &
                                              2.8E-3_RP, &
                                              0.0E-3_RP, &
                                             -2.8E-3_RP, &
                                             -2.0E-3_RP, &
                                                 0.0_RP  /)
    real(RP) :: temp_isa(nref)
    real(RP) :: pres_isa(nref)

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)

    real(RP) :: gmr ! grav / Rdry
    integer  :: i, k
    !---------------------------------------------------------------------------

    gmr = GRAV / Rdry

    !--- ISA profile
    temp_isa(1) = temp_sfc
    pres_isa(1) = pres_sfc

    do i = 2, nref
       temp_isa(i) = temp_isa(i-1) + GAMMA(i-1) * ( CZ_isa(i)-CZ_isa(i-1) )

       if ( GAMMA(i-1) == 0.0_RP ) then
          pres_isa(i) = pres_isa(i-1) * exp( -gmr / temp_isa(i) * ( CZ_isa(i)-CZ_isa(i-1) ) )
       else
          pres_isa(i) = pres_isa(i-1) * ( temp_isa(i)/temp_isa(i-1) ) ** ( -gmr/GAMMA(i-1) )
       endif
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### ICAO International Standard Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:  lapse rate:    pressure: temperature'
    do i = 1, nref
       if( IO_L ) write(IO_FID_LOG,'(4(f13.5))') CZ_isa(i), GAMMA(i), pres_isa(i), temp_isa(i)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    !--- make reference state
    do k = KS, KE
       do i = 2, nref
          if ( CZ(k) > CZ_isa(i-1) .AND. CZ(k) <= CZ_isa(i) ) then

             temp(k) = temp_isa(i-1) + GAMMA(i-1) * ( CZ(k)-CZ_isa(i-1) )
             if ( GAMMA(i-1) == 0.0_RP ) then
                pres(k) = pres_isa(i-1) * exp( -gmr/temp_isa(i-1) * ( CZ(k)-CZ_isa(i-1) ) )
             else
                pres(k) = pres_isa(i-1) * ( temp(k)/temp_isa(i-1) ) ** ( -gmr/GAMMA(i-1) )
             endif

          elseif ( CZ(k) <= CZ_isa(1)    ) then

             temp(k) = temp_isa(1) + GAMMA(1) * ( CZ(k)-CZ_isa(1) )
             pres(k) = pres_isa(1) * ( temp(k)/temp_isa(1) ) ** ( -gmr/GAMMA(1) )

          elseif ( CZ(k)  > CZ_isa(nref) ) then

             temp(k) = temp(k-1)
             pres(k) = pres_isa(i-1) * exp( -gmr/temp_isa(i-1) * ( CZ(k)-CZ_isa(i-1) ) )

          endif
       enddo

       pott(k) = temp(k) * ( P00/pres(k) )**RovCP
    enddo

    return
  end subroutine ATMOS_PROFILE_isa

end module mod_atmos_profile
