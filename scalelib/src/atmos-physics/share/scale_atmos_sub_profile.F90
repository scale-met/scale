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
module scale_atmos_profile
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_const, only: &
     GRAV  => CONST_GRAV,  &
     Rdry  => CONST_Rdry,  &
     RovCP => CONST_RovCP, &
     P00   => CONST_PRE00
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PROFILE_isa

  interface ATMOS_PROFILE_isa
     module procedure ATMOS_PROFILE_isa_1D
     module procedure ATMOS_PROFILE_isa_3D
  end interface ATMOS_PROFILE_isa

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
  ! ISA profile
  integer,  private, parameter :: nref = 8
  real(RP), private, parameter :: z_isa(nref) = (/     0.0_RP, &
                                                   11000.0_RP, &
                                                   20000.0_RP, &
                                                   32000.0_RP, &
                                                   47000.0_RP, &
                                                   51000.0_RP, &
                                                   71000.0_RP, &
                                                   84852.0_RP  /)
  real(RP), private, parameter :: GAMMA(nref)  = (/ -6.5E-3_RP, &
                                                        0.0_RP, &
                                                     1.0E-3_RP, &
                                                     2.8E-3_RP, &
                                                     0.0E-3_RP, &
                                                    -2.8E-3_RP, &
                                                    -2.0E-3_RP, &
                                                        0.0_RP  /)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Generate profile (International Standard Atmosphere)
  subroutine ATMOS_PROFILE_isa_1D( &
       KA, KS, KE, &
       temp_sfc,   &
       pres_sfc,   &
       z,          &
       pott        )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: temp_sfc !< surface temperature   [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure      [Pa]
    real(RP), intent(in)  :: z   (KA) !< potential temperature [K]
    real(RP), intent(out) :: pott(KA) !< potential temperature [K]

    real(RP) :: temp_isa(nref)
    real(RP) :: pres_isa(nref)
    real(RP) :: temp(KA)
    real(RP) :: pres(KA)

    real(RP) :: gmr ! grav / Rdry
    integer  :: k, i, j, n
    !---------------------------------------------------------------------------

    gmr = GRAV / Rdry

    !--- ISA profile
    temp_isa(1) = temp_sfc
    pres_isa(1) = pres_sfc

    do j = JS, JE
    do i = IS, IE
    do n = 2, nref
       temp_isa(n) = temp_isa(n-1) + GAMMA(n-1) * ( z_isa(n)-z_isa(n-1) )

       if ( GAMMA(n-1) == 0.0_RP ) then
          pres_isa(n) = pres_isa(n-1) * exp( -gmr / temp_isa(n) * ( z_isa(n)-z_isa(n-1) ) )
       else
          pres_isa(n) = pres_isa(n-1) * ( temp_isa(n)/temp_isa(n-1) ) ** ( -gmr/GAMMA(n-1) )
       endif
    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### ICAO International Standard Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:  lapse rate:    pressure: temperature'
    do n = 1, nref
       if( IO_L ) write(IO_FID_LOG,'(4(f13.5))') z_isa(n), GAMMA(n), pres_isa(n), temp_isa(n)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    !--- make reference state
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       do n = 2, nref
          if ( z(k) > z_isa(n-1) .AND. z(k) <= z_isa(n) ) then

             temp(k) = temp_isa(n-1) + GAMMA(n-1) * ( z(k)-z_isa(n-1) )
             if ( GAMMA(n-1) == 0.0_RP ) then
                pres(k) = pres_isa(n-1) * exp( -gmr/temp_isa(n-1) * ( z(k)-z_isa(n-1) ) )
             else
                pres(k) = pres_isa(n-1) * ( temp(k)/temp_isa(n-1) ) ** ( -gmr/GAMMA(n-1) )
             endif

          elseif ( z(k) <= z_isa(1)    ) then

             temp(k) = temp_isa(1) + GAMMA(1) * ( z(k)-z_isa(1) )
             pres(k) = pres_isa(1) * ( temp(k)/temp_isa(1) ) ** ( -gmr/GAMMA(1) )

          elseif ( z(k)  > z_isa(nref) ) then

             temp(k) = temp(k-1)
             pres(k) = pres_isa(n-1) * exp( -gmr/temp_isa(n-1) * ( z(k)-z_isa(n-1) ) )

          endif
       enddo

       pott(k) = temp(k) * ( P00/pres(k) )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PROFILE_isa_1D

  !-----------------------------------------------------------------------------
  !> Generate profile (International Standard Atmosphere)
  subroutine ATMOS_PROFILE_isa_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       temp_sfc,   &
       pres_sfc,   &
       z,          &
       pott        )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: temp_sfc(IA,JA)    !< surface temperature   [K]
    real(RP), intent(in)  :: pres_sfc(IA,JA)    !< surface pressure      [Pa]
    real(RP), intent(in)  :: z       (KA,IA,JA) !< potential temperature [K]
    real(RP), intent(out) :: pott    (KA,IA,JA) !< potential temperature [K]

    real(RP) :: temp_isa(nref,IA,JA)
    real(RP) :: pres_isa(nref,IA,JA)
    real(RP) :: temp(KA)
    real(RP) :: pres(KA)

    real(RP) :: gmr ! grav / Rdry
    integer  :: k, i, j, n
    !---------------------------------------------------------------------------

    gmr = GRAV / Rdry

    !--- ISA profile
    do j = JS, JE
    do i = IS, IE
       temp_isa(1,i,j) = temp_sfc(i,j)
       pres_isa(1,i,j) = pres_sfc(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do n = 2, nref
       temp_isa(n,i,j) = temp_isa(n-1,i,j) + GAMMA(n-1) * ( z_isa(n)-z_isa(n-1) )

       if ( GAMMA(n-1) == 0.0_RP ) then
          pres_isa(n,i,j) = pres_isa(n-1,i,j) * exp( -gmr / temp_isa(n,i,j) * ( z_isa(n)-z_isa(n-1) ) )
       else
          pres_isa(n,i,j) = pres_isa(n-1,i,j) * ( temp_isa(n,i,j)/temp_isa(n-1,i,j) ) ** ( -gmr/GAMMA(n-1) )
       endif
    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### ICAO International Standard Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:  lapse rate:    pressure: temperature'
    do n = 1, nref
       if( IO_L ) write(IO_FID_LOG,'(4(f13.5))') z_isa(n), GAMMA(n), pres_isa(n,IS,JS), temp_isa(n,IS,JS)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    !--- make reference state
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       do n = 2, nref
          if ( z(k,i,j) > z_isa(n-1) .AND. z(k,i,j) <= z_isa(n) ) then

             temp(k) = temp_isa(n-1,i,j) + GAMMA(n-1) * ( z(k,i,j)-z_isa(n-1) )
             if ( GAMMA(n-1) == 0.0_RP ) then
                pres(k) = pres_isa(n-1,i,j) * exp( -gmr/temp_isa(n-1,i,j) * ( z(k,i,j)-z_isa(n-1) ) )
             else
                pres(k) = pres_isa(n-1,i,j) * ( temp(k)/temp_isa(n-1,i,j) ) ** ( -gmr/GAMMA(n-1) )
             endif

          elseif ( z(k,i,j) <= z_isa(1)    ) then

             temp(k) = temp_isa(1,i,j) + GAMMA(1) * ( z(k,i,j)-z_isa(1) )
             pres(k) = pres_isa(1,i,j) * ( temp(k)/temp_isa(1,i,j) ) ** ( -gmr/GAMMA(1) )

          elseif ( z(k,i,j)  > z_isa(nref) ) then

             temp(k) = temp(k-1)
             pres(k) = pres_isa(n-1,i,j) * exp( -gmr/temp_isa(n-1,i,j) * ( z(k,i,j)-z_isa(n-1) ) )

          endif
       enddo

       pott(k,i,j) = temp(k) * ( P00/pres(k) )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PROFILE_isa_3D

end module scale_atmos_profile
