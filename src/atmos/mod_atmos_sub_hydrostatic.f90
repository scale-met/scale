!-------------------------------------------------------------------------------
!> module hydrostatic
!!
!! @par Description
!!          make hydrostatic profile in the model
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-02-20 (H.Yashiro)  [new] Extract from the tool of Y.Miyamoto
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_hydro_buildrho

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'

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
  logical, parameter :: use_rapserate = .true.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho( &
      dens,     &
      temp,     &
      pres,     &
      pott,     &
      qv,       &
      temp_sfc, &
      pres_sfc, &
      pott_sfc, &
      qv_sfc    )
    use mod_const, only : &
       GRAV    => CONST_GRAV,    &
       Rdry    => CONST_Rdry,    &
       RovCP   => CONST_RovCP,   &
       CPovR   => CONST_CPovR,   &
       RovCV   => CONST_RovCV,   &
       CVovCP  => CONST_CVovCP,  &
       CPovCV  => CONST_CPovCV,  &
       LASPdry => CONST_LASPdry, &
       EPSTvap => CONST_EPSTvap, &
       P00     => CONST_PRE00
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
    implicit none

    real(8), intent(out) :: dens(KA,IA,JA) !< density [kg/m3]
    real(8), intent(out) :: temp(KA,IA,JA) !< temperature [K]
    real(8), intent(out) :: pres(KA,IA,JA) !< pressure [Pa]
    real(8), intent(in)  :: pott(KA,IA,JA) !< potential temperature [K]
    real(8), intent(in)  :: qv  (KA,IA,JA) !< water vapor [kg/kg]
    real(8), intent(out) :: temp_sfc(1,IA,JA) !< surface temperature [K]
    real(8), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure [Pa]
    real(8), intent(in)  :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(8), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor [kg/kg]

    real(8) :: Rmoist_sfc(1,IA,JA)
    real(8) :: dens_sfc  (1,IA,JA)

    real(8) :: Rmoist(KA,IA,JA)

    real(8) :: dens_s, dhyd, dgrd

    real(8), parameter :: criteria = 1.D-15
    integer, parameter :: itelim = 100

    integer :: k, i, j, ite
    !---------------------------------------------------------------------------

    ! make density at surface
    do j = JS, JE
    do i = IS, IE
       Rmoist_sfc(1,i,j) = Rdry * ( 1.D0 + EPSTvap * qv_sfc(1,i,j) )

       dens_sfc(1,i,j) = P00 / Rmoist_sfc(1,i,j) / pott_sfc(1,i,j) * ( pres_sfc(1,i,j)/P00 )**CVovCP
       temp_sfc(1,i,j) = pres_sfc(1,i,j) / ( dens_sfc(1,i,j) * Rmoist_sfc(1,i,j) )
    enddo
    enddo
!    if( IO_L ) write(IO_FID_LOG,*) &
!    'SFC', dens_sfc(1,IS,JS), temp_sfc(1,IS,JS), pres_sfc(1,IS,JS), pott_sfc(1,IS,JS), qv_sfc(1,IS,JS)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Rmoist(k,i,j) = Rdry * ( 1.D0 + EPSTvap * qv(k,i,j) )
    enddo
    enddo
    enddo

    ! make density at lowermost cell center
    k = KS

    if ( use_rapserate ) then

       do j = JS, JE
       do i = IS, IE
          temp(k,i,j) = pott_sfc(1,i,j) - LASPdry * CZ(k) ! use dry lapse rate
          pres(k,i,j) = P00 * ( temp(k,i,j)/pott(k,i,j) )**CPovR
          dens(k,i,j) = P00 / Rmoist(k,i,j) / pott(k,i,j) * ( pres(k,i,j)/P00 )**CVovCP
       enddo
       enddo

    else ! use itelation

       do j = JS, JE
       do i = IS, IE
          dens_s      = 0.D0
          dens(k,i,j) = dens_sfc(1,i,j) ! first guess

          do ite = 1, itelim
             if( abs(dens(k,i,j)-dens_s) <= criteria ) exit

             dens_s = dens(k,i,j)

             dhyd = + ( P00 * ( dens_sfc(1,i,j) * Rmoist_sfc(1,i,j) * pott_sfc(1,i,j) / P00 )**CPovCV &
                      - P00 * ( dens_s          * Rmoist    (k,i,j) * pott    (k,i,j) / P00 )**CPovCV ) / CZ(k) & ! dp/dz
                    - GRAV * 0.5D0 * ( dens_sfc(1,i,j) + dens_s )                                                 ! rho*g

             dgrd = - P00 * ( Rmoist(k,i,j) * pott(k,i,j) / P00 )**CPovCV / CZ(k) &
                    * CPovCV * dens_s**RovCV                                      &
                    - 0.5D0 * GRAV

             dens(k,i,j) = dens_s - dhyd/dgrd

!             if ( i == IS .AND. j == JS ) then
!                if( IO_L ) write(IO_FID_LOG,*) k, ite, dens(k,i,j)-dens_s, dens(k,i,j), dens_s, dhyd/dgrd, dhyd, dgrd
!             endif

          enddo

          if ( ite > itelim ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k,i,j), dens_s, dhyd, dgrd
          endif
       enddo
       enddo

!       if( IO_L ) write(IO_FID_LOG,*) k, dens(k,IS,JS), pott(k,IS,JS), qv(k,IS,JS), CZ(k)
    endif

    ! make density
    do j = JS, JE
    do i = IS, IE
       do k = KS+1, KE

          dens_s      = 0.D0
          dens(k,i,j) = dens(k-1,i,j)

          do ite = 1, itelim
             if( abs(dens(k,i,j)-dens_s) <= criteria ) exit

             dens_s = dens(k,i,j)

             dhyd = + ( P00 * ( dens(k-1,i,j) * Rmoist(k-1,i,j) * pott(k-1,i,j) / P00 )**CPovCV &
                      - P00 * ( dens_s        * Rmoist(k  ,i,j) * pott(k  ,i,j) / P00 )**CPovCV ) / FDZ(k-1) & ! dp/dz
                    - GRAV * 0.5D0 * ( dens(k-1,i,j) + dens_s )                                                ! rho*g

             dgrd = - P00 * ( Rmoist(k,i,j) * pott(k,i,j) / P00 )**CPovCV / FDZ(k-1) &
                    * CPovCV * dens_s**RovCV                                         &
                    - 0.5D0 * GRAV

             dens(k,i,j) = dens_s - dhyd/dgrd

!             if ( i == IS .AND. j == JS ) then
!                if( IO_L ) write(IO_FID_LOG,*) k, ite, dens(k,i,j)-dens_s, dens(k,i,j), dens_s, dhyd/dgrd, dhyd, dgrd
!             endif

          enddo

!          if ( i == IS .AND. j == JS ) then
!             if( IO_L ) write(IO_FID_LOG,*) k, dens(k,IS,JS), pott(k,IS,JS), qv(k,IS,JS), FDZ(k-1)
!          endif

          if ( ite > itelim ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', &
                  k, ite, dens(k,i,j), dens_s, dhyd, dgrd, Rmoist, FDZ(k-1)
          endif
       enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pres(k,i,j) = P00 * ( dens(k,i,j) * Rmoist(k,i,j) * pott(k,i,j) / P00 )**CPovCV
       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rmoist(k,i,j) )
    enddo
    enddo
    enddo

    ! fill KHALO
    do j  = JS, JE
    do i  = IS, IE
       dens(   1:KS-1,i,j) = dens(KS,i,j)
       dens(KE+1:KA,  i,j) = dens(KE,i,j)
    enddo
    enddo
    ! fill IHALO & JHALO
    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

!    do k = KS+1, KE
!       if( IO_L ) write(IO_FID_LOG,*) k, dens(k,IS,JS), temp(k,IS,JS), pres(k,IS,JS), pott(k,IS,JS), qv(k,IS,JS)
!    enddo

    return
  end subroutine ATMOS_hydro_buildrho

end module mod_atmos_hydrostatic
