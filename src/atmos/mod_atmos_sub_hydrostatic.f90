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
     IO_FID_CONF, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_hydro_buildrho
  public :: ATMOS_hydro_buildrho_1d

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_precision.h'

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
  logical :: HYDROSTATIC_userapserate = .true.
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
      qc,       &
      temp_sfc, &
      pres_sfc, &
      pott_sfc, &
      qc_sfc,&
      qv_sfc    )
    use mod_const, only : &
       GRAV    => CONST_GRAV,    &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CPovR   => CONST_CPovR,   &
       RovCV   => CONST_RovCV,   &
       CVovCP  => CONST_CVovCP,  &
       CPovCV  => CONST_CPovCV,  &
       LASPdry => CONST_LASPdry, &
       P00     => CONST_PRE00
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA,IA,JA) !< density [kg/m3]
    real(RP), intent(out) :: temp(KA,IA,JA) !< temperature [K]
    real(RP), intent(out) :: pres(KA,IA,JA) !< pressure [Pa]
    real(RP), intent(in)  :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA) !< water vapor [kg/kg]

    real(RP), intent(out) :: temp_sfc(1,IA,JA) !< surface temperature [K]
    real(RP), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure [Pa]
    real(RP), intent(in)  :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (1,IA,JA) !< surface water vapor [kg/kg]

    real(RP) :: Rmoist_sfc(1,IA,JA)
    real(RP) :: dens_sfc  (1,IA,JA)

    real(RP) :: Rmoist(KA,IA,JA)

    real(RP) :: dens_s, dhyd, dgrd

    real(RP), parameter :: criteria = 1.0E-15_RP
    integer, parameter :: itelim = 100

    integer :: k, i, j, ite, ierr

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_userapserate

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used!'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif

    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    ! make density at surface
    do j = JS, JE
    do i = IS, IE
       Rmoist_sfc(1,i,j) = Rdry * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                         + Rvap * qv_sfc(1,i,j)

       dens_sfc(1,i,j) = P00 / Rmoist_sfc(1,i,j) / pott_sfc(1,i,j) * ( pres_sfc(1,i,j)/P00 )**CVovCP
       temp_sfc(1,i,j) = pres_sfc(1,i,j) / ( dens_sfc(1,i,j) * Rmoist_sfc(1,i,j) )
    enddo
    enddo
!    if( IO_L ) write(IO_FID_LOG,*) &
!    'SFC', dens_sfc(1,IS,JS), temp_sfc(1,IS,JS), pres_sfc(1,IS,JS), pott_sfc(1,IS,JS), qv_sfc(1,IS,JS)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Rmoist(k,i,j) = Rdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + Rvap * qv(k,i,j)
    enddo
    enddo
    enddo

    ! make density at lowermost cell center
    k = KS

    if ( HYDROSTATIC_userapserate ) then

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
          dens_s      = 0.0_RP
          dens(k,i,j) = dens_sfc(1,i,j) ! first guess

          do ite = 1, itelim
             if( abs(dens(k,i,j)-dens_s) <= criteria ) exit

             dens_s = dens(k,i,j)

             dhyd = + ( P00 * ( dens_sfc(1,i,j) * Rmoist_sfc(1,i,j) * pott_sfc(1,i,j) / P00 )**CPovCV &
                      - P00 * ( dens_s          * Rmoist    (k,i,j) * pott    (k,i,j) / P00 )**CPovCV ) / CZ(k) & ! dp/dz
                    - GRAV * 0.5_RP * ( dens_sfc(1,i,j) + dens_s )                                                 ! rho*g

             dgrd = - P00 * ( Rmoist(k,i,j) * pott(k,i,j) / P00 )**CPovCV / CZ(k) &
                    * CPovCV * dens_s**RovCV                                      &
                    - 0.5_RP * GRAV

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

          dens_s      = 0.0_RP
          dens(k,i,j) = dens(k-1,i,j)

          do ite = 1, itelim
             if( abs(dens(k,i,j)-dens_s) <= criteria ) exit

             dens_s = dens(k,i,j)

             dhyd = + ( P00 * ( dens(k-1,i,j) * Rmoist(k-1,i,j) * pott(k-1,i,j) / P00 )**CPovCV &
                      - P00 * ( dens_s        * Rmoist(k  ,i,j) * pott(k  ,i,j) / P00 )**CPovCV ) / FDZ(k-1) & ! dp/dz
                    - GRAV * 0.5_RP * ( dens(k-1,i,j) + dens_s )                                                ! rho*g

             dgrd = - P00 * ( Rmoist(k,i,j) * pott(k,i,j) / P00 )**CPovCV / FDZ(k-1) &
                    * CPovCV * dens_s**RovCV                                         &
                    - 0.5_RP * GRAV

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

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho_1d( &
      dens,     &
      temp,     &
      pres,     &
      pott,     &
      qv,       &
      qc,       &
      temp_sfc, &
      pres_sfc, &
      pott_sfc, &
      qc_sfc,&
      qv_sfc    )
    use mod_const, only : &
       GRAV    => CONST_GRAV,    &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CPovR   => CONST_CPovR,   &
       RovCV   => CONST_RovCV,   &
       CVovCP  => CONST_CVovCP,  &
       CPovCV  => CONST_CPovCV,  &
       LASPdry => CONST_LASPdry, &
       P00     => CONST_PRE00
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out) :: temp(KA) !< temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)  :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< water vapor [kg/kg]

    real(RP), intent(out) :: temp_sfc !< surface temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(RP), intent(in)  :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface water vapor [kg/kg]

    real(RP) :: Rmoist_sfc
    real(RP) :: dens_sfc

    real(RP) :: Rmoist(KA)

    real(RP) :: dens_s, dhyd, dgrd

    real(RP), parameter :: criteria = 1.0E-15_RP
    integer, parameter :: itelim = 100

    integer :: k, ite, ierr

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_userapserate

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC_1d]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used!'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif

    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    ! make density at surface
    Rmoist_sfc = Rdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap * qv_sfc

    dens_sfc = P00 / Rmoist_sfc / pott_sfc * ( pres_sfc/P00 )**CVovCP
    temp_sfc = pres_sfc / ( dens_sfc * Rmoist_sfc )

    do k = KS, KE
       Rmoist(k) = Rdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + Rvap * qv(k)
    enddo

    ! make density at lowermost cell center
    k = KS

    if ( HYDROSTATIC_userapserate ) then

       temp(k) = pott_sfc - LASPdry * CZ(k) ! use dry lapse rate
       pres(k) = P00 * ( temp(k)/pott(k) )**CPovR
       dens(k) = P00 / Rmoist(k) / pott(k) * ( pres(k)/P00 )**CVovCP

    else ! use itelation

       dens_s      = 0.0_RP
       dens(k) = dens_sfc ! first guess

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens_sfc * Rmoist_sfc * pott_sfc / P00 )**CPovCV &
                   - P00 * ( dens_s   * Rmoist(k)  * pott(k)  / P00 )**CPovCV ) / CZ(k) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc + dens_s )                                                 ! rho*g

          dgrd = - P00 * ( Rmoist(k) * pott(k) / P00 )**CPovCV / CZ(k) &
                 * CPovCV * dens_s**RovCV                                      &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd
       endif

    endif

    ! make density
    do k = KS+1, KE

       dens_s      = 0.0_RP
       dens(k) = dens(k-1)

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens(k-1) * Rmoist(k-1) * pott(k-1) / P00 )**CPovCV &
                   - P00 * ( dens_s    * Rmoist(k  ) * pott(k  ) / P00 )**CPovCV ) / FDZ(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )                                                ! rho*g

          dgrd = - P00 * ( Rmoist(k) * pott(k) / P00 )**CPovCV / FDZ(k-1) &
                 * CPovCV * dens_s**RovCV                                         &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd, Rmoist, FDZ(k-1)
       endif
    enddo

    do k = KS, KE
       pres(k) = P00 * ( dens(k) * Rmoist(k) * pott(k) / P00 )**CPovCV
       temp(k) = pres(k) / ( dens(k) * Rmoist(k) )
    enddo

    ! fill KHALO
    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)

    return
  end subroutine ATMOS_hydro_buildrho_1d

end module mod_atmos_hydrostatic
