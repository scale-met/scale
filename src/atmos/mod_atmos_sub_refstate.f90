!-------------------------------------------------------------------------------
!> module Atmosphere / reference state
!!
!! @par Description
!!          Reference state of Atmosphere
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] integrate
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_refstate
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_REFSTATE_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, allocatable, save :: REF_dens(:,:,:) ! refernce density [kg/m3]
  real(8), public, allocatable, save :: REF_pott(:,:,:) ! refernce potential temperature [K]
  real(8), public, allocatable, save :: REF_pres(:,:,:) ! refernce pressure [Pa]
  real(8), public, allocatable, save :: REF_temp(:,:,:) ! refernce temperature [K]

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
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_setup
    use mod_stdio, only: &
       IO_get_available_fid, &
       IO_FILECHR,  &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF8,          &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,  &
       RovCP  => CONST_RovCP, &
       Pstd   => CONST_Pstd,  &
       Tstd   => CONST_Tstd
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       CZ   => GRID_CZ
    implicit none

    real(8) :: TEMP_SFC ! surface temperature

    NAMELIST / PARAM_ATMOS_REFSTATE / &
       TEMP_SFC

    integer, parameter :: nref = 8
    real(8), parameter :: CZ_isa(nref) = (/     0.0D0,  &
                                            11000.0D0,  &
                                            20000.0D0,  &
                                            32000.0D0,  &
                                            47000.0D0,  &
                                            51000.0D0,  &
                                            71000.0D0,  &
                                            84852.0D0   /)
    real(8), parameter :: GAMMA(nref)  = (/    -6.5D-3, &
                                                0.0D0 , &
                                                1.0D-3, &
                                                2.8D-3, &
                                                0.0D-3, &
                                               -2.8D-3, &
                                               -2.0D-3, &
                                                0.0D0   /)
    real(8) :: temp_isa(nref)
    real(8) :: pres_isa(nref)

    real(8), allocatable :: dens(:) ! refernce density [kg/m3]
    real(8), allocatable :: pott(:) ! refernce potential temperature [K]
    real(8), allocatable :: pres(:) ! refernce pressure [Pa]
    real(8), allocatable :: temp(:) ! refernce temperature [K]

    real(8) :: gmr !! grav / Rair
    integer :: i, k
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    gmr      = GRAV / Rair
    TEMP_SFC = Tstd

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REFSTATE]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_REFSTATE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_REFSTATE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_REFSTATE)


    !--- ISA profile
    temp_isa(1) = TEMP_SFC
    pres_isa(1) = Pstd

    do i = 2, nref
       temp_isa(i) = temp_isa(i-1) + GAMMA(i-1) * ( CZ_isa(i)-CZ_isa(i-1) )

       if ( GAMMA(i-1) == 0.D0 ) then
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


    allocate( pres(KA) )
    allocate( temp(KA) )
    allocate( dens(KA) )
    allocate( pott(KA) )

    !--- make reference state
    do k = KS, KE
       do i = 2, nref
          if ( CZ(k) > CZ_isa(i-1) .AND. CZ(k) <= CZ_isa(i) ) then

             temp(k) = temp_isa(i-1) + GAMMA(i-1) * ( CZ(k)-CZ_isa(i-1) )
             if ( GAMMA(i-1) == 0.D0 ) then
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

       dens(k) = pres(k) / ( temp(k) * Rair )
       pott(k) = temp(k) * ( Pstd/pres(k) )**RovCP
    enddo


    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(5(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    allocate( REF_pres(IA,JA,KA) ); REF_pres(:,:,:) = CONST_UNDEF8
    allocate( REF_temp(IA,JA,KA) ); REF_temp(:,:,:) = CONST_UNDEF8
    allocate( REF_pott(IA,JA,KA) ); REF_pott(:,:,:) = CONST_UNDEF8
    allocate( REF_dens(IA,JA,KA) ); REF_dens(:,:,:) = CONST_UNDEF8

    do k = KS, KE
       REF_pres(:,:,k) = pres(k)
       REF_temp(:,:,k) = temp(k)
       REF_dens(:,:,k) = dens(k)
       REF_pott(:,:,k) = pott(k)
    enddo

    return
  end subroutine ATMOS_REFSTATE_setup

end module mod_atmos_refstate
