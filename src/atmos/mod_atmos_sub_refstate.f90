!-------------------------------------------------------------------------------
!> module Atmosphere / reference state
!!
!! @par Description
!!          Reference state of Atmosphere
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_refstate
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use gtool_file_h, only : &
     File_HLONG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_REFSTATE_setup
  public :: ATMOS_REFSTATE_read
  public :: ATMOS_REFSTATE_write
  public :: ATMOS_REFSTATE_generate

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: ATMOS_REFSTATE_dens(KA) ! refernce density [kg/m3]
  real(RP), public, save :: ATMOS_REFSTATE_pott(KA) ! refernce potential temperature [K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_IN_BASENAME   = ''
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_OUT_BASENAME  = ''
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_TITLE     = 'SCALE3 Refstate'
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_SOURCE    = 'SCALE-LES ver. 3'
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_INSTITUTE = 'AISC/RIKEN'
  character(len=IO_SYSCHR),  private :: ATMOS_REFSTATE_TYPE          = 'ISA'
  real(RP),                  private :: ATMOS_REFSTATE_TEMP_SFC      = 300.0_RP ! surface temperature
  real(RP),                  private :: ATMOS_REFSTATE_POTT_UNIFORM  = 300.0_RP ! uniform potential temperature

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_REFSTATE / &
       ATMOS_REFSTATE_IN_BASENAME,   &
       ATMOS_REFSTATE_OUT_BASENAME,  &
       ATMOS_REFSTATE_OUT_TITLE,     &
       ATMOS_REFSTATE_OUT_SOURCE,    &
       ATMOS_REFSTATE_OUT_INSTITUTE, &
       ATMOS_REFSTATE_TYPE,          &
       ATMOS_REFSTATE_POTT_UNIFORM,  &
       ATMOS_REFSTATE_TEMP_SFC

    integer :: ierr
    !---------------------------------------------------------------------------

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

    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       call ATMOS_REFSTATE_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found reference state file. Generate!'

       if ( trim(ATMOS_REFSTATE_TYPE) == 'ISA' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: ISA'
       elseif ( trim(ATMOS_REFSTATE_TYPE) == 'UNIFORM' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: UNIFORM POTT'
       else
          write(*,*) 'xxx ATMOS_REFSTATE_TYPE must be "ISA" or "UNIFORM". Check!', trim(ATMOS_REFSTATE_TYPE)
          call PRC_MPIstop
       endif
    
       call ATMOS_REFSTATE_generate
    endif

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then
       call ATMOS_REFSTATE_write
    endif

    return
  end subroutine ATMOS_REFSTATE_setup

  !-----------------------------------------------------------------------------
  !> Read Reference state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_read
    use gtool_file, only: &
       FileRead
    use mod_process, only: &
       PRC_myrank
    implicit none

    character(len=IO_FILECHR) :: bname
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input Refstate file ***'

    write(bname,'(A,A,F15.3)') trim(ATMOS_REFSTATE_IN_BASENAME)

    call FileRead( ATMOS_REFSTATE_dens(:), bname, 'DENS', 1, PRC_myrank, single=.true. )
    call FileRead( ATMOS_REFSTATE_pott(:), bname, 'POTT', 1, PRC_myrank, single=.true. )

    return
  end subroutine ATMOS_REFSTATE_read

  !-----------------------------------------------------------------------------
  !> Write Reference state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_write
    use mod_process, only: &
       PRC_myrank, &
       PRC_master
    use gtool_file_h, only: &
       File_HMID,  &
       File_REAL4, &
       File_REAL8
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FileWrite, &
       FileClose
    use mod_grid, only : &
         GRID_CZ
    implicit none

    character(len=IO_FILECHR) :: bname
    integer :: dtype
    integer :: fid, vid_dens, vid_pott
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output grid file ***'
    if( IO_L ) write(IO_FID_LOG,*) '*** Only at Master node ***'

    if ( RP == 8 ) then
       dtype = File_REAL8
    else if ( RP == 4 ) then
       dtype = File_REAL4
    end if

    if ( PRC_myrank == PRC_master ) then
       write(bname,'(A,A,F15.3)') trim(ATMOS_REFSTATE_OUT_BASENAME)
       call FileCreate( fid,              & ! (out)
            bname,                        & ! (in)
            ATMOS_REFSTATE_OUT_TITLE,     & ! (in)
            ATMOS_REFSTATE_OUT_SOURCE,    & ! (in)
            ATMOS_REFSTATE_OUT_INSTITUTE, & ! (in)
            (/'z'/), (/KMAX/), (/'Z'/),   & ! (in)
            (/'m'/), (/File_REAL4/),      & ! (in)
            PRC_master, PRC_myrank,       & ! (in)
            single = .true.               ) ! (in)

       call FileAddVariable( vid_dens,                      & ! (out)
            fid, 'DENS', 'Reference state of rho', 'kg/m3', & ! (in)
            (/'z'/), File_REAL8                             ) ! (in)
       call FileAddVariable( vid_pott,                      & ! (out)
            fid, 'DENS', 'Reference state of theta', 'K',   & ! (in)
            (/'z'/), File_REAL8                             ) ! (in)

       call FilePutAxis(fid, 'z', GRID_CZ(KS:KE))

       call FileWrite( vid_dens, ATMOS_REFSTATE_dens(:), 0.0_RP, 0.0_RP )
       call FileWrite( vid_pott, ATMOS_REFSTATE_pott(:), 0.0_RP, 0.0_RP )

       call FileClose( fid )
    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate Reference state (International Standard Atmosphere)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_generate
    use mod_const, only : &
       GRAV   => CONST_GRAV,  &
       Rdry   => CONST_Rdry,  &
       RovCP  => CONST_RovCP, &
       Pstd   => CONST_Pstd,  &
       P00    => CONST_PRE00
    use mod_grid, only : &
       CZ   => GRID_CZ
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho_1d => ATMOS_HYDRO_buildrho_1d
    implicit none

    integer, parameter :: nref = 8
    real(RP), parameter :: CZ_isa(nref) = (/     0.0_RP,  &
                                            11000.0_RP,  &
                                            20000.0_RP,  &
                                            32000.0_RP,  &
                                            47000.0_RP,  &
                                            51000.0_RP,  &
                                            71000.0_RP,  &
                                            84852.0_RP   /)
    real(RP), parameter :: GAMMA(nref)  = (/    -6.5E-3_RP, &
                                                0.0_RP , &
                                                1.0E-3_RP, &
                                                2.8E-3_RP, &
                                                0.0E-3_RP, &
                                               -2.8E-3_RP, &
                                               -2.0E-3_RP, &
                                                0.0_RP   /)
    real(RP) :: temp_isa(nref)
    real(RP) :: pres_isa(nref)

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)

    real(RP) :: qv(KA) = 0.0_RP
    real(RP) :: qc(KA) = 0.0_RP
    real(RP) :: temp_sfc
    real(RP) :: qv_sfc = 0.0_RP
    real(RP) :: qc_sfc = 0.0_RP

    real(RP) :: gmr !! grav / Rdry
    integer :: i, k
    !---------------------------------------------------------------------------

    gmr      = GRAV / Rdry

    !--- ISA profile
    temp_isa(1) = ATMOS_REFSTATE_TEMP_SFC
    pres_isa(1) = Pstd

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

!       dens(k) = pres(k) / ( temp(k) * Rdry )
       pott(k) = temp(k) * ( P00/pres(k) )**RovCP
    enddo

!    dens(   1:KS-1) = dens(KS)
!    dens(KE+1:KA  ) = dens(KE)
    pott(   1:KS-1) = pott(KS)
    pott(KE+1:KA  ) = pott(KE)

    if ( trim(ATMOS_REFSTATE_TYPE) == 'UNIFORM' ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** pot.temp. is overwrited by uniform value:', ATMOS_REFSTATE_POTT_UNIFORM
       pott(:) = ATMOS_REFSTATE_POTT_UNIFORM
    endif

    ! make density & pressure profile 
    call hydro_buildrho_1d( dens    (:), & ! [OUT]
                            temp    (:), & ! [OUT]
                            pres    (:), & ! [OUT]
                            pott    (:), &
                            qv      (:), &
                            qc      (:), &
                            temp_sfc,    & ! [OUT]
                            pres_isa(1), &
                            temp_isa(1), &
                            qv_sfc,      &
                            qc_sfc       )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(5(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)

    return
  end subroutine ATMOS_REFSTATE_generate

end module mod_atmos_refstate
