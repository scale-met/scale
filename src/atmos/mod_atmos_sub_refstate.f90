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
  real(RP), public, save :: ATMOS_REFSTATE_dens(KA) ! refernce density [kg/m3]
  real(RP), public, save :: ATMOS_REFSTATE_pott(KA) ! refernce potential temperature [K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: ATMOS_REFSTATE_generate_isa
  public :: ATMOS_REFSTATE_generate_uniform
  public :: ATMOS_REFSTATE_generate_frominit

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_IN_BASENAME   = ''
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_OUT_BASENAME  = ''
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_TITLE     = 'SCALE3 Refstate'
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_SOURCE    = 'SCALE-LES ver. 3'
  character(len=File_HLONG), private :: ATMOS_REFSTATE_OUT_INSTITUTE = 'AISC/RIKEN'
  character(len=IO_SYSCHR),  private :: ATMOS_REFSTATE_TYPE          = 'UNIFORM'
  real(RP),                  private :: ATMOS_REFSTATE_TEMP_SFC      = 300.0_RP ! surface temperature
  real(RP),                  private :: ATMOS_REFSTATE_RH            = 0.0_RP   ! surface & environment RH
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
       ATMOS_REFSTATE_TEMP_SFC,      &
       ATMOS_REFSTATE_RH

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
          call ATMOS_REFSTATE_generate_isa
       elseif ( trim(ATMOS_REFSTATE_TYPE) == 'UNIFORM' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: UNIFORM POTT'
          call ATMOS_REFSTATE_generate_uniform
       elseif ( trim(ATMOS_REFSTATE_TYPE) == 'INIT' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: make from initial data'
          call ATMOS_REFSTATE_generate_frominit
       else
          write(*,*) 'xxx ATMOS_REFSTATE_TYPE must be "ISA" or "UNIFORM". Check!', trim(ATMOS_REFSTATE_TYPE)
          call PRC_MPIstop
       endif
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
    use dc_types, only : &
         DP
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

       call FileWrite( vid_dens, ATMOS_REFSTATE_dens(:), 0.0_DP, 0.0_DP )
       call FileWrite( vid_pott, ATMOS_REFSTATE_pott(:), 0.0_DP, 0.0_DP )

       call FileClose( fid )
    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate Reference state (International Standard Atmosphere)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_generate_isa
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
    use mod_atmos_saturation, only: &
       SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    implicit none

    integer, parameter :: nref = 8
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
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc(1)
    real(RP) :: pres_sfc(1)
    real(RP) :: pott_sfc(1)
    real(RP) :: qv_sfc  (1)
    real(RP) :: qc_sfc  (1)

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc(1)

    real(RP) :: gmr !! grav / Rdry
    integer  :: i, k
    !---------------------------------------------------------------------------

    gmr = GRAV / Rdry

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

       pott(k) = temp(k) * ( P00/pres(k) )**RovCP
    enddo

    qv      (:) = 0.0_RP
    qc      (:) = 0.0_RP
    pott_sfc(1) = temp_isa(1)
    pres_sfc(1) = pres_isa(1)
    qv_sfc  (1) = 0.0_RP
    qc_sfc  (1) = 0.0_RP

    ! make density & pressure profile in dry condition
    call hydro_buildrho_1d( dens    (:), & ! [OUT]
                            temp    (:), & ! [OUT]
                            pres    (:), & ! [OUT]
                            pott    (:), & ! [IN]
                            qv      (:), & ! [IN]
                            qc      (:), & ! [IN]
                            temp_sfc(1), & ! [OUT]
                            pres_sfc(1), & ! [IN]
                            pott_sfc(1), & ! [IN]
                            qv_sfc  (1), & ! [IN]
                            qc_sfc  (1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_liq( qsat_sfc(1), temp_sfc(1), pres_sfc(1) )
    call SATURATION_pres2qsat_liq( qsat(:),  temp(:),  pres(:)  )

    qv_sfc(1) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc(1)
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho_1d( dens    (:), & ! [OUT]
                            temp    (:), & ! [OUT]
                            pres    (:), & ! [OUT]
                            pott    (:), & ! [IN]
                            qv      (:), & ! [IN]
                            qc      (:), & ! [IN]
                            temp_sfc(1), & ! [OUT]
                            pres_sfc(1), & ! [IN]
                            pott_sfc(1), & ! [IN]
                            qv_sfc  (1), & ! [IN]
                            qc_sfc  (1)  ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k), qv(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)

    return
  end subroutine ATMOS_REFSTATE_generate_isa

  !-----------------------------------------------------------------------------
  !> Generate Reference state (Uniform Potential Temperature)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_generate_uniform
    use mod_const, only : &
       Pstd   => CONST_Pstd
    use mod_grid, only : &
       CZ   => GRID_CZ
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho_1d => ATMOS_HYDRO_buildrho_1d
    use mod_atmos_saturation, only: &
       SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    implicit none

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc(1)
    real(RP) :: pres_sfc(1)
    real(RP) :: pott_sfc(1)
    real(RP) :: qv_sfc  (1)
    real(RP) :: qc_sfc  (1)

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc(1)

    integer  :: k
    !---------------------------------------------------------------------------

    pres_sfc(1) = Pstd
    pott_sfc(1) = ATMOS_REFSTATE_TEMP_SFC
    qv_sfc  (1) = 0.0_RP
    qc_sfc  (1) = 0.0_RP

    do k = KS, KE
       pott(k) = ATMOS_REFSTATE_POTT_UNIFORM
       qv  (k) = 0.0_RP
       qc  (k) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho_1d( dens    (:), & ! [OUT]
                            temp    (:), & ! [OUT]
                            pres    (:), & ! [OUT]
                            pott    (:), & ! [IN]
                            qv      (:), & ! [IN]
                            qc      (:), & ! [IN]
                            temp_sfc(1), & ! [OUT]
                            pres_sfc(1), & ! [IN]
                            pott_sfc(1), & ! [IN]
                            qv_sfc  (1), & ! [IN]
                            qc_sfc  (1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_liq( qsat_sfc(1), temp_sfc(1), pres_sfc(1) )
    call SATURATION_pres2qsat_liq( qsat(:),  temp(:),  pres(:)  )

    qv_sfc(1) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc(1)
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho_1d( dens    (:), & ! [OUT]
                            temp    (:), & ! [OUT]
                            pres    (:), & ! [OUT]
                            pott    (:), & ! [IN]
                            qv      (:), & ! [IN]
                            qc      (:), & ! [IN]
                            temp_sfc(1), & ! [OUT]
                            pres_sfc(1), & ! [IN]
                            pott_sfc(1), & ! [IN]
                            qv_sfc  (1), & ! [IN]
                            qc_sfc  (1)  ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k), qv(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)

    return
  end subroutine ATMOS_REFSTATE_generate_uniform

  !-----------------------------------------------------------------------------
  !> Generate Reference state (Horizontal average from initial data)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_REFSTATE_generate_frominit
    use mod_const, only : &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       Rvap   => CONST_Rvap,   &
       P00    => CONST_PRE00
    use mod_grid, only : &
       CZ   => GRID_CZ
    use mod_comm, only: &
       COMM_horizontal_mean
    use mod_atmos_thermodyn, only: &
       CPw => AQ_CP
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho_1d => ATMOS_HYDRO_buildrho_1d
    use mod_atmos_vars, only: &
       DENS_3d => DENS, &
       RHOT_3d => RHOT, &
       QTRC_3d => QTRC
    implicit none

    real(RP) :: PRES_3d(KA,IA,JA)
    real(RP) :: TEMP_3d(KA,IA,JA)
    real(RP) :: POTT_3d(KA,IA,JA)

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)

    real(RP) :: QDRY, RTOT, CPTOT, CPovCV

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       POTT_3d(k,i,j) = RHOT_3d(k,i,j) / DENS_3d(k,i,j)

       QDRY  = 1.0_RP
       CPTOT = 0.0_RP
       do iq = QQS, QQE
          QDRY  = QDRY  - QTRC_3d(k,i,j,iq)
          CPTOT = CPTOT + QTRC_3d(k,i,j,iq) * CPw(iq)
       enddo
       RTOT   = Rdry *QDRY + Rvap*QTRC_3d(k,i,j,I_QV)
       CPTOT  = CPdry*QDRY + CPTOT
       CPovCV = CPTOT / ( CPTOT - RTOT )

       PRES_3d(k,i,j) = P00 * ( RHOT_3d(k,i,j) * RTOT / P00 )**CPovCV
       TEMP_3d(k,i,j) = PRES_3d(k,i,j) / ( DENS_3d(k,i,j) * RTOT )
    enddo
    enddo
    enddo

    call COMM_horizontal_mean( pres(:), PRES_3d(:,:,:) )
    call COMM_horizontal_mean( temp(:), TEMP_3d(:,:,:) )
    call COMM_horizontal_mean( dens(:), DENS_3d(:,:,:) )
    call COMM_horizontal_mean( pott(:), POTT_3d(:,:,:) )

    call COMM_horizontal_mean( qv(:), QTRC_3d(:,:,:,I_QV) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k), qv(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)

    return
  end subroutine ATMOS_REFSTATE_generate_frominit

end module mod_atmos_refstate
