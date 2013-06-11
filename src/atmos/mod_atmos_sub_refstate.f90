!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Reference state
!!
!! @par Description
!!          Reference state of Atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2013-02-25 (H.Yashiro)  [mod] Separate ISA profile to mod_atmos_sub_profile
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_refstate
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file_h, only: &
     File_HLONG
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_REFSTATE_setup
  public :: ATMOS_REFSTATE_read
  public :: ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: ATMOS_REFSTATE_dens(KA) !< refernce density [kg/m3]
  real(RP), public, save :: ATMOS_REFSTATE_pott(KA) !< refernce potential temperature [K]
  real(RP), public, save :: ATMOS_REFSTATE_qv  (KA) !< refernce vapor [kg/kg]

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
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_IN_BASENAME  = ''                !< basename of the input  file
  character(len=IO_FILECHR), private :: ATMOS_REFSTATE_OUT_BASENAME = ''                !< basename of the output file
  character(len=IO_SYSCHR),  private :: ATMOS_REFSTATE_OUT_TITLE    = 'SCALE3 Refstate' !< title    of the output file
  character(len=IO_SYSCHR),  private :: ATMOS_REFSTATE_OUT_DTYPE    = 'DEFAULT'         !< REAL4 or REAL8

  character(len=IO_SYSCHR),  private :: ATMOS_REFSTATE_TYPE         = 'UNIFORM'         !< profile type
  real(RP),                  private :: ATMOS_REFSTATE_TEMP_SFC     = 300.0_RP          !< surface temperature           [K]
  real(RP),                  private :: ATMOS_REFSTATE_RH           =   0.0_RP          !< surface & environment RH      [%]
  real(RP),                  private :: ATMOS_REFSTATE_POTT_UNIFORM = 300.0_RP          !< uniform potential temperature [K]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_REFSTATE_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_REFSTATE / &
       ATMOS_REFSTATE_IN_BASENAME,  &
       ATMOS_REFSTATE_OUT_BASENAME, &
       ATMOS_REFSTATE_OUT_TITLE,    &
       ATMOS_REFSTATE_OUT_DTYPE,    &
       ATMOS_REFSTATE_TYPE,         &
       ATMOS_REFSTATE_POTT_UNIFORM, &
       ATMOS_REFSTATE_TEMP_SFC,     &
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


    ! input or generate reference profile
    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then
       call ATMOS_REFSTATE_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found reference state file. Generate!'

       if ( ATMOS_REFSTATE_TYPE == 'ISA' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: ISA'
          call ATMOS_REFSTATE_generate_isa
       elseif ( ATMOS_REFSTATE_TYPE == 'UNIFORM' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: UNIFORM POTT'
          call ATMOS_REFSTATE_generate_uniform
       elseif ( ATMOS_REFSTATE_TYPE == 'INIT' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Reference type: make from initial data'
          call ATMOS_REFSTATE_generate_frominit
       else
          write(*,*) 'xxx ATMOS_REFSTATE_TYPE must be "ISA" or "UNIFORM". Check!', trim(ATMOS_REFSTATE_TYPE)
          call PRC_MPIstop
       endif
    endif

    ! output reference profile
    call ATMOS_REFSTATE_write

    return
  end subroutine ATMOS_REFSTATE_setup

  !-----------------------------------------------------------------------------
  !> Read reference state profile
  subroutine ATMOS_REFSTATE_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input reference state profile ***'

    if ( ATMOS_REFSTATE_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_REFSTATE_dens(:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'DENS_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_pott(:),                             & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'POTT_ref', 'Z', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_REFSTATE_qv(:),                               & ! [OUT]
                         ATMOS_REFSTATE_IN_BASENAME, 'QV_ref',   'Z', step=1 ) ! [IN]

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** refstate file is not specified.'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_REFSTATE_read

  !-----------------------------------------------------------------------------
  !> Write reference state profile
  subroutine ATMOS_REFSTATE_write
    use mod_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( ATMOS_REFSTATE_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output reference state profile ***'

       call FILEIO_write( ATMOS_REFSTATE_dens(:), ATMOS_REFSTATE_OUT_BASENAME,  ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'DENS_ref', 'Reference profile of rho', 'kg/m3', 'Z', ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( ATMOS_REFSTATE_pott(:), ATMOS_REFSTATE_OUT_BASENAME,  ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'POTT_ref', 'Reference profile of theta', 'K', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( ATMOS_REFSTATE_qv(:),   ATMOS_REFSTATE_OUT_BASENAME,  ATMOS_REFSTATE_OUT_TITLE, & ! [IN]
                          'QV_ref',   'Reference profile of qv', 'kg/kg', 'Z',   ATMOS_REFSTATE_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_REFSTATE_write

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (International Standard Atmosphere)
  subroutine ATMOS_REFSTATE_generate_isa
    use mod_const, only: &
       Pstd => CONST_Pstd
    use mod_grid, only: &
       CZ => GRID_CZ
    use mod_atmos_profile, only: &
       PROFILE_isa => ATMOS_PROFILE_isa
    use mod_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use mod_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    implicit none

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc
    real(RP) :: pres_sfc
    real(RP) :: pott_sfc
    real(RP) :: qv_sfc
    real(RP) :: qc_sfc

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    pres_sfc = Pstd

    call PROFILE_isa( pott(:),  & ! [OUT]
                      pott_sfc, & ! [IN]
                      pres_sfc  ) ! [IN]

    qv(:)  = 0.0_RP
    qc(:)  = 0.0_RP
    qv_sfc = 0.0_RP
    qc_sfc = 0.0_RP

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc, temp_sfc, pres_sfc )
    call SATURATION_pres2qsat_all( qsat(:),  temp(:),  pres(:)  )

    qv_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k), qv(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)
    ATMOS_REFSTATE_qv  (:) = qv(:)

    return
  end subroutine ATMOS_REFSTATE_generate_isa

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (Uniform Potential Temperature)
  subroutine ATMOS_REFSTATE_generate_uniform
    use mod_const, only: &
       Pstd   => CONST_Pstd
    use mod_grid, only: &
       CZ   => GRID_CZ
    use mod_atmos_hydrostatic, only: &
       HYDROSTATIC_buildrho => ATMOS_HYDROSTATIC_buildrho
    use mod_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    implicit none

    real(RP) :: temp(KA)
    real(RP) :: pres(KA)
    real(RP) :: dens(KA)
    real(RP) :: pott(KA)
    real(RP) :: qv  (KA)
    real(RP) :: qc  (KA)

    real(RP) :: temp_sfc
    real(RP) :: pres_sfc
    real(RP) :: pott_sfc
    real(RP) :: qv_sfc  
    real(RP) :: qc_sfc  

    real(RP) :: qsat(KA)
    real(RP) :: qsat_sfc

    integer  :: k
    !---------------------------------------------------------------------------

    pres_sfc = Pstd
    pott_sfc = ATMOS_REFSTATE_TEMP_SFC
    qv_sfc   = 0.0_RP
    qc_sfc   = 0.0_RP

    do k = KS, KE
       pott(k) = ATMOS_REFSTATE_POTT_UNIFORM
       qv  (k) = 0.0_RP
       qc  (k) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc, temp_sfc, pres_sfc )
    call SATURATION_pres2qsat_all( qsat(:),  temp(:),  pres(:)  )

    qv_sfc = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat_sfc
    do k = KS, KE
       qv(k) = ATMOS_REFSTATE_RH * 1.E-2_RP * qsat(k)
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( dens(:),  & ! [OUT]
                               temp(:),  & ! [OUT]
                               pres(:),  & ! [OUT]
                               pott(:),  & ! [IN]
                               qv  (:),  & ! [IN]
                               qc  (:),  & ! [IN]
                               temp_sfc, & ! [OUT]
                               pres_sfc, & ! [IN]
                               pott_sfc, & ! [IN]
                               qv_sfc,   & ! [IN]
                               qc_sfc    ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### Generated Reference State of Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) '      height:    pressure: temperature:     density:   pot.temp.: water vapor'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(6(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k), qv(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'

    ATMOS_REFSTATE_dens(:) = dens(:)
    ATMOS_REFSTATE_pott(:) = pott(:)
    ATMOS_REFSTATE_qv  (:) = qv(:)

    return
  end subroutine ATMOS_REFSTATE_generate_uniform

  !-----------------------------------------------------------------------------
  !> Generate reference state profile (Horizontal average from initial data)
  subroutine ATMOS_REFSTATE_generate_frominit
    use mod_const, only: &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       Rvap   => CONST_Rvap,   &
       P00    => CONST_PRE00
    use mod_grid, only: &
       CZ   => GRID_CZ
    use mod_comm, only: &
       COMM_horizontal_mean
    use mod_atmos_thermodyn, only: &
       CPw => AQ_CP
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
    ATMOS_REFSTATE_qv  (:) = qv(:)

    return
  end subroutine ATMOS_REFSTATE_generate_frominit

end module mod_atmos_refstate
