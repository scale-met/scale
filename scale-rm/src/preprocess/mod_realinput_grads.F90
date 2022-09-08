!-------------------------------------------------------------------------------
!> module REAL input GrADS
!!
!! @par Description
!!          read data from GrADS file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module mod_realinput_grads
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_tracer
  use scale_prc, only: &
     myrank => PRC_myrank,  &
     PRC_abort
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtmosSetupGrADS
  public :: ParentAtmosInputGrADS
  public :: ParentLandSetupGrADS
  public :: ParentLandInputGrADS
  public :: ParentOceanSetupGrADS
  public :: ParentOceanInputGrADS

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
  integer :: file_id_atm, file_id_ocn, file_id_lnd

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupGrADS( &
      dims,      &
      timelen,   &
      qtrc_flag, &
      LON_all,   &
      LAT_all,   &
      basename_org, &
      basename_num  )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file_grads, only: &
       FILE_GrADS_open,  &
       FILE_GrADS_varid, &
       FILE_GrADS_get_shape
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    integer,           intent(out) :: dims(6)
    integer,           intent(out) :: timelen
    logical,           intent(out) :: qtrc_flag(QA)
    real(RP), pointer, intent(out) :: LON_all(:,:)
    real(RP), pointer, intent(out) :: LAT_all(:,:)

    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: basename_num

!    namelist / PARAM_MKINIT_REAL_GrADS /

    integer :: var_id
    real(RP), allocatable :: lon1d(:), lat1d(:)


    integer :: i, j, iq
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentAtmosSetupGrADS",*) 'Setup'

!!$    !--- read namelist
!!$    rewind(IO_FID_CONF)
!!$    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_GrADS,iostat=ierr)
!!$
!!$    if( ierr > 0 ) then
!!$       LOG_ERROR("ParentAtmosSetupGrADS",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_GrADS. Check!'
!!$       call PRC_abort
!!$    endif
!!$    LOG_NML(PARAM_MKINIT_REAL_GrADS)


    if ( basename_org == "" ) then
       LOG_ERROR("ParentAtmosSetupGrADS",*) '"BASENAME_ORG" is not specified in "PARAM_MKINIT_REAL_ATMOS"!', trim(basename_org)
       call PRC_abort
    endif

    call FILE_GrADS_open( basename_org, & ! (in)
                          file_id_atm   ) ! (out)

    do iq = 1, QA
       if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
       call FILE_GrADS_varid( file_id_atm, TRACER_NAME(iq), var_id )
       qtrc_flag(iq) =  var_id > 0
    end do

    call FILE_GrADS_get_shape( file_id_atm, "U", & ! (in)
                               dims(:3)          ) ! (out)

    ! half level
    dims(4) = dims(1)
    dims(5) = dims(2)
    dims(6) = dims(3)


    allocate( lon_all(dims(2), dims(3)) )
    allocate( lat_all(dims(2), dims(3)) )


    call read2d( (/1,1/), dims(2:3), lon_all(:,:), "lon", file_id_atm, basename_num, oneD=1 )
    lon_all(:,:) = lon_all(:,:) * D2R

    call read2d( (/1,1/), dims(2:3), lat_all(:,:), "lat", file_id_atm, basename_num, oneD=2 )
    lat_all(:,:) = lat_all(:,:) * D2R

    ! tentative
    timelen = 1

    return
  end subroutine ParentAtmosSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputGrADS( &
       KA_org, KS_org, KE_org, &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       QA,     &
       w_org, &
       u_org, &
       v_org, &
       pres_org, &
       dens_org, &
       pt_org, &
       temp_org, &
       qv_org,   &
       rh_org,   &
       qhyd_org, &
       qtrc_org,&
       cz_org,   &
       nopres, nodens, &
       temp2pt, rh2qv, &
       basename_num,  &
       sfc_diagnoses, &
       nt )
    use scale_const, only: &
       UNDEF   => CONST_UNDEF, &
       D2R     => CONST_D2R,   &
       EPS     => CONST_EPS,   &
       EPSvap  => CONST_EPSvap, &
       EPSTvap => CONST_EPSTvap, &
       GRAV    => CONST_GRAV, &
       LAPS    => CONST_LAPS, &
       P00     => CONST_PRE00, &
       Rdry    => CONST_Rdry, &
       CPdry   => CONST_CPdry
    use scale_atmos_hydrometeor, only: &
       HYD_NAME, &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    integer,          intent(in)  :: KA_org
    integer,          intent(in)  :: KS_org
    integer,          intent(in)  :: KE_org
    integer,          intent(in)  :: IA_org
    integer,          intent(in)  :: IS_org
    integer,          intent(in)  :: IE_org
    integer,          intent(in)  :: JA_org
    integer,          intent(in)  :: JS_org
    integer,          intent(in)  :: JE_org
    integer,          intent(in)  :: QA
    real(RP),         intent(out) :: w_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: u_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: v_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: pres_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: dens_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: pt_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: temp_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: qv_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: rh_org(KA_org,IA_org,JA_org)
    real(RP),         intent(out) :: qhyd_org(KA_org,IA_org,JA_org,N_HYD)
    real(RP),         intent(out) :: qtrc_org(KA_org,IA_org,JA_org,QA)
    real(RP),         intent(out) :: cz_org(KA_org,IA_org,JA_org)
    logical,          intent(out) :: nopres
    logical,          intent(out) :: nodens
    logical,          intent(out) :: temp2pt
    logical,          intent(out) :: rh2qv
    character(len=*), intent(in)  :: basename_num
    logical,          intent(in)  :: sfc_diagnoses
    integer,          intent(in)  :: nt

    character(len=H_SHORT) :: item

    integer  :: lm_layer(IA_org,JA_org)

    real(RP) :: work(KA_org-2,IA_org,JA_org)
    real(RP) :: work2d(IA_org,JA_org)
    integer  :: start(3), count(3)

    logical :: exist

    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    start(:) = (/1,IS_org,JS_org/)
    count(:) = (/KA_org-2,IA_org,JA_org/)

    ! pressure
    nopres = .false.
    call read3d( start(:), count(:), work(:,:,:), "pressure", file_id_atm, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call read3d( start(:), count(:), work(:,:,:), "plev", file_id_atm, basename_num, exist=exist, step=nt )
       if ( .not. exist ) nopres = .true.
    end if
    if ( .not. nopres ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          pres_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
    end if


    ! dens
    call read3d( start(:), count(:), work(:,:,:), "DENS", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          dens_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
       nodens = .false.
    else
       nodens = .true.
    end if

    ! W
    call read3d( start(:), count(:), work(:,:,:), "W", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          W_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
    else
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          W_org(k+2,i,j) = 0.0_RP
       enddo
       enddo
       enddo
    end if

    ! U
    call read3d( start(:), count(:), work(:,:,:), "U", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          U_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
    else
       LOG_ERROR("ParentAtmosInputGrADS",*) '"U" is requierd'
       call PRC_abort
    end if

    ! V
    call read3d( start(:), count(:), work(:,:,:), "V", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          V_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
    else
       LOG_ERROR("ParentAtmosInputGrADS",*) '"V" is requierd'
       call PRC_abort
    end if

    ! T
    call read3d( start(:), count(:), work(:,:,:), "PT", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          PT_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
       temp2pt = .false.
    else
       call read3d( start(:), count(:), work(:,:,:), "T", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             temp_org(k+2,i,j) = work(k,i,j)
          enddo
          enddo
          enddo
          temp2pt = .true.
       else
          LOG_ERROR("ParentAtmosInputGrADS",*) '"PT" or "T" is requierd'
          call PRC_abort
       end if
    end if

    ! height
    call read3d( start(:), count(:), work(:,:,:), "height", file_id_atm, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call read3d( start(:), count(:), work(:,:,:), "HGT", file_id_atm, basename_num, exist=exist, step=nt )
       if ( .not. exist ) then
          LOG_ERROR("ParentAtmosInputGrADS",*) '"height" or "HGT" is requierd'
          call PRC_abort
       end if
    end if
    !$omp parallel do collapse(2)
    do j = 1, JA_org
    do i = 1, IA_org
    do k = 1, KA_org-2
       cz_org(k+2,i,j) = work(k,i,j)
    enddo
    enddo
    enddo

    ! QV
    call read3d( start(:), count(:), work(:,:,:), "QV", file_id_atm, basename_num, exist=exist, step=nt )
    if ( exist ) then
       !$omp parallel do collapse(2)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = 1, KA_org-2
          qv_org(k+2,i,j) = work(k,i,j)
       enddo
       enddo
       enddo
       rh2qv = .false.
    else
       call read3d( start(:), count(:), work(:,:,:), "RH", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             rh_org(k+2,i,j) = work(k,i,j)
          enddo
          enddo
          enddo
          rh2qv = .true.
       else
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qv_org(k+2,i,j) = 0.0_RP
          enddo
          enddo
          enddo
          rh2qv = .false.
       end if
    end if
    
    ! QC, QR, QI, QS, QG, QH
    do iq = 1, N_HYD
       call read3d( start(:), count(:), work(:,:,:), HYD_NAME(iq), file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qhyd_org(k+2,i,j,iq) = work(k,i,j)
          enddo
          enddo
          enddo
       else
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qhyd_org(k+2,i,j,iq) = 0.0_RP
          enddo
          enddo
          enddo
       end if
    end do
    
    ! QTRC
    do iq = 1, QA
       if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
       call read3d( start(:), count(:), work(:,:,:), TRACER_NAME(iq), file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qtrc_org(k+2,i,j,iq) = work(k,i,j)
          enddo
          enddo
          enddo
       else
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qtrc_org(k+2,i,j,iq) = UNDEF
          enddo
          enddo
          enddo
       end if
    end do


    if ( sfc_diagnoses ) then

       ! MSLP
       call read2d( start(2:), count(2:), work2d(:,:), "MSLP", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(1,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(1,i,j) = UNDEF
          enddo
          enddo
       end if

       ! SFC_PRES, PSFC
       call read2d( start(2:), count(2:), work2d(:,:), "SFC_PRES", file_id_atm, basename_num, exist=exist, step=nt )
       if ( .not. exist ) then
          call read2d( start(2:), count(2:), work2d(:,:), "PSFC", file_id_atm, basename_num, exist=exist, step=nt )
       end if
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(2,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(2,i,j) = UNDEF
          enddo
          enddo
       end if

       ! U10
       call read2d( start(2:), count(2:), work2d(:,:), "U10", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             U_org(2,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             U_org(2,i,j) = UNDEF
          enddo
          enddo
       end if

       ! V10
       call read2d( start(2:), count(2:), work2d(:,:), "V10", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             V_org(2,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             V_org(2,i,j) = UNDEF
          enddo
          enddo
       end if

       ! T2
       call read2d( start(2:), count(2:), work2d(:,:), "T2", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             temp_org(2,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             temp_org(2,i,j) = UNDEF
          enddo
          enddo
       end if

       ! Q2, RH2
       if ( rh2qv ) then
          call read2d( start(2:), count(2:), work2d(:,:), "RH2", file_id_atm, basename_num, exist=exist, step=nt )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                rh_org(2,i,j) = work2d(i,j)
             enddo
             enddo
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                rh_org(2,i,j) = UNDEF
             enddo
             enddo
          end if
       else
          call read2d( start(2:), count(2:), work2d(:,:), "Q2", file_id_atm, basename_num, exist=exist, step=nt )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                qv_org(2,i,j) = work2d(i,j)
             enddo
             enddo
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                qv_org(2,i,j) = UNDEF
             enddo
             enddo
          end if

       end if

       ! topo
       call read2d( start(2:), count(2:), work2d(:,:), "topo", file_id_atm, basename_num, exist=exist, step=nt )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             cz_org(2,i,j) = work2d(i,j)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             cz_org(2,i,j) = UNDEF
          enddo
          enddo
       end if

    end if

    return
  end subroutine ParentAtmosInputGrADS

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupGrADS( &
       ldims,            &
       timelen,          &
       lon_all, lat_all, &
       basename,         &
       basename_num      )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid, &
       FILE_GrADS_isOneD
    implicit none
    integer,           intent(out) :: ldims(3)
    integer,           intent(out) :: timelen
    real(RP), pointer, intent(out) :: lon_all(:,:)
    real(RP), pointer, intent(out) :: lat_all(:,:)

    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: basename_num

    character(len=7) :: vname
    integer :: vid
    integer :: shape(2)
    logical :: exist
    !---------------------------------------------------------------------------

    LOG_INFO("ParentLandSetupGrADS",*) 'Real Case/Land Input File Type: GrADS format'


    if ( basename == "" ) then
       LOG_ERROR("ParentLandSetupGrADS",*) '"BASENAME_ORG" is not specified in "PARAM_MKINIT_REAL_LAND"!', trim(basename)
       call PRC_abort
    endif

    ! open
    call FILE_GrADS_open( basename,   & ! (in)
                          file_id_lnd ) ! (out)

    ! get shape
    call FILE_GrADS_varid( file_id_lnd, "LAND_SFC_TEMP", & ! (in)
                           vid                           ) ! (out)
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_lnd, "LAND_TEMP", & ! (in)
                              vid                       ) ! (out)
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_lnd, "SFC_TEMP", & ! (in)
                              vid                      ) ! (out)
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_lnd, "STEMP", & ! (in)
                              vid                   ) ! (out)
    end if
    if ( vid < 0 ) then
       LOG_ERROR("ParentLandSetupGrADS",*) '"LAND_SFC_TEMP", "LAND_TEMP", "SFC_TEMP" or "STEMP" is necessary'
       call PRC_abort
    end if
    call FILE_GrADS_get_shape( file_id_lnd, vid, & ! (in)
                               ldims(:)          ) ! (out)


    ! get lon, lat
    allocate( lon_all(ldims(2), ldims(3)) )
    allocate( lat_all(ldims(2), ldims(3)) )


    call FILE_GrADS_varid( file_id_lnd, "lon_sfc", & ! (in)
                           vid                     ) ! (out)
    if ( vid > 0 ) then
       vname = "lon_sfc"
    else
       call FILE_GrADS_varid( file_id_lnd, "lon", & ! (in)
                              vid                 ) ! (out)
       if ( vid < 0 ) then
          LOG_ERROR("ParentLandSetupGrADS",*) '"lon_sfc" or "lon" is necessary'
          call PRC_abort
       end if
       call FILE_GrADS_get_shape( file_id_lnd, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( FILE_GrADS_isOneD( file_id_lnd, vid ) ) then
          if ( ldims(2) .ne. shape(1) .and. shape(1) .ne. -1 ) then
             LOG_ERROR("ParentLandSetupGrADS",*) 'dimension of "lon" is different! ', ldims(2), shape(1)
             call PRC_abort
          end if
       else
          if ( ldims(2) .ne. shape(1) .or. ldims(3) .ne. shape(2) ) then
             LOG_ERROR("ParentLandSetupGrADS",*) 'dimension of "lon" is different! ', ldims(2), shape(1), ldims(3), shape(2)
             call PRC_abort
          end if
       end if
       vname = "lon"
    end if
    call read2d( (/1,1/), ldims(2:3), lon_all(:,:), vname, file_id_lnd, basename_num, oneD=1 )
    lon_all(:,:) = lon_all(:,:) * D2R


    call FILE_GrADS_varid( file_id_lnd, "lat_sfc", & ! (in)
                           vid                    ) ! (out)
    if ( vid > 0 ) then
       vname = "lat_sfc"
    else
       call FILE_GrADS_varid( file_id_lnd, "lat", & ! (in)
                              vid                 ) ! (out)
       if ( vid < 0 ) then
          LOG_ERROR("ParentLandSetupGrADS",*) '"lat_sfc" or "lat" is necessary'
          call PRC_abort
       end if
       call FILE_GrADS_get_shape( file_id_lnd, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( FILE_GrADS_isOneD( file_id_lnd, vid ) ) then
          if ( ldims(3) .ne. shape(1) .and. shape(1) .ne. -1 ) then
             LOG_ERROR("ParentLandSetupGrADS",*) 'dimension of "lat" is different! ', ldims(2), shape(1)
             call PRC_abort
          end if
       else
          if ( ldims(2) .ne. shape(1) .or. ldims(3) .ne. shape(2) ) then
             LOG_ERROR("ParentLandSetupGrADS",*) 'dimension of "lat" is different! ', ldims(2), shape(1), ldims(3), shape(2)
             call PRC_abort
          end if
       end if
       vname = "lat"
    end if
    call read2d( (/1,1/), ldims(2:3), lat_all(:,:), vname, file_id_atm, basename_num, oneD=2 )
    lat_all(:,:) = lat_all(:,:) * D2R


    ! tentative
    timelen = 1


    return
  end subroutine ParentLandSetupGrADS

  subroutine ParentLandInputGrADS( &
       KA_org, KS_org, KE_org, &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       tg_org,             & ! (out)
       strg_org,           & ! (out)
       smds_org,           & ! (out)
       lst_org,            & ! (out)
       lz_org,             & ! (out)
       topo_org,           & ! (out)
       lmask_org,          & ! (out)
       use_waterratio,     & ! (out)
       ldims,              & ! (in)
       basename_num,       & ! (in)
       use_file_landwater, & ! (in)
       nt                  ) ! (in)
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    use scale_file_grads, only: &
       FILE_GrADS_get_shape, &
       FILE_GrADS_read
    implicit none
    integer, intent(in) :: KA_org, KS_org, KE_org
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out) :: tg_org   (KA_org,IA_org,JA_org)
    real(RP), intent(out) :: strg_org (KA_org,IA_org,JA_org)
    real(RP), intent(out) :: smds_org (KA_org,IA_org,JA_org)
    real(RP), intent(out) :: lst_org  (IA_org,JA_org)
    real(RP), intent(out) :: lz_org   (KA_org)
    real(RP), intent(out) :: topo_org (IA_org,JA_org)
    real(RP), intent(out) :: lmask_org(IA_org,JA_org)
    logical,  intent(out) :: use_waterratio

    integer,          intent(in) :: ldims(3)
    character(len=*), intent(in) :: basename_num
    logical,          intent(in) :: use_file_landwater ! use land water data from files
    integer,          intent(in) :: nt

    integer :: start(3), count(3), shape(2)
    logical :: exist


    integer :: i, j, k
    !---------------------------------------------------------------------------

    start(:) = (/KS_org,IS_org,JS_org/)
    count(:) = (/KA_org,IA_org,JA_org/)

    ! lsmask
    call read2d( start(2:), count(2:), lmask_org(:,:), "lsmask", file_id_lnd, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          lmask_org(i,j) = UNDEF
       end do
       end do
    end if

    ! llev
    call FILE_GrADS_get_shape( file_id_lnd, "llev", & ! (in)
                               shape(:)             ) ! (out)
    if ( ldims(1) .ne. shape(1) )then
       LOG_ERROR("ParentLandInputGrADS",*) '"nz" must be equal to nz of "STEMP" for llev. :', shape(1), ldims(1)
       call PRC_abort
    endif
    call FILE_GrADS_read( file_id_lnd, "llev", & ! (in)
                          lz_org(:)            ) ! (out)

    ! LAND_TEMP, STEMP
    call read3d( start(:), count(:), tg_org(:,:,:), "LAND_TEMP", file_id_lnd, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call read3d( start(:), count(:), tg_org(:,:,:), "STEMP", file_id_lnd, basename_num, exist=exist )
    end if
    if ( .not. exist ) then
       LOG_ERROR("ParentAtmosInputGrADS",*) '"LAND_TEMP" or "STEMP" is necessary'
       call PRC_abort
    endif

    if ( use_file_landwater ) then
       ! LAND_WATER, SMOISVC, SMOISDS
       call read3d( start(:), count(:), strg_org(:,:,:), "LAND_WATER", file_id_lnd, basename_num, exist=exist, step=nt )
       if ( .not. exist ) then
          call read3d( start(:), count(:), strg_org(:,:,:), "SMOISVC", file_id_lnd, basename_num, exist=exist, step=nt )
       end if
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 2, KA_org
             if ( strg_org(k,i,j) == UNDEF ) strg_org(k,i,j) = strg_org(k-1,i,j)
          end do
          end do
          end do
          use_waterratio = .false.
       else
          call read3d( start(:), count(:), smds_org(:,:,:), "SMOISDS", file_id_lnd, basename_num, exist=exist, step=nt )
          if ( exist ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 2, KA_org
                if ( smds_org(k,i,j) == UNDEF ) smds_org(k,i,j) = smds_org(k-1,i,j)
             end do
             end do
             end do
             use_waterratio = .true.
          else
             LOG_ERROR("ParentAtmosInputGrADS",*) '"LAND_WATER", "SMOISVC", or "SMOISDS" is necessary'
             call PRC_abort
          end if
       endif

    end if ! use_file_landwater


    ! LAND_SFC_TEMP, SFC_TEMP, SKINT
    call read2d( start(2:), count(2:), lst_org(:,:), "LAND_SFC_TEMP", file_id_lnd, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call read2d( start(2:), count(2:), lst_org(:,:), "SFC_TEMP", file_id_lnd, basename_num, exist=exist )
    end if
    if ( .not. exist ) then
       call read2d( start(2:), count(2:), lst_org(:,:), "SKINT", file_id_lnd, basename_num, exist=exist )
    end if
    if ( .not. exist ) then
       LOG_ERROR("ParentAtmosInputGrADS",*) '"LAND_SFC_TEMP", "SFC_TEMP", or "SKINT" is necessary'
       call PRC_abort
    endif


    ! topo_sfc, topo
    call read2d( start(2:), count(2:), topo_org(:,:), "topo_sfc", file_id_lnd, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call FILE_GrADS_get_shape( file_id_lnd, "topo", & ! (in)
                                  shape(:)             ) ! (out)
       if ( ldims(2).ne.shape(1) .or. ldims(3).ne.shape(2) ) then
          LOG_WARN("ParentLandInputGrADS",*) 'namelist of "topo_sfc" is not found in grads namelist!'
          LOG_WARN_CONT(*) 'dimension of "topo" is different! ', ldims(2), shape(1), ldims(3), shape(2)
       else
          call read2d( start(2:), count(2:), topo_org(:,:), "topo", file_id_lnd, basename_num, exist=exist )
       end if
    end if

    return
  end subroutine ParentLandInputGrADS

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupGrADS( &
       odims,            & ! (out)
       timelen,          & ! (out)
       lon_all, lat_all, & ! (out)
       basename,         & ! (in)
       basename_num      ) ! (in)
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_varid, &
       FILE_GrADS_get_shape, &
       FILE_GrADS_isOneD
    implicit none

    integer,           intent(out) :: odims(2)
    integer,           intent(out) :: timelen
    real(RP), pointer, intent(out) :: lon_all(:,:)
    real(RP), pointer, intent(out) :: lat_all(:,:)

    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: basename_num

    character(len=7) :: vname
    integer :: vid
    integer :: shape(2)
    logical :: exist
    !---------------------------------------------------------------------------

    LOG_INFO("ParentOceanSetupGrADS",*) 'Real Case/Ocean Input File Type: GrADS format'

    if ( basename == "" ) then
       LOG_ERROR("ParentOceanSetupGrADS",*) '"BASENAME_ORG" is not specified in "PARAM_MKINIT_REAL_OCEAN"!', trim(basename)
       call PRC_abort
    endif
    call FILE_GrADS_open( basename,   & ! (in)
                          file_id_ocn ) ! (out)

    call FILE_GrADS_varid( file_id_ocn, "OCEAN_SFC_TEMP", & ! (in)
                           vid                            ) ! (out)
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_ocn, "SST", & ! (in)
                              vid                 ) ! (out)
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_ocn, "SFC_TEMP", & ! (in)
                              vid                      ) ! (out)
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_ocn, "SKINT", & ! (in)
                              vid                      ) ! (out)
    end if
    if ( vid < 0 ) then
       LOG_ERROR("ParentOceanSetupGrADS",*) '"OCEAN_SFC_TEMP", "SST", "SFC_TEMP", or "SKINT" is necessary'
       call PRC_abort
    end if

    call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                               odims(:)          ) ! (out)


    ! get lon, lat

    allocate( lon_all(odims(1), odims(2)) )
    allocate( lat_all(odims(1), odims(2)) )


    call FILE_GrADS_varid( file_id_ocn, "lon_sst", & ! (in)
                           vid                   ) ! (out)
    if ( vid > 0 ) then
       vname = "lon_sst"
    else
       call FILE_GrADS_varid( file_id_ocn, "lon_sfc", & ! (in)
                              vid                     ) ! (out)
       if ( vid > 0 ) then
          call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                     shape(:)          ) ! (out)
          if ( FILE_GrADS_isOneD( file_id_lnd, vid ) ) then
             if ( odims(1) .eq. shape(1) .or. shape(1) .eq. -1 ) then
                vname = "lon_sfc"
             else
                vid = -1
             end if
          else
             if ( odims(1) .eq. shape(1) .and. odims(2) .eq. shape(2) ) then
                vname = "lon_sfc"
             else
                vid = -1
             end if
          end if
       end if
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_ocn, "lon", & ! (in)
                              vid                 ) ! (out)
       if ( vid < 0 ) then
          LOG_ERROR("ParentLandSetupGrADS",*) '"lon_sst", "lon_sfc", or "lon" is necessary'
          call PRC_abort
       end if
       call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( FILE_GrADS_isOneD( file_id_ocn, vid ) ) then
          if ( odims(1) .eq. shape(1) .or. shape(1) .eq. -1 ) then
             vname = "lon"
          else
             vid = -1
          end if
       else
          if ( odims(1) .eq. shape(1) .and. odims(2) .eq. shape(2) ) then
             vname = "lon"
          else
             vid = -1
          end if
       end if
       if ( vid < 0 ) then
          LOG_ERROR("ParentOceanSetupGrADS",*) 'dimension of "lon_sfc" and "lon" is different! ', odims(:), shape(:)
          call PRC_abort
       end if
    end if
    call read2d( (/1,1/), odims(:), lon_all(:,:), vname, file_id_ocn, basename_num, oneD=1 )
    lon_all(:,:) = lon_all(:,:) * D2R


    call FILE_GrADS_varid( file_id_ocn, "lat_sst", & ! (in)
                           vid                     ) ! (out)
    if ( vid > 0 ) then
       vname = "lat_sst"
    else
       call FILE_GrADS_varid( file_id_ocn, "lat_sfc", & ! (in)
                              vid                     ) ! (out)
       if ( vid > 0 ) then
          call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                     shape(:)          ) ! (out)
          if ( FILE_GrADS_isOneD( file_id_lnd, vid ) ) then
             if ( odims(2) .eq. shape(1) .or. shape(1) .eq. -1 ) then
                vname = "lat_sfc"
             else
                vid = -1
             end if
          else
             if ( odims(1) .eq. shape(1) .and. odims(2) .eq. shape(2) ) then
                vname = "lat_sfc"
             else
                vid = -1
             end if
          end if
       end if
    end if
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id_lnd, "lat", & ! (in)
                              vid                 ) ! (out)
       if ( vid < 0 ) then
          LOG_ERROR("ParentLandSetupGrADS",*) '"lat_sst", "lat_sfc", or "lat" is necessary'
          call PRC_abort
       end if
       call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( FILE_GrADS_isOneD( file_id_lnd, vid ) ) then
          if ( odims(2) .eq. shape(1) .or. shape(1) .eq. -1 ) then
             vname = "lat"
          else
             vid = -1
          end if
       else
          if ( odims(1) .eq. shape(1) .and. odims(2) .eq. shape(2) ) then
             vname = "lat"
          else
             vid = -1
          end if
       end if
       if ( vid < 0 ) then
          LOG_ERROR("ParentOceanSetupGrADS",*) 'dimension of "lat_sfc" and "lat" is different! ', odims(:), shape(:)
          call PRC_abort
       end if
    end if
    call read2d( (/1,1/), odims(:), lat_all(:,:), vname, file_id_ocn, basename_num, oneD=2 )
    lat_all(:,:) = lat_all(:,:) * D2R


    ! tentative
    timelen = 0

    return
  end subroutine ParentOceanSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputGrADS( &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       tw_org,       &
       sst_org,      &
       omask_org,    &
       basename_num, &
       odims,        &
       nt            )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    use scale_file_grads, only: &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid,     &
       FILE_GrADS_read
    implicit none
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out) :: tw_org   (IA_org,JA_org)
    real(RP), intent(out) :: sst_org  (IA_org,JA_org)
    real(RP), intent(out) :: omask_org(IA_org,JA_org)

    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: odims(2)
    integer,          intent(in)  :: nt

    integer :: start(2), count(2), shape(2)
    integer :: vid
    logical :: exist
    integer :: i, j
    !---------------------------------------------------------------------------


    start(:) = (/IS_org,JS_org/)
    count(:) = (/IA_org,JA_org/)

    ! lsmask
    call read2d( start(:), count(:), omask_org(:,:), "lsmask_sst", file_id_ocn, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call FILE_GrADS_varid( file_id_ocn, "lsmask", & ! (in)
                              vid                     ) ! (out)
       call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( odims(1) .ne. shape(1) .or. odims(2) .ne. shape(2) ) then
          LOG_WARN("ParentOceanInputGrADS",*) 'dimension of lsmask is different. not use'
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             omask_org(i,j) = UNDEF
          end do
          end do
       else
          call read2d( start(:), count(:), omask_org(:,:), "lsmask", file_id_ocn, basename_num, exist=exist, step=nt )
       end if
    end if


    ! OCEAN_SFC_TEMP, SST, SFC_TEMP, SKINT
    call read2d( start(:), count(:), sst_org(:,:), "OCEAN_SFC_TEMP", file_id_ocn, basename_num, exist=exist, step=nt )
    if ( .not. exist ) then
       call read2d( start(:), count(:), sst_org(:,:), "SST", file_id_ocn, basename_num, exist=exist, step=nt )
    end if
    if ( .not. exist ) then
       call FILE_GrADS_varid( file_id_ocn, "SFC_TEMP", & ! (in)
                              vid                     ) ! (out)
       if ( vid > 0 ) then
          call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                     shape(:)          ) ! (out)
          if ( odims(1).eq.shape(1) .and. odims(2).eq.shape(2) ) then
             call read2d( start(:), count(:), sst_org(:,:), "SFC_TEMP", file_id_ocn, basename_num, exist=exist, step=nt )
          else
             exist = .false.
          end if
       else
          exist = .false.
       end if
    end if
    if ( .not. exist ) then
       call FILE_GrADS_varid( file_id_ocn, "SKINT", & ! (in)
                              vid                     ) ! (out)
       call FILE_GrADS_get_shape( file_id_ocn, vid, & ! (in)
                                  shape(:)          ) ! (out)
       if ( odims(1).eq.shape(1) .and. odims(2).eq.shape(2) ) then
          call read2d( start(:), count(:), sst_org(:,:), "SKINT", file_id_ocn, basename_num, step=nt )
       else
          LOG_ERROR("ParentOceanInputGrADS",*) 'dimension of "SFC_TEMP" and/or "SKINT" is different'
          call PRC_abort
       end if
    end if

    tw_org = sst_org

    return
  end subroutine ParentOceanInputGrADS

  ! private

  subroutine read2d( start, count, data, name, fid, postfix, exist, oneD, step )
    use scale_file_grads, only: &
       FILE_GrADS_varid, &
       FILE_GrADS_read, &
       FILE_GrADS_isOneD
    integer,  intent(in)  :: start(2)
    integer,  intent(in)  :: count(2)
    real(RP), intent(out) :: data(count(1),count(2))
    character(len=*), intent(in) :: name
    integer,          intent(in) :: fid
    character(len=*), intent(in) :: postfix

    logical, intent(out), optional :: exist
    integer, intent(in),  optional :: oneD
    integer, intent(in),  optional :: step

    integer :: vid
    real(RP), allocatable :: v1d(:)
    integer :: oneD_
    integer :: i, j

    call FILE_GrADS_varid( fid, name, & ! (in)
                           vid        ) ! (out)
    if ( vid < 0 ) then
       if ( present(exist) ) then
          exist = .false.
          return
       end if
       LOG_ERROR("read2d",*) '"', trim(name), '" is required'
       call PRC_abort
    end if
    if ( present(exist) ) exist = .true.

    if( FILE_GrADS_isOneD( fid, vid ) ) then
       if ( present(oneD) ) then
          oneD_ = oneD
       else
          oneD_ = 1
       end if
       allocate( v1d(count(oneD_)) )
       call FILE_GrADS_read( fid, vid,                 & ! (in)
                             v1d(:),                   & ! (out)
                             step=step,                & ! (in)
                             start=start(oneD_:oneD_), & ! (in)
                             count=count(oneD_:oneD_)  ) ! (in)
       if ( oneD_ == 1 ) then
          !$omp parallel do
          do j = 1, count(2)
          do i = 1, count(1)
             data(i,j) = v1d(i)
          enddo
          enddo
       else
          !$omp parallel do
          do j = 1, count(2)
          do i = 1, count(1)
             data(i,j) = v1d(j)
          enddo
          enddo
       end if
       deallocate( v1d )
    else
       call FILE_GrADS_read( fid, vid,       & ! (in)
                             data(:,:),      & ! (out)
                             step=step,      & ! (in)
                             start=start(:), & ! (in)
                             count=count(:), & ! (in)
                             postfix=postfix ) ! (in)
    end if

    return
  end subroutine read2d

  subroutine read3d( start, count, data, name, fid, postfix, exist, step )
    use scale_file_grads, only: &
       FILE_GrADS_varid, &
       FILE_GrADS_read, &
       FILE_GrADS_isOneD
    integer,  intent(in)  :: start(3)
    integer,  intent(in)  :: count(3)
    real(RP), intent(out) :: data(count(1),count(2),count(3))
    character(len=*), intent(in) :: name
    integer,          intent(in) :: fid
    character(len=*), intent(in) :: postfix

    logical, intent(out), optional :: exist
    integer, intent(in), optional :: step

    integer :: vid
    real(RP), allocatable :: v1d(:)
    integer :: k, i, j

    call FILE_GrADS_varid( fid, name, & ! (in)
                           vid        ) ! (out)
    if ( vid < 0 ) then
       if ( present(exist) ) then
          exist = .false.
          return
       end if
       LOG_ERROR("read3d",*) '"', trim(name), '" is required'
       call PRC_abort
    end if
    if ( present(exist) ) exist = .true.

    if( FILE_GrADS_isOneD( fid, vid ) ) then
       allocate( v1d(count(1)) )
       call FILE_GrADS_read( fid, vid,         & ! (in)
                             v1d(:),           & ! (out)
                             step=step,        &
                             start=start(1:1), & ! (in)
                             count=count(1:1)  ) ! (in)
       !$omp parallel do collapse(2)
       do j = 1, count(3)
       do i = 1, count(2)
       do k = 1, count(1)
          data(k,i,j) = v1d(k)
       enddo
       enddo
       enddo
       deallocate( v1d )
    else
       call FILE_GrADS_read( fid, vid,       & ! (in)
                             data(:,:,:),    & ! (out)
                             step=step,      & ! (in)
                             start=start(:), & ! (in)
                             count=count(:), & ! (in)
                             postfix=postfix ) ! (in)
    end if

    return
  end subroutine read3d

end module mod_realinput_grads
