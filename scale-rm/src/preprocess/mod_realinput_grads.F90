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
  public :: ParentAtmosOpenGrADS
  public :: ParentAtmosInputGrADS
  public :: ParentLandSetupGrADS
  public :: ParentLandInputGrADS
  public :: ParentOceanSetupGrADS
  public :: ParentOceanOpenGrADS
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
  integer,  parameter    :: num_item_list_atom  = 25
  integer,  parameter    :: num_item_list_land  = 12
  integer,  parameter    :: num_item_list_ocean = 10
  character(len=H_SHORT) :: item_list_atom (num_item_list_atom)
  character(len=H_SHORT) :: item_list_land (num_item_list_land)
  character(len=H_SHORT) :: item_list_ocean(num_item_list_ocean)
  data item_list_atom  /'lon','lat','plev','DENS','U','V','W','T','HGT','QV','QC','QR','QI','QS','QG','RH', &
                        'MSLP','PSFC','U10','V10','T2','Q2','RH2','topo','RN222' /
  data item_list_land  /'lsmask','lon','lat','lon_sfc','lat_sfc','llev', &
                        'STEMP','SMOISVC','SMOISDS','SKINT','topo','topo_sfc' /
  data item_list_ocean /'lsmask','lsmask_sst','lon','lat','lon_sfc','lat_sfc','lon_sst','lat_sst','SKINT','SST'/

  integer,  parameter    :: num_item_list = 25 ! max of num_item_list_(atom|land|ocean)
  integer                :: var_id(num_item_list,3) ! 1:atom, 2:land, 3:ocean

  integer,  parameter   :: Ia_lon    = 1
  integer,  parameter   :: Ia_lat    = 2
  integer,  parameter   :: Ia_p      = 3  ! Pressure (Pa)
  integer,  parameter   :: Ia_dens   = 4
  integer,  parameter   :: Ia_u      = 5
  integer,  parameter   :: Ia_v      = 6
  integer,  parameter   :: Ia_w      = 7
  integer,  parameter   :: Ia_t      = 8
  integer,  parameter   :: Ia_hgt    = 9  ! Geopotential height (m)
  integer,  parameter   :: Ia_qv     = 10
  integer,  parameter   :: Ia_qc     = 11
  integer,  parameter   :: Ia_qr     = 12
  integer,  parameter   :: Ia_qi     = 13
  integer,  parameter   :: Ia_qs     = 14
  integer,  parameter   :: Ia_qg     = 15
  integer,  parameter   :: Ia_rh     = 16 ! Percentage (%)
  integer,  parameter   :: Ia_slp    = 17 ! Sea level pressure (Pa)
  integer,  parameter   :: Ia_ps     = 18 ! Surface pressure (Pa)
  integer,  parameter   :: Ia_u10    = 19
  integer,  parameter   :: Ia_v10    = 20
  integer,  parameter   :: Ia_t2     = 21
  integer,  parameter   :: Ia_q2     = 22
  integer,  parameter   :: Ia_rh2    = 23 ! Percentage (%)
  integer,  parameter   :: Ia_topo   = 24
  integer,  parameter   :: Ia_rn222  = 25

  integer,  parameter   :: Il_lsmask  = 1
  integer,  parameter   :: Il_lon     = 2
  integer,  parameter   :: Il_lat     = 3
  integer,  parameter   :: Il_lon_sfc = 4
  integer,  parameter   :: Il_lat_sfc = 5
  integer,  parameter   :: Il_lz      = 6  ! Level(depth) of stemp & smois (m)
  integer,  parameter   :: Il_stemp   = 7
  integer,  parameter   :: Il_smoisvc = 8  ! soil moisture (vormetric water content)
  integer,  parameter   :: Il_smoisds = 9  ! soil moisture (degree of saturation)
  integer,  parameter   :: Il_skint   = 10
  integer,  parameter   :: Il_topo    = 11
  integer,  parameter   :: Il_topo_sfc= 12

  integer,  parameter   :: Io_lsmask     = 1
  integer,  parameter   :: Io_lsmask_sst = 2
  integer,  parameter   :: Io_lon        = 3
  integer,  parameter   :: Io_lat        = 4
  integer,  parameter   :: Io_lon_sfc    = 5
  integer,  parameter   :: Io_lat_sfc    = 6
  integer,  parameter   :: Io_lon_sst    = 7
  integer,  parameter   :: Io_lat_sst    = 8
  integer,  parameter   :: Io_skint      = 9
  integer,  parameter   :: Io_sst        = 10

  character(len=H_SHORT) :: upper_qv_type = "ZERO" !< how qv is given at higher level than outer model
                                                   !< "ZERO": 0
                                                   !< "COPY": copy values from the highest level of outer model


  ! atmos data
  integer :: outer_nx     = -1
  integer :: outer_ny     = -1
  integer :: outer_nz     = -1 ! number of atmos layers
  integer :: outer_nl     = -1 ! number of land layers
  ! surface data
  integer :: outer_nx_sfc = -1
  integer :: outer_ny_sfc = -1
  ! sst data
  integer :: outer_nx_sst = -1
  integer :: outer_ny_sst = -1

  integer :: file_id

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupGrADS( &
      dims,    &
      basename )
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid
    implicit none
    integer,          intent(out) :: dims(6)
    character(len=*), intent(in)  :: basename

    namelist / PARAM_MKINIT_REAL_GrADS / &
        upper_qv_type

    integer                :: shape(3)
    integer                :: ielem
    character(len=H_SHORT) :: item

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentAtmosSetupGrADS",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_GrADS,iostat=ierr)

    if( ierr > 0 ) then
       LOG_ERROR("ParentAtmosSetupGrADS",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_GrADS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_GrADS)


    if ( basename == "" ) then
       LOG_ERROR("ParentAtmosSetupGrADS",*) '"BASENAME_ORG" is not specified in "PARAM_MKINIT_ATMOS_GRID_CARTESC_REAL_ATMOS"!', trim(basename)
       call PRC_abort
    endif

    call FILE_GrADS_open( basename, & ! (in)
                          file_id   ) ! (out)

    call FILE_GrADS_get_shape( file_id, "U", & ! (in)
                               shape(:)      ) ! (out)

    ! full level
    dims(1) = shape(1) ! bottom_top
    dims(2) = shape(2) ! west_east
    dims(3) = shape(3) ! south_north
    ! half level
    dims(4) = shape(1) ! bottom_top_stag
    dims(5) = shape(2) ! west_east for 2dim data
    dims(6) = shape(3) ! south_north for 2dim data


    ! check existence
    ! var_id > 0 : exist
    ! var_id < 0 : not exist
    do ielem = 1, num_item_list_atom
       item  = item_list_atom(ielem)
       call FILE_GrADS_varid( file_id, item,  & ! (in)
                              var_id(ielem,1) ) ! (out)
    end do

    ! check necessary data
    do ielem = 1, num_item_list_atom
       item  = item_list_atom(ielem)
       !--- check data
       select case(item)
       case('DENS','W','QC','QR','QI','QS','QG','MSLP','PSFC','U10','V10','T2','topo','RN222')
          if ( var_id(ielem,1) < 0 ) then
             LOG_WARN("ParentAtmosSetupGrADS",*) trim(item),' is not found & will be estimated.'
          endif
          cycle
       case('QV', 'RH')
          if ( var_id(Ia_qv,1) < 0 ) then
             if( var_id(Ia_rh,1) > 0 ) then
                if ( var_id(Ia_t,1) < 0 .or. var_id(Ia_p,1) < 0 ) then
                   LOG_ERROR("ParentAtmosSetupGrADS",*) 'Temperature and pressure are required to convert from RH to QV ! '
                   call PRC_abort
                endif
             else
                LOG_ERROR("ParentAtmosSetupGrADS",*) 'Not found in grads namelist! : QV and RH'
                call PRC_abort
             endif
          else
             var_id(Ia_rh,1) = -1
          endif
          cycle
       case('Q2', 'RH2')
          if ( var_id(Ia_q2,1) < 0 ) then
             if ( var_id(Ia_rh2,1) > 0 ) then
                if ( var_id(Ia_t2,1) < 0 .or.  var_id(Ia_ps,1) < 0 ) then
                   LOG_WARN("ParentAtmosSetupGrADS",*) 'T2 and PSFC are required to convert from RH2 to Q2 !'
                   LOG_INFO_CONT(*)                    'Q2 will be copied from data at above level.'
                   var_id(Ia_rh2,1) = -1
                endif
             else
                LOG_WARN("ParentAtmosSetupGrADS",*) 'Q2 and RH2 are not found, Q2 will be estimated.'
             endif
          else
             var_id(Ia_rh2,1) = -1
          end if
          cycle
       case default ! lon, lat, plev, U, V, T, HGT
          if ( var_id(ielem,1) < 0 ) then
             LOG_ERROR("ParentAtmosSetupGrADS",*) 'Not found in grads namelist! : ',trim(item_list_atom(ielem))
             call PRC_abort
          endif
       end select

    end do

    return
  end subroutine ParentAtmosSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosOpenGrADS( &
       lon_org,  &
       lat_org,  &
       basename_num,  &
       dims )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file_grads, only: &
       FILE_GrADS_isOneD, &
       FILE_GrADS_read
    implicit none

    real(RP),         intent(out) :: lon_org(:,:)
    real(RP),         intent(out) :: lat_org(:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(6)

    real(RP) :: lon1D(dims(2)), lat1D(dims(3))

    character(len=H_SHORT) :: item

    integer  :: i, j, ielem
    !---------------------------------------------------------------------------

    !--- read grads data
    loop_InputAtmosGrADS : do ielem = 1, num_item_list_atom

       if ( var_id(ielem,1) < 0 ) cycle

       item = item_list_atom(ielem)

       ! read data
       select case(item)
       case("lon")

          if( FILE_GrADS_isOneD( file_id, var_id(ielem,1) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   lon1d(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, dims(3)
             do i = 1, dims(2)
                lon_org(i,j) = lon1d(i) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   lon_org(:,:),             & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, dims(3)
             do i = 1, dims(2)
                lon_org(i,j) = lon_org(i,j) * D2R
             enddo
             enddo
          end if

       case("lat")

          if( FILE_GrADS_isOneD( file_id, var_id(ielem,1) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   lat1d(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, dims(3)
             do i = 1, dims(2)
                lat_org(i,j) = lat1d(j) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   lat_org(:,:),             & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, dims(3)
             do i = 1, dims(2)
                lat_org(i,j) = lat_org(i,j) * D2R
             enddo
             enddo
          end if

       end select
    enddo loop_InputAtmosGrADS

    return
  end subroutine ParentAtmosOpenGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputGrADS( &
       velz_org, &
       velx_org, &
       vely_org, &
       pres_org, &
       dens_org, &
       temp_org, &
       qv_org,   &
       qhyd_org, &
       RN222_org,&
       cz_org,   &
       basename_num,  &
       sfc_diagnoses, &
       under_sfc,     &
       KA_org, &
       KS_org, &
       KE_org, &
       IA_org, &
       IS_org, &
       IE_org, &
       JA_org, &
       JS_org, &
       JE_org, &
       dims, &
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
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG
    use scale_atmos_saturation, only: &
       psat => ATMOS_SATURATION_psat_liq
    use scale_file_grads, only: &
       FILE_GrADS_isOneD,    &
       FILE_GrADS_get_shape, &
       FILE_GrADS_read
    implicit none


    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: dens_org(:,:,:)
    real(RP),         intent(out) :: temp_org(:,:,:)
    real(RP),         intent(out) :: qv_org  (:,:,:)
    real(RP),         intent(out) :: qhyd_org(:,:,:,:)
    real(RP),         intent(out) :: RN222_org(:,:,:)
    real(RP),         intent(out) :: cz_org(:,:,:)
    character(len=*), intent(in)  :: basename_num
    logical,          intent(in)  :: sfc_diagnoses
    logical,          intent(in)  :: under_sfc
    integer,          intent(in)  :: KA_org
    integer,          intent(in)  :: KS_org
    integer,          intent(in)  :: KE_org
    integer,          intent(in)  :: IA_org
    integer,          intent(in)  :: IS_org
    integer,          intent(in)  :: IE_org
    integer,          intent(in)  :: JA_org
    integer,          intent(in)  :: JS_org
    integer,          intent(in)  :: JE_org
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: nt

    character(len=H_SHORT) :: item

    integer  :: lm_layer( IA_org, JA_org )

    integer  :: dummy = 1
    real(RP) :: work( dims(1), dims(2), dims(3) )

    logical  :: pressure_coordinates
    real(RP) :: p_sat, qm, dz
    real(RP) :: rh( IA_org, JA_org )
    real(RP) :: Rtot

    integer  :: shape(3)
    integer  :: i, j, k, iq, ielem
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(3)
    do j = 1, JA_org
    do i = 1, IA_org
    do k = 1, KA_org
       dens_org (k,i,j)   = UNDEF
       pres_org (k,i,j)   = UNDEF
       velz_org (k,i,j)   = 0.0_RP
       qv_org   (k,i,j)   = 0.0_RP
       qhyd_org (k,i,j,:) = 0.0_RP
       RN222_org(k,i,j )  = 0.0_RP
    end do
    end do
    end do

    !--- read grads data
    loop_InputAtmosGrADS : do ielem = 1, num_item_list_atom

       if ( var_id(ielem,1) < 0 ) cycle

       item = item_list_atom(ielem)

       ! read data
       select case(item)
       case("plev")

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item), shape(1), KA_org-2
             call PRC_abort
          endif

          if( FILE_GrADS_isOneD( file_id, var_id(ielem,1) ) ) then
             pressure_coordinates = .true. ! use pressure coordinate in the input data
             call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                        shape(:)                  ) ! (out)
             if ( KA_org-2 .ne. shape(1) ) then
                LOG_ERROR("ParentAtmosInputGrADS",*) 'lnum must be same as the nz for plev! ',shape(1), KA_org-2
                call PRC_abort
             endif
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(:,dummy,dummy)       ) ! (out)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                pres_org(k+2,i,j) = work(k,dummy,dummy)
             enddo
             enddo
             enddo
          else
             pressure_coordinates = .false.
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(:,:,:),              & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                pres_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
             enddo
          endif
       case('DENS')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:,:,:),              & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             dens_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

       case('U')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:,:,:),              & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             velx_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
             velx_org(1:2,i,j) = 0.0_RP
          enddo
          enddo

       case('V')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:,:,:),              & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             vely_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
             vely_org(1:2,i,j) = 0.0_RP
          enddo
          enddo

       case('W')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:,:,:),              & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             velz_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
             velz_org(1:2,i,j) = 0.0_RP
          enddo
          enddo

       case('T')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:,:,:),              & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             temp_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

       case('HGT')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( KA_org-2 .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to the default nz for ',trim(item),'. nz:',shape(1),'> outer_nz:',KA_org-2
             call PRC_abort
          endif

          if( FILE_GrADS_isOneD( file_id, var_id(ielem,1) ) ) then
             call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                        shape(:)                  ) ! (out)
             if ( KA_org-2 .ne. shape(1) ) then
                LOG_ERROR("ParentAtmosInputGrADS",*) 'lnum must be same as the nz for HGT! ',KA_org-2, shape(1)
                call PRC_abort
             endif
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(:,dummy,dummy)       ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
                cz_org(1,i,j) = 0.0_RP
                do k = 1, KA_org-2
                   cz_org(k+2,i,j) = work(k,dummy,dummy)
                enddo
             enddo
             enddo
          else
             pressure_coordinates = .false.
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(:,:,:),              & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                cz_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
             enddo
                cz_org(1,i,j) = 0.0_RP
             enddo
             enddo
          endif

       case('QV')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qv_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          if( KA_org-2 > shape(1) ) then
             select case( upper_qv_type )
             case("COPY")
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = shape(1)+1, KA_org-2
                   qv_org(k+2,i,j) = qv_org(shape(1)+2,i,j)
                enddo
                enddo
                enddo
             case("ZERO")
                ! do nothing
             case default
                LOG_ERROR("ParentAtmosInputGrADS",*) 'upper_qv_type in PARAM_MKINIT_REAL_GrADS is invalid! ', upper_qv_type
                call PRC_abort
             end select
          endif

       case('QC')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qhyd_org(k+2,i,j,I_HC) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          ! if shape(1)>knum, QC is assumed to be zero.

       case('QR')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qhyd_org(k+2,i,j,I_HR) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          ! if shape(1)>knum, QR is assumed to be zero.

       case('QI')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qhyd_org(k+2,i,j,I_HI) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          ! if shape(1)>knum, QI is assumed to be zero.

       case('QS')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qhyd_org(k+2,i,j,I_HS) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          ! if shape(1)>knum, QS is assumed to be zero.

       case('QG')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qhyd_org(k+2,i,j,I_HG) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          ! if shape(1)>knum, QG is assumed to be zero.

       case('RH')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(:shape(1),:,:),      & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             qv_org(k+2,i,j) = work(k,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

          !$omp parallel do collapse(2) &
          !$omp private(qm,p_sat)
          do j = 1, JA_org
          do i = 1, IA_org
             do k = 1, shape(1)
                if( qv_org(k+1,i,j) .ne. UNDEF ) then
                   rh(i,j) = qv_org(k+2,i,j) / 100.0_RP         ! relative humidity
                   call psat( temp_org(k+2,i,j), p_sat )   ! satulation pressure
                   qm = EPSvap * rh(i,j) * p_sat &
                      / ( pres_org(k+2,i,j) - rh(i,j) * p_sat ) ! mixing ratio
                   qv_org(k+2,i,j) = qm / ( 1.0_RP + qm )  ! specific humidity
                end if
             enddo
          enddo
          enddo
          if( KA_org-2 > shape(1) ) then
             select case( upper_qv_type )
             case("COPY")
                !$omp parallel do &
                !$omp private(qm,p_sat)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = shape(1)+1, KA_org-2
                   call psat( temp_org(k+2,i,j), p_sat )   ! satulated specific humidity
                   qm = EPSvap * rh(i,j) * p_sat &
                      / ( pres_org(k+2,i,j) - rh(i,j) * p_sat ) ! mixing ratio
                   qv_org(k+2,i,j) = qm / ( 1.0_RP + qm )  ! specific humidity
                   qv_org(k+2,i,j) = min(qv_org(k+2,i,j),qv_org(k+1,i,j))
                enddo
                enddo
                enddo
             case("ZERO")
                ! do nothing
             case default
                LOG_ERROR("ParentAtmosInputGrADS",*) 'upper_qv_type in PARAM_MKINIT_REAL_GrADS is invalid! ', upper_qv_type
                call PRC_abort
             end select
          endif

       case('MSLP')

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(dummy,:,:),          & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(1,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
          enddo
          enddo

       case('PSFC')

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(dummy,:,:),          & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
          enddo
          enddo

       case('U10')

          if ( sfc_diagnoses ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(dummy,:,:),          & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
                velx_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
          end if

       case('V10')

          if ( sfc_diagnoses ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(dummy,:,:),          & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
                vely_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
          end if

       case('T2')

          if ( sfc_diagnoses ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(dummy,:,:),          & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
                temp_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
          end if

       case('Q2')

          if ( sfc_diagnoses ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(dummy,:,:),          & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
                qv_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
          end if

       case('RH2')

          if ( sfc_diagnoses ) then
             call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                   work(dummy,:,:),          & ! (out)
                                   step = nt,                & ! (in)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(3)
             do j = 1, JA_org
             do i = 1, IA_org
                qv_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
             enddo
             enddo
             !$omp parallel do collapse(2) &
             !$omp private (qm,p_sat)
             do j = 1, JA_org
             do i = 1, IA_org
                rh(i,j) = qv_org(2,i,j) / 100.0_RP
                call psat( temp_org(2,i,j), p_sat )   ! satulation pressure
                qm = EPSvap * rh(i,j) * p_sat &
                   / ( pres_org(2,i,j) - rh(i,j) * p_sat ) ! mixing ratio
                qv_org(2,i,j) = qm / ( 1.0_RP + qm )  ! specific humidity
             enddo
             enddo
          end if

       case('topo')

          call FILE_GrADS_read( file_id, var_id(ielem,1), & ! (in)
                                work(dummy,:,:),          & ! (out)
                                postfix = basename_num    ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
             cz_org(2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
          enddo
          enddo

       case('RN222')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,1), & ! (in)
                                     shape(:)                  ) ! (out)

          call FILE_GrADS_read( file_id, var_id(ielem,1),    & ! (in)
                                work(:shape(1),:,:),         & ! (out)
                                step = nt,                   & ! (in)
                                postfix = basename_num       ) ! (in)
          !$omp parallel do collapse(3)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, shape(1)
             RN222_org(k+2,i,j) = work(dummy,i-1+IS_org,j-1+JS_org)
          enddo
          enddo
          enddo

       end select
    enddo loop_InputAtmosGrADS

    lm_layer(:,:) = 3

    !$omp parallel do
    do j = 1, JA_org
    do i = 1, IA_org
       do k = 3, KA_org
          ! search the lowermost layer excluding UNDEF
          if( abs( pres_org(k,i,j) - UNDEF ) < EPS ) then
             lm_layer(i,j) = k + 1
          else
             exit
          end if
       end do
    end do
    end do

    ! density
    if ( var_id(Ia_dens,1) < 0 ) then
       !$omp parallel do &
       !$omp private (Rtot)
       do j = 1, JA_org
       do i = 1, IA_org
       do k = lm_layer(i,j), KA_org
          Rtot = Rdry * ( 1.0_RP + EPSTvap * qv_org(k,i,j) )
          dens_org(k,i,j) = pres_org(k,i,j) / ( Rtot * temp_org(k,i,j) )
       end do
       end do
       end do
    end if

    if ( sfc_diagnoses ) then
       ! surface
       if ( var_id(Ia_topo,1) > 0 ) then
          if ( .not. under_sfc ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                do k = lm_layer(i,j), KA_org
                   if ( cz_org(k,i,j) > cz_org(2,i,j) ) then
                      lm_layer(i,j) = k
                      exit
                   end if
                end do
             end do
             end do
          end if
          if ( var_id(Ia_t2,1) > 0 .and. var_id(Ia_ps,1) > 0 ) then
             !$omp parallel do &
             !$omp private (Rtot)
             do j = 1, JA_org
             do i = 1, IA_org
                Rtot = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                dens_org(2,i,j) = pres_org(2,i,j) / ( Rtot * temp_org(2,i,j) )
             end do
             end do
          else if ( var_id(Ia_ps,1) > 0 ) then
             !$omp parallel do &
             !$omp private (k,dz,Rtot)
             do j = 1, JA_org
             do i = 1, IA_org
                k = lm_layer(i,j)
                dz = cz_org(k,i,j) - cz_org(2,i,j)
                dens_org(2,i,j) = - ( pres_org(k,i,j) - pres_org(2,i,j) ) * 2.0_RP / ( GRAV * dz ) &
                                  - dens_org(k,i,j)
                Rtot = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                temp_org(2,i,j) = pres_org(2,i,j) / ( Rtot * dens_org(2,i,j) )
             end do
             end do
          else if ( var_id(Ia_t2,1) > 0 ) then
             !$omp parallel do &
             !$omp private (k,dz,Rtot)
             do j = 1, JA_org
             do i = 1, IA_org
                k = lm_layer(i,j)
                dz = cz_org(k,i,j) - cz_org(2,i,j)
                Rtot = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                dens_org(2,i,j) = ( pres_org(k,i,j) + GRAV * dens_org(k,i,j) * dz * 0.5_RP ) &
                                / ( Rtot * temp_org(2,i,j) - GRAV * dz * 0.5_RP )
                pres_org(2,i,j) = dens_org(2,i,j) * Rtot * temp_org(2,i,j)
             end do
             end do
          else
             !$omp parallel do &
             !$omp private(k,dz,Rtot)
             do j = 1, JA_org
             do i = 1, IA_org
                k = lm_layer(i,j)
                dz = cz_org(k,i,j) - cz_org(2,i,j)
                temp_org(2,i,j) = temp_org(k,i,j) + LAPS * dz
                Rtot = Rdry * ( 1.0_RP + EPSTvap * qv_org(2,i,j) )
                dens_org(2,i,j) = ( pres_org(k,i,j) + GRAV * dens_org(k,i,j) * dz * 0.5_RP ) &
                                / ( Rtot * temp_org(2,i,j) - GRAV * dz * 0.5_RP )
                pres_org(2,i,j) = dens_org(2,i,j) * Rtot * temp_org(2,i,j)
             end do
             end do
          end if
       else
          !$omp parallel do &
          !$omp private(k)
          do j = 1, JA_org
          do i = 1, IA_org
             k = lm_layer(i,j)
             ! ignore surface variables
             cz_org  (2,i,j)   = cz_org  (k,i,j)
             velz_org(2,i,j)   = velz_org(k,i,j)
             velx_org(2,i,j)   = velx_org(k,i,j)
             vely_org(2,i,j)   = vely_org(k,i,j)
             pres_org(2,i,j)   = pres_org(k,i,j)
             temp_org(2,i,j)   = temp_org(k,i,j)
             dens_org(2,i,j)   = dens_org(k,i,j)
             qv_org  (2,i,j)   = qv_org  (k,i,j)
             qhyd_org(2,i,j,:) = qhyd_org(k,i,j,:)
             RN222_org(2,i,j)  = RN222_org(k,i,j)
!!$          ! guess surface height (elevation)
!!$          if ( pres_org(2,i,j) < pres_org(1,i,j) ) then
!!$             lp2 = log( pres_org(2,i,j) / pres_org(1,i,j) )
!!$          else
!!$             lp2 = -1.0_RP
!!$          end if
!!$          if ( pres_org(k,i,j) < pres_org(1,i,j) ) then
!!$             lp3 = log( pres_org(k,i,j) / pres_org(1,i,j) )
!!$          else
!!$             lp3 = -1.0_RP
!!$          end if
!!$          cz_org(2,i,j) = cz_org(k,i,j) * lp2 / lp3
!!$          if ( cz_org(2,i,j) < 0.0_RP ) cz_org(2,i,j) = cz_org(k,i,j)
          end do
          end do
       end if

       if ( var_id(Ia_q2,1) < 0 .and. var_id(Ia_rh2,1) < 0 ) then
          !$omp parallel do private(k)
          do j = 1, JA_org
          do i = 1, IA_org
             k = lm_layer(i,j)
             qv_org(2,i,j) = qv_org(k,i,j)
          end do
          end do
       end if
       !$omp parallel do private(k)
       do j = 1, JA_org
       do i = 1, IA_org
          k = lm_layer(i,j)
          qv_org   (1,i,j)   = qv_org   (k,i,j)
          qhyd_org (1,i,j,:) = qhyd_org (k,i,j,:)
          qhyd_org (2,i,j,:) = qhyd_org (k,i,j,:)
          RN222_org(1,i,j)   = RN222_org(k,i,j)
          RN222_org(2,i,j)   = RN222_org(k,i,j)
       end do
       end do

       ! sea level
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          temp_org(1,i,j) = temp_org(2,i,j) + LAPS * cz_org(2,i,j)
       end do
       end do
       if ( var_id(Ia_slp,1) > 0 ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             dens_org(1,i,j) = pres_org(1,i,j) / ( Rdry * temp_org(1,i,j) )
          end do
          end do
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             dens_org(1,i,j) = ( pres_org(2,i,j) + GRAV * dens_org(2,i,j) * cz_org(2,i,j) * 0.5_RP ) &
                             / ( Rdry * temp_org(1,i,j) - GRAV * cz_org(2,i,j) * 0.5_RP )
             pres_org(1,i,j) = dens_org(1,i,j) * Rdry * temp_org(1,i,j)
          end do
          end do
       end if

    else
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          velz_org(1:2,i,j)   = UNDEF
          velx_org(1:2,i,j)   = UNDEF
          vely_org(1:2,i,j)   = UNDEF
          dens_org(1:2,i,j)   = UNDEF
          temp_org(1:2,i,j)   = UNDEF
          qv_org  (1:2,i,j)   = UNDEF
          qhyd_org(1:2,i,j,:) = UNDEF
          RN222_org(1:2,i,j)  = UNDEF
          pres_org(1  ,i,j)   = UNDEF
          cz_org  (1  ,i,j)   = UNDEF
       end do
       end do
    end if

    ! check verticaly extrapolated data in outer model
    if( pressure_coordinates .and. var_id(Ia_ps,1) > 0 ) then
       !$omp parallel do private(k)
       do j = 1, JA_org
       do i = 1, IA_org
          if ( under_sfc ) then
             k = lm_layer(i,j)
             if ( pres_org(1,i,j) < pres_org(k,i,j) ) then
                pres_org(1,i,j) = UNDEF
                cz_org  (1,i,j) = UNDEF
             end if
             if ( pres_org(2,i,j) < pres_org(k,i,j) ) then
                pres_org(2,i,j) = pres_org(1,i,j)
                cz_org  (2,i,j) = cz_org(1,i,j)
             end if
             cycle
          end if
          do k = 3, KA_org
             if( pres_org(k,i,j) > pres_org(2,i,j) ) then ! if Pressure is larger than Surface pressure
                if ( sfc_diagnoses ) then
                   velz_org(k,i,j)   = velz_org(2,i,j)
                   velx_org(k,i,j)   = velx_org(2,i,j)
                   vely_org(k,i,j)   = vely_org(2,i,j)
                   pres_org(k,i,j)   = pres_org(2,i,j)
                   dens_org(k,i,j)   = dens_org(2,i,j)
                   temp_org(k,i,j)   = temp_org(2,i,j)
                   qv_org  (k,i,j)   = qv_org  (2,i,j)
                   qhyd_org(k,i,j,:) = qhyd_org(2,i,j,:)
                   cz_org  (k,i,j)   = cz_org  (2,i,j)
                   RN222_org(k,i,j)  = RN222_org(2,i,j)
                else
                   velz_org(k,i,j)   = UNDEF
                   velx_org(k,i,j)   = UNDEF
                   vely_org(k,i,j)   = UNDEF
                   pres_org(k,i,j)   = UNDEF
                   dens_org(k,i,j)   = UNDEF
                   temp_org(k,i,j)   = UNDEF
                   qv_org  (k,i,j)   = UNDEF
                   qhyd_org(k,i,j,:) = UNDEF
                   cz_org  (k,i,j)   = UNDEF
                   RN222_org(k,i,j)  = UNDEF
                end if
             end if
          enddo
       enddo
       enddo
    else if ( var_id(Ia_topo,1) > 0 ) then
       !$omp parallel do private(k)
       do j = 1, JA_org
       do i = 1, IA_org
          if ( under_sfc ) then
             k = lm_layer(i,j)
             if ( cz_org(1,i,j) < cz_org(k,i,j) ) then
                pres_org(1,i,j) = UNDEF
                cz_org  (1,i,j) = UNDEF
             end if
             if ( cz_org(1,i,j) < cz_org(k,i,j) ) then
                pres_org(2,i,j) = pres_org(1,i,j)
                cz_org  (2,i,j) = cz_org(1,i,j)
             end if
             cycle
          end if

       do k = 3, KA_org
          if( cz_org(k,i,j) < cz_org(2,i,j) ) then
             if ( sfc_diagnoses ) then
                velz_org(k,i,j)   = velz_org(2,i,j)
                velx_org(k,i,j)   = velx_org(2,i,j)
                vely_org(k,i,j)   = vely_org(2,i,j)
                pres_org(k,i,j)   = pres_org(2,i,j)
                dens_org(k,i,j)   = dens_org(2,i,j)
                temp_org(k,i,j)   = temp_org(2,i,j)
                qv_org  (k,i,j)   = qv_org  (2,i,j)
                qhyd_org(k,i,j,:) = qhyd_org(2,i,j,:)
                cz_org  (k,i,j)   = cz_org  (2,i,j)
                RN222_org(k,i,j)  = 0.0_RP
             else
                velz_org(k,i,j)   = UNDEF
                velx_org(k,i,j)   = UNDEF
                vely_org(k,i,j)   = UNDEF
                pres_org(k,i,j)   = UNDEF
                dens_org(k,i,j)   = UNDEF
                temp_org(k,i,j)   = UNDEF
                qv_org  (k,i,j)   = UNDEF
                qhyd_org(k,i,j,:) = UNDEF
                cz_org  (k,i,j)   = UNDEF
                RN222_org(k,i,j)  = UNDEF
             end if
          endif
       enddo
       enddo
       enddo
    end if

    return
  end subroutine ParentAtmosInputGrADS

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupGrADS( &
       ldims,                       & ! (out)
       use_waterratio,              & ! (out)
       use_file_landwater,          & ! (in)
       basename                     )
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid
    implicit none
    integer,          intent(out) :: ldims(3)
    logical,          intent(out) :: use_waterratio
    logical,          intent(in)  :: use_file_landwater ! use landwater data from files
    character(len=*), intent(in)  :: basename

    character(len=H_SHORT) :: item
    integer                :: ielem
    !---------------------------------------------------------------------------

    LOG_INFO("ParentLandSetupGrADS",*) 'Real Case/Land Input File Type: GrADS format'

    !--- initialization
    use_waterratio = .false.

    if ( basename == "" ) then
       LOG_ERROR("ParentLandSetupGrADS",*) '"BASEMAAME" is not specified in "PARAM_MKINIT_ATMOS_GRID_CARTESC_REAL_ATOMS"!', trim(basename)
       call PRC_abort
    endif

    call FILE_GrADS_open( basename, & ! (in)
                          file_id   ) ! (out)

    call FILE_GrADS_get_shape( file_id, "STEMP", & ! (in)
                               ldims(:)          ) ! (out)


    ! check existence
    do ielem = 1, num_item_list_land
       item  = item_list_land(ielem)
       call FILE_GrADS_varid( file_id, item,  & ! (in)
                              var_id(ielem,2) ) ! (out)
    end do

    ! check necessary data
    do ielem = 1, num_item_list_land
       item  = item_list_land(ielem)

       select case(item)
       case('lsmask')
          if ( var_id(ielem,2) < 0 ) then
             LOG_WARN("ParentLandSetupGrADS",*) trim(item),' is not found & not used.'
          endif
          cycle
       case('topo','topo_sfc')
          if ( var_id(Il_topo_sfc,2) < 0 ) then
             if ( var_id(Il_topo,2) < 0 ) then
                LOG_WARN("ParentLandSetupGrADS",*) '"topo" and "topo_sfc" are not found & not used.'
             end if
          else
             var_id(Il_topo,2) = -1
          end if
          cycle
       case('lon', 'lat', 'lon_sfc', 'lat_sfc')
          if ( var_id(Il_lon_sfc,2) < 0 ) then
             if ( var_id(Il_lon,2) < 0 ) then
                LOG_ERROR("ParentLandSetupGrADS",*) 'either lon or lon_sfc is required'
                call PRC_abort
             end if
          else
             var_id(Il_lon,2) = -1
          end if
          if ( var_id(Il_lat_sfc,2) < 0 ) then
             if ( var_id(Il_lat,2) < 0 ) then
                LOG_ERROR("ParentLandSetupGrADS",*) 'either lat or lat_sfc is required'
                call PRC_abort
             end if
          else
             var_id(Il_lat,2) = -1
          end if
          cycle
       case('SMOISVC', 'SMOISDS')
          if ( use_file_landwater ) then
             if ( var_id(Il_smoisvc,2) < 0 .and. var_id(Il_smoisds,2) < 0 ) then
                LOG_ERROR("ParentLandSetupGrADS",*) 'Not found in grads namelist! : ',trim(item_list_land(ielem))
                call PRC_abort
             end if
             if ( var_id(Il_smoisds,2) > 0 ) then
                use_waterratio = .true.
                var_id(Il_smoisvc,2) = -1
             end if
          else
             var_id(Il_smoisvc,2) = -1
             var_id(Il_smoisds,2) = -1
          end if
          cycle
       case default ! llev, SKINT, STEMP
          if ( var_id(ielem,2) < 0 ) then
             LOG_ERROR("ParentLandSetupGrADS",*) 'Not found in grads namelist! : ',trim(item_list_land(ielem))
             call PRC_abort
          endif
       end select

    end do

    return
  end subroutine ParentLandSetupGrADS

  subroutine ParentLandInputGrADS( &
      tg_org,             & ! (out)
      strg_org,           & ! (out)
      smds_org,           & ! (out)
      lst_org,            & ! (out)
      llon_org,           & ! (out)
      llat_org,           & ! (out)
      lz_org,             & ! (out)
      topo_org,           & ! (out)
      lmask_org,          & ! (out)
      basename_num,       & ! (in)
      ldims,              & ! (in)
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
       FILE_GrADS_isOneD,    &
       FILE_GrADS_get_shape, &
       FILE_GrADS_read
    implicit none

    real(RP),         intent(out) :: tg_org   (:,:,:)
    real(RP),         intent(out) :: strg_org (:,:,:)
    real(RP),         intent(out) :: smds_org (:,:,:)
    real(RP),         intent(out) :: lst_org  (:,:)
    real(RP),         intent(out) :: llon_org (:,:)
    real(RP),         intent(out) :: llat_org (:,:)
    real(RP),         intent(out) :: lz_org   (:)
    real(RP),         intent(out) :: topo_org (:,:)
    real(RP),         intent(out) :: lmask_org(:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: ldims(3)
    logical,          intent(in)  :: use_file_landwater ! use land water data from files
    integer,          intent(in)  :: nt

    real(RP) :: lon1D(ldims(2)), lat1D(ldims(3))
    integer  :: shape(3)

    character(len=H_SHORT) :: item
    integer                :: ielem

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(2)
    do j = 1, ldims(3)
    do i = 1, ldims(2)
       lmask_org(i,j) = UNDEF
       topo_org(i,j)  = UNDEF
    end do
    end do

    loop_InputLandGrADS : do ielem = 1, num_item_list_land

       item  = item_list_land(ielem)

       if ( var_id(ielem,2) < 0 ) cycle

       ! read data
       select case(item)
       case("lsmask")

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                lmask_org(:,:),           & ! (out)
                                postfix = basename_num    ) ! (in)

       case("lon", "lon_sfc")

          if ( item == "lon" ) then
             if ( FILE_GrADS_isOneD( file_id, var_id(ielem,2) ) ) then
                call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                           shape(1:1)                ) ! (out)
                if ( ldims(2).ne.shape(1) .and. shape(1).ne.-1 ) then
                   LOG_ERROR("ParentLandInputGrADS",*) 'dimension of "lon" is different! ', ldims(2), shape(1)
                   call PRC_abort
                end if
             else
                call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                           shape(1:2)                ) ! (out)
                if ( ldims(2).ne.shape(1) .or. ldims(3).ne.shape(2) ) then
                   LOG_ERROR("ParentLandInputGrADS",*) 'dimension of "lon" is different! ', ldims(2), shape(1), ldims(3), shape(2)
                   call PRC_abort
                end if
             end if
          end if

          if ( FILE_GrADS_isOneD( file_id, var_id(ielem,2) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                   lon1D(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llon_org(i,j) = lon1D(i) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                   llon_org(:,:),            & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llon_org(i,j) = llon_org(i,j) * D2R
             enddo
             enddo
          end if

       case("lat", "lat_sfc")

          if ( item == "lat" ) then
             if ( FILE_GrADS_isOneD( file_id, var_id(ielem,2) ) ) then
                call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                           shape(1:1)                ) ! (out)
                if ( ldims(3).ne.shape(1) .and. shape(1).ne.-1 ) then
                   LOG_ERROR("ParentLandInputGrADS",*) 'dimension of "lat" is different! ', ldims(3), shape(1)
                   call PRC_abort
                end if
             else
                call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                           shape(1:2)                ) ! (out)
                if ( ldims(2).ne.shape(1) .or. ldims(3).ne.shape(2) ) then
                   LOG_ERROR("ParentLandInputGrADS",*) 'dimension of "lat" is different! ', ldims(2), shape(1), ldims(3), shape(2)
                   call PRC_abort
                end if
             end if
          end if

          if ( FILE_GrADS_isOneD( file_id, var_id(ielem,2) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                   lat1D(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llat_org(i,j) = lat1D(j) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                   llat_org(:,:),            & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llat_org(i,j) = llat_org(i,j) * D2R
             enddo
             enddo
          end if

       case("llev")

          call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                     shape(:)                  ) ! (out)
          if( ldims(1) .ne. shape(1) )then
             LOG_ERROR("ParentLandInputGrADS",*) '"nz" must be equal to nz of "STEMP" for llev. :', shape(1), ldims(1)
             call PRC_abort
          endif
          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                lz_org(:)                 ) ! (out)

       case('STEMP')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( ldims(1) .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to nz of "STEMP" for ',trim(item),'. :', shape(1), ldims(1)
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                tg_org(:,:,:),            & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)

       case('SMOISVC')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( ldims(1) .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to nz of "STEMP" for ',trim(item),'. :', shape(1), ldims(1)
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                strg_org(:,:,:),          & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)

       case('SMOISDS')

          call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                     shape(:)                  ) ! (out)
          if ( ldims(1) .ne. shape(1) ) then
             LOG_ERROR("ParentAtmosInputGrADS",*) '"nz" must be equal to nz of "STEMP" for ',trim(item),'. :', shape(1), ldims(1)
             call PRC_abort
          endif

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                smds_org(:,:,:),          & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)

       case('SKINT')

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                lst_org(:,:),             & ! (out)
                                step = nt,                & ! (in)
                                postfix = basename_num    ) ! (in)

       case('topo', 'topo_sfc')

          if ( item == "topo" ) then
             call FILE_GrADS_get_shape( file_id, var_id(ielem,2), & ! (in)
                                        shape(1:2)                ) ! (out)
             if ( ldims(2).ne.shape(1) .or. ldims(3).ne.shape(2) ) then
                LOG_WARN("ParentLandInputGrADS",*) 'namelist of "topo_sfc" is not found in grads namelist!'
                LOG_WARN_CONT(*) 'dimension of "topo" is different! ', ldims(2), shape(1), ldims(3), shape(2)
                cycle
             end if
          end if

          call FILE_GrADS_read( file_id, var_id(ielem,2), & ! (in)
                                topo_org(:,:),            & ! (out)
                                postfix = basename_num    ) ! (in)

       end select
    enddo loop_InputLandGrADS

    return
  end subroutine ParentLandInputGrADS

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupGrADS( &
       odims,   & ! (out)
       timelen, & ! (out)
       basename ) ! (in)
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_varid, &
       FILE_GrADS_get_shape
    implicit none

    integer,          intent(out) :: odims(2)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename

    character(len=H_SHORT) :: item
    integer                :: ielem
    integer                :: vid
    !---------------------------------------------------------------------------

    LOG_INFO("ParentOceanSetupGrADS",*) 'Real Case/Ocean Input File Type: GrADS format'

    !--- read namelist

    call FILE_GrADS_open( basename, & ! (in)
                          file_id   ) ! (out)

    call FILE_GrADS_varid( file_id, "SST", & ! (in)
                           vid             ) ! (out)
    if ( vid < 0 ) then
       call FILE_GrADS_varid( file_id, "SKINT", & ! (in)
                              vid               ) ! (out)
    end if

    if ( vid < 0 ) then
       LOG_ERROR("ParentOceanSetupGrADS",*) 'SST and SKINT are found in grads namelist!'
       call PRC_abort
    end if

    call FILE_GrADS_get_shape( file_id, vid, & ! (in)
                               odims(:)      ) ! (out)


    timelen = 0        ! will be replaced later


    ! check existance
    do ielem = 1, num_item_list_ocean
       item  = item_list_ocean(ielem)
       call FILE_GrADS_varid( file_id, item,  & ! (in)
                              var_id(ielem,3) ) ! (out)
    end do

    ! check necessary datea
    do ielem = 1, num_item_list_ocean
       item  = item_list_ocean(ielem)

       select case(item)
       case('lsmask','lsmask_sst')
          if ( var_id(Io_lsmask_sst,3) < 3 ) then
             if ( var_id(Io_lsmask,3) < 3 ) then
                LOG_WARN("ParentOceanSetupGrADS",*) trim(item),' is not found & not used.'
             end if
          else
             var_id(Io_lsmask,3) = -1
          endif
       case('lon', 'lat', 'lon_sfc', 'lat_sfc', 'lon_sst', 'lat_sst')
          if ( var_id(Io_lon_sst,3) < 0 ) then
             if ( var_id(Io_lon_sfc,3) < 0 ) then
                if ( var_id(Io_lon,3) < 0 ) then
                   LOG_ERROR("ParentOceanSetupGrADS",*) 'either lon_sst, lon_sfc, or lon is necessary!'
                   call PRC_abort
                end if
             else
                var_id(Io_lon,3) = -1
             end if
          else
             var_id(Io_lon_sfc,3) = -1
             var_id(Io_lon,    3) = -1
          end if
          if ( var_id(Io_lat_sst,3) < 0 ) then
             if ( var_id(Io_lat_sfc,3) < 0 ) then
                if ( var_id(Io_lat,3) < 0 ) then
                   LOG_ERROR("ParentOceanSetupGrADS",*) 'either lat_sst, lat_sfc, or lat is necessary!'
                   call PRC_abort
                end if
             else
                var_id(Io_lat,3) = -1
             end if
          else
             var_id(Io_lat_sfc,3) = -1
             var_id(Io_lat,    3) = -1
          end if
       case('SST','SKINT')
          if ( var_id(Io_sst,3) < 0 ) then
             if ( var_id(Io_skint,3) < 0 ) then
                LOG_ERROR("ParentOceanSetupGrADS",*) 'SST and SKINT are found in grads namelist!'
                call PRC_abort
             else
                LOG_WARN("ParentOceanSetupGrADS",*) 'SST is found in grads namelist. SKINT is used in place of SST.'
             end if
          else
             var_id(Io_skint,3) = -1
          end if
       case default !
          if ( var_id(ielem,3) < 0 ) then
             LOG_ERROR("ParentOceanSetupGrADS",*) 'Not found in grads namelist! : ', &
                        trim(item_list_ocean(ielem))
             call PRC_abort
          endif
       end select

    end do

    return
  end subroutine ParentOceanSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenGrADS
    implicit none

    return
  end subroutine ParentOceanOpenGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputGrADS( &
      tw_org,       & ! (out)
      sst_org,      & ! (out)
      omask_org,    & ! (out)
      olon_org,     & ! (out)
      olat_org,     & ! (out)
      basename_num, & ! (in)
      odims,        & ! (in)
      nt            ) ! (in)
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    use scale_file_grads, only: &
       FILE_GrADS_isOneD,    &
       FILE_GrADS_get_shape, &
       FILE_GrADS_read
    implicit none

    real(RP),         intent(out) :: tw_org   (:,:)
    real(RP),         intent(out) :: sst_org  (:,:)
    real(RP),         intent(out) :: omask_org(:,:)
    real(RP),         intent(out) :: olon_org (:,:)
    real(RP),         intent(out) :: olat_org (:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: odims(2)
    integer,          intent(in)  :: nt

    character(len=H_SHORT) :: item
    integer                :: ielem

    real(RP) :: lon1D(odims(1)), lat1D(odims(2))

    integer :: shape(2)
    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(2)
    do j = 1, odims(2)
    do i = 1, odims(1)
       omask_org(i,j) = UNDEF
    end do
    end do

    loop_InputOceanGrADS : do ielem = 1, num_item_list_ocean

       if ( var_id(ielem,3) < 0 ) cycle

       item  = item_list_ocean(ielem)

       ! read data
       select case(item)
       case("lsmask","lsmask_sst")

          if ( item == "lsmask" ) then
             call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                        shape(:)                  ) ! (out)
             if ( odims(1) .ne. shape(1) .or. odims(2) .ne. shape(2) ) then
                LOG_WARN("ParentOceanInputGrADS",*) 'dimension of lsmask is different. not use'
                cycle
             end if
          end if

          call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                omask_org(:,:),           & ! (out)
                                postfix = basename_num    ) ! (in)

       case("lon","lon_sfc","lon_sst")

          if ( item .ne. "lon_sst" ) then
             if ( FILE_GrADS_isOneD( file_id, var_id(ielem,3) ) ) then
                call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                           shape(1:1)                ) ! (out)
                if ( odims(1).ne.shape(1) .and. shape(1).ne.-1 ) then
                   LOG_ERROR("ParentOceanInputGrADS",*) 'dimension of "',trim(item),'" is different! ', odims(1), shape(1)
                   call PRC_abort
                end if
             else
                call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                           shape(:)                  ) ! (out)
                if ( odims(1).ne.shape(1) .or. odims(2).ne.shape(2) ) then
                   LOG_ERROR("ParentOceanInputGrADS",*) 'dimension of "',trim(item),'" is different', odims(1), shape(1), odims(2), shape(2)
                   call PRC_abort
                end if
             end if
          end if

          if ( FILE_GrADS_isOneD( file_id, var_id(ielem,3) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                   lon1D(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, odims(2)
             do i = 1, odims(1)
                olon_org(i,j) = lon1D(i) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                   olon_org(:,:),            & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, odims(2)
             do i = 1, odims(1)
                olon_org(i,j) = olon_org(i,j) * D2R
             enddo
             enddo
          end if

       case("lat","lat_sfc","lat_sst")

          if ( item .ne. "lat_sst" ) then
             if ( FILE_GrADS_isOneD( file_id, var_id(ielem,3) ) ) then
                call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                           shape(1:1)                ) ! (out)
                if ( odims(2).ne.shape(1) .and. shape(1).ne.-1 ) then
                   LOG_ERROR("ParentOceanInputGrADS",*) 'dimension of "',trim(item),'" is different! ', odims(2), shape(1)
                   call PRC_abort
                end if
             else
                call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                           shape(:)                  ) ! (out)
                if ( odims(1).ne.shape(1) .or. odims(2).ne.shape(2) ) then
                   LOG_ERROR("ParentOceanInputGrADS",*) 'dimension of "',trim(item),'" is different', odims(1), shape(1), odims(2), shape(2)
                   call PRC_abort
                end if
             end if
          end if

          if ( FILE_GrADS_isOneD( file_id, var_id(ielem,3) ) ) then
             call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                   lat1D(:)                  ) ! (out)
             !$omp parallel do collapse(2)
             do j = 1, odims(2)
             do i = 1, odims(1)
                olat_org(i,j) = lat1D(j) * D2R
             enddo
             enddo
          else
             call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                   olat_org(:,:),            & ! (out)
                                   postfix = basename_num    ) ! (in)
             !$omp parallel do collapse(2)
             do j = 1, odims(2)
             do i = 1, odims(1)
                olat_org(i,j) = olat_org(i,j) * D2R
             enddo
             enddo
          end if

       case("SKINT","SST")

          if ( item == "SKINT" ) then
             call FILE_GrADS_get_shape( file_id, var_id(ielem,3), & ! (in)
                                        shape(:)                  ) ! (out)
             if ( odims(1).ne.shape(1) .or. odims(2).ne.shape(2) ) then
                LOG_ERROR("ParentLandOceanGrADS",*) 'dimension of "',trim(item),'" is different', odims(1), shape(1), odims(2), shape(2)
                call PRC_abort
             end if
          end if

          call FILE_GrADS_read( file_id, var_id(ielem,3), & ! (in)
                                sst_org(:,:),             & ! (out)
                                postfix = basename_num    ) ! (in)

       end select
    enddo loop_InputOceanGrADS

    tw_org = sst_org

    return
  end subroutine ParentOceanInputGrADS


end module mod_realinput_grads
