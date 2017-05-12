!-------------------------------------------------------------------------------
!> module REAL input GrADS
!!
!! @par Description
!!          read data from GrADS file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-05-24 (S.Nishizawa)   [new] split from mod_realinput.f90
!!
!<
module mod_realinput_grads
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_tracer
  use scale_process, only: &
     myrank => PRC_myrank,  &
     PRC_MPIstop
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetupGrADS
  public :: ParentAtomOpenGrADS
  public :: ParentAtomInputGrADS
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
  integer,  parameter    :: grads_vars_limit = 1000 !> limit of number of values
  integer,  parameter    :: num_item_list = 22
  integer,  parameter    :: num_item_list_atom  = 22
  integer,  parameter    :: num_item_list_land  = 11
  integer,  parameter    :: num_item_list_ocean = 9
  logical                :: data_available(num_item_list_atom,3) ! 1:atom, 2:land, 3:ocean
  character(len=H_SHORT) :: item_list_atom (num_item_list_atom)
  character(len=H_SHORT) :: item_list_land (num_item_list_land)
  character(len=H_SHORT) :: item_list_ocean(num_item_list_ocean)
  data item_list_atom  /'lon','lat','plev','U','V','T','HGT','QV','QC','QR','QI','QS','QG','RH', &
                        'MSLP','PSFC','U10','V10','T2','Q2','RH2','TOPO' /
  data item_list_land  /'lsmask','lon','lat','lon_sfc','lat_sfc','llev', &
                        'STEMP','SMOISVC','SMOISDS','SKINT','TOPO','TOPO_sfc' /
  data item_list_ocean /'lsmask','lsmask_sst','lon','lat','lon_sfc','lat_sfc','lon_sst','lat_sst','SKINT','SST'/

  integer,  parameter   :: Ia_lon    = 1
  integer,  parameter   :: Ia_lat    = 2
  integer,  parameter   :: Ia_p      = 3  ! Pressure (Pa)
  integer,  parameter   :: Ia_u      = 4
  integer,  parameter   :: Ia_v      = 5
  integer,  parameter   :: Ia_t      = 6
  integer,  parameter   :: Ia_hgt    = 7  ! Geopotential height (m)
  integer,  parameter   :: Ia_qv     = 8
  integer,  parameter   :: Ia_qc     = 9 
  integer,  parameter   :: Ia_qr     = 10
  integer,  parameter   :: Ia_qi     = 11
  integer,  parameter   :: Ia_qs     = 12
  integer,  parameter   :: Ia_qg     = 13
  integer,  parameter   :: Ia_rh     = 14 ! Percentage (%)
  integer,  parameter   :: Ia_slp    = 15 ! Sea level pressure (Pa)
  integer,  parameter   :: Ia_ps     = 16 ! Surface pressure (Pa)
  integer,  parameter   :: Ia_u10    = 17
  integer,  parameter   :: Ia_v10    = 18
  integer,  parameter   :: Ia_t2     = 19
  integer,  parameter   :: Ia_q2     = 20
  integer,  parameter   :: Ia_rh2    = 21 ! Percentage (%)
  integer,  parameter   :: Ia_topo   = 22

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


  integer,  parameter   :: lvars_limit = 1000 ! limit of values for levels data
  real(RP), parameter   :: large_number_one = 9.999E+15_RP


  character(len=H_SHORT) :: upper_qv_type = "ZERO" !< how qv is given at higher level than outer model
                                                   !< "ZERO": 0
                                                   !< "COPY": copy values from the highest level of outer model

  character(len=H_SHORT) :: grads_item    (num_item_list,3)
  character(len=H_LONG)  :: grads_dtype   (num_item_list,3)
  character(len=H_LONG)  :: grads_fname   (num_item_list,3)
  character(len=H_SHORT) :: grads_fendian (num_item_list,3)
  character(len=H_SHORT) :: grads_yrev    (num_item_list,3)
  real(RP)               :: grads_swpoint (num_item_list,3)
  real(RP)               :: grads_dd      (num_item_list,3)
  integer                :: grads_lnum    (num_item_list,3)
  real(RP)               :: grads_lvars   (lvars_limit,num_item_list,3)
  integer                :: grads_startrec(num_item_list,3)
  integer                :: grads_totalrec(num_item_list,3)
  integer                :: grads_knum    (num_item_list,3)
  real(SP)               :: grads_missval (num_item_list,3)

  real(SP), allocatable :: gdata2D(:,:)
  real(SP), allocatable :: gdata3D(:,:,:)
  real(SP), allocatable :: gland2D(:,:)
  real(SP), allocatable :: gland3D(:,:,:)
  real(SP), allocatable :: gsst2D (:,:)

  integer :: io_fid_grads_nml  = -1
  integer :: io_fid_grads_data = -1


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

  NAMELIST / nml_grads_grid / &
       outer_nx,     &
       outer_ny,     &
       outer_nz,     &
       outer_nl,     &
       outer_nx_sfc, &
       outer_ny_sfc, &
       outer_nx_sst, &
       outer_ny_sst

  character(len=H_SHORT) :: item                                      ! up to 16 characters
  integer                :: knum                                      ! optional: vertical level
  character(len=H_SHORT) :: dtype                                     ! 'linear','levels','map'
  character(len=H_LONG)  :: fname                                     ! head of file name
  real(RP)               :: swpoint                                   ! start point (south-west point) for "linear"
  real(RP)               :: dd                                        ! dlon,dlat for "linear"
  integer                :: lnum                                      ! number of data
  real(RP)               :: lvars(lvars_limit) = large_number_one     ! values for "levels"
  integer                :: startrec                                  ! record position
  integer                :: totalrec                                  ! total record number per one time
  real(SP)               :: missval                                   ! missing value
  character(len=H_SHORT) :: fendian='big'                             ! option for "map"
  character(len=H_SHORT) :: yrev='off'                                ! option for "map", if yrev=on, order of data is NW to SE.


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtomSetupGrADS( &
      dims,                        & ! (out)
      basename                     ) ! (in)
    implicit none

    integer,          intent(out) :: dims(6)
    character(len=*), intent(in)  :: basename


    NAMELIST / PARAM_MKINIT_REAL_GrADS / &
        upper_qv_type

    integer :: ielem

    integer :: k, n
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: GrADS format'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_GrADS,iostat=ierr)

    if( ierr > 0 ) then
       write(*,*) 'xxx [realinput_grads] Not appropriate names in namelist PARAM_MKINIT_REAL_GrADS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_MKINIT_REAL_GrADS)


    if ( len_trim(basename) == 0 ) then
       write(*,*) 'xxx [realinput_grads] "BASENAME_ORG" is not specified in "PARAM_MKINIT_REAL_ATMOS"!', trim(basename)
       call PRC_MPIstop
    endif

    !--- read namelist
    io_fid_grads_nml = IO_get_available_fid()
    open( io_fid_grads_nml,       &
         file   = trim(basename), &
         form   = 'formatted',    &
         status = 'old',          &
         action = 'read',         &
         iostat = ierr            )
    if ( ierr /= 0 ) then
       write(*,*) 'xxx [realinput_grads] Input file is not found! ', trim(basename)
       call PRC_MPIstop
    endif

    read(io_fid_grads_nml,nml=nml_grads_grid,iostat=ierr)
    if( ierr /= 0 ) then !--- missing or fatal error
       write(*,*) 'xxx [realinput_grads] Not appropriate names in nml_grads_grid in ', trim(basename),'. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=nml_grads_grid)

    ! full level
    dims(1) = outer_nz ! bottom_top
    dims(2) = outer_nx ! west_east
    dims(3) = outer_ny ! south_north
    ! half level
    dims(4) = outer_nz ! bottom_top_stag
    dims(5) = outer_nx ! west_east for 2dim data
    dims(6) = outer_ny ! south_north for 2dim data

    allocate( gdata2D( dims(2), dims(3)          ) )
    allocate( gdata3D( dims(2), dims(3), dims(1) ) )

    call read_namelist( &
         grads_item(:,1),     & ! (out)
         grads_fname(:,1),    & ! (out)
         grads_dtype(:,1),    & ! (out)
         grads_swpoint(:,1),  & ! (out)
         grads_dd(:,1),       & ! (out)
         grads_lnum(:,1),     & ! (out)
         grads_lvars(:,:,1),  & ! (out)
         grads_startrec(:,1), & ! (out)
         grads_totalrec(:,1), & ! (out)
         grads_knum(:,1),     & ! (out)
         grads_yrev(:,1),     & ! (out)
         grads_fendian(:,1),  & ! (out)
         grads_missval(:,1),  & ! (out)
         data_available(:,1), & ! (out)
         item_list_atom,      & ! (in)
         num_item_list_atom,  & ! (in)
         basename,            & ! (in)
         io_fid_grads_nml     ) ! (in)

    close( io_fid_grads_nml )

    do ielem = 1, num_item_list_atom
       item  = item_list_atom(ielem)
       !--- check data
       select case(trim(item))
       case('QV')
          if (.not. data_available(Ia_qv,1)) then
             if (.not.data_available(Ia_rh,1)) then
                write(*,*) 'xxx [realinput_grads] Not found in grads namelist! : QV and RH'
                call PRC_MPIstop
             else ! will read RH
                cycle
             endif
          endif
       case('RH')
          if (.not. data_available(Ia_qv,1))then
             if(data_available(Ia_rh,1)) then
                if ((.not. data_available(Ia_t,1)).or.(.not. data_available(Ia_p,1))) then
                   write(*,*) 'xxx [realinput_grads] Temperature and pressure are required to convert from RH to QV ! '
                   call PRC_MPIstop
                else
                   cycle ! read RH and estimate QV
                endif
             else
                write(*,*) 'xxx [realinput_grads] Not found in grads namelist! : QV and RH'
                call PRC_MPIstop
             endif
          endif
       case('QC','QR','QI','QS','QG')
          if (.not. data_available(ielem,1)) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: ',trim(item),' is not found & will be estimated.'
             cycle
          endif
       case('MSLP','PSFC','U10','V10','T2')
          if (.not. data_available(ielem,1)) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: ',trim(item),' is not found & will be estimated.'
             cycle
          endif
       case('Q2')
          if ( .not. data_available(Ia_q2,1) ) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: Q2 is not found & will be estimated.'
             cycle
          endif
       case('RH2')
          if ( data_available(Ia_q2,1) ) then
             cycle
          else
             if ( data_available(Ia_rh2,1) ) then
                if ((.not. data_available(Ia_t2,1)).or.(.not. data_available(Ia_ps,1))) then
                   if( IO_L ) write(IO_FID_LOG,*) 'warning: T2 and PSFC are required to convert from RH2 to Q2 !'
                   if( IO_L ) write(IO_FID_LOG,*) '         Q2 will be copied from data at above level.'
                   data_available(Ia_rh2,1) = .false.
                   cycle
                endif
             else
                if( IO_L ) write(IO_FID_LOG,*) 'warning: Q2 and RH2 are not found, Q2 will be estimated.'
                cycle
             endif
          endif
       case('TOPO')
          if ( .not. data_available(ielem,1) ) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: ',trim(item),' is not found & not used.'
             cycle
          endif
       case default ! lon, lat, plev, U, V, T, HGT
          if ( .not. data_available(ielem,1) ) then
             write(*,*) 'xxx [realinput_grads] Not found in grads namelist! : ',trim(item_list_atom(ielem))
             call PRC_MPIstop
          endif
       end select

    end do

    return
  end subroutine ParentAtomSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtomOpenGrADS
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtomOpenGrADS]'

    return
  end subroutine ParentAtomOpenGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtomInputGrADS( &
       velz_org, &
       velx_org, &
       vely_org, &
       pres_org, &
       temp_org, &
       qtrc_org, &
       lon_org,  &
       lat_org,  &
       cz_org,   &
       basename_num, &
       dims, &
       nt )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R => CONST_D2R,   &
       EPS => CONST_EPS,   &
       EPSvap => CONST_EPSvap, &
       GRAV => CONST_GRAV, &
       LAPS => CONST_LAPS, &
       P00 => CONST_PRE00, &
       Rdry => CONST_Rdry, &
       CPdry => CONST_CPdry
    use scale_atmos_saturation, only: &
       psat => ATMOS_SATURATION_psat_liq
    implicit none


    real(RP),         intent(out) :: velz_org(:,:,:)
    real(RP),         intent(out) :: velx_org(:,:,:)
    real(RP),         intent(out) :: vely_org(:,:,:)
    real(RP),         intent(out) :: pres_org(:,:,:)
    real(RP),         intent(out) :: temp_org(:,:,:)
    real(RP),         intent(out) :: qtrc_org(:,:,:,:)
    real(RP),         intent(out) :: lon_org(:,:)
    real(RP),         intent(out) :: lat_org(:,:)
    real(RP),         intent(out) :: cz_org(:,:,:)
    character(len=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(6)
    integer,          intent(in)  :: nt

    real(RP) :: rhprs_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: pott(dims(2),dims(3))

    real(RP) :: RovCP
    real(RP) :: CPovR

    integer  :: lm_layer(dims(2),dims(3))

    ! data
    character(len=H_LONG) :: gfile

    real(RP) :: p_sat, qm, rhsfc
    real(RP) :: lp2, lp3

    integer  :: i, j, k, iq, ielem

    logical  :: pressure_coordinates
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[AtomInputGrADS]'

    qtrc_org(:,:,:,:) = 0.0_RP

    !--- read grads data
    loop_InputAtomGrADS : do ielem = 1, num_item_list_atom

       if ( .not. data_available(ielem,1) ) cycle

       item     = grads_item    (ielem,1)
       dtype    = grads_dtype   (ielem,1)
       fname    = grads_fname   (ielem,1)
       lnum     = grads_lnum    (ielem,1)
       missval  = grads_missval (ielem,1)

       if ( dims(1) < grads_knum(ielem,1) ) then
          write(*,*) 'xxx "knum" must be less than or equal to outer_nz. knum:',knum,'> outer_nz:',dims(1),trim(item)
          call PRC_MPIstop
       else if ( grads_knum(ielem,1) > 0 )then
          knum = grads_knum(ielem,1)  ! not missing
       else
          knum = dims(1)
       endif

       select case(trim(dtype))
       case("linear")
          swpoint = grads_swpoint (ielem,1)
          dd      = grads_dd      (ielem,1)
          if( (abs(swpoint-large_number_one)<EPS).or.(abs(dd-large_number_one)<EPS) )then
             write(*,*) 'xxx "swpoint" is required in grads namelist! ',swpoint
             write(*,*) 'xxx "dd"      is required in grads namelist! ',dd
             call PRC_MPIstop
          endif
       case("levels")
          if ( lnum < 0 )then
             write(*,*) 'xxx "lnum" is required in grads namelist for levels data! '
             call PRC_MPIstop
          endif
          do k=1, lnum
             lvars(k)=grads_lvars(k,ielem,1)
          enddo
          if(abs(lvars(1)-large_number_one)<EPS)then
             write(*,*) 'xxx "lvars" must be specified in grads namelist for levels data! '
             call PRC_MPIstop
          endif
       case("map")
          startrec = grads_startrec(ielem,1)
          totalrec = grads_totalrec(ielem,1)
          fendian  = grads_fendian (ielem,1)
          yrev     = grads_yrev (ielem,1)
          if( (startrec<0).or.(totalrec<0) )then
             write(*,*) 'xxx "startrec" is required in grads namelist! ',startrec
             write(*,*) 'xxx "totalrec" is required in grads namelist! ',totalrec
             call PRC_MPIstop
          endif
          ! get file_id
          if(io_fid_grads_data < 0)then
             io_fid_grads_data = IO_get_available_fid()
          endif
          gfile=trim(fname)//trim(basename_num)//'.grd'
          if( len_trim(fname)==0 )then
             write(*,*) 'xxx "fname" is required in grads namelist for map data! ',trim(fname)
             call PRC_MPIstop
          endif
       end select

       ! read data
       select case(trim(item))
       case("lon")
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
                lon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,1,item,startrec,totalrec,yrev,gdata2D)
             lon_org(:,:) = real(gdata2D(:,:), kind=RP) * D2R
          endif
       case("lat")
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
                lat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,1,item,startrec,totalrec,yrev,gdata2D)
             lat_org(:,:) = real(gdata2D(:,:), kind=RP) * D2R
          endif
       case("plev")
          if(dims(1)/=knum)then
             write(*,*) 'xxx "knum" must be equal to outer_nz for plev. knum:',knum,'> outer_nz:',dims(1)
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "levels" ) then
             pressure_coordinates = .true. ! use pressure coordinate in the input data
             if(dims(1)/=lnum)then
                write(*,*) 'xxx lnum must be same as the outer_nz for plev! ',dims(1),lnum
                call PRC_MPIstop
             endif
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)
                pres_org(k+2,i,j) = real(lvars(k), kind=RP)
             enddo
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             pressure_coordinates = .false.
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),dims(1),nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)
                pres_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                ! replace missval with UNDEF
                if( abs( pres_org(k+2,i,j) - missval ) < EPS ) then
                   pres_org(k+2,i,j) = UNDEF
                end if
             enddo
             enddo
             enddo
          endif
       case('U')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                velx_org(1:2,i,j) = 0.0_RP
                do k = 1, knum
                   velx_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( velx_org(k+2,i,j) - missval ) < EPS ) then
                      velx_org(k+2,i,j) = UNDEF
                   end if
                enddo
                if(dims(1)>knum)then
                   do k = knum+1, dims(1)
                      velx_org(k+2,i,j) = velx_org(knum+2,i,j)
                   enddo
                endif
             enddo
             enddo
          endif
       case('V')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                vely_org(1:2,i,j) = 0.0_RP
                do k = 1, knum
                   vely_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( vely_org(k+2,i,j) - missval ) < EPS ) then
                      vely_org(k+2,i,j) = UNDEF
                   end if
                enddo
                if(dims(1)>knum)then
                   do k = knum+1, dims(1)
                      vely_org(k+2,i,j) = vely_org(knum+2,i,j)
                   enddo
                endif
             enddo
             enddo
          endif
       case('T')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   temp_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( temp_org(k+2,i,j) - missval ) < EPS ) then
                      temp_org(k+2,i,j) = UNDEF
                   end if
                enddo
                if(dims(1)>knum)then
                   do k = knum+1, dims(1)
                      temp_org(k+2,i,j) = temp_org(knum+2,i,j)
                   enddo
                endif
             enddo
             enddo
          endif
       case('HGT')
          if(dims(1)/=knum)then
             write(*,*) 'xxx The number of levels for HGT must be same as plevs! knum:',knum,'> outer_nz:',dims(1)
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "levels" ) then
             if(dims(1)/=lnum)then
                write(*,*) 'xxx lnum must be same as the outer_nz for HGT! ',dims(1),lnum
                call PRC_MPIstop
             endif
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, dims(1)
                   cz_org(k+2,i,j) = real(lvars(k), kind=RP)
                enddo
                cz_org(1,i,j) = 0.0_RP
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),dims(1),nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, dims(1)
                   cz_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( cz_org(k+2,i,j) - missval ) < EPS ) then
                      cz_org(k+2,i,j) = UNDEF
                   end if
                enddo
                cz_org(1,i,j) = 0.0_RP
             enddo
             enddo
          endif
       case('QV')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QV) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QV) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QV) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QV) = qtrc_org(3,i,j,I_QV)
             enddo
             enddo
             if( dims(1)>knum ) then
                select case( upper_qv_type )
                case("COPY")
                   do j = 1, dims(3)
                   do i = 1, dims(2)
                   do k = knum+1, dims(1)
                      qtrc_org(k+2,i,j,I_QV) = qtrc_org(knum+2,i,j,I_QV)
                   enddo
                   enddo
                   enddo
                case("ZERO")
                   ! do nothing
                case default
                   write(*,*) 'xxx upper_qv_type in PARAM_MKINIT_REAL_GrADS is invalid! ', upper_qv_type
                   call PRC_MPIstop
                end select
             endif
          endif
       case('QC')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QC) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QC) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QC) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QC) = qtrc_org(3,i,j,I_QC)
             enddo
             enddo
          endif
       case('QR')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QR) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QR) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QR) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QR) = qtrc_org(3,i,j,I_QR)
             enddo
             enddo
          endif
       case('QI')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QI) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QI) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QI) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QI) = qtrc_org(3,i,j,I_QI)
             enddo
             enddo
          endif
       case('QS')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QS) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QS) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QS) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QS) = qtrc_org(3,i,j,I_QS)
             enddo
             enddo
          endif
       case('QG')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QG) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QG) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QG) = UNDEF
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QG) = qtrc_org(3,i,j,I_QG)
             enddo
             enddo
          endif
       case('RH')
          if (data_available(Ia_qv,1)) cycle  ! use QV
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,yrev,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   qtrc_org(k+2,i,j,I_QV) = real(gdata3D(i,j,k), kind=RP)
                   ! replace missval with UNDEF
                   if( abs( qtrc_org(k+2,i,j,I_QV) - missval ) < EPS ) then
                      qtrc_org(k+2,i,j,I_QV) = UNDEF
                   else
                      rhprs_org(k+2,i,j) = qtrc_org(k+2,i,j,I_QV) / 100.0_RP   ! relative humidity
                      call psat( p_sat, temp_org(k+2,i,j) )                    ! satulation pressure
                      qm = EPSvap * rhprs_org(k+2,i,j) * p_sat &
                         / ( pres_org(k+2,i,j) - rhprs_org(k+2,i,j) * p_sat )  ! mixing ratio
                      qtrc_org(k+2,i,j,I_QV) = qm / ( 1.0_RP + qm )            ! specific humidity
                   end if
                enddo
                qtrc_org(1:2,i,j,I_QV) = qtrc_org(3,i,j,I_QV)
             enddo
             enddo
             if( dims(3)>knum ) then
                select case( upper_qv_type )
                case("COPY")
                   do j = 1, dims(3)
                   do i = 1, dims(2)
                   do k = knum+1, dims(1)
                      rhprs_org(k+2,i,j) = rhprs_org(knum+2,i,j)              ! relative humidity
                      call psat( p_sat, temp_org(k+2,i,j) )                   ! satulated specific humidity
                      qm = EPSvap * rhprs_org(k+2,i,j) * p_sat &
                         / ( pres_org(k+2,i,j) - rhprs_org(k+2,i,j) * p_sat ) ! mixing ratio
                      qtrc_org(k+2,i,j,I_QV) = qm / ( 1.0_RP + qm )           ! specific humidity
                      qtrc_org(k+2,i,j,I_QV) = min(qtrc_org(k+2,i,j,I_QV),qtrc_org(k+1,i,j,I_QV))
                   enddo
                   enddo
                   enddo
                case("ZERO")
                   ! do nothing
                case default
                   write(*,*) 'xxx upper_qv_type in PARAM_MKINIT_REAL_GrADS is invalid! ', upper_qv_type
                   call PRC_MPIstop
                end select
             endif
          endif
       case('MSLP')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                pres_org(1,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( pres_org(1,i,j) - missval ) < EPS ) then
                   pres_org(1,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       case('PSFC')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                pres_org(2,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( pres_org(2,i,j) - missval ) < EPS ) then
                   pres_org(2,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       case('U10')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                velx_org(2,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( velx_org(2,i,j) - missval ) < EPS ) then
                   velx_org(2,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       case('V10')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                vely_org(2,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( vely_org(2,i,j) - missval ) < EPS ) then
                   vely_org(2,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       case('T2')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                temp_org(2,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( temp_org(2,i,j) - missval ) < EPS ) then
                   temp_org(2,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       case('Q2')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                qtrc_org(2,i,j,I_QV) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( qtrc_org(2,i,j,I_QV) - missval ) < EPS ) then
                   qtrc_org(2,i,j,I_QV) = UNDEF
                end if
             enddo
             enddo
          endif
       case('RH2')
          if (data_available(Ia_q2,1)) cycle  ! use QV
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                qtrc_org(2,i,j,I_QV) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( qtrc_org(2,i,j,I_QV) - missval ) < EPS ) then
                   qtrc_org(2,i,j,I_QV) = UNDEF
                else
                   rhsfc = qtrc_org(2,i,j,I_QV) / 100.0_RP
                   call psat( p_sat, temp_org(2,i,j) )         ! satulation pressure
                   qm = EPSvap * rhsfc * p_sat &
                      / ( pres_org(2,i,j) - rhsfc * p_sat )    ! mixing ratio
                   qtrc_org(2,i,j,I_QV) = qm / ( 1.0_RP + qm ) ! specific humidity
                end if
             enddo
             enddo
          endif
       case('TOPO')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,yrev,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                cz_org(2,i,j) = real(gdata2D(i,j), kind=RP)
                ! replace missval with UNDEF
                if( abs( cz_org(2,i,j) - missval ) < EPS ) then
                   cz_org(2,i,j) = UNDEF
                end if
             enddo
             enddo
          endif
       end select
    enddo loop_InputAtomGrADS

    lm_layer(:,:) = 3

    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 3, dims(1)+2
      ! search the lowermost layer excluding UNDEF
      if( abs( pres_org(k,i,j) - UNDEF ) < EPS ) then
        lm_layer(i,j) = k + 1
      else
        exit
      end if
    end do
    end do
    end do

    RovCP = Rdry / CPdry
    CPovR = CPdry / Rdry

    if ( data_available(Ia_t2,1) .and. data_available(Ia_ps,1) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          pott(i,j) = temp_org(2,i,j) * (P00/pres_org(2,i,j))**RovCP
       end do
       end do
    else
       do j = 1, dims(3)
       do i = 1, dims(2)
          k = lm_layer(i,j)
          pott(i,j) = temp_org(k,i,j) * (P00/pres_org(k,i,j))**RovCP
       end do
       end do
    end if

    if ( .not. data_available(Ia_t2,1) ) then
       if ( data_available(Ia_ps,1) ) then
          do j = 1, dims(3)
          do i = 1, dims(2)
             temp_org(2,i,j) = pott(i,j) * (pres_org(2,i,j)/P00)**RovCP
          end do
          end do
       else
          if ( data_available(Ia_topo,1) ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
                k = lm_layer(i,j)
                temp_org(2,i,j) = temp_org(k,i,j) &
                                + LAPS * (cz_org(k,i,j)-cz_org(2,i,j))
             end do
             end do
          else
             do j = 1, dims(3)
             do i = 1, dims(2)
                k = lm_layer(i,j)
                temp_org(2,i,j) = temp_org(k,i,j)
             end do
             end do
          end if
       end if
    end if

    if ( .not. data_available(Ia_ps,1) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          pres_org(2,i,j) = P00 * (temp_org(2,i,j)/pott(i,j))**CPovR
       end do
       end do
    end if

    if ( data_available(Ia_slp,1) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          temp_org(1,i,j) = pott(i,j) * (pres_org(1,i,j)/P00)**RovCP
       end do
       end do
    else
       if ( data_available(Ia_t2,1) .and. data_available(Ia_topo,1) ) then
          do j = 1, dims(3)
          do i = 1, dims(2)
             temp_org(1,i,j) = temp_org(2,i,j) + LAPS * cz_org(2,i,j)
          end do
          end do
       else
          do j = 1, dims(3)
          do i = 1, dims(2)
             k = lm_layer(i,j)
             temp_org(1,i,j) = temp_org(k,i,j) + LAPS * cz_org(k,i,j)
          end do
          end do
       end if
       do j = 1, dims(3)
       do i = 1, dims(2)
          pres_org(1,i,j) = P00 * (temp_org(1,i,j)/pott(i,j))**CPovR
       end do
       end do
    end if

    if ( .not. data_available(Ia_topo,1) ) then
       ! guess surface height (elevation)
       do j = 1, dims(3)
       do i = 1, dims(2)
          k = lm_layer(i,j)
          if ( pres_org(2,i,j) < pres_org(1,i,j) ) then
             lp2 = log( pres_org(2,i,j) / pres_org(1,i,j) )
          else
             lp2 = 1.0_RP
          end if
          if ( pres_org(k,i,j) < pres_org(1,i,j) ) then
             lp3 = log( pres_org(k,i,j) / pres_org(1,i,j) )
          else
             lp3 = 1.0_RP
          end if
          cz_org(2,i,j) = max( 0.0_RP, cz_org(k,i,j) * lp2 / lp3 )
       end do
       end do
    end if

    velz_org = 0.0_RP


    ! check verticaly extrapolated data in outer model
    if( pressure_coordinates ) then
      do j = 1, dims(3)
      do i = 1, dims(2)
      do k = 3, dims(1)+2
        if( pres_org(k,i,j) > pres_org(2,i,j) ) then ! if Pressure is larger than Surface pressure
          velx_org(k,i,j)   = velx_org(2,i,j)
          vely_org(k,i,j)   = vely_org(2,i,j)
          temp_org(k,i,j)   = temp_org(2,i,j)
          qtrc_org(k,i,j,:) = qtrc_org(2,i,j,:)
          cz_org  (k,i,j)   = cz_org  (2,i,j)
        end if
      enddo
      enddo
      enddo
    else
      do j = 1, dims(3)
      do i = 1, dims(2)
      do k = 3, dims(1)+2
        if( abs( pres_org(k,i,j) - UNDEF ) < EPS ) pres_org(k,i,j) = pres_org(2,i,j)
        if( abs( velx_org(k,i,j) - UNDEF ) < EPS ) velx_org(k,i,j) = velx_org(2,i,j)
        if( abs( vely_org(k,i,j) - UNDEF ) < EPS ) vely_org(k,i,j) = vely_org(2,i,j)
        if( abs( temp_org(k,i,j) - UNDEF ) < EPS ) temp_org(k,i,j) = temp_org(2,i,j)
        do iq = 1, QA
          if( abs( qtrc_org(k,i,j,iq) - UNDEF ) < EPS ) qtrc_org(k,i,j,iq) = 0.0_RP
        end do
      enddo
      enddo
      enddo
    end if


    !do it = 1, nt
    !   k=1 ; j=int(dims(3)/2) ; i=int(dims(2)/2) ; iq = 1
    !   write(*,*) "read 3D grads data",i,j,k
    !   write(*,*) "lon_org    ",lon_org   (i,j)/D2R
    !   write(*,*) "lat_org    ",lat_org   (i,j)/D2R
    !   write(*,*) "pres_org   ",pres_org  (k,i,j)
    !   write(*,*) "usfc_org   ",usfc_org  (i,j)
    !   write(*,*) "vsfc_org   ",vsfc_org  (i,j)
    !   write(*,*) "tsfc_org   ",tsfc_org  (i,j)
    !   write(*,*) "qsfc_org   ",qsfc_org  (i,j,iq)
    !   write(*,*) "rhsfc_org  ",rhsfc_org (i,j)
    !   write(*,*) "velx_org   ",velx_org  (k,i,j)
    !   write(*,*) "vely_org   ",vely_org  (k,i,j)
    !   write(*,*) "temp_org   ",temp_org  (k,i,j)
    !   write(*,*) "hgt_org    ",hgt_org   (k,i,j)
    !   write(*,*) "qtrc_org   ",qtrc_org  (k,i,j,iq)
    !   write(*,*) "rhprs_org  ",rhprs_org (k,i,j)
    !enddo

    return
  end subroutine ParentAtomInputGrADS

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupGrADS( &
       ldims,                       & ! (out)
       use_waterratio,              & ! (out)
       use_file_landwater,          & ! (in)
       basename                     )
    implicit none

    integer,          intent(out) :: ldims(3)
    logical,          intent(out) :: use_waterratio
    logical,          intent(in)  :: use_file_landwater ! use landwater data from files
    character(len=*), intent(in)  :: basename

    integer             :: ielem
    integer             :: k, n

    integer             :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Land Input File Type: GrADS format'

    !--- initialization
    use_waterratio = .false.

    if ( len_trim(basename) == 0 ) then
       write(*,*) 'xxx [realinput_grads] "BASEMAAME" is not specified in "PARAM_MKINIT_REAL_ATOMS"!', trim(basename)
       call PRC_MPIstop
    endif

    !--- read namelist
    io_fid_grads_nml = IO_get_available_fid()
    open( io_fid_grads_nml,       &
         file   = trim(basename), &
         form   = 'formatted',    &
         status = 'old',          &
         action = 'read',         &
         iostat = ierr            )
    if ( ierr /= 0 ) then
       write(*,*) 'xxx [realinput_grads] Input file is not found! ', trim(basename)
       call PRC_MPIstop
    endif

    read(io_fid_grads_nml,nml=nml_grads_grid,iostat=ierr)
    if( ierr /= 0 ) then !--- missing or fatal error
       write(*,*) 'xxx [realinput_grads] Not appropriate names in nml_grads_grid in ', trim(basename),'. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=nml_grads_grid)

    ! land
    ldims(1) = outer_nl ! soil_layers_stag
    if(outer_nx_sfc > 0)then
       ldims(2) = outer_nx_sfc
    else
       ldims(2) = outer_nx
       outer_nx_sfc = outer_nx
    endif
    if(outer_ny_sfc > 0)then
       ldims(3) = outer_ny_sfc
    else
       ldims(3) = outer_ny
       outer_ny_sfc = outer_ny
    endif

    allocate( gland2D( ldims(2), ldims(3)           ) )
    allocate( gland3D( ldims(2), ldims(3), ldims(1) ) )

    call read_namelist( &
         grads_item(:,2),     & ! (out)
         grads_fname(:,2),    & ! (out)
         grads_dtype(:,2),    & ! (out)
         grads_swpoint(:,2),  & ! (out)
         grads_dd(:,2),       & ! (out)
         grads_lnum(:,2),     & ! (out)
         grads_lvars(:,:,2),  & ! (out)
         grads_startrec(:,2), & ! (out)
         grads_totalrec(:,2), & ! (out)
         grads_knum(:,2),     & ! (out)
         grads_yrev(:,2),     & ! (out)
         grads_fendian(:,2),  & ! (out)
         grads_missval(:,2),  & ! (out)
         data_available(:,2), & ! (out)
         item_list_land,      & ! (in)
         num_item_list_land,  & ! (in)
         basename,            & ! (in)
         io_fid_grads_nml     ) ! (in)

    close( io_fid_grads_nml )

    do ielem = 1, num_item_list_land
       item  = item_list_land(ielem)
       !--- check data
       select case(trim(item))
       case('TOPO','TOPO_sfc', 'lsmask')
          if ( .not. data_available(ielem,2) ) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: ',trim(item),' is not found & not used.'
             cycle
          endif
       case('lon', 'lat', 'lon_sfc', 'lat_sfc')
          cycle
       case('SMOISVC', 'SMOISDS')
          if ( use_file_landwater ) then
             if (.not. data_available(Il_smoisvc,2) .and. .not. data_available(Il_smoisds,2)) then
                write(*,*) 'xxx [realinput_grads] Not found in grads namelist! : ',trim(item_list_land(ielem))
                call PRC_MPIstop
             end if
             use_waterratio =  data_available(Il_smoisds,2)
          else
             cycle
          end if
       case default ! llev, SKINT, STEMP
          if ( .not. data_available(ielem,2) ) then
             write(*,*) 'xxx [realinput_grads] Not found in grads namelist! : ',trim(item_list_land(ielem))
             call PRC_MPIstop
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
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R,   &
       TEM00 => CONST_TEM00, &
       EPS   => CONST_EPS
    implicit none

    real(RP), intent(out) :: tg_org    (:,:,:)
    real(RP), intent(out) :: strg_org  (:,:,:)
    real(RP), intent(out) :: smds_org  (:,:,:)
    real(RP), intent(out) :: lst_org   (:,:)
    real(RP), intent(out) :: llon_org  (:,:)
    real(RP), intent(out) :: llat_org  (:,:)
    real(RP), intent(out) :: lz_org    (:)
    real(RP), intent(out) :: topo_org(:,:)
    real(RP), intent(out) :: lmask_org(:,:)

    character(len=*), intent(in) :: basename_num
    integer,          intent(in) :: ldims(3)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: nt
    ! ----------------------------------------------------------------

    !> grads data
    character(len=H_LONG) :: gfile

    real(RP) :: qvsat, qm

    integer :: i, j, k, ielem, n
    integer :: ierr

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[LandInputGrADS]'


    loop_InputLandGrADS : do ielem = 1, num_item_list_land

       item  = item_list_land(ielem)

       dtype    = grads_dtype   (ielem,2)
       fname    = grads_fname   (ielem,2)
       lnum     = grads_lnum    (ielem,2)
       missval  = grads_missval (ielem,2)

       if ( grads_knum(ielem,2) > 0 )then
          knum = grads_knum(ielem,2)
       else
          knum = ldims(1)
       endif

       select case(trim(dtype))
       case("linear")
          swpoint = grads_swpoint (ielem,2)
          dd      = grads_dd      (ielem,2)
          if( (abs(swpoint-large_number_one)<EPS).or.(abs(dd-large_number_one)<EPS) )then
             write(*,*) 'xxx "swpoint" is required in grads namelist! ',swpoint
             write(*,*) 'xxx "dd"      is required in grads namelist! ',dd
             call PRC_MPIstop
          endif
       case("levels")
          if ( lnum < 0 )then
             write(*,*) 'xxx "lnum" in grads namelist is required for levels data! '
             call PRC_MPIstop
          endif
          do k=1, lnum
             lvars(k)=grads_lvars(k,ielem,2)
          enddo
          if(abs(lvars(1)-large_number_one)<EPS)then
             write(*,*) 'xxx "lvars" must be specified in grads namelist for levels data!',(lvars(k),k=1,lnum)
             call PRC_MPIstop
          endif
       case("map")
          startrec = grads_startrec(ielem,2)
          totalrec = grads_totalrec(ielem,2)
          fendian  = grads_fendian (ielem,2)
          yrev     = grads_yrev (ielem,2)
          if( (startrec<0).or.(totalrec<0) )then
             write(*,*) 'xxx "startrec" is required in grads namelist! ',startrec
             write(*,*) 'xxx "totalrec" is required in grads namelist! ',totalrec
             call PRC_MPIstop
          endif
          ! get file_io
          if(io_fid_grads_data < 0)then
             io_fid_grads_data = IO_get_available_fid()
          endif
          gfile=trim(fname)//trim(basename_num)//'.grd'
          if( len_trim(fname)==0 )then
             write(*,*) 'xxx "fname" is required in grads namelist for map data! ',trim(fname)
             call PRC_MPIstop
          endif
       end select

       ! read data
       select case(trim(item))
       case("lsmask")
          if ( data_available(Il_lsmask,2) ) then
             if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,1,item,startrec,totalrec,yrev,gland2D)
                lmask_org(:,:) = real(gland2D(:,:), kind=RP)
             endif
          else
             lmask_org = UNDEF
          end if
       case("lon")
          if ( .not. data_available(Il_lon_sfc,2) ) then
             if ( ldims(2).ne.outer_nx .or. ldims(3).ne.outer_ny ) then
                write(*,*) 'xxx namelist of "lon_sfc" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sfc! ', outer_nx, ldims(2)
                write(*,*) '                          : outer_ny and outer_ny_sfc! ', outer_ny, ldims(3)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, ldims(3)
                do i = 1, ldims(2)
                   llon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,1,item,startrec,totalrec,yrev,gland2D)
                llon_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lon_sfc")
          if ( .not. data_available(Il_lon_sfc,2) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,1,item,startrec,totalrec,yrev,gland2D)
             llon_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
          endif
       case("lat")
          if ( .not. data_available(Il_lat_sfc,2) ) then
             if ( ldims(2).ne.outer_nx .or. ldims(3).ne.outer_ny ) then
                write(*,*) 'xxx namelist of "lat_sfc" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sfc! ', outer_nx, ldims(2)
                write(*,*) '                          : outer_ny and outer_ny_sfc! ', outer_nx, ldims(3)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, ldims(3)
                do i = 1, ldims(2)
                   llat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,1,item,startrec,totalrec,yrev,gland2D)
                llat_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lat_sfc")
          if ( .not. data_available(Il_lat_sfc,2) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                llat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,1,item,startrec,totalrec,yrev,gland2D)
             llat_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
          endif
       case("llev")
          if(ldims(1)/=knum)then
             write(*,*) 'xxx "knum" must be equal to outer_nl for llev. knum:',knum,'> outer_nl:',ldims(1)
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "levels" ) then
             if(ldims(1)/=lnum)then
                write(*,*) 'xxx lnum must be same as the outer_nl for llev! ',ldims(1),lnum
                call PRC_MPIstop
             endif
             do k = 1, ldims(1)
                lz_org(k) = real(lvars(k), kind=RP)
             enddo
!          else if ( trim(dtype) == "map" ) then
!             call read_grads_file_3d(io_fid_grads_data,gfile,ldims(2),ldims(3),ldims(1),nt,item,startrec,totalrec,yrev,gland)
!             do j = 1, ldims(3)
!             do i = 1, ldims(2)
!             do k = 1, ldims(1)
!                lz_org(k,i,j) = real(gland(i,j,k), kind=RP)
!             enddo
!             enddo
!             enddo
          endif
       case('STEMP')
          if(ldims(1)/=knum)then
             write(*,*) 'xxx The number of levels for STEMP must be same as llevs! ',ldims(1),knum
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,ldims(2),ldims(3),ldims(1),nt,item,startrec,totalrec,yrev,gland3D)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
             do k = 1, ldims(1)
                if ( abs(gland3D(i,j,k)-missval) < EPS ) then
                   tg_org(k,i,j) = UNDEF
                else
                   tg_org(k,i,j) = real(gland3D(i,j,k), kind=RP)
                end if
             enddo
             enddo
             enddo
          endif
       case('SMOISVC')
          if ( use_file_landwater ) then
             if(ldims(1)/=knum)then
                write(*,*) 'xxx The number of levels for SMOISVC must be same as llevs! ',ldims(1),knum
                call PRC_MPIstop
             endif
             if ( trim(dtype) == "map" ) then
                call read_grads_file_3d(io_fid_grads_data,gfile,ldims(2),ldims(3),ldims(1),nt,item,startrec,totalrec,yrev,gland3D)
                do j = 1, ldims(3)
                do i = 1, ldims(2)
                do k = 1, ldims(1)
                   if ( abs(gland3D(i,j,k)-missval) < EPS ) then
                      strg_org(k,i,j) = UNDEF
                   else
                      strg_org(k,i,j) = real(gland3D(i,j,k), kind=RP)
                   end if
                enddo
                enddo
                enddo
             endif
          endif
       case('SMOISDS')
          if ( use_file_landwater ) then
             if(ldims(1)/=knum)then
                write(*,*) 'xxx The number of levels for SMOISDS must be same as llevs! ',ldims(1),knum
                call PRC_MPIstop
             endif
             if ( trim(dtype) == "map" ) then
                call read_grads_file_3d(io_fid_grads_data,gfile,ldims(2),ldims(3),ldims(1),nt,item,startrec,totalrec,yrev,gland3D)
                do j = 1, ldims(3)
                do i = 1, ldims(2)
                do k = 1, ldims(1)
                   if ( abs(gland3D(i,j,k)-missval) < EPS ) then
                      smds_org(k,i,j) = UNDEF
                   else
                      smds_org(k,i,j) = real(gland3D(i,j,k), kind=RP)
                   end if
                enddo
                enddo
                enddo
             endif
          endif
       case('SKINT')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,nt,item,startrec,totalrec,yrev,gland2D)
             do j = 1, ldims(3)
             do i = 1, ldims(2)
                if ( abs(gland2D(i,j)-missval) < EPS ) then
                   lst_org(i,j) = UNDEF
                else
                   lst_org(i,j) = real(gland2D(i,j), kind=RP)
                end if
             enddo
             enddo
          endif
       case('TOPO')
          if ( .not. data_available(Il_topo_sfc,2) ) then
             if ( ldims(2)==outer_nx .or. ldims(3)==outer_ny ) then
                if ( trim(dtype) == "map" ) then
                   call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,nt,item,startrec,totalrec,yrev,gland2D)
                   do j = 1, ldims(3)
                   do i = 1, ldims(2)
                      if ( abs(gland2D(i,j)-missval) < EPS ) then
                         topo_org(i,j) = UNDEF
                      else
                         topo_org(i,j) = real(gland2D(i,j), kind=RP)
                      end if
                   enddo
                   enddo
                end if
             else
                topo_org = UNDEF
             endif
          end if
       case('TOPO_sfc')
          if ( data_available(Il_topo_sfc,2) ) then
             if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,ldims(2),ldims(3),1,nt,item,startrec,totalrec,yrev,gland2D)
                do j = 1, ldims(3)
                do i = 1, ldims(2)
                   if ( abs(gland2D(i,j)-missval) < EPS ) then
                      topo_org(i,j) = UNDEF
                   else
                      topo_org(i,j) = real(gland2D(i,j), kind=RP)
                   end if
                enddo
                enddo
             endif
          else if ( .not. data_available(Il_topo,2) ) then
             topo_org = UNDEF
          endif
       end select
    enddo loop_InputLandGrADS

    !do it = 1, nt
    !   i=int(ldims(2)/2) ; j=int(ldims(3)/2)
    !   write(*,*) "read 2D grads data",ldims(2),ldims(3),i,j,it
    !   write(*,*) "lon_org    ",lon_org  (i,j)
    !   write(*,*) "lat_org    ",lat_org  (i,j)
    !   write(*,*) "lst_org  ",lst_org(i,j)
    !   do k=1,dims(7)
    !      write(*,*) "tg_org    ",tg_org   (k,i,j)," k= ",k
    !      write(*,*) "strg_org  ",strg_org (k,i,j)," k= ",k
    !   enddo
    !enddo

    return
  end subroutine ParentLandInputGrADS

  !-----------------------------------------------------------------------------
  !> Ocean Setup
  subroutine ParentOceanSetupGrADS( &
       odims,   & ! (out)
       timelen, & ! (out)
       basename ) ! (in)
    implicit none

    integer,          intent(out) :: odims(2)
    integer,          intent(out) :: timelen
    character(len=*), intent(in)  :: basename

    character(len=H_LONG) :: grads_ctl
    integer             :: ielem
    integer             :: n

    integer             :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Ocean Input File Type: GrADS format'

    !--- read namelist

    if ( len_trim(basename) == 0 ) then
       grads_ctl = "namelist.grads_boundary"
    else
       grads_ctl = basename
    endif

    !--- read namelist
    io_fid_grads_nml = IO_get_available_fid()
    open( io_fid_grads_nml,                    &
         file   = trim(grads_ctl), &
         form   = 'formatted',                   &
         status = 'old',                         &
         action = 'read',                        &
         iostat = ierr                           )
    if ( ierr /= 0 ) then
       write(*,*) 'xxx [realinput_grads] Input file is not found! ', trim(grads_ctl)
       call PRC_MPIstop
    endif

    read(io_fid_grads_nml,nml=nml_grads_grid,iostat=ierr)
    if( ierr /= 0 ) then !--- missing or fatal error
       write(*,*) 'xxx [realinput_grads] Not appropriate names in nml_grads_grid in ', trim(grads_ctl),'. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=nml_grads_grid)

    timelen = 0        ! will be replaced later

    ! sst
    if(outer_nx_sst > 0)then
       odims(1) = outer_nx_sst
    else if (outer_nx_sfc > 0) then
       odims(1) = outer_nx_sfc
       outer_nx_sst = outer_nx_sfc
    else
       odims(1) = outer_nx
       outer_nx_sst = outer_nx
    endif
    if(outer_ny_sst > 0)then
       odims(2) = outer_ny_sst
    else if(outer_ny_sfc > 0)then
       odims(2) = outer_ny_sfc
       outer_ny_sst = outer_ny_sfc
    else
       odims(2) = outer_ny
       outer_ny_sst = outer_ny
    endif

    allocate( gsst2D ( odims(1), odims(2)        ) )


    call read_namelist( &
         grads_item(:,3),     & ! (out)
         grads_fname(:,3),    & ! (out)
         grads_dtype(:,3),    & ! (out)
         grads_swpoint(:,3),  & ! (out)
         grads_dd(:,3),       & ! (out)
         grads_lnum(:,3),     & ! (out)
         grads_lvars(:,:,3),  & ! (out)
         grads_startrec(:,3), & ! (out)
         grads_totalrec(:,3), & ! (out)
         grads_knum(:,3),     & ! (out)
         grads_yrev(:,3),     & ! (out)
         grads_fendian(:,3),  & ! (out)
         grads_missval(:,3),  & ! (out)
         data_available(:,3), & ! (out)
         item_list_ocean,     & ! (in)
         num_item_list_ocean, & ! (in)
         grads_ctl,           & ! (in)
         io_fid_grads_nml     ) ! (in)

    close( io_fid_grads_nml )

    do ielem = 1, num_item_list_ocean
       item  = item_list_ocean(ielem)
       !--- check data
       select case(trim(item))
       case('lsmask','lsmask_sst')
          if ( .not. data_available(Io_lsmask,3) .and. .not. data_available(Io_lsmask_sst,3) ) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: ',trim(item),' is not found & not used.'
             cycle
          endif
       case('lon', 'lat', 'lon_sfc', 'lat_sfc', 'lon_sst', 'lat_sst')
          cycle
       case('SST')
          if (.not. data_available(Io_sst,3) .and. .not. data_available(Io_skint,3) ) then
             write(*,*) 'xxx [realinput_grads] SST and SKINT are found in grads namelist!'
             call PRC_MPIstop
          endif
          if (.not. data_available(Io_sst,3)) then
             if( IO_L ) write(IO_FID_LOG,*) 'warning: SST is found in grads namelist. SKINT is used in place of SST.'
             cycle
          endif
       case('SKINT')
          cycle
       case default !
          if ( .not. data_available(ielem,3) ) then
             write(*,*) 'xxx [realinput_grads/ParentOceanSetupGrADS] Not found in grads namelist! : ', &
                        trim(item_list_ocean(ielem))
             call PRC_MPIstop
          endif
       end select

    end do

    return
  end subroutine ParentOceanSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentOceanOpenGrADS
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanOpenGrADS]'

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
    implicit none

    real(RP), intent(out) :: tw_org   (:,:)
    real(RP), intent(out) :: sst_org  (:,:)
    real(RP), intent(out) :: omask_org(:,:)
    real(RP), intent(out) :: olon_org (:,:)
    real(RP), intent(out) :: olat_org (:,:)

    character(len=*), intent(in) :: basename_num
    integer,          intent(in) :: odims(2)
    integer,          intent(in) :: nt
    ! ----------------------------------------------------------------

    !> grads data
    character(len=H_LONG) :: gfile

    real(RP) :: qvsat, qm

    integer :: i, j, ielem, n
    integer :: ierr

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[OceanInputGrADS]'


    loop_InputOceanGrADS : do ielem = 1, num_item_list_ocean

       item  = item_list_ocean(ielem)

       dtype    = grads_dtype   (ielem,3)
       fname    = grads_fname   (ielem,3)
       lnum     = grads_lnum    (ielem,3)
       missval  = grads_missval (ielem,3)

       select case(trim(dtype))
       case("linear")
          swpoint = grads_swpoint (ielem,3)
          dd      = grads_dd      (ielem,3)
          if( (abs(swpoint-large_number_one)<EPS).or.(abs(dd-large_number_one)<EPS) )then
             write(*,*) 'xxx "swpoint" is required in grads namelist! ',swpoint
             write(*,*) 'xxx "dd"      is required in grads namelist! ',dd
             call PRC_MPIstop
          endif
       case("levels")
          write(*,*) 'xxx "lnum" in grads namelist is invalid for ocean data'
          call PRC_MPIstop
       case("map")
          startrec = grads_startrec(ielem,3)
          totalrec = grads_totalrec(ielem,3)
          fendian  = grads_fendian (ielem,3)
          yrev     = grads_yrev (ielem,3)
          if( (startrec<0).or.(totalrec<0) )then
             write(*,*) 'xxx "startrec" is required in grads namelist! ',startrec
             write(*,*) 'xxx "totalrec" is required in grads namelist! ',totalrec
             call PRC_MPIstop
          endif
          ! get file_io
          if(io_fid_grads_data < 0)then
             io_fid_grads_data = IO_get_available_fid()
          endif
          gfile=trim(fname)//trim(basename_num)//'.grd'
          if( len_trim(fname)==0 )then
             write(*,*) 'xxx "fname" is required in grads namelist for map data! ',trim(fname)
             call PRC_MPIstop
          endif
       end select

       ! read data
       select case(trim(item))
       case("lsmask")
          if ( .not. data_available(Io_lsmask_sst,3) .and. data_available(Io_lsmask,3) ) then
             if ( odims(1)==outer_nx_sfc .and. odims(2)==outer_ny_sfc ) then
                if ( trim(dtype) == "map" ) then
                   call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                   omask_org(:,:) = real(gsst2D(:,:), kind=RP)
                endif
             else
                omask_org = UNDEF
             end if
          end if
       case("lsmask_sst")
          if ( data_available(Io_lsmask_sst,3) ) then
             if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                omask_org(:,:) = real(gsst2D(:,:), kind=RP)
             endif
          else if ( .not. data_available(Io_lsmask,3) ) then
             omask_org = UNDEF
          end if
       case("lon")
          if ( .not. data_available(Io_lon_sst,3) .and. .not. data_available(Io_lon_sfc,3) ) then
             if ( odims(1).ne.outer_nx .or. odims(2).ne.outer_ny ) then
                write(*,*) 'xxx namelist of "lon_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sst! ', outer_nx, odims(1)
                write(*,*) '                          : outer_ny and outer_ny_sst! ', outer_ny, odims(2)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, odims(2)
                do i = 1, odims(1)
                   olon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                olon_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lon_sfc")
          if ( .not. data_available(Io_lon_sst,3) .and. data_available(Io_lon_sfc,3) ) then
             if ( odims(1).ne.outer_nx_sfc .or. odims(2).ne.outer_ny_sfc ) then
                write(*,*) 'xxx namelist of "lon_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx_sfc and outer_nx_sst! ', outer_nx_sfc, odims(1)
                write(*,*) '                          : outer_ny_sfc and outer_ny_sst! ', outer_ny_sfc, odims(2)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, odims(2)
                do i = 1, odims(1)
                   olon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                olon_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lon_sst")
          if ( .not. data_available(Io_lon_sst,3) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, odims(2)
             do i = 1, odims(1)
                olon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
             olon_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
          endif
       case("lat")
          if ( .not. data_available(Io_lat_sfc,3) .and. .not. data_available(Io_lat_sst,3) ) then
             if ( odims(1).ne.outer_nx .or. odims(2).ne.outer_ny ) then
                write(*,*) 'xxx namelist of "lat_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sst! ', outer_nx, odims(1)
                write(*,*) '                          : outer_ny and outer_ny_sst! ', outer_ny, odims(2)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, odims(2)
                do i = 1, odims(1)
                   olat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                olat_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lat_sfc")
          if ( .not. data_available(Io_lat_sst,3) .and. data_available(Io_lat_sfc,3) ) then
             if ( odims(1).ne.outer_nx .or. odims(1).ne.outer_ny ) then
                write(*,*) 'xxx namelist of "lat_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx_sfc and outer_nx_sst! ', outer_nx_sfc, odims(1)
                write(*,*) '                          : outer_ny_sfc and outer_ny_sst! ', outer_ny_sfc, odims(2)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, odims(2)
                do i = 1, odims(1)
                   olat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
                olat_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
             endif
          end if
       case("lat_sst")
          if ( .not. data_available(Io_lat_sst,3) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, odims(2)
             do i = 1, odims(1)
                olat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,1,item,startrec,totalrec,yrev,gsst2D)
             olat_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
          endif
       case('SKINT')
          if ( .not. data_available(Io_sst,3) ) then
             if ( odims(1).ne.outer_nx_sfc .or. odims(2).ne.outer_ny_sfc ) then
                write(*,*) 'xxx dimsntion is different: outer_nx_sst/outer_nx_sfc and outer_nx_sst! ', odims(1), outer_nx_sfc
                write(*,*) '                          : outer_ny_sst/outer_ny_sfc and outer_ny_sst! ', odims(2), outer_ny_sfc
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,nt,item,startrec,totalrec,yrev,gsst2D)
                do j = 1, odims(2)
                do i = 1, odims(1)
                   if ( abs(gsst2D(i,j)-missval) < EPS ) then
                      sst_org(i,j) = UNDEF
                   else
                      sst_org(i,j) = real(gsst2D(i,j), kind=RP)
                   end if
                enddo
                enddo
             end if
          endif
       case('SST')
          if ( .not. data_available(Io_sst,3) ) cycle
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,odims(1),odims(2),1,nt,item,startrec,totalrec,yrev,gsst2D)
             do j = 1, odims(2)
             do i = 1, odims(1)
                if ( abs(gsst2D(i,j)-missval) < EPS ) then
                   sst_org(i,j) = UNDEF
                else
                   sst_org(i,j) = real(gsst2D(i,j), kind=RP)
                end if
             enddo
             enddo
          end if
       end select
    enddo loop_InputOceanGrADS

    tw_org = sst_org

    !do it = 1, nt
    !   i=int(dims(8)/2) ; j=int(dims(9)/2)
    !   write(*,*) "read 2D grads data",dims(8),dims(9),i,j,it
    !   write(*,*) "lon_org    ",lon_org  (i,j)
    !   write(*,*) "lat_org    ",lat_org  (i,j)
    !   write(*,*) "sst_org    ",sst_org  (i,j)
    !   write(*,*) "lst_org  ",lst_org(i,j)
    !   do k=1,dims(7)
    !      write(*,*) "tg_org    ",tg_org   (k,i,j)," k= ",k
    !      write(*,*) "strg_org  ",strg_org (k,i,j)," k= ",k
    !   enddo
    !enddo

    return
  end subroutine ParentOceanInputGrADS

  subroutine read_namelist( &
       grads_item,      &
       grads_fname,     &
       grads_dtype,     &
       grads_swpoint,   &
       grads_dd,        &
       grads_lnum,      &
       grads_lvars,     &
       grads_startrec,  &
       grads_totalrec,  &
       grads_knum,      &
       grads_yrev,      &
       grads_fendian,   &
       grads_missval,   &
       data_available,  &
       item_list,       &
       num_item_list,   &
       basename,        &
       io_fid_grads_nml )
    implicit none
    character(len=H_SHORT), intent(out) :: grads_item    (:)
    character(len=H_LONG),  intent(out) :: grads_fname   (:)
    character(len=H_LONG),  intent(out) :: grads_dtype   (:)
    real(RP),               intent(out) :: grads_swpoint (:)
    real(RP),               intent(out) :: grads_dd      (:)
    integer,                intent(out) :: grads_lnum    (:)
    real(RP),               intent(out) :: grads_lvars   (:,:)
    integer,                intent(out) :: grads_startrec(:)
    integer,                intent(out) :: grads_totalrec(:)
    integer,                intent(out) :: grads_knum    (:)
    character(len=H_SHORT), intent(out) :: grads_yrev    (:)
    character(len=H_SHORT), intent(out) :: grads_fendian (:)
    real(SP),               intent(out) :: grads_missval (:)
    logical,                intent(out) :: data_available(:)
    character(len=*),       intent(in)  :: item_list     (:)
    integer,                intent(in)  :: num_item_list
    character(len=*),       intent(in)  :: basename
    integer,                intent(in)  :: io_fid_grads_nml

    integer :: grads_vars_nmax
    integer :: k, n, ielem, ierr

    namelist /grdvar/ &
         item,      &  ! necessary
         dtype,     &  ! necessary
         fname,     &  ! necessary except for linear data
         swpoint,   &  ! for linear data
         dd,        &  ! for linear data
         lnum,      &  ! for levels data
         lvars,     &  ! for levels data
         startrec,  &  ! for map data
         totalrec,  &  ! for map data
         missval,   &  ! option
         knum,      &  ! option
         yrev,      &  ! option
         fendian       ! option

    ! listup variables
    if ( io_fid_grads_nml > 0 ) then
       rewind( io_fid_grads_nml )
       grads_vars_nmax = 0
       do n = 1, grads_vars_limit
          read(io_fid_grads_nml, nml=grdvar, iostat=ierr)
          if( ierr > 0 )then
             write(*,*) 'xxx [realinput_grads/read_namelist] Not appropriate names in grdvar in ', &
                        trim(basename),'. Check!'
             call PRC_MPIstop
          else if( ierr < 0 )then
             exit
          endif
          grads_vars_nmax = grads_vars_nmax + 1
       enddo
    else
       write(*,*) 'xxx [realinput_grads/read_namelist] namelist file is not open! ', trim(basename)
       call PRC_MPIstop
    endif

    if ( grads_vars_nmax > grads_vars_limit ) then
       write(*,*) 'xxx [realinput_grads/read_namelist] The number of grads vars exceeds grads_vars_limit! ', &
                  grads_vars_nmax, ' > ', grads_vars_limit
       call PRC_MPIstop
    endif

    ! check data availability
    data_available(:) = .false.
    do ielem = 1, num_item_list
       if ( io_fid_grads_nml > 0 ) rewind( io_fid_grads_nml )
       do n = 1, grads_vars_nmax

          ! set default
          item     = ''
          dtype    = ''
          fname    = ''
          swpoint  = large_number_one
          dd       = large_number_one
          lnum     = -99
          lvars    = large_number_one
          startrec = -99
          totalrec = -99
          knum     = -99
          yrev     = 'off'
          fendian  = 'big'
          missval  = large_number_one

          ! read namelist
          if ( io_fid_grads_nml > 0 ) then
             read(io_fid_grads_nml, nml=grdvar, iostat=ierr)
             if( ierr /= 0 ) exit
          endif

          if(item == item_list(ielem))then
             grads_item    (ielem) = item
             grads_fname   (ielem) = fname
             grads_dtype   (ielem) = dtype
             grads_swpoint (ielem) = swpoint
             grads_dd      (ielem) = dd
             grads_lnum    (ielem) = lnum
             do k = 1, lvars_limit
                grads_lvars(k,ielem) = lvars(k)
             enddo
             grads_startrec(ielem) = startrec
             grads_totalrec(ielem) = totalrec
             grads_knum    (ielem) = knum
             grads_yrev    (ielem) = yrev
             grads_fendian (ielem) = fendian
             grads_missval (ielem) = missval
             data_available(ielem) = .true.

             exit
          endif
       enddo ! n
       if( IO_L ) write(IO_FID_LOG,*) 'GrADS data availability ',trim(item_list(ielem)),data_available(ielem)
    enddo ! ielem

  end subroutine read_namelist

  !-----------------------------------------------------------------------------
  subroutine open_grads_file(io_fid,filename,irecl)
    implicit none

    integer,          intent(in) :: io_fid
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: irecl

    integer  :: ierr

    open(io_fid,                   &
         file   = trim(filename),  &
         form   = 'unformatted',   &
         access = 'direct',        &
         recl   = irecl,           &
         status = 'old',           &
         iostat = ierr             )
    if ( ierr /= 0 ) then
       write(*,*) 'xxx grads file does not found! ', trim(filename)
       call PRC_MPIstop
    endif

    return
  end subroutine open_grads_file

  !-----------------------------------------------------------------------------
  subroutine read_grads_file_2d(  &
       io_fid,                    &
       gfile,                     &
       nx,ny,nz,it,               &
       item,                      &
       startrec,                  &
       totalrec,                  &
       yrev,                      &
       gdata                      )
    implicit none

    integer,          intent(in)  :: io_fid
    character(len=*), intent(in)  :: gfile
    integer,          intent(in)  :: nx,ny,nz,it
    character(len=*), intent(in)  :: item
    integer,          intent(in)  :: startrec
    integer,          intent(in)  :: totalrec
    character(len=*), intent(in)  :: yrev
    real(SP),         intent(out) :: gdata(nx,ny)

    real(SP) :: work(nx,ny)

    integer  :: ierr
    integer  :: irec, irecl
    integer  :: i,j,k

    irecl=nx*ny*4
    call open_grads_file(io_fid, gfile, irecl)
    irec = totalrec * (it-1) + startrec
    read(io_fid, rec=irec, iostat=ierr) gdata(:,:)
    if ( ierr /= 0 ) then
       write(*,*) 'xxx grads data is not found! ',trim(item),it
       write(*,*) 'xxx namelist or grads data might be wrong.'
       call PRC_MPIstop
    endif

    if( trim(yrev) == "on" )then
       work(:,:)=gdata(:,:)
       do j=1,ny
       do i=1,nx
          gdata(i,j)=work(i,ny-j+1)
       enddo
       enddo
    endif

    call close_grads_file(io_fid,gfile)

    return
  end subroutine read_grads_file_2d

  !-----------------------------------------------------------------------------
  subroutine read_grads_file_3d(  &
       io_fid,                    &
       gfile,                     &
       nx,ny,nz,it,               &
       item,                      &
       startrec,                  &
       totalrec,                  &
       yrev,                      &
       gdata                      )
    implicit none

    integer,          intent(in)  :: io_fid
    character(len=*), intent(in)  :: gfile
    integer,          intent(in)  :: nx,ny,nz,it
    character(len=*), intent(in)  :: item
    integer,          intent(in)  :: startrec
    integer,          intent(in)  :: totalrec
    character(len=*), intent(in)  :: yrev
    real(SP),         intent(out) :: gdata(nx,ny,nz)

    real(SP) :: work(nx,ny,nz)

    integer  :: ierr
    integer  :: irec,irecl
    integer  :: i,j,k

    irecl=nx*ny*4
    call open_grads_file(io_fid, gfile, irecl)
    do k = 1, nz
       irec = totalrec * (it-1) + startrec + (k-1)
       read(io_fid, rec=irec, iostat=ierr) gdata(:,:,k)
       if ( ierr /= 0 ) then
          write(*,*) 'xxx grads data does not found! ',trim(item),', k=',k,', it=',it,' in ', trim(gfile)
          call PRC_MPIstop
       endif
    enddo

    if( trim(yrev) == "on" )then
       work(:,:,:)=gdata(:,:,:)
       do k=1,nz
       do j=1,ny
       do i=1,nx
          gdata(i,j,k)=work(i,ny-j+1,k)
       enddo
       enddo
       enddo
    endif

    call close_grads_file(io_fid,gfile)

    return
  end subroutine read_grads_file_3d

  !-----------------------------------------------------------------------------
  subroutine close_grads_file(io_fid,filename)
    implicit none

    integer,          intent(in) :: io_fid
    character(len=*), intent(in) :: filename
    integer                      :: ierr

    close(io_fid, iostat=ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'xxx grads file was not closed peacefully! ',trim(filename)
       call PRC_MPIstop
    endif

    return
  end subroutine close_grads_file

end module mod_realinput_grads
