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
     PRC_master,            &
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
  public :: ParentSurfaceInputGrADS

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
  integer,  parameter   :: num_item_list = 26
  integer,  parameter   :: num_item_list_atm = 16
  logical               :: data_available(num_item_list)
  character(H_SHORT)    :: item_list     (num_item_list)
  data item_list /'lon','lat','plev','U','V','T','HGT','QV','RH','MSLP','PSFC','U10','V10','T2','Q2','RH2', &
                  'lsmask','lon_sfc','lat_sfc','lon_sst','lat_sst','llev','STEMP','SMOIS','SKINT','SST'/

  integer,  parameter   :: Ig_lon    = 1
  integer,  parameter   :: Ig_lat    = 2
  integer,  parameter   :: Ig_p      = 3  ! Pressure (Pa)
  integer,  parameter   :: Ig_u      = 4
  integer,  parameter   :: Ig_v      = 5
  integer,  parameter   :: Ig_t      = 6
  integer,  parameter   :: Ig_hgt    = 7  ! Geopotential height (m)
  integer,  parameter   :: Ig_qv     = 8
  integer,  parameter   :: Ig_rh     = 9  ! Percentage (%)
  integer,  parameter   :: Ig_slp    = 10 ! Sea level pressure (Pa)
  integer,  parameter   :: Ig_ps     = 11 ! Surface pressure (Pa)
  integer,  parameter   :: Ig_u10    = 12
  integer,  parameter   :: Ig_v10    = 13
  integer,  parameter   :: Ig_t2     = 14
  integer,  parameter   :: Ig_q2     = 15
  integer,  parameter   :: Ig_rh2    = 16 ! Percentage (%)

  integer,  parameter   :: Ig_lsmask  = 17
  integer,  parameter   :: Ig_lon_sfc = 18
  integer,  parameter   :: Ig_lat_sfc = 19
  integer,  parameter   :: Ig_lon_sst = 20
  integer,  parameter   :: Ig_lat_sst = 21
  integer,  parameter   :: Ig_lz      = 22  ! Level(depth) of stemp & smois (m)
  integer,  parameter   :: Ig_stemp   = 23
  integer,  parameter   :: Ig_smois   = 24
  integer,  parameter   :: Ig_skint   = 25
  integer,  parameter   :: Ig_sst     = 26

 
  integer,  parameter   :: lvars_limit = 1000 ! limit of values for levels data
  real(RP), parameter   :: large_number_one = 9.999E+15_RP


  character(len=H_SHORT) :: upper_qv_type = "ZERO" !< how qv is given at higher level than outer model
                                                     !< "ZERO": 0
                                                     !< "COPY": copy values from the highest level of outer model

  
  character(H_SHORT)    :: grads_item    (num_item_list)
  character(H_SHORT)    :: grads_kytpe   (num_item_list)
  character(H_LONG)     :: grads_dtype   (num_item_list)
  character(H_LONG)     :: grads_fname   (num_item_list)
  character(H_SHORT)    :: grads_fendian (num_item_list)
  real(RP)              :: grads_swpoint (num_item_list)
  real(RP)              :: grads_dd      (num_item_list)
  integer               :: grads_lnum    (num_item_list)
  real(RP)              :: grads_lvars   (num_item_list,lvars_limit)
  integer               :: grads_startrec(num_item_list)
  integer               :: grads_totalrec(num_item_list)
  integer               :: grads_knum    (num_item_list)
  real(SP)              :: grads_missval (num_item_list)

  real(SP), allocatable :: gdata2D(:,:)
  real(SP), allocatable :: gdata3D(:,:,:)
  real(SP), allocatable :: gland2D(:,:)
  real(SP), allocatable :: gland3D(:,:,:)
  real(SP), allocatable :: gsst2D (:,:)

  integer :: io_fid_grads_nml  = -1
  integer :: io_fid_grads_data = -1
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetupGrADS( &
      dims,   &
      timelen )
    implicit none

    integer,          intent(out) :: dims(11)
    integer,          intent(out) :: timelen

    character(len=H_LONG) :: grads_boundary_namelist = "namelist.grads_boundary"

    NAMELIST / PARAM_MKINIT_REAL_GrADS / &
        grads_boundary_namelist, &
        upper_qv_type


    integer :: outer_nx     = 0
    integer :: outer_ny     = 0
    integer :: outer_nz     = 0 ! number of atmos layers
    integer :: outer_nl     = 0 ! number of land layers
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

    integer                       :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: GrADS format'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_GrADS,iostat=ierr)

    if( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_GrADS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_REAL_GrADS)


    if ( len_trim(GrADS_BOUNDARY_namelist) == 0 ) then
       write(*,*) 'xxx "GrADS_BOUNDARY_namelist" is not specified in "PARAM_MKINIT_REAL"!',trim(GrADS_BOUNDARY_namelist)
       call PRC_MPIstop
    endif

    !--- read namelist
    io_fid_grads_nml = IO_get_available_fid()
    open( io_fid_grads_nml,                    &
         file   = trim(grads_boundary_namelist), &
         form   = 'formatted',                   &
         status = 'old',                         &
         action = 'read',                        &
         iostat = ierr                           )
    if ( ierr /= 0 ) then
       write(*,*) 'xxx Input data file does not found! ', trim(grads_boundary_namelist)
       stop
    endif

    read(io_fid_grads_nml,nml=nml_grads_grid,iostat=ierr)
    if( ierr /= 0 ) then !--- missing or fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_REAL_GrADS_grid. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=nml_grads_grid)

    timelen = 0        ! will be replaced later

    ! full level
    dims(1) = outer_nz ! bottom_top
    dims(2) = outer_nx ! west_east
    dims(3) = outer_ny ! south_north
    ! half level
    dims(4) = outer_nz ! bottom_top_stag
    dims(5) = outer_nx ! west_east for 2dim data
    dims(6) = outer_ny ! south_north for 2dim data
    ! land
    dims(7) = outer_nl ! soil_layers_stag
    if(outer_nx_sfc > 0)then
       dims(8) = outer_nx_sfc
    else
       dims(8) = outer_nx
    endif
    if(outer_ny_sfc > 0)then
       dims(9) = outer_ny_sfc
    else
       dims(9) = outer_ny
    endif
    ! sst
    if(outer_nx_sst > 0)then
       dims(10) = outer_nx_sst
    else
       dims(10) = outer_nx
    endif
    if(outer_ny_sst > 0)then
       dims(11) = outer_ny_sst
    else
       dims(11) = outer_ny
    endif

    allocate( gdata2D( dims(2), dims(3)          ) )
    allocate( gdata3D( dims(2), dims(3), dims(1) ) )

    allocate( gland2D( dims(8), dims(9)          ) )
    allocate( gland3D( dims(8), dims(9), dims(7) ) )

    allocate( gsst2D ( dims(10), dims(11)        ) )

    return
  end subroutine ParentAtomSetupGrADS

  !-----------------------------------------------------------------------------
  subroutine ParentAtomOpenGrADS
    implicit none

    ! namelist.grads_boundary
    integer, parameter   :: grads_vars_limit = 1000 !> limit of number of values
    integer              :: grads_vars_nmax = 0     !> number of variables in grads file

    character(H_SHORT)   :: item                                      ! up to 16 characters
    integer              :: knum                                      ! optional: vertical level
    character(H_SHORT)   :: dtype                                     ! 'linear','levels','map'
    character(H_LONG)    :: fname                                     ! head of file name
    real(RP)             :: swpoint                                   ! start point (south-west point) for linear
    real(RP)             :: dd                                        ! dlon,dlat for linear
    integer              :: lnum                                      ! number of data
    real(RP)             :: lvars(lvars_limit) = large_number_one     ! values for levels
    integer              :: startrec=1                                ! record position
    integer              :: totalrec=1                                ! total record number per one time
    real(SP)             :: missval                                   ! missing value
    character(H_SHORT)   :: fendian                                   ! option

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
        fendian       ! option


    integer :: ielem
    integer :: ierr
    integer :: k, n

    !---------------------------------------------------------------------------

    ! listup variables
    if ( io_fid_grads_nml > 0 ) then
       rewind( io_fid_grads_nml )
       grads_vars_nmax = 0
       do n = 1, grads_vars_limit
          read(io_fid_grads_nml, nml=grdvar, iostat=ierr)
          if( ierr > 0 )then
             write(IO_FID_LOG,*) '*** cannot read namelist successfully! '
             call PRC_MPIstop
          else if( ierr < 0 )then
             exit
          endif
          grads_vars_nmax = grads_vars_nmax + 1
       enddo
    else
       write(IO_FID_LOG,*) '*** namelist file is not open! '
       call PRC_MPIstop
    endif

    if ( grads_vars_nmax > grads_vars_limit ) then
       write(IO_FID_LOG,*) '*** number of grads vars is exceed! ',grads_vars_nmax,' >', grads_vars_limit
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
          startrec = 0
          totalrec = 0
          knum     = -99
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
                grads_lvars(ielem,k) = lvars(k)
             enddo
             grads_startrec(ielem) = startrec
             grads_totalrec(ielem) = totalrec
             grads_knum    (ielem) = knum
             grads_missval (ielem) = missval
             grads_fendian (ielem) = fendian
             data_available(ielem) = .true.

             exit
          endif
       enddo ! n
       if( IO_L ) write(IO_FID_LOG,*) 'GrADS data availability ',trim(item_list(ielem)),data_available(ielem)
    enddo ! ielem


    do ielem = 1, num_item_list
       item  = item_list(ielem)
       !--- check data
       select case(trim(item))
       case('QV')
          if (.not. data_available(Ig_qv)) then
             if (.not.data_available(Ig_rh)) then
                write(IO_FID_LOG,*) '*** cannot found QV and RH in grads namelist! '
                call PRC_MPIstop
             else ! will read RH
                cycle
             endif
          endif
       case('RH')
          if ( data_available(Ig_rh)) then
             if ((.not. data_available(Ig_t)).or.(.not. data_available(Ig_p))) then
                write(IO_FID_LOG,*) '*** Temperature and pressure are required to convert from RH to QV ! '
                call PRC_MPIstop
             endif
          else
             if (.not.data_available(Ig_qv)) then
                write(IO_FID_LOG,*) '*** cannot found QV and RH in grads namelist! '
                call PRC_MPIstop
             else ! will read QV
                cycle
             endif
          endif
       case('Q2')
          if (.not. data_available(Ig_q2)) then
             if (.not.data_available(Ig_rh2)) then
                write(IO_FID_LOG,*) '*** cannot found Q2 and RH2 in grads namelist! '
                call PRC_MPIstop
             else
                cycle ! will read RH2
             endif
          endif
       case('RH2')
          if ( data_available(Ig_rh2)) then
             if ((.not. data_available(Ig_t2)).or.(.not. data_available(Ig_ps))) then
                write(IO_FID_LOG,*) '*** T2 and PSFC are required to convert from RH2 to Q2 ! '
                stop
             endif
          else
             if (.not.data_available(Ig_q2)) then
                write(IO_FID_LOG,*) '*** cannot found Q2 and RH2 in grads namelist! ',trim(item)
                call PRC_MPIstop
             else
                cycle ! will read Q2
             endif
          endif
       case('lsmask')
          if ( .not. data_available(Ig_lsmask) ) then
             write(IO_FID_LOG,*) '*** not use ',trim(item),' data!'
             cycle
          endif
       case('lon_sfc', 'lat_sfc', 'lon_sst', 'lat_sst')
          cycle
       case('SST')
          if (.not. data_available(Ig_sst)) then
             write(IO_FID_LOG,*) 'warning: cannot found SST. SKINT is used for SST! '
             cycle
          endif
       case('SMOIS')
          if (.not. data_available(Ig_smois)) then
             cycle
          end if
       case default
          if ( .not. data_available(ielem) ) then
             write(*,*) 'xxx cannot found data in grads namelist! ',trim(item), ", ", trim(item_list(ielem))
             call PRC_MPIstop
          endif
       end select

    end do

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
    character(LEN=*), intent(in)  :: basename_num
    integer,          intent(in)  :: dims(11)
    integer,          intent(in)  :: nt

    real(RP) :: rhprs_org(dims(1)+2,dims(2),dims(3))
    real(RP) :: pott(dims(2),dims(3))

    real(RP) :: RovCP
    real(RP) :: CPovR

    ! for namelist
    character(H_SHORT)   :: item                                      ! up to 16 characters
    integer              :: knum                                      ! optional: vertical level
    character(H_SHORT)   :: dtype                                     ! 'linear','levels','map'
    character(H_LONG)    :: fname                                     ! head of file name
    real(RP)             :: swpoint                                   ! start point (south-west point) for linear
    real(RP)             :: dd                                        ! dlon,dlat for linear
    integer              :: lnum                                      ! number of data
    real(RP)             :: lvars(lvars_limit) = large_number_one     ! values for levels
    integer              :: startrec=1                                ! record position
    integer              :: totalrec=1                                ! total record number per one time
    character(H_SHORT)   :: fendian='big_endian'                      ! option

    ! data
    character(len=H_LONG) :: gfile

    integer  :: QA_outer = 1
    real(RP) :: p_sat, qm, rhsfc

    integer  :: i, j, k, ielem

    !---------------------------------------------------------------------------



    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputGrADS reanalysis]'

    qtrc_org = 0.0_RP

    !--- read grads data
    loop_InputAtomGrADS : do ielem = 1, num_item_list_atm

       item     = grads_item    (ielem)
       dtype    = grads_dtype   (ielem)
       fname    = grads_fname   (ielem)
       lnum     = grads_lnum    (ielem)

       if ( dims(1) < grads_knum(ielem) ) then
          write(IO_FID_LOG,*) '*** please check plev data! knum must be less than or equal to outer_nz',knum,dims(1)
          call PRC_MPIstop
       else if ( grads_knum(ielem) > 0 )then
          knum = grads_knum(ielem)  ! not missing
       else
          knum = dims(1)
       endif

       select case (trim(dtype))
       case("linear")
          swpoint = grads_swpoint (ielem)
          dd      = grads_dd      (ielem)
          if( (abs(swpoint-large_number_one)<EPS).or.(abs(dd-large_number_one)<EPS) )then
             write(IO_FID_LOG,*) '*** "swpoint" is required in grads namelist! ',swpoint
             write(IO_FID_LOG,*) '*** "dd"      is required in grads namelist! ',dd
             call PRC_MPIstop
          endif
       case("levels")
          if ( lnum < 0 )then
             write(IO_FID_LOG,*) '*** "lnum" in grads namelist is required for levels data! '
             call PRC_MPIstop
          endif
          do k=1, lnum
             lvars(k)=grads_lvars(ielem,k)
          enddo
          if(abs(lvars(1)-large_number_one)<EPS)then
             write(IO_FID_LOG,*) '*** "lvars" must be specified in grads namelist for levels data! '
             call PRC_MPIstop
          endif
       case("map")
          startrec = grads_startrec(ielem)
          totalrec = grads_totalrec(ielem)
          fendian  = grads_fendian (ielem)
          if( (startrec==0).or.(totalrec==0) )then
             write(IO_FID_LOG,*) '*** "startrec" is required in grads namelist! ',startrec
             write(IO_FID_LOG,*) '*** "totalrec" is required in grads namelist! ',totalrec
             call PRC_MPIstop
          endif
          ! get file_id
          if(io_fid_grads_data < 0)then
             io_fid_grads_data = IO_get_available_fid()
          endif
          gfile=trim(fname)//trim(basename_num)//'.grd'
       end select

       ! read data
       select case (trim(item))
       case("lon")
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
                lon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,1,item,startrec,totalrec,gdata2D)
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
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,1,item,startrec,totalrec,gdata2D)
             lat_org(:,:) = real(gdata2D(:,:), kind=RP) * D2R
          endif
       case("plev")
          if(dims(1)/=knum)then
             write(IO_FID_LOG,*) '*** please check plev data! ',dims(1),knum
             stop
          endif
          if ( trim(dtype) == "levels" ) then
             if(dims(1)/=lnum)then
                write(IO_FID_LOG,*) '*** lnum must be same as the outer_nz! ',dims(1),lnum
                stop
             endif
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)
                pres_org(k+2,i,j) = real(lvars(k), kind=RP)
             enddo
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),dims(1),nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)
                pres_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
             enddo
             enddo
             enddo
          endif
       case('U')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                velx_org(1:2,i,j) = 0.0_RP
                do k = 1, knum
                   velx_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
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
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                vely_org(1:2,i,j) = 0.0_RP
                do k = 1, knum
                   vely_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
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
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, knum
                   temp_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
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
             write(IO_FID_LOG,*) '*** The number of levels for HGT must be same as plevs! ',dims(1),knum
             stop
          endif
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),dims(1),nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                do k = 1, dims(1)
                   cz_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP)
                enddo
                cz_org(1,i,j) = 0.0_RP
             enddo
             enddo
          endif
       case('QV')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, knum
                qtrc_org(k+2,i,j,QA_outer) = real(gdata3D(i,j,k), kind=RP)
             enddo
             enddo
             enddo
             if( dims(1)>knum ) then
                select case ( upper_qv_type )
                case ("COPY")
                   do j = 1, dims(3)
                   do i = 1, dims(2)
                   do k = knum+1, dims(1)
                      qtrc_org(k+2,i,j,QA_outer) = qtrc_org(knum+2,i,j,QA_outer)
                   enddo
                   enddo
                   enddo
                case ("ZERO")
                   ! do nothing
                case default
                   write(*,*) 'xxx upper_qv_type is invalid! ', upper_qv_type
                   call PRC_MPIstop
                end select
             endif
          endif
       case('RH')
          if (data_available(Ig_qv)) cycle  ! use QV
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(2),dims(3),knum,nt,item,startrec,totalrec,gdata3D)
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, knum
                rhprs_org(k+2,i,j) = real(gdata3D(i,j,k), kind=RP) / 100.0_RP  ! relative humidity
                call psat( p_sat, temp_org(k+2,i,j) )                          ! satulation pressure
                qm = EPSvap * rhprs_org(k+2,i,j) * p_sat &
                   / ( pres_org(k+2,i,j) - rhprs_org(k+2,i,j) * p_sat )          ! mixing ratio
                qtrc_org(k+2,i,j,QA_outer) = qm / ( 1.0_RP + qm )              ! specific humidity
             enddo
             enddo
             enddo
             if( dims(3)>knum ) then
                select case ( upper_qv_type )
                case ("COPY")
                   do j = 1, dims(3)
                   do i = 1, dims(2)
                   do k = knum+1, dims(1)
                      rhprs_org(k+2,i,j) = rhprs_org(knum+2,i,j)              ! relative humidity
                      call psat( p_sat, temp_org(k+2,i,j) )                 ! satulated specific humidity
                      qm = EPSvap * rhprs_org(k+2,i,j) * p_sat &
                         / ( pres_org(k+2,i,j) - rhprs_org(k+2,i,j) * p_sat ) ! mixing ratio
                      qtrc_org(k+2,i,j,QA_outer) = qm / ( 1.0_RP + qm )     ! specific humidity
                      qtrc_org(k+2,i,j,QA_outer) = min(qtrc_org(k+2,i,j,QA_outer),qtrc_org(k+1,i,j,QA_outer))
                   enddo
                   enddo
                   enddo
                case ("ZERO")
                   ! do nothing
                case default
                   write(*,*) 'xxx upper_qv_type is invalid! ', upper_qv_type
                   call PRC_MPIstop
                end select
             endif
          endif
       case('MSLP')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                pres_org(1,i,j) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('PSFC')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                pres_org(2,i,j) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('U10')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                velx_org(2,i,j) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('V10')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                vely_org(2,i,j) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('T2')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                temp_org(2,i,j) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('Q2')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                qtrc_org(2,i,j,QA_outer) = real(gdata2D(i,j), kind=RP)
             enddo
             enddo
          endif
       case('RH2')
          if (data_available(Ig_q2)) cycle  ! use QV
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(2),dims(3),1,nt,item,startrec,totalrec,gdata2D)
             do j = 1, dims(3)
             do i = 1, dims(2)
                rhsfc = real(gdata2D(i,j), kind=RP) / 100.0_RP   ! relative humidity
                call psat( p_sat, temp_org(2,i,j) )                         ! satulation pressure
                qm = EPSvap * rhsfc * p_sat &
                   / ( pres_org(2,i,j) - rhsfc * p_sat )           ! mixing ratio
                qtrc_org(2,i,j,QA_outer) = qm / ( 1.0_RP + qm )             ! specific humidity
             enddo
             enddo
          endif
       end select
    enddo loop_InputAtomGrADS


    RovCP = Rdry / CPdry
    CPovR = CPdry / Rdry

    if ( .not. data_available(Ig_t2) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          temp_org(2,i,j) = temp_org(3,i,j)
       end do
       end do
    end if
    if ( data_available(Ig_t2) .and. data_available(Ig_ps) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          pott(i,j) = temp_org(2,i,j) * (P00/pres_org(2,i,j))**RovCP
       end do
       end do
    else
       do j = 1, dims(3)
       do i = 1, dims(2)
          pott(i,j) = temp_org(3,i,j) * (P00/pres_org(3,i,j))**RovCP
       end do
       end do
    end if
    if ( .not. data_available(Ig_ps) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          pres_org(2,i,j) = P00 * (temp_org(2,i,j)/pott(i,j))**CPovR
       end do
       end do
    end if
    if ( data_available(Ig_slp) ) then
       do j = 1, dims(3)
       do i = 1, dims(2)
          temp_org(1,i,j) = pott(i,j) * (pres_org(1,i,j)/P00)**RovCP
       end do
       end do
    else
       do j = 1, dims(3)
       do i = 1, dims(2)
          temp_org(1,i,j) = temp_org(3,i,j) + LAPS * cz_org(3,i,j)
          pres_org(1,i,j) = P00 * (temp_org(1,i,j)/pott(i,j))**CPovR
       end do
       end do
    end if

    ! guess surface height (elevation)
    do j = 1, dims(3)
    do i = 1, dims(2)
       cz_org(2,i,j) = max( 0.0_RP, &
                            cz_org(3,i,j) &
                            * ( log(pres_org(2,i,j)/pres_org(1,i,j)) ) &
                            / ( log(pres_org(3,i,j)/pres_org(1,i,j)) ) )
    end do
    end do


    velz_org = 0.0_RP


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
  subroutine ParentSurfaceInputGrADS( &
      tg_org,             & ! (out)
      strg_org,           & ! (out)
      tw_org,             & ! (out)
      lst_org,            & ! (out)
      sst_org,            & ! (out)
      lsmask_org,         & ! (out)
      lz_org,             & ! (out)
      llon_org,           & ! (out)
      llat_org,           & ! (out)
      olon_org,           & ! (out)
      olat_org,           & ! (out)
      basename_num,       & ! (in)
      dims,               & ! (in)
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
    real(RP), intent(out) :: tw_org    (:,:)
    real(RP), intent(out) :: lst_org   (:,:)
    real(RP), intent(out) :: sst_org   (:,:)
    real(RP), intent(out) :: lsmask_org(:,:)
    real(RP), intent(out) :: lz_org    (:)
    real(RP), intent(out) :: llon_org  (:,:)
    real(RP), intent(out) :: llat_org  (:,:)
    real(RP), intent(out) :: olon_org  (:,:)
    real(RP), intent(out) :: olat_org  (:,:)

    character(LEN=*), intent(in) :: basename_num
    integer,          intent(in) :: dims(11)
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    integer,          intent(in) :: nt
    ! ----------------------------------------------------------------

    ! namelist.grads_boundary
    integer              :: grads_vars_nmax = 0     !> number of variables in grads file

    character(H_SHORT)   :: item                                      ! up to 16 characters
    integer              :: knum                                      ! optional: vertical level
    character(H_SHORT)   :: dtype                                     ! 'linear','levels','map'
    character(H_LONG)    :: fname                                     ! head of file name
    real(RP)             :: swpoint                                   ! start point (south-west point) for linear
    real(RP)             :: dd                                        ! dlon,dlat for linear
    integer, parameter   :: lvars_limit = 1000                        ! limit of values for levels data
    integer              :: lnum                                      ! number of data
    real(RP)             :: lvars(lvars_limit) = large_number_one     ! values for levels
    integer              :: startrec=1                                ! record position
    integer              :: totalrec=1                                ! total record number per one time
    real(SP)             :: missval                                   ! option
    character(H_SHORT)   :: fendian                                   ! option

    !> grads data
    character(len=H_LONG) :: gfile

    real(RP) :: qvsat, qm

    integer :: i, j, k, ielem, n
    integer :: ierr

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputGrADS-Surface]'


    loop_InputSurfaceGrADS : do ielem = 1, num_item_list

       item  = item_list(ielem)

       dtype    = grads_dtype   (ielem)
       fname    = grads_fname   (ielem)
       lnum     = grads_lnum    (ielem)
       missval  = grads_missval (ielem)

       if ( grads_knum(ielem) > 0 )then
          knum = grads_knum(ielem)
       else
          knum = dims(7)
       endif

       select case (trim(dtype))
       case("linear")
          swpoint = grads_swpoint (ielem)
          dd      = grads_dd      (ielem)
          if( (abs(swpoint-large_number_one)<EPS).or.(abs(dd-large_number_one)<EPS) )then
             write(IO_FID_LOG,*) '*** "swpoint" is required in grads namelist! ',swpoint
             write(IO_FID_LOG,*) '*** "dd"      is required in grads namelist! ',dd
             call PRC_MPIstop
          endif
       case("levels")
          if ( lnum < 0 )then
             write(IO_FID_LOG,*) '*** "lnum" in grads namelist is required for levels data! '
             call PRC_MPIstop
          endif
          do k=1, lnum
             lvars(k)=grads_lvars(ielem,k)
          enddo
          if(abs(lvars(1)-large_number_one)<EPS)then
             write(IO_FID_LOG,*) '*** "lvars" is required in grads namelist! ',(lvars(k),k=1,lnum)
             call PRC_MPIstop
          endif
       case("map")
          startrec = grads_startrec(ielem)
          totalrec = grads_totalrec(ielem)
          fendian  = grads_fendian (ielem)
          if( (startrec==0).or.(totalrec==0) )then
             write(IO_FID_LOG,*) '*** "startrec" is required in grads namelist! ',startrec
             write(IO_FID_LOG,*) '*** "totalrec" is required in grads namelist! ',totalrec
             call PRC_MPIstop
          endif
          ! get file_io
          if(io_fid_grads_data < 0)then
             io_fid_grads_data = IO_get_available_fid()
          endif
          gfile=trim(fname)//trim(basename_num)//'.grd'
       end select

       ! read data
       select case (trim(item))
       case("lsmask")
          if ( data_available(Ig_lsmask) ) then
             if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,1,item,startrec,totalrec,gland2D)
                lsmask_org(:,:) = real(gland2D(:,:), kind=RP)
             endif
          else
             lsmask_org = UNDEF
          end if
       case("lon")
          if ( .not. data_available(Ig_lon_sfc) ) then
             if ( dims(2).ne.dims(8) .or. dims(3).ne.dims(9) ) then
                write(*,*) 'xxx namelist of "lon_sfc" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sfc! ',dims(2), dims(8)
                write(*,*) '                          : outer_ny and outer_ny_sfc! ',dims(3), dims(9)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, dims(9)
                do i = 1, dims(8)
                   llon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,1,item,startrec,totalrec,gland2D)
                llon_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
             endif
          end if
          if ( .not. data_available(Ig_lon_sfc) .and. .not. data_available(Ig_lon_sst) ) then
             if ( dims(2).ne.dims(10) .or. dims(3).ne.dims(11) ) then
                write(*,*) 'xxx namelist of "lon_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sst! ',dims(2), dims(10)
                write(*,*) '                          : outer_ny and outer_ny_sst! ',dims(3), dims(11)
                call PRC_MPIstop
             end if
             olon_org = llon_org
          end if
       case("lon_sfc")
          if ( .not. data_available(Ig_lon_sfc) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(9)
             do i = 1, dims(8)
                llon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,1,item,startrec,totalrec,gland2D)
             llon_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
          endif
          if ( .not. data_available(Ig_lon_sst) ) then
             if ( dims(8)==dims(10) .and. dims(9)==dims(11) ) then
                olon_org = llon_org
             else
                write(*,*) 'xxx namelist of "lon_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx/outer_nx_sfc and outer_nx_sst! ',dims(8), dims(10)
                write(*,*) '                          : outer_ny/outer_ny_sfc and outer_ny_sst! ',dims(9), dims(11)
                call PRC_MPIstop
             end if
          end if
       case("lon_sst")
          if ( .not. data_available(Ig_lon_sst) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(11)
             do i = 1, dims(10)
                llon_org(i,j) = real(swpoint+real(i-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(10),dims(11),1,1,item,startrec,totalrec,gsst2D)
             llon_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
          endif
       case("lat")
          if ( .not. data_available(Ig_lat_sfc) ) then
             if ( dims(2).ne.dims(8) .or. dims(3).ne.dims(9) ) then
                write(*,*) 'xxx namelist of "lat_sfc" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sfc! ',dims(2), dims(8)
                write(*,*) '                          : outer_ny and outer_ny_sfc! ',dims(3), dims(9)
                call PRC_MPIstop
             end if
             if ( trim(dtype) == "linear" ) then
                do j = 1, dims(9)
                do i = 1, dims(8)
                   llat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
                enddo
                enddo
             else if ( trim(dtype) == "map" ) then
                call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,1,item,startrec,totalrec,gland2D)
                llat_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
             endif
             if ( .not. data_available(Ig_lat_sst) ) then
                if ( dims(2).ne.dims(10) .or. dims(3).ne.dims(11) ) then
                   write(*,*) 'xxx namelist of "lat_sst" is not found in grads namelist!'
                   write(*,*) 'xxx dimension is different: outer_nx and outer_nx_sfc! ',dims(2), dims(10)
                   write(*,*) '                          : outer_ny and outer_ny_sfc! ',dims(3), dims(11)
                   call PRC_MPIstop
                end if
                olat_org = llat_org
             end if
          end if
       case("lat_sfc")
          if ( .not. data_available(Ig_lat_sfc) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(9)
             do i = 1, dims(8)
                llat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,1,item,startrec,totalrec,gland2D)
             llat_org(:,:) = real(gland2D(:,:), kind=RP) * D2R
          endif
          if ( .not. data_available(Ig_lat_sst) ) then
             if ( dims(8)==dims(10) .and. dims(9)==dims(11) ) then
                olat_org = llat_org
             else
                write(*,*) 'xxx namelist of "lat_sst" is not found in grads namelist!'
                write(*,*) 'xxx dimension is different: outer_nx/outer_nx_sfc and outer_nx_sst! ',dims(8), dims(10)
                write(*,*) '                          : outer_ny/outer_ny_sfc and outer_ny_sst! ',dims(9), dims(11)
                call PRC_MPIstop
             end if
          end if
       case("lat_sst")
          if ( .not. data_available(Ig_lat_sst) ) cycle
          if ( trim(dtype) == "linear" ) then
             do j = 1, dims(11)
             do i = 1, dims(10)
                olat_org(i,j) = real(swpoint+real(j-1)*dd, kind=RP) * D2R
             enddo
             enddo
          else if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(10),dims(11),1,1,item,startrec,totalrec,gsst2D)
             olat_org(:,:) = real(gsst2D(:,:), kind=RP) * D2R
          endif
       case("llev")
          if(dims(7)/=knum)then
             write(IO_FID_LOG,*) '*** please check llev data! ',dims(7),knum
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "levels" ) then
             if(dims(7)/=lnum)then
                write(IO_FID_LOG,*) '*** lnum must be same as the outer_nz! ',dims(7),lnum
                call PRC_MPIstop
             endif
             do k = 1, dims(7)
                lz_org(k) = real(lvars(k), kind=RP)
             enddo
!          else if ( trim(dtype) == "map" ) then
!             call read_grads_file_3d(io_fid_grads_data,gfile,dims(8),dims(9),dims(7),nt,item,startrec,totalrec,gland)
!             do j = 1, dims(9)
!             do i = 1, dims(8)
!             do k = 1, dims(7)
!                lz_org(k,i,j) = real(gland(i,j,k), kind=RP)
!             enddo
!             enddo
!             enddo
          endif
       case('STEMP')
          if(dims(7)/=knum)then
             write(IO_FID_LOG,*) '*** The number of levels for STEMP must be same as llevs! ',dims(7),knum
             call PRC_MPIstop
          endif
          if ( trim(dtype) == "map" ) then
             call read_grads_file_3d(io_fid_grads_data,gfile,dims(8),dims(9),dims(7),nt,item,startrec,totalrec,gland3D)
             do j = 1, dims(9)
             do i = 1, dims(8)
             do k = 1, dims(7)
                tg_org(k,i,j) = real(gland3D(i,j,k), kind=RP)
                if ( abs(tg_org(k,i,j)-missval)<EPS ) tg_org(k,i,j) = UNDEF
             enddo
             enddo
             enddo
          endif
       case('SMOIS')
          if ( use_file_landwater ) then
             if ( .not. data_available(Ig_smois) ) then
                write(IO_FID_LOG,*) '*** cannot found SMOIS in grads namelist! ',trim(item),trim(item_list(ielem))
                call PRC_MPIstop
             end if
             if(dims(7)/=knum)then
                write(IO_FID_LOG,*) '*** The number of levels for SMOIS must be same as llevs! ',dims(7),knum
                call PRC_MPIstop
             endif
             if ( trim(dtype) == "map" ) then
                call read_grads_file_3d(io_fid_grads_data,gfile,dims(8),dims(9),dims(7),nt,item,startrec,totalrec,gland3D)
                do j = 1, dims(9)
                do i = 1, dims(8)
                do k = 1, dims(7)
                   strg_org(k,i,j) = real(gland3D(i,j,k), kind=RP)
                   if ( abs(strg_org(k,i,j)-missval)<EPS ) strg_org(k,i,j) = UNDEF
                enddo
                enddo
                enddo
             endif
          endif
       case('SKINT')
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(8),dims(9),1,nt,item,startrec,totalrec,gland2D)
             do j = 1, dims(9)
             do i = 1, dims(8)
                lst_org(i,j) = real(gland2D(i,j), kind=RP)
                if ( abs(lst_org(i,j)-missval)<EPS ) lst_org(i,j) = UNDEF
             enddo
             enddo
          endif
       case('SST')
          if ( .not. data_available(Ig_sst) ) cycle
          if ( trim(dtype) == "map" ) then
             call read_grads_file_2d(io_fid_grads_data,gfile,dims(10),dims(11),1,nt,item,startrec,totalrec,gsst2D)
             do j = 1, dims(11)
             do i = 1, dims(10)
                sst_org(i,j) = real(gsst2D(i,j), kind=RP)
                if ( abs(sst_org(i,j)-missval)<EPS ) sst_org(i,j) = UNDEF
             enddo
             enddo
          end if
       end select
    enddo loop_InputSurfaceGrADS

    !--SST: use skin temp
    if (.not. data_available(Ig_sst)) then
       if ( dims(8)==dims(10) .and. dims(9)==dims(11) ) then
          ! use skin temp
          sst_org(:,:) = lst_org(:,:)
       else
          write(*,*) 'xxx dimsntion is different: outer_nx/outer_nx_sfc and outer_nx_sst! ', dims(8), dims(10)
          write(*,*) '                          : outer_ny/outer_ny_sfc and outer_ny_sst! ', dims(9), dims(11)
          call PRC_MPIstop
       end if
    endif

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
  end subroutine ParentSurfaceInputGrADS

  !-----------------------------------------------------------------------------
  subroutine open_grads_file(io_fid,filename,irecl)
    implicit none
    integer,intent(in)      :: io_fid
    character(*),intent(in) :: filename
    integer,intent(in)      :: irecl
    integer  :: ierr

    open(io_fid,                   &
         file   = trim(filename),  &
         form   = 'unformatted',   &
         access = 'direct',        &
         recl   = irecl,           &
         status = 'old',           &
         iostat = ierr             )
    if ( ierr /= 0 ) then
       if (IO_L) write(IO_FID_LOG,*) 'xxx grads file does not found! ', trim(filename)
       stop
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
       gdata                      )
    implicit none
    integer,     intent(in)  :: io_fid
    character(*),intent(in)  :: gfile
    integer,     intent(in)  :: nx,ny,nz,it
    character(*),intent(in)  :: item
    integer,     intent(in)  :: startrec
    integer,     intent(in)  :: totalrec
    real(SP),    intent(out) :: gdata(nx,ny)

    integer  :: ierr
    integer  :: irec, irecl
    integer  :: i,j,k

    irecl=nx*ny*4
    call open_grads_file(io_fid, gfile, irecl)
    irec = totalrec * (it-1) + startrec
    read(io_fid, rec=irec, iostat=ierr) gdata(:,:)
    if ( ierr /= 0 ) then
       if (IO_L) write(IO_FID_LOG,*) 'xxx grads data does not found! ',trim(item),it
       stop
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
       gdata                      )
    implicit none
    integer,intent(in)      :: io_fid
    character(*),intent(in) :: gfile
    integer,     intent(in) :: nx,ny,nz,it
    character(*),intent(in) :: item
    integer,intent(in)      :: startrec
    integer,intent(in)      :: totalrec
    real(SP),intent(out)    :: gdata(nx,ny,nz)

    integer  :: ierr
    integer  :: irec,irecl
    integer  :: i,j,k

    irecl=nx*ny*4
    call open_grads_file(io_fid, gfile, irecl)
    do k = 1, nz
       irec = totalrec * (it-1) + startrec + (k-1)
       read(io_fid, rec=irec, iostat=ierr) gdata(:,:,k)
       if ( ierr /= 0 ) then
          if (IO_L) write(IO_FID_LOG,*) 'xxx grads data does not found! ',trim(item),k,it
          stop
       endif
    enddo
    call close_grads_file(io_fid,gfile)

    return
  end subroutine read_grads_file_3d

  !-----------------------------------------------------------------------------
  subroutine close_grads_file(io_fid,filename)
    implicit none
    integer,intent(in)      :: io_fid
    character(*),intent(in) :: filename
    integer                 :: ierr

    close(io_fid, iostat=ierr)
    if ( ierr /= 0 ) then
       if (IO_L) write(IO_FID_LOG,*) 'xxx grads file was not closed peacefully! ',trim(filename)
       stop
    endif

    return
  end subroutine close_grads_file

end module mod_realinput_grads
!-------------------------------------------------------------------------------
