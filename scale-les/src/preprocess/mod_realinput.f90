!-------------------------------------------------------------------------------
!> module REAL input
!!
!! @par Description
!!          read data from file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-04-28 (R.Yoshida)   [new]
!!
!<
module mod_realinput
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_land_grid_index
  use scale_grid_real, only: &
     LON  => REAL_LON,  &
     LAT  => REAL_LAT,  &
     LONX => REAL_LONX, &
     LATY => REAL_LATY, &
     CZ   => REAL_CZ,   &
     FZ   => REAL_FZ
  use scale_grid_nest, only: &
     NEST_INTERP_LEVEL
  use scale_index
  use scale_tracer
  use scale_process, only: &
     myrank => PRC_myrank,  &
     PRC_master,            &
     PRC_MPIstop
  use scale_external_io, only: &
     iSCALE, &
     iWRFARW, &
     iNICAM, &
     iGrADS
  use scale_comm, only: &
     COMM_bcast
  use scale_interpolation_nest, only: &
     INTRPNEST_domain_compatibility, &
     INTRPNEST_interp_fact_llz,      &
     INTRPNEST_interp_fact_latlon,   &
     INTRPNEST_interp_2d,            &
     INTRPNEST_interp_3d

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtomSetup
  public :: ParentAtomInput
  public :: ParentAtomBoundary
  public :: ParentSurfaceInput

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: ndims = 11 !> 1-3: normal, 4-6: staggerd, 7-9: soil-layer, 10-11: sst

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: interp_OceanLand_data
  private :: diagnose_number_concentration

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), allocatable :: lon_org (:,:)
  real(RP), allocatable :: lat_org (:,:)
  real(RP), allocatable :: cz_org(:,:,:)

  real(RP), allocatable :: dens_org(:,:,:)
  real(RP), allocatable :: qtrc_org(:,:,:,:)

  real(RP), allocatable :: velz_org(:,:,:)
  real(RP), allocatable :: velx_org(:,:,:)
  real(RP), allocatable :: vely_org(:,:,:)
  real(RP), allocatable :: pott_org(:,:,:)
  real(RP), allocatable :: temp_org(:,:,:)
  real(RP), allocatable :: pres_org(:,:,:)

  real(RP), allocatable :: hfact(:,:,:)
  real(RP), allocatable :: vfact(:,:,:,:,:)
  real(RP), allocatable :: vfactl(:,:,:,:,:)
  integer,  allocatable :: igrd (:,:,:)
  integer,  allocatable :: jgrd (:,:,:)
  integer,  allocatable :: kgrd (:,:,:,:,:)
  integer,  allocatable :: kgrdl(:,:,:,:,:)
  integer,  allocatable :: ncopy(:,:,:)

  integer, private :: itp_nh = 4
  integer, private :: itp_nv = 2

  integer, private :: interp_search_divnum = 10
  logical, private :: wrfout = .false.  ! file type switch (wrfout or wrfrst)

  integer, parameter :: cosin = 1
  integer, parameter :: sine  = 2

  integer, private   :: io_fid_grads_nml  = -1
  integer, private   :: io_fid_grads_data = -1

  logical, private   :: do_read
  logical, private   :: rotate
  logical, private   :: use_waterratio
  logical, private   :: use_file_density
  logical, private   :: update_coord
  logical, private   :: use_temp
  logical, private   :: serial
  logical, private   :: first = .true.

  ! replace missing value
  real(RP), private, parameter :: maskval_tg   = 298.0_RP !> mask value 50K => 298K
  real(RP), private, parameter :: maskval_strg = 0.02_RP  !> mask value 0.0 => 0.02
                                ! default value 0.02: set as value of forest at 40% of evapolation rate.
                                ! forest is considered as a typical landuse over Japan area.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ParentAtomSetup( &
      dims,                  &
      timelen,               &
      mdlid,                 &
      basename_org,          &
      filetype,              &
      use_file_density_in,   &
      use_file_landwater_in, &
      serial_in,             &
      search_divnum_in       )
    use scale_external_io, only: &
         iSCALE, &
         iWRFARW, &
         iNICAM, &
         iGrADS
    use mod_realinput_scale, only: &
         ParentAtomSetupSCALE
    use mod_realinput_wrfarw, only: &
         ParentAtomSetupWRFARW
    use mod_realinput_nicam, only: &
         ParentAtomSetupNICAM
    use mod_realinput_grads, only: &
         ParentAtomSetupGrADS
    implicit none

    integer,          intent(out) :: dims(ndims)
    integer,          intent(out) :: timelen
    integer,          intent(out) :: mdlid
    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: filetype
    logical,          intent(in)  :: serial_in             ! read by a serial process
    logical,          intent(in)  :: use_file_density_in   ! use density data from files
    logical,          intent(in)  :: use_file_landwater_in ! use landwater data from files
    integer,          intent(in)  :: search_divnum_in

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Setup]'

    serial = serial_in
    if( serial ) then
       if( myrank == PRC_master ) then
          do_read = .true.
       else
          do_read = .false.
       endif
    else
       do_read = .true.
    endif

    select case(trim(filetype))
    case('SCALE-LES')

       mdlid = iSCALE
       do_read = .true.
       serial = .false.
       if ( do_read ) call ParentAtomSetupSCALE( dims, timelen ) ! (out)
       update_coord = .false.
       use_file_density = use_file_density_in
       use_temp = .false.
       rotate = .false.
       use_waterratio = .false.

    case('WRF-ARW')

       mdlid = iWRFARW
       if ( do_read ) call ParentAtomSetupWRFARW( dims, timelen, & ! (out)
                                                  basename_org   ) ! (in)
       update_coord = .true.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.
       use_waterratio = .true.

    case('NICAM-NETCDF')

       mdlid = iNICAM
       if ( do_read ) call ParentAtomSetupNICAM( dims, timelen, & ! (out)
                                                 basename_org   ) ! (in)
       update_coord = .false.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.
       use_waterratio = .false.

    case('GrADS')

       mdlid = iGrADS
       if ( do_read ) call ParentAtomSetupGrADS( dims, timelen,        & ! (out)
                                                 use_waterratio,       & ! (out)
                                                 use_file_landwater_in ) ! (in)
       update_coord = .true.
       use_file_density = .false.
       use_temp = .true.
       rotate = .true.

    case default

       write(*,*) ' xxx Unsupported FILE TYPE:', trim(filetype)
       call PRC_MPIstop

    endselect

    if( serial ) then
       call COMM_bcast( dims(:), ndims )
       call COMM_bcast( timelen )
       call COMM_bcast( use_waterratio )
    endif

    interp_search_divnum = search_divnum_in
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolation Search Block Dividing Num:', &
                                    interp_search_divnum
    if( IO_L ) write(IO_FID_LOG,*) '+++ Horizontal Interpolation Level:', &
                                    NEST_INTERP_LEVEL
    itp_nh = int( NEST_INTERP_LEVEL )
    itp_nv = 2

    allocate( hfact (        IA, JA, itp_nh         ) )
    allocate( vfact ( KA,    IA, JA, itp_nh, itp_nv ) )
    allocate( vfactl( LKMAX, IA, JA, itp_nh, itp_nv ) )
    allocate( igrd  (        IA, JA, itp_nh         ) )
    allocate( jgrd  (        IA, JA, itp_nh         ) )
    allocate( kgrd  ( KA,    IA, JA, itp_nh, itp_nv ) )
    allocate( kgrdl ( LKMAX, IA, JA, itp_nh, itp_nv ) )
    allocate( ncopy (        IA, JA, itp_nh         ) )

    allocate( lon_org (            dims(2), dims(3) ) )
    allocate( lat_org (            dims(2), dims(3) ) )
    allocate( cz_org  ( dims(1)+2, dims(2), dims(3) ) )

    allocate( velz_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( velx_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( vely_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( pott_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( temp_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( pres_org( dims(1)+2, dims(2), dims(3) ) )
    allocate( qtrc_org( dims(1)+2, dims(2), dims(3), QA ) )
    allocate( dens_org( dims(1)+2, dims(2), dims(3) ) )

    return
  end subroutine ParentAtomSetup

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Read
  subroutine ParentAtomInput( &
      dens,             &
      momz,             &
      momx,             &
      momy,             &
      rhot,             &
      qtrc,             &
      basename_org,     &
      dims,             &
      mdlid,            &
      mptype_parent,    &
      timelen           )
    use scale_comm, only: &
         COMM_vars8, &
         COMM_wait
    use scale_gridtrans, only: &
         rotc => GTRANS_ROTC
    use scale_atmos_thermodyn, only: &
         THERMODYN_pott => ATMOS_THERMODYN_pott
    use scale_atmos_hydrostatic, only: &
         HYDROSTATIC_buildrho_real => ATMOS_HYDROSTATIC_buildrho_real
    use scale_interpolation_nest, only: &
         INTRPNEST_domain_compatibility, &
         INTRPNEST_interp_fact_llz, &
         INTRPNEST_interp_fact_latlon, &
         INTRPNEST_interp_3d
    use mod_realinput_scale, only: &
         ParentAtomOpenSCALE, &
         ParentAtomInputSCALE
    use mod_realinput_wrfarw, only: &
         ParentAtomOpenWRFARW, &
         ParentAtomInputWRFARW
    use mod_realinput_nicam, only: &
         ParentAtomOpenNICAM, &
         ParentAtomInputNICAM
    use mod_realinput_grads, only: &
         ParentAtomOpenGrADS, &
         ParentAtomInputGrADS
    implicit none

    real(RP),         intent(out) :: dens(:,:,:,:)
    real(RP),         intent(out) :: momz(:,:,:,:)
    real(RP),         intent(out) :: momx(:,:,:,:)
    real(RP),         intent(out) :: momy(:,:,:,:)
    real(RP),         intent(out) :: rhot(:,:,:,:)
    real(RP),         intent(out) :: qtrc(:,:,:,:,:)
    character(LEN=*), intent(in)  :: basename_org
    integer,          intent(in)  :: dims(ndims)
    integer,          intent(in)  :: mdlid            ! model type id
    integer,          intent(in)  :: mptype_parent    ! microphysics type of the parent model (number of classes)
    integer,          intent(in)  :: timelen          ! time steps in one file

    real(RP) :: velz  (KA,IA,JA)
    real(RP) :: velx  (KA,IA,JA)
    real(RP) :: vely  (KA,IA,JA)
    real(RP) :: llvelx(KA,IA,JA)
    real(RP) :: llvely(KA,IA,JA)
    real(RP) :: work  (KA,IA,JA)
    real(RP) :: pott  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: qc(KA,IA,JA)

    integer :: k, i, j, iq
    integer :: n
    character(len=H_SHORT)  :: mptype_run
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input]'

    select case (TRACER_TYPE)
    case ("DRY")
       mptype_run = 'dry'
    case ("KESSLER")
       mptype_run = 'single'
    case ("TOMITA08")
       mptype_run = 'single'
    case ("SN14")
       mptype_run = 'double'
    case ("SUZUKI10")
       mptype_run = 'double'
    case default
       write(*,*) 'xxx Unsupported TRACER_TYPE (', trim(TRACER_TYPE), '). Check!'
       call PRC_MPIstop
    end select


    if ( do_read ) then

       select case ( mdlid )
       case ( iSCALE ) ! TYPE: SCALE-LES

          call ParentAtomOpenSCALE( lon_org, lat_org, & ! (out)
                                    cz_org,           & ! (out)
                                    basename_org,     & ! (in)
                                    dims              ) ! (in)

       case ( iWRFARW ) ! TYPE: WRF-ARW

          call ParentAtomOpenWRFARW

       case ( iNICAM ) ! TYPE: NICAM-NETCDF

          call ParentAtomOpenNICAM( lon_org, lat_org, & ! (out)
                                    cz_org,           & ! (out)
                                    basename_org,     & ! (in)
                                    dims              ) ! (in)

       case ( iGrADS ) ! TYPE: GrADS format

          call ParentAtomOpenGrADS

       end select

    end if

    do n = 1, timelen

       if ( do_read ) then

          select case ( mdlid )
          case ( iSCALE ) ! TYPE: SCALE-LES

             call ParentAtomInputSCALE( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, dens_org, pott_org, & ! (out)
                                        qtrc_org,                     & ! (out)
                                        basename_org, dims, n         ) ! (in)

          case ( iWRFARW ) ! TYPE: WRF-ARW

             call ParentAtomInputWRFARW( velz_org, velx_org, vely_org, & ! (out)
                                         pres_org, temp_org, qtrc_org, & ! (out)
                                         lon_org, lat_org, cz_org,     & ! (out)
                                         basename_org, mptype_parent,  & ! (in)
                                         dims, n                       ) ! (in)

          case ( iNICAM ) ! TYPE: NICAM-NETCDF

             call ParentAtomInputNICAM( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, temp_org, qtrc_org, & ! (out)
                                        basename_org, dims, n         ) ! (in)

          case ( iGrADS ) ! TYPE: GrADS format

             call ParentAtomInputGrADS( velz_org, velx_org, vely_org, & ! (out)
                                        pres_org, temp_org, qtrc_org, & ! (out)
                                        lon_org, lat_org, cz_org,     & ! (out)
                                        basename_org, dims, n         ) ! (in)

          end select

          if ( use_temp ) then
             do j = 1, dims(3)
             do i = 1, dims(2)
             do k = 1, dims(1)+2
                call THERMODYN_pott( pott_org(k,i,j),  & ! [OUT]
                                     temp_org(k,i,j),  & ! [IN]
                                     pres_org(k,i,j),  & ! [IN]
                                     qtrc_org(k,i,j,:) ) ! [IN]
             end do
             end do
             end do
          end if

       end if

       if ( serial ) then
          call COMM_bcast( velz_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( velx_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( vely_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( pott_org, dims(1)+2, dims(2), dims(3) )
          call COMM_bcast( qtrc_org, dims(1)+2, dims(2), dims(3), QA )
          if ( use_file_density ) then
             call COMM_bcast( dens_org, dims(1)+2, dims(2), dims(3) )
          else
             call COMM_bcast( pres_org, dims(1)+2, dims(2), dims(3) )
          end if

          if ( first .or. update_coord ) then

             call COMM_bcast( lon_org, dims(2), dims(3) )
             call COMM_bcast( lat_org, dims(2), dims(3) )
             call COMM_bcast( cz_org,  dims(1)+2, dims(2), dims(3) )

          end if

       end if

       if ( first .or. update_coord ) then

          call INTRPNEST_domain_compatibility( lon_org(:,:), lat_org(:,:), cz_org(:,:,:), &
                                               LON(:,:), LAT(:,:), CZ(KS:KE,:,:) )

          ! full level
          call INTRPNEST_interp_fact_llz( hfact, vfact,               & ! [OUT]
                                          kgrd, igrd, jgrd,           & ! [OUT]
                                          ncopy,                      & ! [OUT]
                                          CZ, LAT, LON,               & ! [IN]
                                          KS, KE, IA, JA,             & ! [IN]
                                          cz_org, lat_org, lon_org,   & ! [IN]
                                          dims(1)+2, dims(2), dims(3) ) ! [IN]

       end if

       call INTRPNEST_interp_3d( velz(:,:,:),      &
                                 velz_org(:,:,:),  &
                                 hfact(:,:,:),     &
                                 vfact(:,:,:,:,:), &
                                 kgrd(:,:,:,:,:),  &
                                 igrd(:,:,:),      &
                                 jgrd(:,:,:),      &
                                 IA, JA, KS, KE-1  )

       call INTRPNEST_interp_3d( llvelx  (:,:,:),     &
                                 velx_org(:,:,:),     &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE       )

       call INTRPNEST_interp_3d( llvely  (:,:,:),     &
                                 vely_org(:,:,:),     &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE       )

       if ( rotate ) then
          ! convert from latlon coordinate to local mapping (x)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             work(k,i,j) = llvelx(k,i,j) * rotc(i,j,cosin) + llvely(k,i,j) * rotc(i,j,sine )
          end do
          end do
          end do

          ! from scalar point to staggered point
          do j = 1, JA
          do i = 1, IA-1
          do k = KS, KE
             velx(k,i,j) = ( work(k,i+1,j) + work(k,i,j) ) * 0.5_RP
          end do
          end do
          end do
          do j = 1, JA
          do k = KS, KE
             velx(k,IA,j) = work(k,IA,j)
          end do
          end do
          velx(KS-1,:,:) = 0.0_RP
          velx(KS-2,:,:) = 0.0_RP
          call COMM_vars8( velx(:,:,:), 1 )
          call COMM_wait ( velx(:,:,:), 1, .false. )

          ! convert from latlon coordinate to local mapping (y)
          do j = 1, JA
          do i = 1, IA
          do k = KS, KE
             work(k,i,j) = - llvelx(k,i,j) * rotc(i,j,sine ) + llvely(k,i,j) * rotc(i,j,cosin)
          end do
          end do
          end do

          do j = 1, JA-1
          do i = 1, IA
          do k = KS, KE
             vely(k,i,j) = ( work(k,i,j+1) + work(k,i,j) ) * 0.5_RP
          end do
          end do
          end do
          do i = 1, IA
          do k = KS, KE
             vely(k,i,JA) = work(k,i,JA)
          end do
          end do
          vely(KS-1,:,:) = 0.0_RP
          vely(KS-2,:,:) = 0.0_RP
          call COMM_vars8( vely(:,:,:), 1 )
          call COMM_wait ( vely(:,:,:), 1, .false. )

       else
          velx = llvelx
          vely = llvely
       end if

       if( trim(mptype_run)=='double' .and. mptype_parent <= 6 )then
          if( IO_L ) write(IO_FID_LOG,*) '--- Diagnose Number Concentration from Mixing Ratio'
          call diagnose_number_concentration( qtrc_org(:,:,:,:) ) ! [inout]
       endif

       do j = 1, dims(3)
       do i = 1, dims(2)
       do k = 1, dims(1)+2
          do iq = 1, QA
             qtrc_org(k,i,j,iq) = max( qtrc_org(k,i,j,iq), 0.0_RP )
          end do
       end do
       end do
       end do

       call INTRPNEST_interp_3d( pott    (:,:,:),       &
                                 pott_org(:,:,:),       &
                                 hfact   (:,:,:),     &
                                 vfact   (:,:,:,:,:), &
                                 kgrd    (:,:,:,:,:), &
                                 igrd    (:,:,:),     &
                                 jgrd    (:,:,:),     &
                                 IA, JA, KS, KE         )

       do iq = 1, QA
          call INTRPNEST_interp_3d( qtrc    (:,:,:,iq,n),  &
                                    qtrc_org(:,:,:,iq),    &
                                    hfact   (:,:,:),     &
                                    vfact   (:,:,:,:,:), &
                                    kgrd    (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, KS, KE         )
       end do

       if( use_file_density ) then
          ! use logarithmic density to interpolate more accurately

          dens_org = log( dens_org )
          call INTRPNEST_interp_3d( dens    (:,:,:,n),   &
                                    dens_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfact   (:,:,:,:,:), &
                                    kgrd    (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, KS, KE,      &
                                    logwegt=.true.       )
       else

          pres_org = log( pres_org )
          call INTRPNEST_interp_3d( pres    (:,:,:),     &
                                    pres_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfact   (:,:,:,:,:), &
                                    kgrd    (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, KS, KE,      &
                                    logwegt=.true.       )

#ifdef DRY
          qc = 0.0_RP
#else
          if ( I_QC > 0 ) then
             qc(:,:,:) = QTRC(:,:,:,I_QC,n)
          else
             qc = 0.0_RP
          end if
#endif
          ! make density & pressure profile in moist condition
          call HYDROSTATIC_buildrho_real( dens    (:,:,:,n),      & ! [OUT]
                                          temp    (:,:,:),        & ! [OUT]
                                          pres    (:,:,:),        & ! [INOUT]
                                          pott    (:,:,:),        & ! [IN]
                                          qtrc    (:,:,:,I_QV,n), & ! [IN]
                                          qc      (:,:,:)        ) ! [IN]

          call COMM_vars8( dens(:,:,:,n), 1 )
          call COMM_wait ( dens(:,:,:,n), 1 )

       end if


       do j = 1, JA
       do i = 1, IA
       do k = KS, KE-1
          momz(k,i,j,n) = velz(k,i,j) * ( dens(k+1,i  ,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA
          momz(KE,i,j,n) = 0.0_RP
       end do
       end do
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          rhot(k,i,j,n) = pott(k,i,j) * dens(k,i,j,n)
       end do
       end do
       end do
       do j = 1, JA
       do i = 1, IA-1
       do k = KS, KE
          momx(k,i,j,n) = velx(k,i,j) * ( dens(k  ,i+1,j  ,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JA
       do k = KS, KE
          momx(k,IA,j,n) = velx(k,IA,j) * dens(k,IA,j,n)
       end do
       end do
       call COMM_vars8( momx(:,:,:,n), 1 )

       do j = 1, JA-1
       do i = 1, IA
       do k = KS, KE
          momy(k,i,j,n) = vely(k,i,j) * ( dens(k  ,i  ,j+1,n) + dens(k,i,j,n) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          momy(k,i,JA,n) = vely(k,i,JA) * dens(k,i,JA,n)
       end do
       end do
       call COMM_vars8( momy(:,:,:,n), 2 )

       call COMM_wait ( momx(:,:,:,n), 1, .false. )
       call COMM_wait ( momy(:,:,:,n), 2, .false. )

    end do

    first = .false.

    return
  end subroutine ParentAtomInput

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  subroutine ParentAtomBoundary( &
      dens,      &
      momz,      &
      momx,      &
      momy,      &
      rhot,      &
      qtrc,      &
      numsteps,  &
      update_dt, &
      basename,  &
      title      )
    use scale_comm, only: &
       COMM_vars, &
       COMM_wait
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    real(RP),         intent(in)   :: dens(:,:,:,:)
    real(RP),         intent(in)   :: momz(:,:,:,:)
    real(RP),         intent(in)   :: momx(:,:,:,:)
    real(RP),         intent(in)   :: momy(:,:,:,:)
    real(RP),         intent(in)   :: rhot(:,:,:,:)
    real(RP),         intent(in)   :: qtrc(:,:,:,:,:)
    real(RP),         intent(in)   :: update_dt
    character(LEN=*), intent(in)   :: basename
    character(LEN=*), intent(in)   :: title
    integer,          intent(in)   :: numsteps ! total time steps

    character(len=H_MID)  :: atmos_boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    real(RP), allocatable :: buffer(:,:,:,:)

    integer :: k, i, j, n, iq
    integer :: ts, te
    !---------------------------------------------------------------------------

    ts = 1
    te = numsteps

    allocate( buffer(KA,IA,JA,te-ts+1) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Boundary]'

    call FILEIO_write( DENS(:,:,:,ts:te), basename, title,           &
                       'DENS', 'Reference Density', 'kg/m3', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt           )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMZ(k,i,j,n) / ( DENS(k+1,i,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELZ', 'Reference VELZ', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMX(k,i,j,n) / ( DENS(k,i+1,j,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = ts, te
       buffer(:,IA,:,n-ts+1) = buffer(:,IA-1,:,n-ts+1)
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELX', 'Reference VELX', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = 2.0_RP * MOMY(k,i,j,n) / ( DENS(k,i,j+1,n) + DENS(k,i,j,n) )
    end do
    end do
    end do
    end do
    do n = ts, te
       buffer(:,:,JA,n-ts+1) = buffer(:,:,JA-1,n-ts+1)
    end do
    call FILEIO_write( buffer, basename, title,                 &
                       'VELY', 'Reference VELY', 'm/s', 'ZXYT', &
                       atmos_boundary_out_dtype, update_dt      )

    do n = ts, te
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       buffer(k,i,j,n-ts+1) = RHOT(k,i,j,n) / DENS(k,i,j,n)
    end do
    end do
    end do
    end do
    call FILEIO_write( buffer, basename, title,                       &
                       'POTT', 'Reference PT', 'K', 'ZXYT',           &
                       atmos_boundary_out_dtype, update_dt            )

    do iq = 1, QA
       do n = ts, te
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          buffer(k,i,j,n-ts+1) = QTRC(k,i,j,iq,n)
       end do
       end do
       end do
       end do
       call FILEIO_write( buffer, basename, title,                                           &
                          AQ_NAME(iq), 'Reference '//trim(AQ_NAME(iq)), AQ_UNIT(iq), 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt                                )
    end do

    deallocate( buffer )

    return
  end subroutine ParentAtomBoundary

  !-----------------------------------------------------------------------------
  !> Land&Ocean Data Read/Write
  !-----------------------------------------------------------------------------
  subroutine ParentSurfaceInput( &
      dens,                 &
      momz,                 &
      momx,                 &
      momy,                 &
      rhot,                 &
      qtrc,                 &
      basename_org,         &
      dims,                 &
      it,                   &
      use_file_landwater,   &
      init_landwater_ratio, &
      intrp_land_temp,      &
      intrp_land_water,     &
      intrp_land_sfc_temp,  &
      intrp_ocean_temp,     &
      intrp_ocean_sfc_temp, &
      soilwater_DS2VC_flag, &
      mdlid                 )
    use scale_land_grid, only: &
         LCZ  => GRID_LCZ
    use scale_const, only: &
         EPS => CONST_EPS, &
         UNDEF => CONST_UNDEF, &
         I_SW => CONST_I_SW, &
         I_LW => CONST_I_LW
    use scale_comm, only: &
         COMM_bcast, &
         COMM_vars8, &
         COMM_wait
    use scale_landuse, only: &
         lsmask_nest => LANDUSE_frac_land, &
         fact_ocean  => LANDUSE_fact_ocean, &
         fact_land   => LANDUSE_fact_land, &
         fact_urban  => LANDUSE_fact_urban
    use scale_interpolation_nest, only: &
         INTRPNEST_interp_3d, &
         INTRPNEST_interp_2d
    use mod_ocean_vars, only: &
         OCEAN_vars_restart_write, &
         OCEAN_vars_external_in
    use mod_land_vars, only: &
         LAND_vars_restart_write, &
         LAND_vars_external_in, &
         convert_WS2VWC
    use mod_urban_vars, only: &
         URBAN_vars_restart_write, &
         URBAN_vars_external_in
    use mod_atmos_phy_sf_vars, only: &
         ATMOS_PHY_SF_vars_external_in
    use scale_atmos_thermodyn, only: &
         THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_realinput_scale, only: &
         ParentSurfaceInputSCALE
    use mod_realinput_wrfarw, only: &
         ParentSurfaceInputWRFARW
    use mod_realinput_nicam, only: &
         ParentSurfaceInputNICAM
    use mod_realinput_grads, only: &
         ParentSurfaceInputGrADS
    implicit none

    real(RP),         intent(in) :: dens(KA,IA,JA)
    real(RP),         intent(in) :: momz(KA,IA,JA)
    real(RP),         intent(in) :: momx(KA,IA,JA)
    real(RP),         intent(in) :: momy(KA,IA,JA)
    real(RP),         intent(in) :: rhot(KA,IA,JA)
    real(RP),         intent(in) :: qtrc(KA,IA,JA,QA)
    character(LEN=*), intent(in) :: basename_org
    integer,          intent(in) :: dims(ndims)
    integer,          intent(in) :: it
    logical,          intent(in) :: use_file_landwater   ! use land water data from files
    real(RP),         intent(in) :: init_landwater_ratio ! Ratio of land water to storage is constant,
                                                         !          if use_file_landwater is ".false."
    character(LEN=*), intent(in) :: intrp_land_temp
    character(LEN=*), intent(in) :: intrp_land_water
    character(LEN=*), intent(in) :: intrp_land_sfc_temp
    character(LEN=*), intent(in) :: intrp_ocean_temp
    character(LEN=*), intent(in) :: intrp_ocean_sfc_temp
    logical,          intent(in) :: soilwater_DS2VC_flag
    integer,          intent(in) :: mdlid                ! model type id

    integer, parameter :: i_intrp_off  = 0
    integer, parameter :: i_intrp_mask = 1
    integer, parameter :: i_intrp_fill = 2

    integer :: i_intrp_land_temp
    integer :: i_intrp_land_water
    integer :: i_intrp_land_sfc_temp
    integer :: i_intrp_ocean_temp
    integer :: i_intrp_ocean_sfc_temp

    real(RP) :: temp
    real(RP) :: pres

    ! land
    real(RP) :: tg_org    (dims(7),dims(8),dims(9))
    real(RP) :: strg_org  (dims(7),dims(8),dims(9))
    real(RP) :: smds_org  (dims(7),dims(8),dims(9))
    real(RP) :: skint_org (        dims(8),dims(9))
    real(RP) :: lst_org   (        dims(8),dims(9))
    real(RP) :: ust_org   (        dims(8),dims(9))
    real(RP) :: albg_org  (        dims(8),dims(9),2)
    real(RP) :: z0w_org   (        dims(8),dims(9))
    real(RP) :: lsmask_org(        dims(8),dims(9))
    real(RP) :: lz_org    (dims(7)                )
    real(RP) :: lz3d_org  (dims(7),dims(8),dims(9))
    real(RP) :: llon_org  (        dims(8),dims(9))
    real(RP) :: llat_org  (        dims(8),dims(9))
    real(RP) :: work      (        dims(8),dims(9))
    real(RP) :: lmask     (        dims(8),dims(9))

    ! ocean
    real(RP) :: tw_org    (        dims(10),dims(11))
    real(RP) :: sst_org   (        dims(10),dims(11))
    real(RP) :: albw_org  (        dims(10),dims(11),2)
    real(RP) :: olon_org  (        dims(10),dims(11))
    real(RP) :: olat_org  (        dims(10),dims(11))
    real(RP) :: omask     (        dims(10),dims(11))

    ! for missing value
    real(RP) :: hfact_l(dims(8), dims(9), itp_nh)
    integer  :: igrd_l (dims(8), dims(9), itp_nh)
    integer  :: jgrd_l (dims(8), dims(9), itp_nh)
    real(RP) :: hfact_o(dims(10),dims(11),itp_nh)
    integer  :: igrd_o (dims(10),dims(11),itp_nh)
    integer  :: jgrd_o (dims(10),dims(11),itp_nh)

    real(RP) :: tg   (LKMAX,IA,JA)
    real(RP) :: strg (LKMAX,IA,JA)
    real(RP) :: smds (LKMAX,IA,JA)
!    real(RP) :: roff (IA,JA)
!    real(RP) :: qvef (IA,JA)
    real(RP) :: tw   (IA,JA)
    real(RP) :: lst  (IA,JA)
    real(RP) :: ust  (IA,JA)
    real(RP) :: sst  (IA,JA)
    real(RP) :: albw (IA,JA,2)  ! albedo for ocean
    real(RP) :: albg (IA,JA,2)  ! albedo for land
    real(RP) :: albu (IA,JA,2)  ! albedo for urban, this is not used.
    real(RP) :: z0w  (IA,JA)
    real(RP) :: skint(IA,JA)
!    real(RP) :: skinw(IA,JA)
!    real(RP) :: snowq(IA,JA)
!    real(RP) :: snowt(IA,JA)
    real(RP) :: tc_urb(IA,JA)
    real(RP) :: qc_urb(IA,JA)
    real(RP) :: uc_urb(IA,JA)
    real(RP) :: lcz_3d(LKMAX,IA,JA)

    integer :: k, i, j, n

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input-Surface]'

    if( LKMAX < 4 )then
       write(*,*) 'xxx LKMAX less than 4: ', LKMAX
       write(*,*) 'xxx in Real Case, LKMAX should be set more than 4'
       call PRC_MPIstop
    endif


    if ( do_read ) then

       ! Read Data
       select case ( mdlid )
       case ( iSCALE ) ! TYPE: SCALE-LES

          call ParentSurfaceInputSCALE( &
               tg_org, strg_org, tw_org,    & ! (out)
               lst_org, ust_org, sst_org,   & ! (out)
               albw_org, albg_org, z0w_org, & ! (out)
               lsmask_org, lz_org,          & ! (out)
               basename_org, dims,          & ! (in)
               use_file_landwater, it       ) ! (in)
          llon_org = lon_org
          llat_org = lat_org
          olon_org = lon_org
          olat_org = lat_org

       case ( iWRFARW ) ! TYPE: WRF-ARW

          call ParentSurfaceInputWRFARW( &
               tg_org, smds_org, tw_org,     & ! (out)
               lst_org, ust_org, sst_org,    & ! (out)
               albw_org, albg_org, z0w_org,  & ! (out)
               lsmask_org, lz_org,           & ! (out)
               basename_org, dims,           & ! (in)
               use_file_landwater, it        ) ! (in)
          llon_org = lon_org
          llat_org = lat_org
          olon_org = lon_org
          olat_org = lat_org

       case ( iNICAM ) ! TYPE: NICAM-NETCDF

          call ParentSurfaceInputNICAM( &
               tg_org, strg_org,      & ! (out)
               lst_org, sst_org,      & ! (out)
               lsmask_org, lz_org,    & ! (out)
               basename_org, dims,    & ! (in)
               use_file_landwater, it ) ! (in)
          llon_org = lon_org
          llat_org = lat_org
          olon_org = lon_org
          olat_org = lat_org
          tw_org = sst_org
          ust_org = UNDEF
          albw_org = UNDEF
          albg_org = UNDEF
          z0w_org = UNDEF

       case ( iGrADS ) ! TYPE: GrADS

          call ParentSurfaceInputGrADS( &
               tg_org, strg_org, smds_org, & ! (out)
               tw_org, lst_org, sst_org,   & ! (out)
               lsmask_org,                 & ! (out)
               lz_org,                     & ! (out)
               llon_org, llat_org,         & ! (out)
               olon_org, olat_org,         & ! (out)
               basename_org, dims,         & ! (in)
               use_file_landwater, it      ) ! (in)
          ust_org = UNDEF
          albw_org = UNDEF
          albg_org = UNDEF
          z0w_org = UNDEF

       end select

    end if

    if ( serial ) then
       call COMM_bcast( tg_org, dims(7), dims(8), dims(9) )
       if ( use_waterratio ) then
          call COMM_bcast( smds_org, dims(7), dims(8), dims(9) )
       else
          call COMM_bcast( strg_org, dims(7), dims(8), dims(9) )
       end if
       call COMM_bcast( lst_org,         dims(8), dims(9) )
       call COMM_bcast( ust_org,         dims(8), dims(9) )
       call COMM_bcast( albg_org(:,:,1), dims(8), dims(9) )
       call COMM_bcast( albg_org(:,:,2), dims(8), dims(9) )
       call COMM_bcast( z0w_org(:,:),    dims(8), dims(9) )
       call COMM_bcast( lsmask_org(:,:), dims(8), dims(9) )
       call COMM_bcast( lz_org,          dims(7) )
       call COMM_bcast( llon_org(:,:),   dims(8), dims(9) )
       call COMM_bcast( llat_org(:,:),   dims(8), dims(9) )

       call COMM_bcast( tw_org,          dims(10), dims(11) )
       call COMM_bcast( sst_org,         dims(10), dims(11) )
       call COMM_bcast( albw_org(:,:,1), dims(10), dims(11) )
       call COMM_bcast( albw_org(:,:,2), dims(10), dims(11) )
       call COMM_bcast( olon_org(:,:),   dims(10), dims(11) )
       call COMM_bcast( olat_org(:,:),   dims(10), dims(11) )
    end if

    select case ( INTRP_LAND_TEMP )
    case ( 'off' )
       i_intrp_land_temp = i_intrp_off
    case ( 'mask' )
       i_intrp_land_temp = i_intrp_mask
    case ( 'fill' )
       i_intrp_land_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_TEMP is invalid. ', INTRP_LAND_TEMP
       call PRC_MPIstop
    end select
    select case ( INTRP_LAND_SFC_TEMP )
    case ( 'off' )
       i_intrp_land_sfc_temp = i_intrp_off
    case ( 'mask' )
       i_intrp_land_sfc_temp = i_intrp_mask
    case ( 'fill' )
       i_intrp_land_sfc_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_SFC_TEMP is invalid. ', INTRP_LAND_SFC_TEMP
       call PRC_MPIstop
    end select
    select case ( INTRP_LAND_WATER )
    case ( 'off' )
       i_intrp_land_water = i_intrp_off
    case ( 'mask' )
       i_intrp_land_water = i_intrp_mask
    case ( 'fill' )
       i_intrp_land_water = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_LAND_WATER is invalid. ', INTRP_LAND_WATER
       call PRC_MPIstop
    end select
    select case ( INTRP_OCEAN_TEMP )
    case ( 'off' )
       i_intrp_ocean_temp = i_intrp_off
    case ( 'mask' )
       i_intrp_ocean_temp = i_intrp_mask
    case ( 'fill' )
       i_intrp_ocean_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_OCEAN_TEMP is invalid. ', INTRP_OCEAN_TEMP
       call PRC_MPIstop
    end select
    select case ( INTRP_OCEAN_SFC_TEMP )
    case ( 'off' )
       i_intrp_ocean_sfc_temp = i_intrp_off
    case ( 'mask' )
       i_intrp_ocean_sfc_temp = i_intrp_mask
    case ( 'fill' )
       i_intrp_ocean_sfc_temp = i_intrp_fill
    case default
       write(*,*) 'xxx INTRP_OCEAN_SFC_TEMP is invalid. ', INTRP_OCEAN_SFC_TEMP
       call PRC_MPIstop
    end select

    select case ( mdlid )
    case ( iSCALE, iWRFARW, iNICAM )
       i_intrp_land_temp      = i_intrp_mask
       i_intrp_land_sfc_temp  = i_intrp_mask
       i_intrp_land_water     = i_intrp_mask
       i_intrp_ocean_temp     = i_intrp_mask
       i_intrp_ocean_sfc_temp = i_intrp_mask
    end select


    ! fill grid data
    do j = 1, dims(9)
    do i = 1, dims(8)
      lz3d_org(:,i,j) = lz_org(:)
    end do
    end do

    do j = 1, JA
    do i = 1, IA
      lcz_3D(:,i,j) = LCZ(:)
    enddo
    enddo

    ! interpolation facter between outer land grid and ocean grid
    call INTRPNEST_interp_fact_latlon( hfact_l(:,:,:),               & ! [OUT]
                                       igrd_l(:,:,:), jgrd_l(:,:,:), & ! [OUT]
                                       llat_org(:,:), llon_org(:,:), & ! [IN]
                                       dims(8), dims(9),             & ! [IN]
                                       olat_org(:,:), olon_org(:,:), & ! [IN]
                                       dims(10), dims(11)            ) ! [IN]

    call INTRPNEST_interp_fact_latlon( hfact_o(:,:,:),               & ! [OUT]
                                       igrd_o(:,:,:), jgrd_o(:,:,:), & ! [OUT]
                                       olat_org(:,:), olon_org(:,:), & ! [IN]
                                       dims(10), dims(11),           & ! [IN]
                                       llat_org(:,:), llon_org(:,:), & ! [IN]
                                       dims(8), dims(9)              ) ! [IN]

    ! Ocean temp: interpolate over the land
    if ( i_INTRP_OCEAN_TEMP .ne. i_intrp_off ) then
       select case ( i_INTRP_OCEAN_TEMP )
       case ( i_intrp_mask )
          if ( dims(8)==dims(10) .and. dims(9)==dims(11) ) then
             omask = lsmask_org
          else
             write(*,*) 'grid size for land data and sst data is different'
             write(*,*) '"mask" option is not available for INTRP_OCEAN_TEMP'
             call PRC_MPIstop
          end if
       case ( i_intrp_fill )
          call make_mask( omask, tw_org, dims(10), dims(11), landdata=.false.)
       end select
       call interp_OceanLand_data(tw_org, omask, dims(10), dims(11), landdata=.false.)
    end if

    ! SST: interpolate over the land
    if ( i_INTRP_OCEAN_SFC_TEMP .ne. i_intrp_off ) then
       select case ( i_INTRP_OCEAN_SFC_TEMP )
       case ( i_intrp_mask )
          if ( dims(8)==dims(10) .and. dims(9)==dims(11) ) then
             omask = lsmask_org
          else
             write(*,*) 'grid size for land data and sst data is different'
             write(*,*) '"mask" option is not available for INTRP_OCEAN_SFC_TEMP'
             call PRC_MPIstop
          end if
       case ( i_intrp_fill )
          call make_mask( omask, sst_org, dims(10), dims(11), landdata=.false.)
       end select
       call interp_OceanLand_data(sst_org, omask, dims(10), dims(11), landdata=.false.)
    end if

    ! Surface skin temp: interpolate over the ocean
    if ( i_INTRP_LAND_SFC_TEMP .ne. i_intrp_off ) then
       select case ( i_INTRP_LAND_SFC_TEMP )
       case ( i_intrp_mask )
          lmask = lsmask_org
       case ( i_intrp_fill )
          call make_mask( lmask, lst_org, dims(8), dims(9), landdata=.true.)
       case default
          write(*,*) 'xxx INTRP_LAND_SFC_TEMP is invalid. ', INTRP_LAND_SFC_TEMP
          call PRC_MPIstop
       end select
       call interp_OceanLand_data(lst_org, lmask, dims(8), dims(9), landdata=.true.)
    end if

    ! Urban surface temp: interpolate over the ocean
    ! if ( i_INTRP_URB_SFC_TEMP .ne. i_intrp_off ) then
    !   select case ( i_INTRP_URB_SFC_TEMP )
    !   case ( i_intrp_mask )
    !      lmask = lsmask_org
    !   case ( i_intrp_fill )
    !      call make_mask( lmask, ust_org, dims(8), dims(9), landdata=.true.)
    !   case default
    !      write(*,*) 'xxx INTRP_URB_SFC_TEMP is invalid. ', INTRP_URB_SFC_TEMP
    !      call PRC_MPIstop
    !   end select
    !   call interp_OceanLand_data(ust_org, lmask, dims(8), dims(9), landdata=.true.)
    !end if


    !replace lmask using SST and omask using SKINT
    call INTRPNEST_interp_2d( lmask(:,:), sst_org(:,:), hfact_l(:,:,:),      &
                              igrd_l(:,:,:), jgrd_l(:,:,:), dims(8), dims(9) )
    call INTRPNEST_interp_2d( omask(:,:), lst_org(:,:), hfact_o(:,:,:),        &
                              igrd_o(:,:,:), jgrd_o(:,:,:), dims(10), dims(11) )

    call replace_misval_map( sst_org, omask, dims(10), dims(11), "SST")
    call replace_misval_map( tw_org,  omask, dims(10), dims(11), "OCEAN_TEMP")
    call replace_misval_map( lst_org, lmask, dims(8),  dims(9),  "SKINT")

    !!! replace missing value
    do j = 1, dims(9)
    do i = 1, dims(8)
       if ( ust_org(i,j) == UNDEF ) ust_org(i,j) = lst_org(i,j)
!       if ( skinw_org(i,j) == UNDEF ) skinw_org(i,j) = 0.0_RP
!       if ( snowq_org(i,j) == UNDEF ) snowq_org(i,j) = 0.0_RP
!       if ( snowt_org(i,j) == UNDEF ) snowt_org(i,j) = TEM00
       if ( albg_org(i,j,I_LW) == UNDEF ) albg_org(i,j,I_LW) = 0.03_RP  ! emissivity of general ground surface : 0.95-0.98
       if ( albg_org(i,j,I_SW) == UNDEF ) albg_org(i,j,I_SW) = 0.22_RP
       if ( z0w_org(i,j) == UNDEF ) z0w_org(i,j) = 0.001_RP
    end do
    end do
    do j = 1, dims(11)
    do i = 1, dims(10)
       if ( albw_org(i,j,I_LW) == UNDEF ) albw_org(i,j,I_LW) = 0.04_RP  ! emissivity of water surface : 0.96
       if ( albw_org(i,j,I_SW) == UNDEF ) albw_org(i,j,I_SW) = 0.10_RP
    end do
    end do

    ! Land temp: interpolate over the ocean
    if ( i_INTRP_LAND_TEMP .ne. i_intrp_off ) then
       do k = 1, dims(7)
          work(:,:) = tg_org(k,:,:)
          select case ( i_INTRP_LAND_TEMP )
          case ( i_intrp_mask )
             lmask = lsmask_org
          case ( i_intrp_fill )
             call make_mask( lmask, work, dims(8), dims(9), landdata=.true.)
          end select
          call interp_OceanLand_data( work, lmask, dims(8), dims(9), landdata=.true. )
          !replace land temp using skin temp
          call replace_misval_map( work, lst_org, dims(8),  dims(9),  "STEMP")
          tg_org(k,:,:) = work(:,:)
       end do
    end if



    call INTRPNEST_interp_fact_llz( hfact  (:,:,:),               & ! [OUT]
                                    vfactl (:,:,:,:,:),           & ! [OUT]
                                    kgrdl  (:,:,:,:,:),           & ! [OUT]
                                    igrd   (:,:,:), jgrd(:,:,:),  & ! [OUT]
                                    ncopy  (:,:,:),               & ! [OUT]
                                    lcz_3D (:,:,:),               & ! [IN]
                                    LAT    (:,:), LON    (:,:),   & ! [IN]
                                    1, LKMAX, IA, JA,             & ! [IN]
                                    lz3d_org(:,:,:),              & ! [IN]
                                    llat_org(:,:), llon_org(:,:), & ! [IN]
                                    dims(7), dims(8), dims(9),    & ! [IN]
                                    landgrid=.true.               ) ! [IN]

    call INTRPNEST_interp_2d( z0w(:,:),   z0w_org(:,:),   hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( lst(:,:),   lst_org(:,:),   hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( ust(:,:),   ust_org(:,:),   hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( skinw(:,:), skinw_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( snowq(:,:), snowq_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
!    call INTRPNEST_interp_2d( snowt(:,:), snowt_org(:,:), hfact(:,:,:), &
!                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    do n = 1, 2 ! 1:LW, 2:SW
       call INTRPNEST_interp_2d( albg(:,:,n), albg_org(:,:,n), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
    enddo


    call INTRPNEST_interp_3d( tg    (:,:,:),     &
                              tg_org(:,:,:),     &
                              hfact (:,:,:),     &
                              vfactl(:,:,:,:,:), &
                              kgrdl (:,:,:,:,:), &
                              igrd  (:,:,:),     &
                              jgrd  (:,:,:),     &
                              IA, JA, 1, LKMAX-1 )

    ! interpolation
    do j = 1, JA
    do i = 1, IA
       tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo ! i
    enddo ! j

    ! replace values over the ocean
    do k = 1, LKMAX
       call replace_misval_const( tg(k,:,:), maskval_tg, lsmask_nest )
    enddo

    ! Land water: interpolate over the ocean
    if( use_file_landwater )then

       if ( use_waterratio ) then

          if ( i_INTRP_LAND_WATER .ne. i_intrp_off ) then
             do k = 1, dims(7)
                work(:,:) = smds_org(k,:,:)
                select case ( i_INTRP_LAND_WATER )
                case ( i_intrp_mask )
                   lmask = lsmask_org
                case ( i_intrp_fill )
                   call make_mask( lmask, work, dims(8), dims(9), landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, dims(8), dims(9), landdata=.true.)
                lmask(:,:) = init_landwater_ratio
                !replace missing value to init_landwater_ratio
                call replace_misval_map( work, lmask, dims(8),  dims(9),  "SMOISDS")
                smds_org(k,:,:) = work(:,:)
             enddo
          end if

          call INTRPNEST_interp_3d( smds    (:,:,:),     &
                                    smds_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfactl  (:,:,:,:,:), &
                                    kgrdl   (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, 1, LKMAX-1   )
          do k = 1, LKMAX
             strg(k,:,:) = convert_WS2VWC( smds(k,:,:), critical=soilwater_DS2VC_flag )
          end do

       else

          if ( i_INTRP_LAND_WATER .ne. i_intrp_off ) then
             do k = 1, dims(7)
                work(:,:) = strg_org(k,:,:)
                select case ( i_INTRP_LAND_WATER )
                case ( i_intrp_mask )
                   lmask = lsmask_org
                case ( i_intrp_fill )
                   call make_mask( lmask, work, dims(8), dims(9), landdata=.true.)
                end select
                call interp_OceanLand_data(work, lmask, dims(8),dims(9), landdata=.true.)
                lmask(:,:) = maskval_strg
                !replace missing value to init_landwater_ratio
                call replace_misval_map( work, lmask, dims(8),  dims(9),  "SMOIS")
                strg_org(k,:,:) = work(:,:)
             enddo
          end if

          call INTRPNEST_interp_3d( strg    (:,:,:),     &
                                    strg_org(:,:,:),     &
                                    hfact   (:,:,:),     &
                                    vfactl  (:,:,:,:,:), &
                                    kgrdl   (:,:,:,:,:), &
                                    igrd    (:,:,:),     &
                                    jgrd    (:,:,:),     &
                                    IA, JA, 1, LKMAX-1   )
          ! interpolation
          do j = 1, JA
          do i = 1, IA
             strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
          enddo
          enddo

       end if

       ! replace values over the ocean
       do k = 1, LKMAX
          call replace_misval_const( strg(k,:,:), maskval_strg, lsmask_nest )
       enddo

    else  ! not read from boundary file

       smds(:,:,:) = init_landwater_ratio
       ! conversion from water saturation [fraction] to volumetric water content [m3/m3]
       do k = 1, LKMAX
          strg(k,:,:) = convert_WS2VWC( smds(k,:,:), critical=.true. )
       end do

    endif


    ! interporation for ocean variables
    call INTRPNEST_interp_fact_latlon( hfact(:,:,:),                 & ! [OUT]
                                       igrd(:,:,:), jgrd(:,:,:),     & ! [OUT]
                                       LAT(:,:), LON(:,:),           & ! [IN]
                                       IA, JA,                       & ! [IN]
                                       olat_org(:,:), olon_org(:,:), & ! [IN]
                                       dims(10), dims(11)            ) ! [IN]

    call INTRPNEST_interp_2d( tw(:,:), tw_org(:,:), hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    call INTRPNEST_interp_2d( sst(:,:), sst_org(:,:), hfact(:,:,:), &
                              igrd(:,:,:), jgrd(:,:,:), IA, JA )
    do n = 1, 2 ! 1:LW, 2:SW
       call INTRPNEST_interp_2d( albw(:,:,n), albw_org(:,:,n), hfact(:,:,:), &
                                 igrd(:,:,:), jgrd(:,:,:), IA, JA )
    enddo



    ! replace values over the ocean ####
     do j = 1, JA
     do i = 1, IA
        if( abs(lsmask_nest(i,j)-0.0_RP) < EPS )then ! ocean grid
           lst(i,j)   = sst(i,j)
           ust(i,j)   = sst(i,j)
        endif
        skint(i,j) = fact_ocean(i,j) * sst(i,j) &
                   + fact_land (i,j) * lst(i,j) &
                   + fact_urban(i,j) * ust(i,j)
     enddo
     enddo


    do j = 1, JA
    do i = 1, IA
       call THERMODYN_temp_pres( temp,          & ! [OUT]
                                 pres,          & ! [OUT] not used
                                 dens(KS,i,j),  & ! [IN]
                                 rhot(KS,i,j),  & ! [IN]
                                 qtrc(KS,i,j,:) ) ! [IN]

       tc_urb(i,j) = temp
#ifdef DRY
       qc_urb(i,j) = 0.0_RP
#else
       qc_urb(i,j) = qtrc(KS,i,j,I_QV)
#endif
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
       uc_urb(i,j) = max(sqrt( ( momx(KS,i,j) / (dens(KS,i+1,  j)+dens(KS,i,j)) * 2.0_RP )**2.0_RP &
                             + ( momy(KS,i,j) / (dens(KS,  i,j+1)+dens(KS,i,j)) * 2.0_RP )**2.0_RP ), &
                         0.01_RP)
    enddo
    enddo
    do j = 1, JA-1
       uc_urb(IA,j) = max(sqrt( ( momx(KS,IA,j) /  dens(KS,IA,j  ) )**2.0_RP &
                              + ( momy(KS,IA,j) / (dens(KS,IA,j+1)+dens(KS,IA,j)) * 2.0_RP )**2.0_RP ), &
                          0.01_RP)
    enddo
    do i = 1, IA-1
       uc_urb(i,JA) = max(sqrt( ( momx(KS,i,JA) / (dens(KS,i+1,JA)+dens(KS,i,JA)) * 2.0_RP )**2.0_RP &
                              + ( momy(KS,i,JA) /  dens(KS,i  ,JA) )**2.0_RP ), 0.01_RP)
    enddo
    uc_urb(IA,JA) = max(sqrt( ( momx(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP &
                            + ( momy(KS,IA,JA) / dens(KS,IA,JA) )**2.0_RP ), 0.01_RP)

    ! copy albedo of land to urban
    do j = 1, JA
    do i = 1, IA
       albu(i,j,:) = albg(i,j,:)
    enddo
    enddo

    call COMM_vars8( uc_urb, 1 )
    call COMM_wait ( uc_urb, 1, .false. )

    ! Write out: Ocean
    call OCEAN_vars_external_in( tw, sst, albw, z0w, z0w, z0w )

    ! Write out: Land
    call LAND_vars_external_in( tg, strg, lst, albg )

    ! Write out: Urban
    call URBAN_vars_external_in( tc_urb, qc_urb, uc_urb, ust, albu )

    ! Input to PHY_SF Container
    call ATMOS_PHY_SF_vars_external_in( skint, albw(:,:,I_LW), albw(:,:,I_SW), z0w, z0w, z0w )

    return
  end subroutine ParentSurfaceInput


!-------------------------------
  subroutine make_mask( &
      gmask,     & ! (out)
      data,      & ! (in)
      nx,        & ! (in)
      ny,        & ! (in)
      landdata   ) ! (in)
    use scale_const, only: &
       EPS => CONST_EPS, &
       UNDEF => CONST_UNDEF
    implicit none
    real(RP), intent(out)  :: gmask(:,:)
    real(RP), intent(in)   :: data(:,:)
    integer,  intent(in)   :: nx
    integer,  intent(in)   :: ny
    logical,  intent(in)   :: landdata   ! .true. => land data , .false. => ocean data

    real(RP)               :: dd
    integer                :: i,j

    if( landdata )then
       gmask(:,:) = 1.0_RP  ! gmask=1 will be skip in "interp_OceanLand_data"
       dd         = 0.0_RP
    else
       gmask(:,:) = 0.0_RP  ! gmask=0 will be skip in "interp_OceanLand_data"
       dd         = 1.0_RP
    endif

    do j = 1, ny
    do i = 1, nx
       if( abs(data(i,j) - UNDEF) < sqrt(EPS) )then
          gmask(i,j) = dd
       endif
    enddo
    enddo

    return
  end subroutine make_mask
  !-----------------------------------------------------------------------------
  subroutine interp_OceanLand_data( &
      data,      & ! (inout)
      lsmask,    & ! (in)
      nx,        & ! (in)
      ny,        & ! (in)
      landdata,  & ! (in)
      maskval    & ! (out)
      )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout)         :: data(:,:)
    real(RP), intent(in)            :: lsmask(:,:)
    integer,  intent(in)            :: nx
    integer,  intent(in)            :: ny
    logical,  intent(in)            :: landdata   ! .true. => land data , .false. => ocean data
    real(RP), intent(out), optional :: maskval

    integer                 :: untarget_mask
    integer, allocatable    :: imask(:,:),imaskr(:,:)
    real(RP),allocatable    :: newdata(:,:)
    real(RP)                :: flag(8)
    real(RP)                :: ydata(8)
    real(RP)                :: dd

    integer :: i, j, kk

    !---------------------------------------------------------------------------
    allocate( imask  (nx,ny) )
    allocate( imaskr (nx,ny) )
    allocate( newdata(nx,ny) )
    newdata = 0.0_RP
    if( present(maskval) ) maskval = 999.99_RP

    ! search target cell for interpolation
     do j = 1, ny
     do i = 1, nx
        if( abs(lsmask(i,j)-1.0_RP) < EPS )then
           imask(i,j) = 1  ! land grid
        else
           imask(i,j) = 0  ! ocean grid
        endif
     enddo
     enddo
     if ( landdata ) then  ! interpolation for land data
       untarget_mask = 1
     else                  ! interpolation for ocean data
       untarget_mask = 0
     endif

    ! start interpolation

    imaskr = imask
    do kk = 1, 20
    do j  = 1, ny
    do i  = 1, nx
       if ( imask(i,j) == untarget_mask ) then  ! not missing value
           newdata(i,j) = data(i,j)
           cycle
       else

         if ( present(maskval) ) then
         if ( abs(maskval-999.99_RP)<EPS ) then
            if (abs(lsmask(i,j)-0.0_RP) < EPS) maskval = data(i,j)
         endif
         endif

        !--------------------------------------
        ! check data of neighbor grid
        !
        !           flag(i,j,8)
        !
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   6   +   7   +   8   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   4   + (i,j) +   5   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !     |       |       |       |
        !     +   1   +   2   +   3   +
        !     |       |       |       |
        !     +---+---+---+---+---+---+
        !---------------------------------------
        flag    = 0.0_RP
        ydata   = 0.0_RP
        if((i==1).and.(j==1))then
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if((i==nx).and.(j==1))then
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
        else if((i==1).and.(j==ny))then
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
        else if((i==nx).and.(j==ny))then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
        else if(i==1)then
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if(i==nx)then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
        else if(j==1)then
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        else if(j==ny)then
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
        else
           if(imask(i-1,j-1)==untarget_mask) flag(1)=1.0_RP
           if(imask(i,  j-1)==untarget_mask) flag(2)=1.0_RP
           if(imask(i+1,j-1)==untarget_mask) flag(3)=1.0_RP
           if(imask(i-1,j  )==untarget_mask) flag(4)=1.0_RP
           if(imask(i+1,j  )==untarget_mask) flag(5)=1.0_RP
           if(imask(i-1,j+1)==untarget_mask) flag(6)=1.0_RP
           if(imask(i,  j+1)==untarget_mask) flag(7)=1.0_RP
           if(imask(i+1,j+1)==untarget_mask) flag(8)=1.0_RP
        endif

        if(int(sum(flag(:)))>=3)then  ! coast grid : interpolate
           dd = 0.0_RP
           newdata(i,j) = 0.0_RP
           if(int(flag(1))==1) newdata(i,j) = newdata(i,j) + data(i-1,j-1)
           if(int(flag(2))==1) newdata(i,j) = newdata(i,j) + data(i,j-1)
           if(int(flag(3))==1) newdata(i,j) = newdata(i,j) + data(i+1,j-1)
           if(int(flag(4))==1) newdata(i,j) = newdata(i,j) + data(i-1,j)
           if(int(flag(5))==1) newdata(i,j) = newdata(i,j) + data(i+1,j)
           if(int(flag(6))==1) newdata(i,j) = newdata(i,j) + data(i-1,j+1)
           if(int(flag(7))==1) newdata(i,j) = newdata(i,j) + data(i,j+1)
           if(int(flag(8))==1) newdata(i,j) = newdata(i,j) + data(i+1,j+1)

           dd = sum(flag(:))
           newdata(i,j) = newdata(i,j) / dd

           imaskr(i,j) = untarget_mask
        else
          newdata(i,j) = data(i,j)
        endif

       endif ! sea/land

    enddo
    enddo

     imask(:,:) = imaskr(:,:)
     data(:,:)  = newdata(:,:)
    enddo ! kk

    deallocate( imask   )
    deallocate( imaskr  )
    deallocate( newdata )

    return
  end subroutine interp_OceanLand_data

  !-----------------------------------------------------------------------------
  subroutine replace_misval_const( data, maskval, frac_land )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none
    real(RP), intent(inout) :: data(:,:)
    real(RP), intent(in)    :: maskval
    real(RP), intent(in)    :: frac_land(:,:)
    integer                 :: i, j

    do j = 1, JA
    do i = 1, IA
       if( abs(frac_land(i,j)-0.0_RP) < EPS )then ! ocean grid
          data(i,j) = maskval
       endif
    enddo
    enddo

  end subroutine replace_misval_const

  !-----------------------------------------------------------------------------
  subroutine replace_misval_map( data, maskval, nx, ny, elem)
    use scale_const, only: &
       EPS => CONST_EPS, &
       UNDEF => CONST_UNDEF
    implicit none
    real(RP),     intent(inout) :: data(:,:)
    real(RP),     intent(in)    :: maskval(:,:)
    integer,      intent(in)    :: nx, ny
    character(*), intent(in)    :: elem
    integer                     :: i, j

    do j = 1, ny
    do i = 1, nx
       if( abs(data(i,j) - UNDEF) < sqrt(EPS) )then
          if( abs(maskval(i,j) - UNDEF) < sqrt(EPS) )then
             write(*,*) "Data for mask has missing value. ",trim(elem),i,j
          else
             data(i,j) = maskval(i,j)
          endif
       endif
    enddo
    enddo

  end subroutine replace_misval_map

  !-----------------------------------------------------------------------------
  subroutine diagnose_number_concentration( &
      qvars      & ! (inout)
      )
    use scale_const, only: &
         PI => CONST_PI
    implicit none
    real(RP), intent(inout) :: qvars(:,:,:,:)

    real(RP), parameter :: Dc   =  20.D-6  ! typical particle diameter for cloud  [m]
    real(RP), parameter :: Dr   = 200.D-6  ! typical particle diameter for rain   [m]
    real(RP), parameter :: Di   =  80.D-6  ! typical particle diameter for ice    [m]
    real(RP), parameter :: Ds   =  80.D-6  ! typical particle diameter for snow   [m]
    real(RP), parameter :: Dg   = 200.D-6  ! typical particle diameter for grapel [m]
    real(RP), parameter :: RHOw = 1000.D0  ! typical density for warm particles   [kg/m3]
    real(RP), parameter :: RHOf =  100.D0  ! typical density for frozen particles [kg/m3]
    real(RP), parameter :: RHOg =  400.D0  ! typical density for grapel particles [kg/m3]
    real(RP), parameter :: b    =  3.D0    ! assume spherical form

    real(RP) :: piov6
    !---------------------------------------------------------------------------

#ifndef DRY
    piov6 = pi / 6.0_RP

    qvars(:,:,:,I_NC) = qvars(:,:,:,I_QC) / ( (piov6*RHOw) * Dc**b )
    qvars(:,:,:,I_NR) = qvars(:,:,:,I_QR) / ( (piov6*RHOw) * Dr**b )
    qvars(:,:,:,I_NI) = qvars(:,:,:,I_QI) / ( (piov6*RHOf) * Di**b )
    qvars(:,:,:,I_NS) = qvars(:,:,:,I_QS) / ( (piov6*RHOf) * Ds**b )
    qvars(:,:,:,I_NG) = qvars(:,:,:,I_QG) / ( (piov6*RHOg) * Dg**b )
#endif

    return
  end subroutine diagnose_number_concentration

end module mod_realinput
!-------------------------------------------------------------------------------
