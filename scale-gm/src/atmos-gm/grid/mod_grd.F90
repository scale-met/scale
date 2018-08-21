!-------------------------------------------------------------------------------
!> Module grid
!!
!! @par Description
!!         This module is for the management of the icosahedral grid system
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_grd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: GRD_setup
  public :: GRD_input_hgrid
  public :: GRD_output_hgrid
  public :: GRD_input_vgrid
  public :: GRD_output_vgrid
  public :: GRD_makelatlon
  public :: GRD_scaling

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Indentifiers for the directions in the Cartesian coordinate
  integer,  public, parameter :: GRD_XDIR = 1
  integer,  public, parameter :: GRD_YDIR = 2
  integer,  public, parameter :: GRD_ZDIR = 3

  ! Indentifiers for the directions in the spherical coordinate
  integer,  public, parameter :: I_LAT = 1
  integer,  public, parameter :: I_LON = 2

  !====== Horizontal Grid ======
  !
  ! Grid points ( X: CELL CENTER )
  !          .___.
  !         /     \
  !        .   p   .
  !         \ ___ /
  !          '   '
  !
  ! Grid points ( Xt: CELL VERTEX )
  !          p___p
  !         /     \
  !        p       p
  !         \ ___ /
  !          p   p
  !
  ! Grid points ( Xr: CELL ARC )
  !          ._p_.
  !         p     p
  !        .       .
  !         p _ _ p
  !          ' p '
  !
  real(RP), public              :: GRD_rscale ! scaling factor for the radius of the sphere

  real(RP), public, allocatable :: GRD_x     (:,:,:,:)
  real(RP), public, allocatable :: GRD_x_pl  (:,:,:,:)
  real(RP), public, allocatable :: GRD_xt    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_xt_pl (:,:,:,:)
  real(RP), public, allocatable :: GRD_xr    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_xr_pl (:,:,:,:)

  real(RP), public, allocatable :: GRD_s     (:,:,:,:)
  real(RP), public, allocatable :: GRD_s_pl  (:,:,:,:)
  real(RP), public, allocatable :: GRD_st    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_st_pl (:,:,:,:)

  real(RP), public, allocatable :: GRD_LAT   (:,:)
  real(RP), public, allocatable :: GRD_LAT_pl(:,:)
  real(RP), public, allocatable :: GRD_LON   (:,:)
  real(RP), public, allocatable :: GRD_LON_pl(:,:)



  !====== Topography ======
  integer,  public, parameter   :: GRD_ZSFC = 1
  integer,  public, parameter   :: GRD_ZSD  = 2

  real(RP), public, allocatable :: GRD_zs   (:,:,:,:)
  real(RP), public, allocatable :: GRD_zs_pl(:,:,:,:)



  !====== Vertical Grid ======
  !
  ! z coordinate ( actual height )
  !
  !         GRD_vz(1:ADM_gall,  &
  !                1:ADM_kall,  &
  !                1:ADM_lall,  &
  !                GRD_Z:GRD_ZH )
  !
  !         GRD_vz_pl(1:ADM_gall_pl, &
  !                   1:ADM_kall,    &
  !                   1:ADM_lall_pl, &
  !                   GRD_Z:GRD_ZH   )
  !

  integer,  public, parameter   :: GRD_Z  = 1
  integer,  public, parameter   :: GRD_ZH = 2

  real(RP), public              :: GRD_htop ! model top height [m]

  real(RP), public, allocatable :: GRD_gz   (:) ! gsi (z-star) coordinate
  real(RP), public, allocatable :: GRD_gzh  (:) ! gsi (z-star) coordinate at the half point
  real(RP), public, allocatable :: GRD_dgz  (:) ! d(gsi)
  real(RP), public, allocatable :: GRD_dgzh (:) ! d(gsi) at the half point
  real(RP), public, allocatable :: GRD_rdgz (:)
  real(RP), public, allocatable :: GRD_rdgzh(:)

  real(RP), public, allocatable :: GRD_afac (:) ! From the cell center value to the cell wall value
  real(RP), public, allocatable :: GRD_bfac (:) !  A(k-1/2) = ( afact(k) A(k) + bfact(k) * A(k-1) )
  real(RP), public, allocatable :: GRD_cfac (:) ! From the cell wall value to the cell center value
  real(RP), public, allocatable :: GRD_dfac (:) !  A(k) = ( cfact(k) A(k+1/2) + dfact(k) * A(k-1/2) )

  real(RP), public, allocatable :: GRD_vz   (:,:,:,:)
  real(RP), public, allocatable :: GRD_vz_pl(:,:,:,:)

  integer,  public, parameter   :: GRD_grid_type_on_sphere = 1
  integer,  public, parameter   :: GRD_grid_type_on_plane  = 2
  integer,  public              :: GRD_grid_type           = GRD_grid_type_on_sphere

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: GRD_gen_plgrid
  private :: GRD_makearc

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: hgrid_io_mode  = 'ADVANCED'
  character(len=H_SHORT), private :: topo_io_mode   = 'ADVANCED'
  character(len=H_LONG),  private :: hgrid_fname    = ''         ! horizontal grid file
  character(len=H_LONG),  private :: topo_fname     = ''         ! topography file
  character(len=H_LONG),  private :: toposd_fname   = ''         ! topography file

  character(len=H_LONG),  private :: vgrid_fname    = ''         ! vertical grid file
  character(len=H_SHORT), private :: vgrid_scheme   = 'LINEAR'   ! vertical coordinate scheme
  real(RP),               private :: h_efold        = 10000.0_RP ! e-folding height for hybrid vertical coordinate [m]
  real(RP),               private :: hflat          =  -999.0_RP ! [m]
  logical,                private :: output_vgrid   = .false.    ! output verical grid file?

  logical,                private :: hgrid_comm_flg = .true.     ! communicate GRD_x?          [add] T.Ohno 110722
  real(RP),               private :: triangle_size  = 0.0_RP     ! length of sides of triangle [add] T.Ohno 110722

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GRD_setup
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_have_pl
    use scale_const, only: &
       UNDEF  => CONST_UNDEF,  &
       RADIUS => CONST_RADIUS
    use scale_comm_icoA, only:  &
       COMM_data_transfer
    implicit none

    namelist / GRDPARAM / &
       GRD_grid_type,  &
       hgrid_io_mode,  &
       topo_io_mode,   &
       hgrid_fname,    &
       topo_fname,     &
       toposd_fname,   &
       vgrid_fname,    &
       vgrid_scheme,   &
       h_efold,        &
       hflat,          &
       output_vgrid,   &
       hgrid_comm_flg, &
       triangle_size

    real(RP) :: htop
    integer  :: nstart, nend
    integer  :: kflat

    integer  :: ierr
    integer  :: n, k, l, k0

    integer  :: suf, i, j
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[grd]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=GRDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** GRDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist GRDPARAM. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=GRDPARAM)



    !---< horizontal grid >---
    allocate( GRD_x     (ADM_gall   ,k0,ADM_lall   ,              ADM_nxyz) )
    allocate( GRD_x_pl  (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )
    allocate( GRD_xt    (ADM_gall   ,k0,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz) )
    allocate( GRD_xt_pl (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )
    allocate( GRD_xr    (ADM_gall   ,k0,ADM_lall   ,ADM_AI:ADM_AJ,ADM_nxyz) )
    allocate( GRD_xr_pl (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )

    allocate( GRD_s     (ADM_gall   ,k0,ADM_lall   ,              2) )
    allocate( GRD_s_pl  (ADM_gall_pl,k0,ADM_lall_pl,              2) )
    allocate( GRD_st    (ADM_gall   ,k0,ADM_lall   ,ADM_TI:ADM_TJ,2) )
    allocate( GRD_st_pl (ADM_gall_pl,k0,ADM_lall_pl,              2) )

    allocate( GRD_LAT   (ADM_gall   ,ADM_lall   ) )
    allocate( GRD_LAT_pl(ADM_gall_pl,ADM_lall_pl) )
    allocate( GRD_LON   (ADM_gall   ,ADM_lall   ) )
    allocate( GRD_LON_pl(ADM_gall_pl,ADM_lall_pl) )

    GRD_x     (:,:,:,:)   = UNDEF
    GRD_x_pl  (:,:,:,:)   = UNDEF
    GRD_xt    (:,:,:,:,:) = UNDEF
    GRD_xt_pl (:,:,:,:)   = UNDEF
    GRD_xr    (:,:,:,:,:) = UNDEF
    GRD_xr_pl (:,:,:,:)   = UNDEF

    GRD_s     (:,:,:,:)   = UNDEF
    GRD_s_pl  (:,:,:,:)   = UNDEF
    GRD_st    (:,:,:,:,:) = UNDEF
    GRD_st_pl (:,:,:,:)   = UNDEF

    GRD_LAT   (:,:)       = UNDEF
    GRD_LAT_pl(:,:)       = UNDEF
    GRD_LON   (:,:)       = UNDEF
    GRD_LON_pl(:,:)       = UNDEF

    call GRD_input_hgrid( hgrid_fname,  & ![IN]
                          .true.,       & ![IN]
                          hgrid_io_mode ) ![IN]

    ! data transfer for GRD_x (note: do not communicate GRD_xt)
    if( hgrid_comm_flg ) call COMM_data_transfer( GRD_x, GRD_x_pl )

    ! scaling
    if ( GRD_grid_type == GRD_grid_type_on_plane ) then
       call GRD_scaling(triangle_size)
    else
       call GRD_scaling(RADIUS)
    endif

    ! calc latitude/longitude of each grid point
    call GRD_makelatlon

    ! calc position of cell arc
    call GRD_makearc



    !---< surface height >---
    allocate( GRD_zs   (ADM_gall,   k0,ADM_lall,   GRD_ZSFC:GRD_ZSD) )
    allocate( GRD_zs_pl(ADM_gall_pl,k0,ADM_lall_pl,GRD_ZSFC:GRD_ZSD) )
    GRD_zs   (:,:,:,:) = 0.0_RP
    GRD_zs_pl(:,:,:,:) = 0.0_RP

    call GRD_input_topograph( topo_fname,   & ![IN]
                              toposd_fname, & ![IN]
                              topo_io_mode  ) ![IN]



    !---< vertical coordinate >---
    if ( ADM_kall /= ADM_KNONE ) then
       allocate( GRD_gz   (ADM_kall) )
       allocate( GRD_gzh  (ADM_kall) )
       allocate( GRD_dgz  (ADM_kall) )
       allocate( GRD_dgzh (ADM_kall) )
       allocate( GRD_rdgz (ADM_kall) )
       allocate( GRD_rdgzh(ADM_kall) )

       allocate( GRD_afac (ADM_kall) )
       allocate( GRD_bfac (ADM_kall) )
       allocate( GRD_cfac (ADM_kall) )
       allocate( GRD_dfac (ADM_kall) )

       allocate( GRD_vz   (ADM_gall   ,ADM_kall,ADM_lall   ,GRD_Z:GRD_ZH) )
       allocate( GRD_vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_Z:GRD_ZH) )

       call GRD_input_vgrid(vgrid_fname)

       ! calculation of grid intervals ( cell center )
       do k = ADM_kmin-1, ADM_kmax
          GRD_dgz(k) = GRD_gzh(k+1) - GRD_gzh(k)
       enddo
       GRD_dgz(ADM_kmax+1) = GRD_dgz(ADM_kmax)

       ! calculation of grid intervals ( cell wall )
       do k = ADM_kmin, ADM_kmax+1
          GRD_dgzh(k) = GRD_gz(k) - GRD_gz(k-1)
       enddo
       GRD_dgzh(ADM_kmin-1) = GRD_dgzh(ADM_kmin)

       do k = 1, ADM_kall
          GRD_rdgz (k) = 1.0_RP / grd_dgz (k)
          GRD_rdgzh(k) = 1.0_RP / grd_dgzh(k)
       enddo

       ! hight top
       GRD_htop = GRD_gzh(ADM_kmax+1) - GRD_gzh(ADM_kmin)

       ! vertical interpolation factor
       do k = ADM_kmin, ADM_kmax+1
          GRD_afac(k) = ( GRD_gzh(k) - GRD_gz(k-1) ) &
                      / ( GRD_gz (k) - GRD_gz(k-1) )
       enddo
       GRD_afac(ADM_kmin-1) = 1.0_RP

       GRD_bfac(:) = 1.0_RP - GRD_afac(:)

       do k = ADM_kmin, ADM_kmax
          GRD_cfac(k) = ( GRD_gz (k  ) - GRD_gzh(k) ) &
                      / ( GRD_gzh(k+1) - GRD_gzh(k) )
       enddo
       GRD_cfac(ADM_kmin-1) = 1.0_RP
       GRD_cfac(ADM_kmax+1) = 0.0_RP

       GRD_dfac(:) = 1.0_RP - GRD_cfac(:)

       ! setup z-coordinate
       nstart = suf(ADM_gmin,ADM_gmin)
       nend   = suf(ADM_gmax,ADM_gmax)

       select case(vgrid_scheme)
       case('LINEAR')

          ! linear transfromation : (Gal-Chen & Sommerville(1975)
          !    gz = H(z-zs)/(H-zs) -> z = (H-zs)/H * gz + zs

          kflat = -1
          if ( hflat > 0.0_RP ) then ! default = -999.0
             do k = ADM_kmin+1, ADM_kmax+1
                if ( hflat < GRD_gzh(k) ) then
                   kflat = k
                   exit
                endif
             enddo
          endif

          if ( kflat == -1 ) then
             kflat = ADM_kmax + 1
             htop  = GRD_htop
          else
             htop  = GRD_gzh(kflat) - GRD_gzh(ADM_kmin)
          endif

          do l = 1, ADM_lall
             do k = ADM_kmin-1, kflat
             do n = nstart, nend
                GRD_vz(n,k,l,GRD_Z ) = GRD_zs(n,k0,l,GRD_ZSFC) &
                                     + ( htop - GRD_zs(n,k0,l,GRD_ZSFC) ) / htop * GRD_gz(k)
                GRD_vz(n,k,l,GRD_ZH) = GRD_zs(n,k0,l,GRD_ZSFC) &
                                     + ( htop - GRD_zs(n,k0,l,GRD_ZSFC) ) / htop * GRD_gzh(k)
             enddo
             enddo

             if ( kflat < ADM_kmax+1 ) then
                do k = kflat+1, ADM_kmax+1
                do n = nstart, nend
                   GRD_vz(n,k,l,GRD_Z ) = GRD_gz (k)
                   GRD_vz(n,k,l,GRD_ZH) = GRD_gzh(k)
                enddo
                enddo
             endif
          enddo

          if ( PRC_have_pl ) then
             n = ADM_gslf_pl
             do l = 1, ADM_lall_pl
                do k = ADM_kmin-1, kflat
                   GRD_vz_pl(n,k,l,GRD_Z)  = GRD_zs_pl(n,k0,l,GRD_ZSFC) &
                                           + ( htop - GRD_zs_pl(n,k0,l,GRD_ZSFC) ) / htop * GRD_gz(k)
                   GRD_vz_pl(n,k,l,GRD_ZH) = GRD_zs_pl(n,k0,l,GRD_ZSFC) &
                                           + ( htop - GRD_zs_pl(n,k0,l,GRD_ZSFC) ) / htop * GRD_gzh(k)
                enddo

                if ( kflat < ADM_kmax+1 ) then
                   do k = kflat+1, ADM_kmax+1
                      GRD_vz_pl(n,k,l,GRD_Z ) = GRD_gz (k)
                      GRD_vz_pl(n,k,l,GRD_ZH) = GRD_gzh(k)
                   enddo
                endif
             enddo
          endif

       case('HYBRID')

          ! Hybrid transformation : like as Simmons & Buridge(1981)

          do l = 1, ADM_lall
          do k = ADM_kmin-1, ADM_kmax+1
          do n = nstart, nend
             GRD_vz(n,k,l,GRD_Z)  = GRD_gz(k)                               &
                                  + GRD_zs(n,k0,l,GRD_ZSFC)                 &
                                  * sinh( (GRD_htop-GRD_gz (k)) / h_efold ) &
                                  / sinh(  GRD_htop             / h_efold )
             GRD_vz(n,k,l,GRD_ZH) = GRD_gzh(k)                              &
                                  + GRD_zs(n,k0,l,GRD_ZSFC)                 &
                                  * sinh( (GRD_htop-GRD_gzh(k)) / h_efold ) &
                                  / sinh(  GRD_htop             / h_efold )
          enddo
          enddo
          enddo

          if ( PRC_have_pl ) then
             n = ADM_gslf_pl
             do l = 1, ADM_lall_pl
             do k = ADM_kmin-1, ADM_kmax+1
                GRD_vz_pl(n,k,l,GRD_Z)  = GRD_gz(k)                               &
                                        + GRD_zs_pl(n,k0,l,GRD_ZSFC)              &
                                        * sinh( (GRD_htop-GRD_gz (k)) / h_efold ) &
                                        / sinh(  GRD_htop             / h_efold )
                GRD_vz_pl(n,k,l,GRD_ZH) = GRD_gzh(k)                              &
                                        + GRD_zs_pl(n,k0,l,GRD_ZSFC)              &
                                        * sinh( (GRD_htop-GRD_gzh(k)) / h_efold ) &
                                        / sinh(  GRD_htop             / h_efold )
             enddo
             enddo
          endif

       endselect

       ! fill HALO
       call COMM_data_transfer( GRD_vz, GRD_vz_pl )
    else
       allocate( GRD_gz (ADM_KNONE) )
       allocate( GRD_gzh(ADM_KNONE) )
       GRD_gz (:) = 1.0_RP
       GRD_gzh(:) = 1.0_RP
    endif

    !--- output information about grid.
    if ( ADM_kall /= ADM_KNONE ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|======      Vertical Coordinate [m]      ======|'
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|                                               |'
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|          -GRID CENTER-       -GRID INTERFACE- |'
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|  k        gz     d(gz)      gzh    d(gzh)   k |'
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|                                               |'
       k = ADM_kmax + 1
       if( IO_L ) write(IO_FID_LOG,'(5x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
       if( IO_L ) write(IO_FID_LOG,'(5x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | TOA'
       k = ADM_kmax
       if( IO_L ) write(IO_FID_LOG,'(5x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmax'
       if( IO_L ) write(IO_FID_LOG,'(5x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
       do k = ADM_kmax-1, ADM_kmin+1, -1
       if( IO_L ) write(IO_FID_LOG,'(5x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        |'
       if( IO_L ) write(IO_FID_LOG,'(5x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
       enddo
       k = ADM_kmin
       if( IO_L ) write(IO_FID_LOG,'(5x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmin'
       if( IO_L ) write(IO_FID_LOG,'(5x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | ground'
       k = ADM_kmin-1
       if( IO_L ) write(IO_FID_LOG,'(5x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
       if( IO_L ) write(IO_FID_LOG,'(5x,A)')             '|===============================================|'

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- Vertical layer scheme = ', trim(vgrid_scheme)
       if ( vgrid_scheme == 'HYBRID' ) then
          if( IO_L ) write(IO_FID_LOG,*) '--- e-folding height = ', h_efold
       endif

       if ( output_vgrid ) then
          if ( PRC_IsMaster ) then
             call GRD_output_vgrid('./vgrid_used.dat')
          endif
       endif
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '--- vartical layer = 1'
    endif

    return
  end subroutine GRD_setup

  !-----------------------------------------------------------------------------
  !> Input horizontal grid
  subroutine GRD_input_hgrid( &
       basename,     &
       input_vertex, &
       io_mode       )
    use scale_prc, only: &
       PRC_abort
    use mod_fio, only: &
       FIO_input
    implicit none

    character(len=*), intent(in) :: basename     ! input basename
    logical,          intent(in) :: input_vertex ! flag of B-grid input
    character(len=*), intent(in) :: io_mode      ! io_mode
    !---------------------------------------------------------------------------

    if ( io_mode == 'ADVANCED' ) then

       call FIO_input( GRD_x(:,:,:,GRD_XDIR),basename,'grd_x_x','ZSSFC1',1,1,1 )
       call FIO_input( GRD_x(:,:,:,GRD_YDIR),basename,'grd_x_y','ZSSFC1',1,1,1 )
       call FIO_input( GRD_x(:,:,:,GRD_ZDIR),basename,'grd_x_z','ZSSFC1',1,1,1 )
       if ( input_vertex ) then
          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),basename,'grd_xt_ix','ZSSFC1',1,1,1 )
          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),basename,'grd_xt_jx','ZSSFC1',1,1,1 )
          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),basename,'grd_xt_iy','ZSSFC1',1,1,1 )
          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),basename,'grd_xt_jy','ZSSFC1',1,1,1 )
          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),basename,'grd_xt_iz','ZSSFC1',1,1,1 )
          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),basename,'grd_xt_jz','ZSSFC1',1,1,1 )
       endif

    else
       if( IO_L ) write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_abort
    endif

    call GRD_gen_plgrid

    return
  end subroutine GRD_input_hgrid

  !-----------------------------------------------------------------------------
  !> Output horizontal grid
  subroutine GRD_output_hgrid( &
       basename,      &
       output_vertex, &
       io_mode        )
    use mod_io_param, only: &
       IO_REAL8
    use scale_prc, only: &
       PRC_abort
    use mod_fio, only: &
       FIO_output
    implicit none

    character(len=*), intent(in) :: basename      ! output basename
    logical,          intent(in) :: output_vertex ! output flag of B-grid
    character(len=*), intent(in) :: io_mode       ! io_mode

    character(len=H_MID) :: desc = 'HORIZONTAL GRID FILE'
    !---------------------------------------------------------------------------

    if ( io_mode == 'ADVANCED' ) then

       call FIO_output( GRD_x(:,:,:,GRD_XDIR),                            & ! [IN]
                        basename, desc, '',                               & ! [IN]
                       'grd_x_x', 'GRD_x (X_DIR)', '',                    & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( GRD_x(:,:,:,GRD_YDIR),                            & ! [IN]
                        basename, desc, '',                               & ! [IN]
                       'grd_x_y', 'GRD_x (Y_DIR)', '',                    & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       call FIO_output( GRD_x(:,:,:,GRD_ZDIR),                            & ! [IN]
                        basename, desc, '',                               & ! [IN]
                       'grd_x_z', 'GRD_x (Z_DIR)', '',                    & ! [IN]
                       'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]

       if ( output_vertex ) then
          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_ix', 'GRD_xt (TI,X_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_jx', 'GRD_xt (TJ,X_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_iy', 'GRD_xt (TI,Y_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_jy', 'GRD_xt (TJ,Y_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_iz', 'GRD_xt (TI,Z_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),                    & ! [IN]
                           basename, desc, '',                               & ! [IN]
                          'grd_xt_jz', 'GRD_xt (TJ,Z_DIR)', '',              & ! [IN]
                          'NIL', IO_REAL8, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
       endif

    else
       if( IO_L ) write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_abort
    endif

    return
  end subroutine GRD_output_hgrid

  !-----------------------------------------------------------------------------
  !> Input vertical grid
  subroutine GRD_input_vgrid( fname )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: fname ! vertical grid file name

    integer               :: num_of_layer
    real(DP), allocatable :: gz (:)
    real(DP), allocatable :: gzh(:)

    integer :: fid, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Read vertical grid file: ', trim(fname)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'old',         &
          form   = 'unformatted', &
          iostat = ierr           )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx [GRD_input_vgrid] No vertical grid file.'
          call PRC_abort
       endif

       read(fid) num_of_layer

       allocate( gz (1+num_of_layer+1) )
       allocate( gzh(1+num_of_layer+1) )

       read(fid) gz (:)
       read(fid) gzh(:)

       if ( num_of_layer /= ADM_vlayer ) then
          write(*,*) 'xxx [GRD_input_vgrid] inconsistency in number of vertical layers.'
          call PRC_abort
       endif

       GRD_gz (:) = real(gz ,kind=RP)
       GRD_gzh(:) = real(gzh,kind=RP)

    close(fid)

    return
  end subroutine GRD_input_vgrid

  !-----------------------------------------------------------------------------
  !> Output vertical grid
  subroutine GRD_output_vgrid( fname )
    implicit none

    character(len=*), intent(in) :: fname

    integer :: fid, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Write vertical grid file: ', trim(fname)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'new',         &
          form   = 'unformatted', &
          iostat = ierr           )

       write(fid) ADM_vlayer
       write(fid) GRD_gz
       write(fid) GRD_gzh

    close(fid)

    return
  end subroutine GRD_output_vgrid

  !-----------------------------------------------------------------------------
  !> Input topography data
  subroutine GRD_input_topograph( &
       topo_basename,   &
       toposd_basename, &
       io_mode          )
    use scale_prc, only: &
       PRC_abort
    use scale_comm_icoA, only: &
       COMM_var
    use mod_fio, only: &
       FIO_input
    use mod_ideal_topo, only: &
       IDEAL_topo
    implicit none

    character(len=*), intent(in) :: topo_basename   ! input basename (topography)
    character(len=*), intent(in) :: toposd_basename ! input basename (std.dev. of topography)
    character(len=*), intent(in) :: io_mode         ! io_mode
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** topography data input'

    if ( io_mode == 'ADVANCED' ) then

       if ( topo_basename /= 'NONE' ) then
          call FIO_input(GRD_zs(:,:,:,GRD_ZSFC),topo_basename,  'topo'       ,'ZSSFC1',1,1,1)
       endif
       if ( toposd_basename /= 'NONE' ) then
          call FIO_input(GRD_zs(:,:,:,GRD_ZSD ),toposd_basename,'topo_stddev','ZSSFC1',1,1,1)
       endif

    elseif( io_mode == 'IDEAL' ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** make ideal topography'

       call IDEAL_topo( GRD_s (:,:,:,I_LAT),   & ! [IN]
                        GRD_s (:,:,:,I_LON),   & ! [IN]
                        GRD_zs(:,:,:,GRD_ZSFC) ) ! [OUT]

    else
       write(*,*) 'xxx [grd/GRD_input_topograph] Invalid io_mode!', trim(io_mode)
       call PRC_abort
    endif ! io_mode

    call COMM_var( GRD_zs, GRD_zs_pl, ADM_KNONE, 2 )

    return
  end subroutine GRD_input_topograph

  !-----------------------------------------------------------------------------
  !> Communicate grid data for pole region: This routine is NOT same as COMM_var
  subroutine GRD_gen_plgrid
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_myrank,     &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_have_pl,      &
       PRC_RGN_vlink,    &
       PRC_RGN_total,    &
       PRC_RGN_local,    &
       PRC_RGN_vert_num, &
       PRC_RGN_vert_tab, &
       PRC_RGN_l2r,      &
       PRC_RGN_r2lp,     &
       PRC_RGN_r2p_pl,   &
       I_prc,            &
       I_RGNID,          &
       I_N,              &
       I_S,              &
       I_NPL,            &
       I_SPL
    use scale_comm_icoA, only: &
       COMM_var
    implicit none

    integer :: prctab   (PRC_RGN_vlink)
    integer :: rgntab   (PRC_RGN_vlink)
    integer :: sreq     (PRC_RGN_vlink)
    integer :: rreq     (PRC_RGN_vlink)
    logical :: send_flag(PRC_RGN_vlink)

    real(RP) :: vsend_pl (ADM_nxyz,PRC_RGN_vlink)
    real(RP) :: vrecv_pl (ADM_nxyz,PRC_RGN_vlink)

    integer  :: datatype

    integer  :: istat(MPI_STATUS_SIZE)
    integer  :: l, n, nn
    integer  :: ierr

    integer  :: suf, i, j
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    if ( RP == DP ) then
       datatype = MPI_DOUBLE_PRECISION
    elseif( RP == SP ) then
       datatype = MPI_REAL
    else
       write(*,*) 'xxx precision is not supportd'
       call PRC_abort
    endif

    !--- send information of grid around north pole from regular region

    ! find region which has the north pole
    do l = PRC_RGN_total, 1, -1
       if ( PRC_RGN_vert_num(I_N,l) == PRC_RGN_vlink ) then
          do n = 1, PRC_RGN_vlink
             nn = n + 1
             if( nn > ADM_vlink ) nn = 1
             rgntab(n) = PRC_RGN_vert_tab(I_RGNID,I_N,l,nn)
             prctab(n) = PRC_RGN_r2lp(I_prc,rgntab(n))
          enddo
          exit
       endif
    enddo

    send_flag(:) = .false.

    ! send grid position from regular region
    do n = 1, PRC_RGN_vlink
       do l = 1, PRC_RGN_local
          if ( PRC_RGN_l2r(l) == rgntab(n) ) then
             vsend_pl(:,n) = GRD_xt(suf(ADM_gmin,ADM_gmax),ADM_KNONE,l,ADM_TJ,:) ! [mod] H.Yashiro 20120525

             call MPI_ISEND( vsend_pl(:,n),         &
                             3,                     &
                             datatype,              &
                             PRC_RGN_r2p_pl(I_NPL), &
                             rgntab(n),             &
                             PRC_LOCAL_COMM_WORLD,  &
                             sreq(n),               &
                             ierr                   )

             send_flag(n) = .true.
          endif
       enddo
    enddo

    if ( PRC_myrank == PRC_RGN_r2p_pl(I_NPL) ) then
       do n = 1, PRC_RGN_vlink
          call MPI_IRECV( vrecv_pl(:,n),        &
                          3,                    &
                          datatype,             &
                          prctab(n),            &
                          rgntab(n),            &
                          PRC_LOCAL_COMM_WORLD, &
                          rreq(n),              &
                          ierr                  )
       enddo
    endif

    ! wait and store
    do n = 1, PRC_RGN_vlink
       if ( send_flag(n) ) then
          call MPI_WAIT(sreq(n),istat,ierr)
       endif
    enddo

    if ( PRC_myrank == PRC_RGN_r2p_pl(I_NPL) ) then
       do n = 1, PRC_RGN_vlink
          call MPI_WAIT(rreq(n),istat,ierr)
          GRD_xt_pl(n+1,ADM_KNONE,I_NPL,:) = vrecv_pl(:,n) ! [mod] H.Yashiro 20120525
       enddo
    endif

    !--- send information of grid around south pole from regular region

    ! find region which has the south pole
    do l = PRC_RGN_total, 1, -1
       if ( PRC_RGN_vert_num(I_S,l) == PRC_RGN_vlink ) then
          do n = 1, PRC_RGN_vlink
             rgntab(n) = PRC_RGN_vert_tab(I_RGNID,I_S,l,n)
             prctab(n) = PRC_RGN_r2lp(I_prc,rgntab(n))
          enddo
          exit
       endif
    enddo

    call MPI_Barrier(PRC_LOCAL_COMM_WORLD,ierr)

    send_flag(:) = .false.

    do n = 1, PRC_RGN_vlink
       do l =1, PRC_RGN_local
          if (PRC_RGN_l2r(l) == rgntab(n) ) then
             vsend_pl(:,n) = GRD_xt(suf(ADM_gmax,ADM_gmin),ADM_KNONE,l,ADM_TI,:) ! [mod] H.Yashiro 20120525

             call MPI_ISEND( vsend_pl(:,n),         &
                             3,                     &
                             datatype,              &
                             PRC_RGN_r2p_pl(I_SPL), &
                             rgntab(n),             &
                             PRC_LOCAL_COMM_WORLD,  &
                             sreq(n),               &
                             ierr                   )

             send_flag(n) = .true.
          endif
       enddo
    enddo

    if ( PRC_myrank == PRC_RGN_r2p_pl(I_SPL) ) then
       do n = 1, PRC_RGN_vlink
          call MPI_IRECV( vrecv_pl(:,n),        &
                          3,                    &
                          datatype,             &
                          prctab(n),            &
                          rgntab(n),            &
                          PRC_LOCAL_COMM_WORLD, &
                          rreq(n),              &
                          ierr                  )
       enddo
    endif

    ! wait and store
    do n = 1, PRC_RGN_vlink
       if ( send_flag(n) ) then
          call MPI_WAIT(sreq(n),istat,ierr)
       endif
    enddo

    if ( PRC_myrank == PRC_RGN_r2p_pl(I_SPL) ) then
       do n = 1, PRC_RGN_vlink
          call MPI_WAIT(rreq(n),istat,ierr)
          GRD_xt_pl(n+1,ADM_KNONE,I_SPL,:) = vrecv_pl(:,n) ! [mod] H.Yashiro 20120525
       enddo
    endif

    ! grid point communication
    call COMM_var( GRD_x, GRD_x_pl, ADM_KNONE, 3 )

    if (      PRC_myrank == PRC_RGN_r2p_pl(I_NPL) &
         .OR. PRC_myrank == PRC_RGN_r2p_pl(I_SPL) ) then
       GRD_xt_pl(ADM_gslf_pl,:,:,:) = GRD_x_pl(ADM_gslf_pl,:,:,:)
    endif

    return
  end subroutine GRD_gen_plgrid

  !-----------------------------------------------------------------------------
  !> scaling grid position on the plane/sphere
  subroutine GRD_scaling( fact )
    use scale_prc_icoA, only: &
       PRC_have_pl
    implicit none

    real(RP), intent(in) :: fact ! scaling factor
    !---------------------------------------------------------------------------

    GRD_x (:,:,:,:)   = GRD_x (:,:,:,:)   * fact
    GRD_xt(:,:,:,:,:) = GRD_xt(:,:,:,:,:) * fact

    if ( PRC_have_pl ) then
       GRD_x_pl (:,:,:,:)   = GRD_x_pl (:,:,:,:) * fact
       GRD_xt_pl(:,:,:,:)   = GRD_xt_pl(:,:,:,:) * fact
    endif

    if ( GRD_grid_type == GRD_grid_type_on_plane ) then
       ! do nothing
    else
       ! setting the sphere radius
       GRD_rscale = fact
    endif

    return
  end subroutine GRD_scaling

  !-----------------------------------------------------------------------------
  !> calculate longitude and latitude
  subroutine GRD_makelatlon
    use scale_prc_icoA, only: &
       PRC_have_pl
    use scale_vector, only: &
       VECTR_xyz2latlon
    implicit none

    integer :: g, k0, l
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       call VECTR_xyz2latlon( GRD_x(g,k0,l,GRD_XDIR), & ! [IN]
                              GRD_x(g,k0,l,GRD_YDIR), & ! [IN]
                              GRD_x(g,k0,l,GRD_ZDIR), & ! [IN]
                              GRD_s(g,k0,l,I_LAT),    & ! [OUT]
                              GRD_s(g,k0,l,I_LON)     ) ! [OUT]

       call VECTR_xyz2latlon( GRD_xt(g,k0,l,ADM_TI,GRD_XDIR), & ! [IN]
                              GRD_xt(g,k0,l,ADM_TI,GRD_YDIR), & ! [IN]
                              GRD_xt(g,k0,l,ADM_TI,GRD_ZDIR), & ! [IN]
                              GRD_st(g,k0,l,ADM_TI,I_LAT),    & ! [OUT]
                              GRD_st(g,k0,l,ADM_TI,I_LON)     ) ! [OUT]

       call VECTR_xyz2latlon( GRD_xt(g,k0,l,ADM_TJ,GRD_XDIR), & ! [IN]
                              GRD_xt(g,k0,l,ADM_TJ,GRD_YDIR), & ! [IN]
                              GRD_xt(g,k0,l,ADM_TJ,GRD_ZDIR), & ! [IN]
                              GRD_st(g,k0,l,ADM_TJ,I_LAT),    & ! [OUT]
                              GRD_st(g,k0,l,ADM_TJ,I_LON)     ) ! [OUT]

       GRD_LAT(g,l) = GRD_s(g,k0,l,I_LAT)
       GRD_LON(g,l) = GRD_s(g,k0,l,I_LON)
    enddo
    enddo

    if ( PRC_have_pl ) then
       do l = 1,ADM_lall_pl
       do g = 1,ADM_gall_pl
          call VECTR_xyz2latlon( GRD_x_pl(g,k0,l,GRD_XDIR), & ! [IN]
                                 GRD_x_pl(g,k0,l,GRD_YDIR), & ! [IN]
                                 GRD_x_pl(g,k0,l,GRD_ZDIR), & ! [IN]
                                 GRD_s_pl(g,k0,l,I_LAT),    & ! [OUT]
                                 GRD_s_pl(g,k0,l,I_LON)     ) ! [OUT]

          call VECTR_xyz2latlon( GRD_xt_pl(g,k0,l,GRD_XDIR), & ! [IN]
                                 GRD_xt_pl(g,k0,l,GRD_YDIR), & ! [IN]
                                 GRD_xt_pl(g,k0,l,GRD_ZDIR), & ! [IN]
                                 GRD_st_pl(g,k0,l,I_LAT),    & ! [OUT]
                                 GRD_st_pl(g,k0,l,I_LON)     ) ! [OUT]

          GRD_LAT_pl(g,l) = GRD_s_pl(g,k0,l,I_LAT)
          GRD_LON_pl(g,l) = GRD_s_pl(g,k0,l,I_LON)
       enddo
       enddo
    endif

    return
  end subroutine GRD_makelatlon

  !-----------------------------------------------------------------------------
  !> calculate location of the mid-point of cell arc
  subroutine GRD_makearc
    use scale_prc_icoA, only: &
       PRC_have_pl
    implicit none

    integer  :: ij
    integer  :: im1j, ijm1

    integer  :: nstart,nend
    integer  :: n, l, v, k0

    integer  :: suf, i, j
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d

          GRD_xr(n,k0,l,ADM_AI ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_XDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_XDIR) )
          GRD_xr(n,k0,l,ADM_AI ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_YDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_YDIR) )
          GRD_xr(n,k0,l,ADM_AI ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_ZDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_ZDIR) )
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n

          GRD_xr(n,k0,l,ADM_AIJ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_XDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_XDIR) )
          GRD_xr(n,k0,l,ADM_AIJ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_YDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_YDIR) )
          GRD_xr(n,k0,l,ADM_AIJ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_ZDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_ZDIR) )
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1

          GRD_xr(n,k0,l,ADM_AJ ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_XDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_XDIR) )
          GRD_xr(n,k0,l,ADM_AJ ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_YDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_YDIR) )
          GRD_xr(n,k0,l,ADM_AJ ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_ZDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_ZDIR) )
       enddo
    enddo

    if ( PRC_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             GRD_xr_pl(v,k0,l,GRD_XDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_XDIR)+GRD_xt_pl(ij,k0,l,GRD_XDIR))
             GRD_xr_pl(v,k0,l,GRD_YDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_YDIR)+GRD_xt_pl(ij,k0,l,GRD_YDIR))
             GRD_xr_pl(v,k0,l,GRD_ZDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_ZDIR)+GRD_xt_pl(ij,k0,l,GRD_ZDIR))
          enddo
       enddo
    endif

    return
  end subroutine GRD_makearc

end module mod_grd
