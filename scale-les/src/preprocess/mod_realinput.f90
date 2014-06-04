!-------------------------------------------------------------------------------
!> module FILE I/O (netcdf)
!!
!! @par Description
!!          general file I/O module
!!          frontend interface of netCDF I/O routine
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
  use scale_grid_real, only: &
     lon   => REAL_LON,   &
     lat   => REAL_LAT,   &
     lon_u => REAL_LONX,  &
     lat_v => REAL_LATY,  &
     cz    => REAL_CZ,    &
     fz    => REAL_FZ
  use scale_tracer
  use gtool_file_h
  use gtool_file, only: &
     FileGetShape, &
     FileRead
  use scale_process, only: &
     myrank => PRC_myrank,  &
     PRC_MPIstop
  use scale_const, only:    &
     pi     => CONST_PI,     &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CVvap  => CONST_CVvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     PRE00  => CONST_PRE00,  &
     TEM00  => CONST_TEM00,  &
     grav   => CONST_GRAV,   &
     r_in_m => CONST_RADIUS, &
     dwatr  => CONST_DWATR,  &
     I_SW   => CONST_I_SW,   &
     I_LW   => CONST_I_LW,   &
     CONST_UNDEF
  use scale_atmos_hydrostatic, only: &
     HYDROSTATIC_buildrho_atmos  => ATMOS_HYDROSTATIC_buildrho_atmos
  use scale_external_io
  use scale_land_grid_index
  use scale_land_grid, only: &
     LCZ  => GRID_LCZ
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
  public :: ParentLndOcnInput

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: InputAtomSCALE
  private :: InputAtomWRF
  !private :: InputLndOcnSCALE
  !private :: InputLndOcnWRF
  private :: LatLonZ_Interpolation_Linear
  private :: latlonz_interporation_fact
  private :: haversine

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: iSCALE  = 1
  integer, private, parameter :: iWRFARW = 2
  integer, private, parameter :: iJMAMSM = 3
  integer, private, parameter :: iNICAM  = 4

  real(RP), private, parameter :: large_number_one   = 9.999E+15_RP
  real(RP), private, parameter :: large_number_two   = 8.888E+15_RP
  real(RP), private, parameter :: large_number_three = 7.777E+15_RP
  integer, parameter :: itp_nh = 3
  integer, parameter :: itp_nv = 2

  integer, private :: interp_search_divnum = 10
  logical, private :: wrfout = .false.  ! file type switch (wrfout or wrfrst)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine ParentAtomSetup( &
      dims,             & ! (out)
      timelen,          & ! (out)
      mdlid,            & ! (out)
      basename_org,     & ! (in)
      filetype,         & ! (in)
      search_divnum_in, & ! (in)
      wrf_file_type_in  & ! (in)
      )
    implicit none
    integer,          intent(out)   :: dims(:)
    integer,          intent(out)   :: timelen
    integer,          intent(out)   :: mdlid
    character(LEN=*), intent( in)  :: basename_org
    character(LEN=*), intent( in)  :: filetype
    integer,          intent( in)   :: search_divnum_in
    logical,          intent( in)   :: wrf_file_type_in

    intrinsic shape
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Setup]'

    select case(trim(filetype))
    case('SCALE-LES')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: SCALE-LES'
       mdlid = iSCALE
       call FileGetShape( dims(:),                            &
                           BASENAME_ORG, "DENS", myrank, single=.true. )
       timelen = 1   !temtative approach

    case('WRF-ARW')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Real Case/Atom Input File Type: WRF-ARW'
       mdlid = iWRFARW
       call ExternalFileGetShape( dims(:), timelen, mdlid,        &
                                   BASENAME_ORG, myrank, single=.true. )
       if ( wrf_file_type_in ) then
          wrfout = .true.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF History Output'
       else
          wrfout = .false.
          if( IO_L ) write(IO_FID_LOG,*) '+++ WRF-ARW FILE-TYPE: WRF Restart'
       endif

    !case('JMA-MSM')

    !case('NCEP-FNL')

    !case('NICAM')

    case default
       write(*,*) ' xxx Unsupported FILE TYPE:', trim(filetype)
       call PRC_MPIstop
    endselect

    interp_search_divnum = search_divnum_in
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolation Search Block Dividing Num:', &
                                    interp_search_divnum

    return
  end subroutine ParentAtomSetup

  !-----------------------------------------------------------------------------
  !> Atmosphere Data Read
  !-----------------------------------------------------------------------------
  subroutine ParentAtomInput( &
      dens,          & ! (out)
      momz,          & ! (out)
      momx,          & ! (out)
      momy,          & ! (out)
      rhot,          & ! (out)
      qtrc,          & ! (out)
      basename_org,  & ! (in)
      dims,          & ! (in)
      timelen,       & ! (in)
      mdlid,         & ! (in)
      step           & ! (in)
      )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP),         intent(out)           :: dens(:,:,:)
    real(RP),         intent(out)           :: momz(:,:,:)
    real(RP),         intent(out)           :: momx(:,:,:)
    real(RP),         intent(out)           :: momy(:,:,:)
    real(RP),         intent(out)           :: rhot(:,:,:)
    real(RP),         intent(out)           :: qtrc(:,:,:,:)
    character(LEN=*), intent( in)          :: basename_org
    integer,          intent( in)           :: dims(:)
    integer,          intent( in)           :: timelen  ! time steps in one file
    integer,          intent( in)           :: mdlid    ! model type id
    integer,          intent( in)           :: step     ! initial file number

    real(RP), allocatable :: W(:,:,:)
    real(RP), allocatable :: U(:,:,:)
    real(RP), allocatable :: V(:,:,:)
    real(RP), allocatable :: POTT(:,:,:)

    real(RP), allocatable :: PRES(:,:,:)
    real(RP), allocatable :: TEMP(:,:,:)
    real(RP), allocatable :: QV(:,:,:)
    real(RP), allocatable :: QC(:,:,:)

    integer :: k, i, j

    intrinsic shape
    !---------------------------------------------------------------------------

    allocate( w(KA,IA,JA) )
    allocate( u(KA,IA,JA) )
    allocate( v(KA,IA,JA) )
    allocate( pott(KA,IA,JA) )
    allocate( pres(KA,IA,JA) )
    allocate( temp(KA,IA,JA) )
    allocate( qv(KA,IA,JA) )
    allocate( qc(KA,IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input]'

    if( mdlid == iSCALE )then      !TYPE: SCALE-LES
       !call InputAtomSCALE( dens(:,:,:),         &
       !                      w(:,:,:),            &
       !                      u(:,:,:),            &
       !                      v(:,:,:),            &
       !                      pott(:,:,:),         &
       !                      qtrc(:,:,:,:),       &
       !                      basename_org,        &
       !                      dims(:),             &
       !                      step )
       write(*,*) 'xxx NOT WORK YET: InputAtomSCALE'
       call PRC_MPIstop

    elseif( mdlid == iWRFARW )then !TYPE: WRF-ARW
       call InputAtomWRF(   dens(:,:,:),         &
                             w(:,:,:),            &
                             u(:,:,:),            &
                             v(:,:,:),            &
                             pott(:,:,:),         &
                             qtrc(:,:,:,:),       &
                             basename_org,        &
                             dims(:),             &
                             mdlid,               &
                             "wsm6", step )
    endif

    if ( I_QV > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qv(k,i,j) = qtrc(k,i,j,I_QV)
       end do
       end do
       end do
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qv(k,i,j) = 0.0_RP
       end do
       end do
       end do
    end if

    if ( I_QC > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qc(k,i,j) = qtrc(k,i,j,I_QC)
       end do
       end do
       end do
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qc(k,i,j) = 0.0_RP
       end do
       end do
       end do
    end if

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho_atmos( dens(:,:,:), & ! [INOUT]
                                     temp(:,:,:), & ! [OUT]
                                     pres(:,:,:), & ! [OUT]
                                     pott(:,:,:), & ! [IN]
                                     qv  (:,:,:), & ! [IN]
                                     qc  (:,:,:)  ) ! [IN]

    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       momz(k,i,j) = 0.5_RP * w(k,i,j) * ( dens(k,i,j) + dens(k+1,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       momx(k,i,j) = 0.5_RP * u(k,i,j) * ( dens(k,i+1,j) + dens(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       momy(k,i,j) = 0.5_RP * v(k,i,j) * ( dens(k,i,j+1) + dens(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rhot(k,i,j) = pott(k,i,j) * dens(k,i,j)
    enddo
    enddo
    enddo

    deallocate( w )
    deallocate( u )
    deallocate( v )
    deallocate( pott )
    deallocate( pres )
    deallocate( temp )
    deallocate( qv )
    deallocate( qc )

    return
  end subroutine ParentAtomInput

  !-----------------------------------------------------------------------------
  !> Boundary Data Write
  !-----------------------------------------------------------------------------
  subroutine ParentAtomBoundary( &
      dens,          & ! (in)
      momz,          & ! (in)
      momx,          & ! (in)
      momy,          & ! (in)
      rhot,          & ! (in)
      qtrc,          & ! (in)
      numsteps,      & ! (in)
      initstep,      & ! (in)
      update_dt,     & ! (in)
      basename,      & ! (in)
      title          & ! (in)
      )
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
    integer,          intent(in)   :: initstep ! initial step

    character(len=H_MID)  :: atmos_boundary_out_dtype = 'DEFAULT'  !< REAL4 or REAL8
    real(RP), allocatable :: atmos_boundary_var(:,:,:,:,:)         !> reference container (with HALO)

    integer :: k, i, j, n, iv
    integer :: ts, te, ta
    integer :: timetarg

    intrinsic shape
    !---------------------------------------------------------------------------

    ts = initstep
    te = numsteps
    ta = numsteps

    allocate( atmos_boundary_var(KA,IA,JA,ta,5) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Boundary]'

    ! copy from "ATMOS_BOUNDARY_setinitval/scale_atmos_boundary" //-------
    atmos_boundary_var(:,:,:,:,I_BND_VELZ) = CONST_UNDEF
    atmos_boundary_var(:,:,:,:,I_BND_VELY) = CONST_UNDEF
    atmos_boundary_var(:,:,:,:,I_BND_VELX) = CONST_UNDEF
    atmos_boundary_var(:,:,:,:,I_BND_POTT) = CONST_UNDEF
    atmos_boundary_var(:,:,:,:,I_BND_QV  ) = CONST_UNDEF

    do n = ts, te
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       atmos_boundary_var(k,i,j,n,I_BND_VELZ) = 2.0_RP * MOMZ(k,i,j,n) / ( DENS(k+1,i,j,n) + DENS(k,i,j,n) )
       atmos_boundary_var(k,i,j,n,I_BND_VELX) = 2.0_RP * MOMX(k,i,j,n) / ( DENS(k,i+1,j,n) + DENS(k,i,j,n) )
       atmos_boundary_var(k,i,j,n,I_BND_VELY) = 2.0_RP * MOMY(k,i,j,n) / ( DENS(k,i,j+1,n) + DENS(k,i,j,n) )
       atmos_boundary_var(k,i,j,n,I_BND_POTT) = RHOT(k,i,j,n) / DENS(k,i,j,n)
       atmos_boundary_var(k,i,j,n,I_BND_QV  ) = QTRC(k,i,j,I_QV,n)
    enddo
    enddo
    enddo
    enddo

    ! fill KHALO
    do n = ts, te
    do iv = I_BND_VELZ, I_BND_QV
    do j  = JS, JE
    do i  = IS, IE
       atmos_boundary_var(   1:KS-1,i,j,n,iv) = atmos_boundary_var(KS,i,j,n,iv)
       atmos_boundary_var(KE+1:KA,  i,j,n,iv) = atmos_boundary_var(KE,i,j,n,iv)
    enddo
    enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    do n = ts, te
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars( atmos_boundary_var(:,:,:,n,iv), iv )
    enddo
    enddo

    do n = ts, te
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait( atmos_boundary_var(:,:,:,n,iv), iv )
    enddo
    enddo
    !--------//

    do n = ts, te
       call FILEIO_write( atmos_boundary_var(:,:,:,:,I_BND_VELZ),        &
                          basename, title,                               &
                          'VELZ', 'Reference Velocity w', 'm/s', 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt, timetarg=n )

       call FILEIO_write( atmos_boundary_var(:,:,:,:,I_BND_VELX),        &
                          basename, title,                               &
                          'VELX', 'Reference Velocity u', 'm/s', 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt, timetarg=n )

       call FILEIO_write( atmos_boundary_var(:,:,:,:,I_BND_VELY),        &
                          basename, title,                               &
                          'VELY', 'Reference Velocity y', 'm/s', 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt, timetarg=n )

       call FILEIO_write( atmos_boundary_var(:,:,:,:,I_BND_POTT),        &
                          basename, title,                               &
                          'POTT', 'Reference PT', 'K', 'ZXYT',           &
                          atmos_boundary_out_dtype, update_dt, timetarg=n )

       call FILEIO_write( atmos_boundary_var(:,:,:,:,I_BND_QV),           &
                          basename, title,                                &
                          'QV', 'Reference water vapor', 'kg/kg', 'ZXYT', &
                          atmos_boundary_out_dtype, update_dt, timetarg=n )
    enddo

    deallocate( atmos_boundary_var )

    return
  end subroutine ParentAtomBoundary

  !-----------------------------------------------------------------------------
  !> Land&Ocean Data Read/Write
  !-----------------------------------------------------------------------------
  subroutine ParentLndOcnInput( &
      basename_org,  & ! (in)
      dims,          & ! (in)
      timelen,       & ! (in)
      mdlid,         & ! (in)
      step           & ! (in)
      )
    use mod_land_vars, only:   &
       LAND_vars_restart_write, &
       LAND_vars_external_in
    use mod_ocean_vars, only:   &
       OCEAN_vars_restart_write, &
       OCEAN_vars_external_in
    use mod_cpl_vars, only:   &
       CPL_vars_restart_write, &
       CPL_vars_external_in
    implicit none

    character(LEN=*), intent( in)          :: basename_org
    integer,          intent( in)           :: dims(:)
    integer,          intent( in)           :: timelen  ! time steps in one file
    integer,          intent( in)           :: mdlid    ! model type id
    integer,          intent( in)           :: step     ! initial file number

    real(RP), allocatable :: tg(:,:,:)
    real(RP), allocatable :: strg(:,:,:)
    real(RP), allocatable :: roff(:,:)
    real(RP), allocatable :: qvef(:,:)
    real(RP), allocatable :: tw(:,:)
    real(RP), allocatable :: lst(:,:)
    real(RP), allocatable :: ust(:,:)
    real(RP), allocatable :: sst(:,:)
    real(RP), allocatable :: albw(:,:,:)
    real(RP), allocatable :: albg(:,:,:)
    real(RP), allocatable :: z0w(:,:)
    real(RP), allocatable :: skint(:,:)
    real(RP), allocatable :: skinw(:,:)
    real(RP), allocatable :: snowq(:,:)
    real(RP), allocatable :: snowt(:,:)

    intrinsic shape
    !---------------------------------------------------------------------------

    allocate( tg(LKMAX,IA,JA)   )
    allocate( strg(LKMAX,IA,JA) )
    allocate( roff(IA,JA)       )
    allocate( qvef(IA,JA)       )
    allocate( tw(IA,JA)         )
    allocate( lst(IA,JA)        )
    allocate( ust(IA,JA)        )
    allocate( sst(IA,JA)        )
    allocate( albw(IA,JA,2)     )
    allocate( albg(IA,JA,2)     )
    allocate( z0w(IA,JA)        )
    allocate( skint(IA,JA)      )
    allocate( skinw(IA,JA)      )
    allocate( snowq(IA,JA)      )
    allocate( snowt(IA,JA)      )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[Input-LndOcn]'

    if( LKMAX < 4 )then
       write(*,*) 'xxx LKMAX less than 4: ', LKMAX
       write(*,*) 'xxx in Real Case, LKMAX should be set more than 4'
       call PRC_MPIstop
    endif

    ! Read Data
    if( mdlid == iSCALE )then      !TYPE: SCALE-LES
       !
       ! under construction
       !
       write(*,*) 'xxx NOT WORK YET: InputLndOcnSCALE'
       call PRC_MPIstop
    elseif( mdlid == iWRFARW )then !TYPE: WRF-ARW
       call InputLndOcnWRF( tg(:,:,:),        &
                             strg(:,:,:),      &
                             roff(:,:),        &
                             qvef(:,:),        &
                             tw(:,:),          &
                             lst(:,:),         &
                             ust(:,:),         &
                             sst(:,:),         &
                             albw(:,:,:),      &
                             albg(:,:,:),      &
                             z0w(:,:),         &
                             skint(:,:),       &
                             skinw(:,:),       &
                             snowq(:,:),       &
                             snowt(:,:),       &
                             basename_org,     &
                             dims(:),          &
                             mdlid,            &
                             step              )
    endif

    ! Write out: Land
    call LAND_vars_external_in( tg, strg, roff, qvef )
    call LAND_vars_restart_write

    ! Write out: Ocean
    call OCEAN_vars_external_in( tw )
    call OCEAN_vars_restart_write

    ! Write out: Coupler
    call CPL_vars_external_in( lst, ust, sst, albw, albg, &
                                z0w, skint, skinw, snowq, snowt )
    call CPL_vars_restart_write

    deallocate( tg    )
    deallocate( strg  )
    deallocate( roff  )
    deallocate( qvef  )
    deallocate( tw    )
    deallocate( lst   )
    deallocate( ust   )
    deallocate( sst   )
    deallocate( albw  )
    deallocate( albg  )
    deallocate( z0w   )
    deallocate( skint )
    deallocate( skinw )
    deallocate( snowq )
    deallocate( snowt )

    return
  end subroutine ParentLndOcnInput

  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------
  subroutine InputAtomScale( &
      dens,          & ! (out)
      w,             & ! (out)
      u,             & ! (out)
      v,             & ! (out)
      pott,          & ! (out)
      qtrc,          & ! (out)
      basename_org,  & ! (in)
      dims,          & ! (in)
      step           & ! (in)
      )
    implicit none

    real(RP),         intent(out)           :: dens(:,:,:)
    real(RP),         intent(out)           :: w(:,:,:)
    real(RP),         intent(out)           :: u(:,:,:)
    real(RP),         intent(out)           :: v(:,:,:)
    real(RP),         intent(out)           :: pott(:,:,:)
    real(RP),         intent(out)           :: qtrc(:,:,:,:)
    character(LEN=*), intent( in)          :: basename_org
    integer,          intent( in)           :: dims(:)
    integer,          intent( in)           :: step

    real(RP), allocatable :: DENS_ORG(:,:,:)
    real(RP), allocatable :: MOMZ_ORG(:,:,:)
    real(RP), allocatable :: MOMX_ORG(:,:,:)
    real(RP), allocatable :: MOMY_ORG(:,:,:)
    real(RP), allocatable :: RHOT_ORG(:,:,:)
    real(RP), allocatable :: POTT_ORG(:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:)

    real(RP), allocatable :: W_ORG(:,:,:)
    real(RP), allocatable :: U_ORG(:,:,:)
    real(RP), allocatable :: V_ORG(:,:,:)

    integer :: k, i, j, iq

    !logical :: single_ = .false.

    intrinsic shape
    !---------------------------------------------------------------------------

    allocate( dens_org(dims(1),dims(2),dims(3)) )
    allocate( momz_org(dims(1),dims(2),dims(3)) )
    allocate( momx_org(dims(1),dims(2),dims(3)) )
    allocate( momy_org(dims(1),dims(2),dims(3)) )
    allocate( rhot_org(dims(1),dims(2),dims(3)) )
    allocate( pott_org(dims(1),dims(2),dims(3)) )
    allocate( qtrc_org(dims(1),dims(2),dims(3),QA) )

    allocate( w_org(dims(1),dims(2),dims(3)) )
    allocate( u_org(dims(1),dims(2),dims(3)) )
    allocate( v_org(dims(1),dims(2),dims(3)) )

    ! Need allocation for grid information of parent model (indicated as ORG)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputSCALE]'

    call FileRead( dens_org(:,:,:),                          &
                   BASENAME_ORG, "DENS", step, myrank, single=.true. )
    call FileRead( momz_org(:,:,:),                          &
                   BASENAME_ORG, "MOMZ", step, myrank, single=.true. )
    call FileRead( momx_org(:,:,:),                          &
                   BASENAME_ORG, "MOMX", step, myrank, single=.true. )
    call FileRead( momy_org(:,:,:),                          &
                   BASENAME_ORG, "MOMY", step, myrank, single=.true. )
    call FileRead( rhot_org(:,:,:),                          &
                   BASENAME_ORG, "RHOT", step, myrank, single=.true. )
    do iq = 1, QA
       call FileRead( qtrc_org(:,:,:,iq),                            &
                      BASENAME_ORG, AQ_NAME(iq), step, myrank, single=.true. )
    end do

    ! convert from momentum to velocity
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)-1
       w_org(k,i,j) = 2.0_RP * momz_org(k,i,j) / ( dens_org(k+1,i,j) + dens_org(k,i,j) )
    end do
    end do
    end do
    do j = 1, dims(3)
    do i = 1, dims(2)
       w_org(dims(1),i,j) = 0.0_RP
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)-1
    do k = 1, dims(1)
       u_org(k,i,j) = 2.0_RP * momx_org(k,i,j) / ( dens_org(k,i+1,j) + dens_org(k,i,j) )
    end do
    end do
    end do
    do j = 1, dims(3)
    do k = 1, dims(1)
       u_org(k,dims(2),j) = 2.0_RP * momx_org(k,dims(2),j) / ( dens_org(k,1,j) + dens_org(k,dims(2),j) )
    end do
    end do

    do j = 1, dims(3)-1
    do i = 1, dims(2)
    do k = 1, dims(1)
       v_org(k,i,j) = 2.0_RP * momy_org(k,i,j) / ( dens_org(k,i,j+1) + dens_org(k,i,j) )
    end do
    end do
    end do
    do i = 1, dims(2)
    do k = 1, dims(1)
       v_org(k,i,dims(3)) = 2.0_RP * momy_org(k,i,dims(3)) / ( dens_org(k,i,1) + dens_org(k,i,dims(3)) )
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       pott_org(k,i,j) = rhot_org(k,i,j) / dens_org(k,i,j)
    end do
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       dens_org(k,i,j) = log( dens_org(k,i,j) )
    end do
    end do
    end do

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Create Interpolation system here
    ! ----------------------------------
    ! tentatively dummy copy
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       dens(k,i,j) = dens_org(k,i,j)
       w(k,i,j) = w_org(k,i,j)
       u(k,i,j) = u_org(k,i,j)
       v(k,i,j) = v_org(k,i,j)
       pott(k,i,j) = pott_org(k,i,j)

       do iq = 1, QA
          qtrc(k,i,j,iq) = qtrc_org(k,i,j,iq)
       end do
    end do
    end do
    end do
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    deallocate( dens_org )
    deallocate( momz_org )
    deallocate( momx_org )
    deallocate( momy_org )
    deallocate( rhot_org )
    deallocate( qtrc_org )
    deallocate( w_org )
    deallocate( u_org )
    deallocate( v_org )
    deallocate( pott_org )

    return
  end subroutine InputAtomScale

  !-----------------------------------------------------------------------------
  subroutine InputAtomWRF( &
      dens,          & ! (out)
      w,             & ! (out)
      u,             & ! (out)
      v,             & ! (out)
      pott,          & ! (out)
      qtrc,          & ! (out)
      basename,      & ! (in)
      dims,          & ! (in)
      mdlid,         & ! (in)
      mptype_org,    & ! (in)
      step           & ! (in)
      )
    implicit none

    real(RP),         intent(out)  :: dens(:,:,:)
    real(RP),         intent(out)  :: w(:,:,:)
    real(RP),         intent(out)  :: u(:,:,:)
    real(RP),         intent(out)  :: v(:,:,:)
    real(RP),         intent(out)  :: pott(:,:,:)
    real(RP),         intent(out)  :: qtrc(:,:,:,:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: mptype_org
    integer,          intent( in)  :: mdlid
    integer,          intent( in)  :: dims(:)
    integer,          intent( in)  :: step

    real(RP), allocatable :: DENS_ORG(:,:,:,:)
    real(RP), allocatable :: P_ORG(:,:,:,:)
    real(RP), allocatable :: PBASE_ORG(:,:,:,:)
    real(RP), allocatable :: W_ORG(:,:,:,:)
    real(RP), allocatable :: U_ORG(:,:,:,:)
    real(RP), allocatable :: V_ORG(:,:,:,:)
    real(RP), allocatable :: POTT_ORG(:,:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:,:)
    real(RP), allocatable :: LAT_ORG(:,:,:)
    real(RP), allocatable :: LON_ORG(:,:,:)
    real(RP), allocatable :: LATU_ORG(:,:,:)
    real(RP), allocatable :: LONU_ORG(:,:,:)
    real(RP), allocatable :: LATV_ORG(:,:,:)
    real(RP), allocatable :: LONV_ORG(:,:,:)
    real(RP), allocatable :: GEOF_ORG(:,:,:,:)
    real(RP), allocatable :: GEOH_ORG(:,:,:,:)

    real(RP), allocatable :: AQ_CV(:) !< CV for each hydrometeors [J/kg/K]

    real(SP), allocatable :: dummy_3D(:,:,:)     ! for WRF restart file
    real(SP), allocatable :: dummy_4D(:,:,:,:)   ! for WRF restart file
    real(SP), allocatable :: dummy_5D(:,:,:,:,:) ! for WRF restart file

    real(RP) :: qdry
    real(RP) :: CVtot
    real(RP) :: Rtot
    real(RP) :: RovCP
    real(RP) :: temp
    real(RP) :: t0 = 300.0_RP  ! Defined in WRF
    real(RP) :: d2r

    integer :: n, k, i, j, iq, iqw
    integer :: tcount = 1

    character(len=H_MID)  :: varname_wrf

    intrinsic shape
    !---------------------------------------------------------------------------

    d2r = pi / 180.0_RP

    allocate( AQ_CV(QQA) )
    AQ_CV(I_QV) = CVvap
    if ( QWS /= 0 ) then
       do n = QWS, QWE
          AQ_CV(n) = CL
       enddo
    endif
    if ( QIS /= 0 ) then
       do n = QIS, QIE
          AQ_CV(n) = CI
       enddo
    endif

    allocate( p_org(dims(1),dims(2),dims(3),tcount) )
    allocate( pbase_org(dims(1),dims(2),dims(3),tcount) )
    allocate( dens_org(dims(1),dims(2),dims(3),tcount) )
    allocate( w_org(dims(4),dims(2),dims(3),tcount) )
    allocate( u_org(dims(1),dims(5),dims(3),tcount) )
    allocate( v_org(dims(1),dims(2),dims(6),tcount) )
    allocate( pott_org(dims(1),dims(2),dims(3),tcount) )
    allocate( qtrc_org(dims(1),dims(2),dims(3),tcount,QA) )
    allocate( lat_org(dims(2),dims(3),tcount) )
    allocate( lon_org(dims(2),dims(3),tcount) )
    allocate( latu_org(dims(5),dims(3),tcount) )
    allocate( lonu_org(dims(5),dims(3),tcount) )
    allocate( latv_org(dims(2),dims(6),tcount) )
    allocate( lonv_org(dims(2),dims(6),tcount) )
    allocate( geoh_org(dims(1),dims(2),dims(3),tcount) )
    allocate( geof_org(dims(4),dims(2),dims(3),tcount) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF]'

    allocate( dummy_4D(dims(1),dims(2),dims(3),tcount) )
    call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                   BASENAME, "P",  step, tcount, myrank, mdlid, single=.true. )
    p_org(:,:,:,:) = dummy_4D(:,:,:,:)
    call ExternalFileRead( dummy_4D(:,:,:,:),                         &
                   BASENAME, "PB", step, tcount, myrank, mdlid, single=.true. )
    pbase_org(:,:,:,:) = dummy_4D(:,:,:,:)
    if ( wrfout ) then
       varname_wrf = "T"
    else
       varname_wrf = "T_1"
    endif
    call ExternalFileRead( dummy_4D(:,:,:,:),                          &
                   BASENAME, varname_wrf,  step, tcount, myrank, mdlid, single=.true. )
    pott_org(:,:,:,:) = dummy_4D(:,:,:,:)
    deallocate( dummy_4D )
    allocate( dummy_4D(dims(4),dims(2),dims(3),tcount) )
    if ( wrfout ) then
       varname_wrf = "W"
    else
       varname_wrf = "W_1"
    endif
    call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                   BASENAME, varname_wrf,  step, tcount, myrank, mdlid, single=.true., zstag=.true. )
    w_org(:,:,:,:) = dummy_4D(:,:,:,:)
    deallocate( dummy_4D )
    allocate( dummy_4D(dims(1),dims(5),dims(3),tcount) )
    if ( wrfout ) then
       varname_wrf = "U"
    else
       varname_wrf = "U_1"
    endif
    call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                   BASENAME, varname_wrf,  step, tcount, myrank, mdlid, single=.true., xstag=.true. )
    u_org(:,:,:,:) = dummy_4D(:,:,:,:)
    deallocate( dummy_4D )
    allocate( dummy_4D(dims(1),dims(2),dims(6),tcount) )
    if ( wrfout ) then
       varname_wrf = "V"
    else
       varname_wrf = "V_1"
    endif
    call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                   BASENAME, varname_wrf,  step, tcount, myrank, mdlid, single=.true., ystag=.true. )
    v_org(:,:,:,:) = dummy_4D(:,:,:,:)
    deallocate( dummy_4D )

    do n = 1, tcount
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       p_org(k,i,j,n) = pbase_org(k,i,j,n) + p_org(k,i,j,n)
       pott_org(k,i,j,n) = pott_org(k,i,j,n) + t0
    end do
    end do
    end do
    end do

    allocate( dummy_5D(dims(1),dims(2),dims(3),tcount,QA) )
    dummy_5D(:,:,:,:,:) = 0.0_SP
    if( trim(mptype_org)=='wsm6' .or. &
         trim(mptype_org)=='wdm6' )then
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QV),                     &
                      BASENAME, "QVAPOR", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QC),                     &
                      BASENAME, "QCLOUD", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QR),                     &
                      BASENAME, "QRAIN",  step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QI),                     &
                      BASENAME, "QICE",   step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QS),                     &
                      BASENAME, "QSNOW",  step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_QG),                     &
                      BASENAME, "QGRAUP", step, tcount, myrank, mdlid, single=.true. )
    endif

    if( trim(mptype_org)=='wdm6' )then
       call ExternalFileRead( dummy_5D(:,:,:,:,I_NC),                     &
                      BASENAME, "NC", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_NR),                     &
                      BASENAME, "NR", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_NI),                     &
                      BASENAME, "NI", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_NS),                     &
                      BASENAME, "NS", step, tcount, myrank, mdlid, single=.true. )
       call ExternalFileRead( dummy_5D(:,:,:,:,I_NG),                     &
                      BASENAME, "NG", step, tcount, myrank, mdlid, single=.true. )
    endif

    qtrc_org(:,:,:,:,:) = dummy_5D(:,:,:,:,:)
    deallocate( dummy_5D )

    allocate( dummy_3D(dims(2),dims(3),tcount) )
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLAT",    step, tcount, myrank, mdlid, single=.true. )
    lat_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLONG",   step, tcount, myrank, mdlid, single=.true. )
    lon_org(:,:,:) = dummy_3D(:,:,:) * d2r
    deallocate( dummy_3D )
    allocate( dummy_3D(dims(5),dims(3),tcount) )
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLAT_U",  step, tcount, myrank, mdlid, single=.true., xstag=.true. )
    latu_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLONG_U", step, tcount, myrank, mdlid, single=.true., xstag=.true. )
    lonu_org(:,:,:) = dummy_3D(:,:,:) * d2r
    deallocate( dummy_3D )
    allocate( dummy_3D(dims(2),dims(6),tcount) )
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLAT_V",  step, tcount, myrank, mdlid, single=.true., ystag=.true. )
    latv_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLONG_V", step, tcount, myrank, mdlid, single=.true., ystag=.true. )
    lonv_org(:,:,:) = dummy_3D(:,:,:) * d2r
    deallocate( dummy_3D )
    allocate( dummy_4D(dims(4),dims(2),dims(3),tcount) )
    call ExternalFileRead( dummy_4D(:,:,:,:),                           &
                   BASENAME, "PHB",     step, tcount, myrank, mdlid, single=.true., zstag=.true. )
    geof_org(:,:,:,:) = dummy_4D(:,:,:,:)
    deallocate( dummy_4D )
    ! convert to geopotential height to use as real height in WRF
    do n = 1, tcount
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       geof_org(k,i,j,n) = geof_org(k,i,j,n) / grav
    end do
    end do
    end do
    end do
    ! make half level of geopotential height from face level
    do n = 1, tcount
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       geoh_org(k,i,j,n) = (geof_org(k,i,j,n) + geof_org(k+1,i,j,n))*0.5_RP
    end do
    end do
    end do
    end do

    ! calc dens
    do n = 1, tcount
    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - qtrc_org(k,i,j,n,iqw)
          CVtot = CVtot + qtrc_org(k,i,j,n,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * qtrc_org(k,i,j,n,I_QV)
       RovCP = Rtot / ( CVtot + Rtot )
       temp = pott_org(k,i,j,n) * ( p_org(k,i,j,n) / PRE00 )**RovCP

       dens_org(k,i,j,n) = p_org(k,i,j,n) / ( Rtot * temp )
    end do
    end do
    end do
    end do

    ! make initial condition by interpolation
    call LatLonZ_Interpolation_Linear( dens,w,u,v,pott,qtrc,                               &
                                       dens_org,w_org,u_org,v_org,pott_org,qtrc_org,        &
                                       lat_org,lon_org,latu_org,lonu_org,latv_org,lonv_org, &
                                       geof_org,geoh_org,dims,tcount,step )

    deallocate( p_org )
    deallocate( pbase_org )
    deallocate( dens_org )
    deallocate( w_org )
    deallocate( u_org )
    deallocate( v_org )
    deallocate( pott_org )
    deallocate( qtrc_org )

    deallocate( lat_org )
    deallocate( lon_org )
    deallocate( latu_org )
    deallocate( lonu_org )
    deallocate( latv_org )
    deallocate( lonv_org )
    deallocate( geoh_org )
    deallocate( geof_org )

    deallocate( AQ_CV )

    return
  end subroutine InputAtomWRF

  !-----------------------------------------------------------------------------
  subroutine InputLndOcnWRF( &
      tg,          & ! (out)
      strg,        & ! (out)
      roff,        & ! (out)
      qvef,        & ! (out)
      tw,          & ! (out)
      lst,         & ! (out)
      ust,         & ! (out)
      sst,         & ! (out)
      albw,        & ! (out)
      albg,        & ! (out)
      z0w,         & ! (out)
      skint,       & ! (out)
      skinw,       & ! (out)
      snowq,       & ! (out)
      snowt,       & ! (out)
      basename,    & ! (in)
      dims,        & ! (in)
      mdlid,       & ! (in)
      step         & ! (in)
      )
    implicit none

    real(RP),         intent(out)  :: tg(:,:,:)
    real(RP),         intent(out)  :: strg(:,:,:)
    real(RP),         intent(out)  :: roff(:,:)
    real(RP),         intent(out)  :: qvef(:,:)
    real(RP),         intent(out)  :: tw(:,:)
    real(RP),         intent(out)  :: lst(:,:)
    real(RP),         intent(out)  :: ust(:,:)
    real(RP),         intent(out)  :: sst(:,:)
    real(RP),         intent(out)  :: albw(:,:,:)
    real(RP),         intent(out)  :: albg(:,:,:)
    real(RP),         intent(out)  :: z0w(:,:)
    real(RP),         intent(out)  :: skint(:,:)
    real(RP),         intent(out)  :: skinw(:,:)
    real(RP),         intent(out)  :: snowq(:,:)
    real(RP),         intent(out)  :: snowt(:,:)
    character(LEN=*), intent( in) :: basename
    integer,          intent( in)  :: mdlid
    integer,          intent( in)  :: dims(:)
    integer,          intent( in)  :: step

    real(RP), allocatable :: org_3D(:,:,:)
    real(RP), allocatable :: org_4D(:,:,:,:)
    real(RP), allocatable :: lat_ORG(:,:,:)
    real(RP), allocatable :: lon_ORG(:,:,:)
    real(RP), allocatable :: zs_org(:,:,:,:)
    real(RP), allocatable :: dzs_org(:,:)

    real(SP), allocatable :: dummy_2D(:,:)       ! for WRF restart file
    real(SP), allocatable :: dummy_3D(:,:,:)     ! for WRF restart file
    real(SP), allocatable :: dummy_4D(:,:,:,:)   ! for WRF restart file

    real(RP), allocatable :: lcz_3D(:,:,:)
    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)
    integer, allocatable :: igrd(:,:,:)
    integer, allocatable :: jgrd(:,:,:)
    integer, allocatable :: kgrd(:,:,:,:,:)

    real(RP) :: d2r
    integer :: n, k, i, j, iq, iqw
    integer :: tcount = 1

    intrinsic shape
    !---------------------------------------------------------------------------
    d2r = pi / 180.0_RP

    allocate( lcz_3D(LKMAX,IA,JA) )
    allocate( hfact(IA,JA,itp_nh) )
    allocate( vfact(LKMAX,IA,JA,itp_nh,itp_nv) )
    allocate( igrd(IA,JA,itp_nh) )
    allocate( jgrd(IA,JA,itp_nh) )
    allocate( kgrd(LKMAX,IA,JA,itp_nh,itp_nv) )

    allocate( org_3D(dims(2),dims(3),tcount) )
    allocate( org_4D(dims(7),dims(2),dims(3),tcount) )
    allocate( lat_org(dims(2),dims(3),tcount) )
    allocate( lon_org(dims(2),dims(3),tcount) )
    allocate( zs_org(dims(7),dims(2),dims(3),tcount) )
    allocate( dzs_org(dims(7),tcount) )
    allocate( dummy_2D(dims(7),tcount) )
    allocate( dummy_3D(dims(2),dims(3),tcount) )
    allocate( dummy_4D(dims(7),dims(2),dims(3),tcount) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ ScaleLib/IO[realinput]/Categ[InputWRF-LndOcn]'

    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLAT",    step, tcount, myrank, mdlid, single=.true. )
    lat_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "XLONG",   step, tcount, myrank, mdlid, single=.true. )
    lon_org(:,:,:) = dummy_3D(:,:,:) * d2r
    call ExternalFileRead( dummy_2D(:,:),                               &
                   BASENAME, "ZS",      step, tcount, myrank, mdlid, dims(7), single=.true. )
    do n = 1, tcount
    do j = 1, dims(3)
    do i = 1, dims(2)
       zs_org(:,i,j,n) = dummy_2D(:,n)
    enddo
    enddo
    enddo
    call ExternalFileRead( dummy_2D(:,:),                               &
                   BASENAME, "DZS",      step, tcount, myrank, mdlid, dims(7), single=.true. )
    dzs_org(:,:) = dummy_2D(:,:)

    ! preparing for interpolation
    do j = 1, JA
    do i = 1, IA
       lcz_3D(:,i,j) = lcz(:)
    enddo
    enddo
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,lcz_3D,lat,lon, &
                                      zs_org,lat_org,lon_org,dims(7),dims(2),dims(3),step,landgrid=.true. )

    ! soil temperature [K]
    call ExternalFileRead( dummy_4D(:,:,:,:),                             &
                   BASENAME, "TSLB",  step, tcount, myrank, mdlid, single=.true., landgrid=.true. )
    org_4D(:,:,:,:) = dummy_4D(:,:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = 1, LKMAX-1  ! interpolation
       tg(k,i,j) =  org_4D(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                  + org_4D(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                  + org_4D(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                  + org_4D(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                  + org_4D(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                  + org_4D(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    tg(LKMAX,i,j) = tg(LKMAX-1,i,j)
    enddo
    enddo

    ! soil liquid water m3 m-3] (no wrfout-default)
    ! --- not available; appropriate convert method is not found
    strg(:,:,:) = 0.0_DP
    ! ---
    !call ExternalFileRead( dummy_4D(:,:,:,:),                             &
    !               BASENAME, "SH2O",  step, tcount, myrank, mdlid, single=.true., landgrid=.true. )
    !do n = 1, tcount 
    !do j = 1, dims(3)
    !do i = 1, dims(2)
    !do k = 1, dims(7)
    !   org_4D(k,i,j,n) = dummy_4D(k,i,j,n) * dwatr * dzs_org(k,n)
    !enddo
    !enddo
    !enddo
    !enddo
    !n = step
    !do j = JS-1, JE+1
    !do i = IS-1, IE+1
    !do k = 1, LKMAX-1  ! interpolation
    !   strg(k,i,j) =  org_4D(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
    !                + org_4D(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
    !                + org_4D(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
    !                + org_4D(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
    !                + org_4D(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
    !                + org_4D(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    !enddo
    !strg(LKMAX,i,j) = strg(LKMAX-1,i,j)
    !enddo
    !enddo

    ! surface runoff [mm]
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "SFROFF",  step, tcount, myrank, mdlid, single=.true. )
    do n = 1, tcount 
    do j = 1, dims(3)
    do i = 1, dims(2)
       org_3D(i,j,n) = dummy_3D(i,j,n) * 1000.0_DP * dwatr
    enddo
    enddo
    enddo
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       roff(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                  + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                  + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! accumulated surface evaporation [kg m-2]
    ! --- not available; convert method is not found
    !call ExternalFileRead( dummy_3D(:,:,:),                             &
    !               BASENAME, "SFCEVP",  step, tcount, myrank, mdlid, single=.true. )
    !org_3D(:,:,:) = dummy_3D(:,:,:)
    !n = step
    qvef(:,:) = 0.0_DP

    ! sea surface temperature [K]
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "SST",  step, tcount, myrank, mdlid, single=.true. )
    org_3D(:,:,:) = dummy_3D(:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       tw(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! surface skin temperature [K]
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "TSK",  step, tcount, myrank, mdlid, single=.true. )
    org_3D(:,:,:) = dummy_3D(:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       lst(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                 + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                 + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! urban skin temperature [K] (no wrfout-default)
    ! --- tentative approach;
    !     coz, appropriate value does not retrieve from default wrfout file
    ust(:,:) = lst(:,:)
    ! ---
    !call ExternalFileRead( dummy_3D(:,:,:),                             &
    !               BASENAME, "TS_URB",  step, tcount, myrank, mdlid, single=.true. )
    !org_3D(:,:,:) = dummy_3D(:,:,:)
    !n = step
    !do j = JS-1, JE+1
    !do i = IS-1, IE+1
    !   ust(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
    !             + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
    !             + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    !enddo
    !enddo

    ! skin sea surface temperature [K] (no wrfout-default)
    ! --- tentative approach;
    !     coz, appropriate value does not retrieve from default wrfout file
    sst(:,:) = lst(:,:)
    ! ---
    !call ExternalFileRead( dummy_3D(:,:,:),                             &
    !               BASENAME, "SSTSK",  step, tcount, myrank, mdlid, single=.true. )
    !org_3D(:,:,:) = dummy_3D(:,:,:)
    !n = step
    !do j = JS-1, JE+1
    !do i = IS-1, IE+1
    !   sst(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
    !             + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
    !             + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    !enddo
    !enddo

    ! ALBEDO [-]
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "ALBEDO",  step, tcount, myrank, mdlid, single=.true. )
    org_3D(:,:,:) = dummy_3D(:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       albw(i,j,I_SW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_SW) = albw(:,:,I_SW)

    ! SURFACE EMISSIVITY [-] (no wrfout-default)
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "EMISS",  step, tcount, myrank, mdlid, single=.true. )
    do n = 1, tcount 
    do j = 1, dims(3)
    do i = 1, dims(2)
       org_3D(i,j,n) = 1.0_DP - dummy_3D(i,j,n)
    enddo
    enddo
    enddo
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       albw(i,j,I_LW) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                       + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                       + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo
    albg(:,:,I_LW) = albw(:,:,I_LW)

    ! TIME-VARYING ROUGHNESS LENGTH [m] (no wrfout-default)
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "ZNT",  step, tcount, myrank, mdlid, single=.true. )
    org_3D(:,:,:) = dummy_3D(:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       z0w(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                 + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                 + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    !tentative approach for skin
    skint(:,:) = lst(:,:)
    skinw(:,:) = 0.0_DP

    ! SNOW WATER EQUIVALENT [kg m-2] (no wrfout-default)
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "SNOW",  step, tcount, myrank, mdlid, single=.true. )
    org_3D(:,:,:) = dummy_3D(:,:,:)
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       snowq(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    ! AVERAGE SNOW TEMPERATURE [C] (no wrfout-default)
    call ExternalFileRead( dummy_3D(:,:,:),                             &
                   BASENAME, "TSNAV",  step, tcount, myrank, mdlid, single=.true. )
    do n = 1, tcount 
    do j = 1, dims(3)
    do i = 1, dims(2)
       org_3D(i,j,n) = dummy_3D(i,j,n) + TEM00
    enddo
    enddo
    enddo
    n = step
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       snowt(i,j) =  org_3D(igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) &
                   + org_3D(igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) &
                   + org_3D(igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3)
    enddo
    enddo

    deallocate( dummy_2D )
    deallocate( dummy_3D )
    deallocate( dummy_4D )
    deallocate( org_3D )
    deallocate( org_4D )

    deallocate( lat_org )
    deallocate( lon_org )
    deallocate( zs_org )
    deallocate( dzs_org )

    deallocate( lcz_3D )
    deallocate( hfact )
    deallocate( vfact )
    deallocate( igrd )
    deallocate( jgrd )
    deallocate( kgrd )

    return
  end subroutine InputLndOcnWRF

  !-----------------------------------------------------------------------------
  subroutine LatLonZ_Interpolation_Linear( &
      dens,          & ! (out)
      w,             & ! (out)
      u,             & ! (out)
      v,             & ! (out)
      pott,          & ! (out)
      qtrc,          & ! (out)
      dens_org,      & ! (in)
      w_org,         & ! (in)
      u_org,         & ! (in)
      v_org,         & ! (in)
      pott_org,      & ! (in)
      qtrc_org,      & ! (in)
      lat_org,       & ! (in)
      lon_org,       & ! (in)
      latu_org,      & ! (in)
      lonu_org,      & ! (in)
      latv_org,      & ! (in)
      lonv_org,      & ! (in)
      geof_org,      & ! (in)
      geoh_org,      & ! (in)
      dims,          & ! (in)
      tcount,        & ! (in)
      step           & ! (in)
      )
    implicit none

    real(RP), intent(out)  :: dens(:,:,:)
    real(RP), intent(out)  :: w(:,:,:)
    real(RP), intent(out)  :: u(:,:,:)
    real(RP), intent(out)  :: v(:,:,:)
    real(RP), intent(out)  :: pott(:,:,:)
    real(RP), intent(out)  :: qtrc(:,:,:,:)
    real(RP), intent( in)  :: dens_org(:,:,:,:)
    real(RP), intent( in)  :: w_org(:,:,:,:)
    real(RP), intent( in)  :: u_org(:,:,:,:)
    real(RP), intent( in)  :: v_org(:,:,:,:)
    real(RP), intent( in)  :: pott_org(:,:,:,:)
    real(RP), intent( in)  :: qtrc_org(:,:,:,:,:)
    real(RP), intent( in)  :: lat_org(:,:,:)
    real(RP), intent( in)  :: lon_org(:,:,:)
    real(RP), intent( in)  :: latu_org(:,:,:)
    real(RP), intent( in)  :: lonu_org(:,:,:)
    real(RP), intent( in)  :: latv_org(:,:,:)
    real(RP), intent( in)  :: lonv_org(:,:,:)
    real(RP), intent( in)  :: geof_org(:,:,:,:)
    real(RP), intent( in)  :: geoh_org(:,:,:,:)
    integer,  intent( in)  :: dims(:)
    integer,  intent( in)  :: tcount
    integer,  intent( in)  :: step

    real(RP), allocatable :: hfact(:,:,:)
    real(RP), allocatable :: vfact(:,:,:,:,:)
    integer, allocatable :: igrd(:,:,:)
    integer, allocatable :: jgrd(:,:,:)
    integer, allocatable :: kgrd(:,:,:,:,:)

    integer :: n, k, i, j, iq

    intrinsic shape
    !---------------------------------------------------------------------------

    allocate( hfact(IA,JA,itp_nh) )
    allocate( vfact(KA,IA,JA,itp_nh,itp_nv) )
    allocate( igrd(IA,JA,itp_nh) )
    allocate( jgrd(IA,JA,itp_nh) )
    allocate( kgrd(KA,IA,JA,itp_nh,itp_nv) )

    ! for scalar points
    n = step
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat,lon, &
                                      geoh_org,lat_org,lon_org,dims(1),dims(2),dims(3),n )

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       dens(k,i,j) =  dens_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                    + dens_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                    + dens_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                    + dens_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                    + dens_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                    + dens_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       pott(k,i,j) =  pott_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                    + pott_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                    + pott_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                    + pott_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                    + pott_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                    + pott_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)

       do iq = 1, QA
          qtrc(k,i,j,iq) =  qtrc_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                          + qtrc_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                          + qtrc_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                          + qtrc_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n,iq) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                          + qtrc_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n,iq) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                          + qtrc_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n,iq) * hfact(i,j,3) * vfact(k,i,j,3,2)
       enddo
    enddo
    enddo
    enddo

    ! for vector (u) points
    n = step
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat,lon_u, &
                                      geoh_org,latu_org,lonu_org,dims(1),dims(2),dims(3),n )
    !                                 dims(5) wasn't used to keep consistency with geoh_org
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       u(k,i,j) =  u_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                 + u_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                 + u_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                 + u_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                 + u_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                 + u_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    ! for vector (v) points
    n = step
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,cz,lat_v,lon, &
                                      geoh_org,latv_org,lonv_org,dims(1),dims(2),dims(3),n )
    !                                 dims(6) wasn't used to keep consistency with geoh_org
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       v(k,i,j) =  v_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                 + v_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                 + v_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                 + v_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                 + v_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                 + v_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    ! for vector (w) points
    n = step
    call latlonz_interporation_fact( hfact,vfact,kgrd,igrd,jgrd,fz,lat,lon, &
                                      geof_org,lat_org,lon_org,dims(4),dims(2),dims(3),n )
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS-1, KE+1
       w(k,i,j) =  w_org(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,1) &
                 + w_org(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,1) &
                 + w_org(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,1) &
                 + w_org(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1),n) * hfact(i,j,1) * vfact(k,i,j,1,2) &
                 + w_org(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2),n) * hfact(i,j,2) * vfact(k,i,j,2,2) &
                 + w_org(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3),n) * hfact(i,j,3) * vfact(k,i,j,3,2)
    enddo
    enddo
    enddo

    deallocate( hfact )
    deallocate( vfact )
    deallocate( kgrd )
    deallocate( igrd )
    deallocate( jgrd )

    return
  end subroutine LatLonZ_Interpolation_Linear

  !-----------------------------------------------------------------------------
  subroutine latlonz_interporation_fact( &
      hfact,      & ! (out)
      vfact,      & ! (out)
      kgrd,       & ! (out)
      igrd,       & ! (out)
      jgrd,       & ! (out)
      myhgt,      & ! (in)
      mylat,      & ! (in)
      mylon,      & ! (in)
      inhgt,      & ! (in)
      inlat,      & ! (in)
      inlon,      & ! (in)
      nz,         & ! (in)
      nx,         & ! (in)
      ny,         & ! (in)
      step,       & ! (in)
      landgrid    & ! (in)
      )
    implicit none

    real(RP),           intent(out)  :: hfact(:,:,:)
    real(RP),           intent(out)  :: vfact(:,:,:,:,:)
    integer,           intent(out)  :: kgrd(:,:,:,:,:)
    integer,           intent(out)  :: igrd(:,:,:)
    integer,           intent(out)  :: jgrd(:,:,:)
    real(RP),           intent( in)  :: myhgt(:,:,:)
    real(RP),           intent( in)  :: mylat(:,:)
    real(RP),           intent( in)  :: mylon(:,:)
    real(RP),           intent( in)  :: inhgt(:,:,:,:)
    real(RP),           intent( in)  :: inlat(:,:,:)
    real(RP),           intent( in)  :: inlon(:,:,:)
    integer,            intent( in)  :: nz
    integer,            intent( in)  :: nx
    integer,            intent( in)  :: ny
    integer,            intent( in)  :: step
    logical, optional, intent( in)  :: landgrid

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(itp_nh)
    integer :: i, j, k, ii, jj, kk
    integer :: idx
    integer :: istart, iend, iinc, blk_i
    integer :: jstart, jend, jinc, blk_j
    integer :: kstart, kend
    logical :: lndgrd
    !---------------------------------------------------------------------------

    lndgrd = .false.
    if ( present(landgrid) .and. landgrid ) then
       lndgrd = .true.
    endif

    hfact(:,:,:) = 0.0_RP
    vfact(:,:,:,:,:) = 0.0_RP

    do j = JS-1, JE+1
    do i = IS-1, IE+1
       ! nearest block search
       iinc = (nx + 1) / interp_search_divnum
       jinc = (ny + 1) / interp_search_divnum
       dist(1) = large_number_one
       jj = 1 + (jinc/2)
       do while (jj <= ny)
          ii = 1 + (iinc/2)
          do while (ii <= nx)
             distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj,step),inlon(ii,jj,step) )

             if( distance < dist(1) )then
                dist(1) = distance
                blk_i = ii
                blk_j = jj
             endif
             ii = ii + iinc
          enddo
          jj = jj + jinc
       enddo
       istart = blk_i - (iinc/2) - 1
       if( istart < 1 ) istart = 1
       iend   = blk_i + (iinc/2) + 1
       if( iend  > nx ) iend   = nx
       jstart = blk_j - (jinc/2) - 1
       if( jstart < 1 ) jstart = 1
       jend   = blk_j + (jinc/2) + 1
       if( jend  > ny ) jend   = ny

       ! main search
       dist(1) = large_number_three
       dist(2) = large_number_two
       dist(3) = large_number_one
       do jj = jstart, jend
       do ii = istart, iend
          distance = haversine( mylat(i,j),mylon(i,j),inlat(ii,jj,step),inlon(ii,jj,step) )
          if ( distance <= dist(1) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = dist(1);     igrd(i,j,2) = igrd(i,j,1);  jgrd(i,j,2) = jgrd(i,j,1)
             dist(1) = distance;    igrd(i,j,1) = ii;           jgrd(i,j,1) = jj
          elseif ( dist(1) < distance .and. distance <= dist(2) ) then
             dist(3) = dist(2);     igrd(i,j,3) = igrd(i,j,2);  jgrd(i,j,3) = jgrd(i,j,2)
             dist(2) = distance;    igrd(i,j,2) = ii;           jgrd(i,j,2) = jj
          elseif ( dist(2) < distance .and. distance <= dist(3) ) then
             dist(3) = distance;    igrd(i,j,3) = ii;           jgrd(i,j,3) = jj
          endif
       enddo
       enddo
       if( dist(1)==0.0_RP )then
          hfact(i,j,1) = 1.0_RP
          hfact(i,j,2) = 0.0_RP
          hfact(i,j,3) = 0.0_RP
       else
          denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) + (1.0_RP/dist(3)) )
          hfact(i,j,1) = ( 1.0_RP/dist(1) ) * denom
          hfact(i,j,2) = ( 1.0_RP/dist(2) ) * denom
          hfact(i,j,3) = ( 1.0_RP/dist(3) ) * denom
       endif

       if ( lndgrd ) then
          kstart = 1;     kend = LKMAX
       else
          kstart = KS-1;  kend = KE+1
       endif
       do idx = 1, itp_nh
          ii = igrd(i,j,idx)
          jj = jgrd(i,j,idx)
          do k = kstart, kend
             dist(1) = large_number_two
             dist(2) = large_number_one
             do kk = 1, nz
                distance = abs( myhgt(k,i,j) - inhgt(kk,ii,jj,step) )
                if ( distance <= dist(1) ) then
                   dist(2) = dist(1);     kgrd(k,i,j,idx,2) = kgrd(k,i,j,idx,1)
                   dist(1) = distance;    kgrd(k,i,j,idx,1) = kk
                elseif ( dist(1) < distance .and. distance <= dist(2) ) then
                   dist(2) = distance;    kgrd(k,i,j,idx,2) = kk
                endif
             enddo
             if( dist(1)==0.0_RP )then
                vfact(k,i,j,idx,1) = 1.0_RP
                vfact(k,i,j,idx,2) = 0.0_RP
             else
                denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
                vfact(k,i,j,idx,1) = ( 1.0_RP/dist(1) ) * denom
                vfact(k,i,j,idx,2) = ( 1.0_RP/dist(2) ) * denom
             endif
          enddo
       enddo
    enddo
    enddo

    return
  end subroutine latlonz_interporation_fact

  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function haversine( &
      la0,       &
      lo0,       &
      la,        &
      lo )       &
      result( d )
    implicit none
    real(RP), intent(in) :: la0, lo0, la, lo   ! la,la0: Lat, lo,lo0: Lon; [rad]
    real(RP) :: d, dlon, dlat, work1, work2
    !---------------------------------------------------------------------------

    ! output unit : [m]
    dlon = lo0 - lo
    dlat = la0 - la
    work1 = (sin(dlat/2.0_RP))**2.0_RP + &
            cos(la0) * cos(la) * (sin(dlon/2.0_RP))**2.0_RP
    work2 = 2.0_RP * asin(min( 1.0_RP, sqrt(work1) ))
    d = r_in_m * work2

  end function haversine


end module mod_realinput
!-------------------------------------------------------------------------------
