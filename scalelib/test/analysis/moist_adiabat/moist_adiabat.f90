program moist_adiabat
  use scale
  use scale_const, only: &
     PRE00 => CONST_PRE00, &
     Rdry  => CONST_Rdry,  &
     Rvap  => CONST_Rvap,  &
     CPdry => CONST_CPdry
  use scale_atmos_hydrometeor, only: &
     CPvap => CP_VAPOR, &
     CL    => CP_WATER, &
     CI    => CP_ICE
  use scale_stdio, only: &
     IO_get_available_fid
  use scale_atmos_saturation, only: &
     ATMOS_SATURATION_setup, &
     ATMOS_SATURATION_tdew_liq, &
     ATMOS_SATURATION_pote
  use scale_atmos_adiabat, only: &
     ATMOS_ADIABAT_setup, &
     ATMOS_ADIABAT_cape
  use netcdf
  implicit none

  integer :: nx, ny, nz
  real(RP), allocatable :: xaxis(:), yaxis(:)
  real(RP), allocatable :: z(:), zh(:)

  real(RP), allocatable :: TEMP(:), PRES(:)
  real(RP), allocatable :: QV(:), QC(:), QR(:), QI(:), QS(:), QG(:)
  real(RP), allocatable :: POTT(:), DENS(:), QDRY(:)
  real(RP), allocatable :: Rtot(:), CPtot(:)
  real(RP), allocatable :: Tdew(:), POTE(:)
  real(RP), allocatable :: DENS_p(:), TEMP_p(:), BUOY_p(:), QV_p(:)
  real(RP) :: CAPE, CIN, LCL, LFC, LNB

  logical :: converged

  ! conf data
  character(len=H_MID) :: basename_in
  integer :: rankid
  character(len=H_MID) :: filename_out
  real(RP) :: x, y
  integer  :: nstep

  NAMELIST / PARAM_MOIST_ADIABAT / &
       basename_in, &
       rankid, &
       filename_out, &
       x, y, &
       nstep

  integer :: fid
  integer :: k
  integer :: ierr


  call SCALE_init

  call ATMOS_SATURATION_setup
  call ATMOS_ADIABAT_setup

  ! default value
  basename_in  = "history"
  filename_out = "output.dat"
  x = 0.0_RP
  y = 0.0_RP
  nstep = 1

  !--- read namelist
  if ( IO_FID_CONF > 0 ) then
     rewind(IO_FID_CONF)
     read(IO_FID_CONF,nml=PARAM_MOIST_ADIABAT,iostat=ierr)
     if( ierr < 0 ) then !--- missing
        if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
     elseif( ierr > 0 ) then !--- fatal error
        write(*,*) 'xxx Not appropriate names in namelist PARAM_MOIST_ADIABAT. Check!'
        call PRC_abort
     endif
  end if
  if( IO_NML ) write(IO_FID_NML,nml=PARAM_MOIST_ADIABAT)
  

  call read_data


  ! approximate estimation (ignore other tracers)
  allocate( POTT(nz), DENS(nz), QDRY(nz), Rtot(nz), CPtot(nz) )
  do k = 1, nz
     QDRY(k) = 1.0_RP - QV(k) - QC(k) - QR(k) - QI(k) - QS(k) - QG(k)
     Rtot(k)  =  Rdry * QDRY(k) +  Rvap * QV(k)
     CPtot(k) = CPdry * QDRY(k) + CPvap * QV(k) + CL * ( QC(k) + QR(k) ) + CI * ( QI(k) + QS(k) + QG(k) )
     POTT(k) = TEMP(k) * ( PRE00 / PRES(k) )**( Rtot(k) / CPtot(k) )
     DENS(k) = PRES(k) / ( Rtot(k) * TEMP(k) )
  end do

  ! dew point
  allocate( Tdew(nz) )
  call ATMOS_SATURATION_tdew_liq( nz, 1, nz, &
                                  DENS(:), TEMP(:), QV(:), & ! [IN]
                                  Tdew(:), converged       ) ! [OUT]

  ! equivalent potential temperature
  allocate( POTE(nz) )
  call ATMOS_SATURATION_pote( nz, 1, nz, &
                              DENS(:), POTT(:), TEMP(:), QV(:), & ! [IN]
                              POTE(:)                           ) ! [OUT]

  ! lift parcel
  allocate( DENS_p(nz), TEMP_p(nz), BUOY_p(nz), QV_p(nz) )
  call ATMOS_ADIABAT_cape( nz, 1, nz, 1, &
                           DENS(:), TEMP(:), PRES(:),                & ! [IN]
                           QV(:), QC(:), QDRY(:), Rtot(:), CPtot(:), & ! [IN]
                           z(:), zh(:),                              & ! [IN]
                           CAPE, CIN, LCL, LFC, LNB,                 & ! [OUT]
                           DENS_p(:), TEMP_p(:), BUOY_p(:), QV_p(:), & ! [OUT]
                           converged                                 ) ! [OUT]
  if ( .not. converged ) then
     write(*,*) 'xxx calculation was not converged!'
     call PRC_abort
  end if
  !write(*,*) buoy_p
  !write(*,*) zh


  write(fid,'(5(a,f10.2))') 'CAPE=', CAPE, ', CIN=', CIN, ', LCL=', LCL, ', LFC=', LFC, ', LNB=', LNB


  fid = IO_get_available_fid()
  open( unit=fid, file=filename_out )
  write(fid,*)nz
  write(fid,'(2a8,5a10)') "z", "pres", "temp", "tdew", "pott", "pote", "temp_p"
  do k = 1, nz
     write(fid,'(2f8.1,5f10.2)') z(k), pres(k)/100.0_RP, temp(k), tdew(k), pott(k), pote(k), temp_p(k)
  end do
  close(fid)

  !write(*,*)temp_p(nz)*(PRE00/pres(nz))**(Rdry/CPdry)


  call SCALE_finalize

  stop

contains

  subroutine read_data
    use scale_file, only: &
         FILE_open, &
         FILE_get_datainfo, &
         FILE_read, &
         FILE_close
    integer :: fid
    integer :: dims(3)
    integer :: start(3), count(3)

    integer :: i, j

    ! file open
    call FILE_open( fid, & ! (out)
                    basename_in, rankid=rankid ) ! (in)

    ! get dimension size
    call FILE_get_datainfo( fid, "height", dim_size=dims(:) )
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    ! get x and y axis
    allocate( xaxis(nx), yaxis(ny) )
    call FILE_read( xaxis(:), fid, "x" )
    call FILE_read( yaxis(:), fid, "y" )

    ! search index
    do i = 1, nx
       if ( x <= xaxis(i) ) exit
    end do
    do j = 1, ny
       if ( y <= yaxis(j) ) exit
    end do
    start(:) = (/i,j,1/)

    ! get z axis
    allocate( z(nz), zh(0:nz) )
    call FILE_read( z(:),  fid, "height",     start=start(:), count=(/1,1,nz/) )
    call FILE_read( zh(:), fid, "height_xyw", start=start(:), count=(/1,1,nz+1/) )

    ! get variables
    allocate( TEMP(nz), PRES(nz), QV(nz), QC(nz), QR(nz), QI(nz), QS(nz), QG(nz) )
    count(:) = (/1,1,nz/)

    call FILE_read( TEMP(:), fid, "T",    start=start(:), count=count(:), step=nstep )
    call FILE_read( PRES(:), fid, "PRES", start=start(:), count=count(:), step=nstep )
    call FILE_read( QV  (:), fid, "QV",   start=start(:), count=count(:), step=nstep )
    call FILE_read( QC  (:), fid, "QC",   start=start(:), count=count(:), step=nstep, allow_missing=.true. )
    call FILE_read( QR  (:), fid, "QR",   start=start(:), count=count(:), step=nstep, allow_missing=.true. )
    call FILE_read( QI  (:), fid, "QI",   start=start(:), count=count(:), step=nstep, allow_missing=.true. )
    call FILE_read( QS  (:), fid, "QS",   start=start(:), count=count(:), step=nstep, allow_missing=.true. )
    call FILE_read( QG  (:), fid, "QG",   start=start(:), count=count(:), step=nstep, allow_missing=.true. )


    ! close
    call FILE_close( fid )

    return
  end subroutine read_data

end program moist_adiabat
