program moist_adiabat
  use scale
  use scale_const, only: &
     PRE00 => CONST_PRE00, &
     Rdry  => CONST_Rdry,  &
     Rvap  => CONST_Rvap,  &
     CPdry => CONST_CPdry, &
     CPvap => CONST_CPvap, &
     CL    => CONST_CL
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

  real(RP), allocatable :: TEMP(:), PRES(:), QV(:), QC(:)
  real(RP), allocatable :: POTT(:), DENS(:), QDRY(:)
  real(RP), allocatable :: Rtot(:), CPtot(:)
  real(RP), allocatable :: Tdew(:), POTE(:)
  real(RP), allocatable :: DENS_p(:), TEMP_p(:), BUOY_p(:), QV_p(:)
  real(RP) :: CAPE, CIN, LCL, LFC, LNB

  logical :: converged

  ! netcdf
  integer :: ncid, vid
  integer :: dimid(1)
  integer :: status
  integer :: start(4), count(4)

  ! conf data
  character(len=H_MID) :: filename_in
  character(len=H_MID) :: filename_out
  real(RP) :: x, y
  integer  :: nstep

  NAMELIST / PARAM_MOIST_ADIABAT / &
       filename_in, &
       filename_out, &
       x, y, &
       nstep

  integer :: fid
  integer :: i, j, k
  integer :: ierr


  call SCALE_init

  call ATMOS_SATURATION_setup
  call ATMOS_ADIABAT_setup

  filename_in  = "history.pe000000.nc"
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
  
  status = nf90_open( filename_in, nf90_nowrite, ncid )
  if ( status /= nf90_noerr) then
     write(*,*) 'xxx failed to open file: ', trim(filename_in)
     call PRC_abort
  end if

  ! get x axis
  status = nf90_inq_varid( ncid, "x", vid )
  status = nf90_inquire_variable( ncid, vid, dimids=dimid(:) )
  status = nf90_inquire_dimension( ncid, dimid(1), len=nx )
  allocate( xaxis(nx) )
  status = nf90_get_var( ncid, vid, xaxis(:) )

  ! get y axis
  status = nf90_inq_varid( ncid, "y", vid )
  status = nf90_inquire_variable( ncid, vid, dimids=dimid(:) )
  status = nf90_inquire_dimension( ncid, dimid(1), len=ny )
  allocate( yaxis(ny) )
  status = nf90_get_var( ncid, vid, yaxis(:) )

  ! get z axis length
  status = nf90_inq_varid( ncid, "z", vid )
  status = nf90_inquire_variable( ncid, vid, dimids=dimid(:) )
  status = nf90_inquire_dimension( ncid, dimid(1), len=nz )


  ! search index
  do i = 1, nx
     if ( x <= xaxis(i) ) exit
  end do
  do j = 1, ny
     if ( y <= yaxis(j) ) exit
  end do

  ! get z axis
  allocate( z(nz), zh(0:nz) )
  status = nf90_inq_varid( ncid, "height", vid )
  status = nf90_get_var( ncid, vid, z(:), start=(/i,j,1/), count=(/1,1,nz/) )
  status = nf90_inq_varid( ncid, "height_xyw", vid )
  status = nf90_get_var( ncid, vid, zh(:), start=(/i,j,1/), count=(/1,1,nz+1/) )

  ! get variables
  allocate( TEMP(nz), PRES(nz), QV(nz), QC(nz) )
  start(:) = (/i,j,1,nstep/)
  count(:) = (/1,1,nz,1/)

  status = nf90_inq_varid( ncid, "T", vid )
  status = nf90_get_var( ncid, vid, TEMP(:), start=start, count=count)

  status = nf90_inq_varid( ncid, "PRES", vid )
  status = nf90_get_var( ncid, vid, PRES(:), start=start, count=count)

  status = nf90_inq_varid( ncid, "QV", vid )
  status = nf90_get_var( ncid, vid, QV(:), start=start, count=count)

  status = nf90_inq_varid( ncid, "QC", vid )
  status = nf90_get_var( ncid, vid, QC(:), start=start, count=count)


  status = nf90_close( ncid )



  ! approximate estimation (ignore other tracers than qv and qc)
  allocate( POTT(nz), DENS(nz), QDRY(nz), Rtot(nz), CPtot(nz) )
  do k = 1, nz
     QDRY(k) = 1.0_RP - QV(k) - QC(k)
     Rtot(k)  =  Rdry * QDRY(k) +  Rvap * QV(k)
     CPtot(k) = CPdry * QDRY(k) + CPvap * QV(k) + CL * QC(k)
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
  write(fid,'(2a8,5a9)') "z", "pres", "temp", "tdew", "pott", "pote", "temp_p"
  do k = 1, nz
     write(fid,'(2f8.1,5f9.2)') z(k), pres(k)/100.0_RP, temp(k), tdew(k), pott(k), pote(k), temp_p(k)
  end do
  close(fid)

  !write(*,*)temp_p(nz)*(PRE00/pres(nz))**(Rdry/CPdry)


  call SCALE_finalize

  stop
end program moist_adiabat
