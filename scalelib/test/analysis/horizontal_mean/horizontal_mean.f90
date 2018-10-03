program horizontal_mean
  use scale
  use scale_atmos_grid_cartesC_index, only: &
       KMAX, IMAX, JMAX, &
       KA, KS, KE
  use scale_prc_cartesC, only: &
       PRC_CARTESC_setup
  use scale_statistics, only: &
       STATISTICS_setup, &
       STATISTICS_horizontal_mean
  use netcdf
  implicit none

  real(RP), allocatable :: U(:,:,:), V(:,:,:)
  real(RP), allocatable :: U_mean(:), V_mean(:)
  real(RP), allocatable :: area(:,:)

  ! conf data
  character(len=H_MID) :: basename_in
  character(len=H_MID) :: basename_out

  NAMELIST / PARAM_HORIZONTAL_MEAN / &
       basename_in, &
       basename_out

  integer :: ierr


  call SCALE_init

  call PRC_CARTESC_setup

  ! default value
  basename_in  = "history"
  basename_out = "history_hmean"

  !--- read namelist
  if ( IO_FID_CONF > 0 ) then
     rewind(IO_FID_CONF)
     read(IO_FID_CONF,nml=PARAM_HORIZONTAL_MEAN,iostat=ierr)
     if( ierr < 0 ) then !--- missing
        if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
     elseif( ierr > 0 ) then !--- fatal error
        write(*,*) 'xxx Not appropriate names in namelist PARAM_HORIZONTAL_MEAN. Check!'
        call PRC_abort
     endif
  end if
  if( IO_NML ) write(IO_FID_NML,nml=PARAM_HORIZONTAL_MEAN)
  

  call read_data

  call STATISTICS_setup
  call STATISTICS_horizontal_mean( KMAX, 1, KMAX, IMAX, 1, IMAX, JMAX, 1, JMAX, &
                                   U(:,:,:), area(:,:), & ! (in)
                                   U_mean(:)            ) ! (out)
  call STATISTICS_horizontal_mean( KMAX, 1, KMAX, IMAX, 1, IMAX, JMAX, 1, JMAX, &
                                   V(:,:,:), area(:,:), & ! (in)
                                   V_mean(:)            ) ! (out)

  call write_data

  call SCALE_finalize

  stop

contains

  subroutine read_data
    use scale_atmos_grid_cartesC_index, only: &
         ATMOS_GRID_CARTESC_INDEX_setup
    use scale_file, only: &
         FILE_get_dimLength
    use scale_file_cartesC, only: &
         FILE_CARTESC_open, &
         FILE_CARTESC_close, &
         FILE_CARTESC_setup, &
         FILE_CARTESC_get_size, &
         FILE_CARTESC_read
    use scale_atmos_grid_cartesC, only: &
         ATMOS_GRID_CARTESC_setup

    integer :: KMAX, OKMAX, LKMAX, UKMAX
    integer :: IMAXG, JMAXG
    integer :: KHALO, IHALO, JHALO
    integer :: nt
    integer :: fid

    ! get grid size
    call FILE_CARTESC_get_size( basename_in,                      & ! (in)
                                KMAX, OKMAX, LKMAX, UKMAX,        & ! (out)
                                IMAXG, JMAXG, KHALO, IHALO, JHALO ) ! (out)

    call ATMOS_GRID_CARTESC_INDEX_setup( KMAX=KMAX, IMAXG=IMAXG, JMAXG=JMAXG, &
                                         KHALO=KHALO, IHALO=IHALO, JHALO=JHALO ) ! (in)

    ! get axis info (used in writing file)
    call ATMOS_GRID_CARTESC_setup( basename = basename_in )
    call FILE_CARTESC_setup

    ! allocate
    allocate( U(KMAX,IMAX,JMAX), V(KMAX,IMAX,JMAX) )
    allocate( area(IMAX,JMAX) )
    allocate( U_mean(KMAX), V_mean(KMAX) )

    ! get data
    call FILE_CARTESC_open( basename_in, & ! (in)
                            fid          ) ! (out)

    call FILE_get_dimLength( fid, "time", & ! (in)
                             nt           ) ! (out)

    call FILE_CARTESC_read( fid, "U", & ! (in)
                            U(:,:,:), & ! (out)
                            step = nt ) ! (in)

    call FILE_CARTESC_read( fid, "V", & ! (in)
                            V(:,:,:), & ! (out)
                            step = nt ) ! (in)

    call FILE_CARTESC_read( fid, "cell_area", & ! (in)
                            area(:,:)         ) ! (out)

    call FILE_CARTESC_close( fid )

    return
  end subroutine read_data

  subroutine write_data
    use scale_file_cartesC, only: &
         FILE_CARTESC_write
    character(len=*), parameter :: title = "Horizontal Mean"
    real(RP) :: buf(KA)

    buf(KS:KE) = U_mean(:)
    call FILE_CARTESC_write( buf(:), basename_out, title, "U_mean", "horizontal mean U", "m/s", "Z", "REAL4" )
    buf(KS:KE) = V_mean(:)
    call FILE_CARTESC_write( buf(:), basename_out, title, "V_mean", "horizontal mean V", "m/s", "Z", "REAL4", append=.true. )

    return
  end subroutine write_data

end program horizontal_mean
