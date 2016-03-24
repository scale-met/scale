!-------------------------------------------------------------------------------
!> module EXTERNAL INPUT
!!
!! @par Description
!!          External file input module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_external_input
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: EXTIN_setup
  public :: EXTIN_update

  interface EXTIN_update
     module procedure EXTIN_update_1D
     module procedure EXTIN_update_2D
     module procedure EXTIN_update_3D
  end interface EXTIN_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: EXTIN_time_advance

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_prev = 1 !< [index] previous
  integer, private, parameter :: I_next = 2 !< [index] next

  integer, private, parameter :: I_periodic_year  = 1
  integer, private, parameter :: I_periodic_month = 2
  integer, private, parameter :: I_periodic_day   = 3

  integer, private, parameter :: EXTIN_item_limit = 1000  !< limit of item
  integer, private, parameter :: EXTIN_step_limit = 10000 !< limit of steps          for each item
  integer, private, parameter :: EXTIN_dim_limit  = 3     !< limit of dimension rank for each item

  type, private :: EXTIN_itemcontainer
     character(len=H_LONG)  :: basename                  !< base file name
     character(len=H_SHORT) :: varname                   !< variable name
     integer                :: dim_size(EXTIN_dim_limit) !< size of dimension (z,x,y)
     integer                :: step_num                  !< size of dimension (t)
     real(DP), allocatable  :: time(:)                   !< time of each step [sec]
     logical                :: fixed_step                !< fix step position?
     integer                :: flag_periodic             !< treat as periodic data? (0:no 1:yearly 2:monthly 3:daily)
     real(RP)               :: offset                    !< offset value (default is set to 0)
     integer                :: data_steppos(2)           !< current step position to read (1:previous,2:next)
     real(RP), allocatable  :: value(:,:,:,:)            !< data value                    (1:previous,2:next)
  end type EXTIN_itemcontainer


  integer,                   private              :: EXTIN_item_count = 0     !< number of item to output
  type(EXTIN_itemcontainer), private, allocatable :: EXTIN_item(:)            !< item to output

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine EXTIN_setup
    use gtool_file, only: &
       FileGetAllDataInfo, &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_calendar, only: &
       CALENDAR_adjust_daysec,  &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec, &
       CALENDAR_CFunits2sec,    &
       I_year,                  &
       I_month,                 &
       I_day
    use scale_time, only: &
       TIME_STARTDAYSEC, &
       TIME_NOWDAYSEC, &
       TIME_OFFSET_year
    implicit none

    character(len=H_LONG)  :: basename
    character(len=H_SHORT) :: varname
    integer                :: step_limit            ! limit number for reading data
    integer                :: step_fixed            ! fixed step position to read
    logical                :: enable_periodic_year  ! treat as yearly               periodic data?
    logical                :: enable_periodic_month ! treat as yearly,monthly       periodic data?
    logical                :: enable_periodic_day   ! treat as yearly,monthly,daily periodic data?
    real(RP)               :: offset
    real(RP)               :: defval

    namelist /EXTITEM/ &
       basename,              &
       varname,               &
       step_limit,            &
       step_fixed,            &
       enable_periodic_year,  &
       enable_periodic_month, &
       enable_periodic_day,   &
       offset,                &
       defval

    integer                :: step_nmax
    character(len=H_LONG)  :: description
    character(len=H_SHORT) :: unit
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (3)
    integer                :: dim_size  (3)
    real(DP)               :: time_start(EXTIN_step_limit)
    real(DP)               :: time_end  (EXTIN_step_limit)
    character(len=H_MID)   :: time_units

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year
    real(DP) :: cftime        !< time in the CF convention

    integer  :: count, nid, t, n
    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[EXTIN] / Categ[IO] / Origin[SCALElib]'

    ! count external data from namelist
    rewind(IO_FID_CONF)
    do count = 1, EXTIN_item_limit
       read(IO_FID_CONF,nml=EXTITEM,iostat=ierr)

       if( ierr < 0 ) then !--- no more items
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist EXTITEM. Check!', count
          call PRC_MPIstop
       endif
    enddo
    EXTIN_item_count = count - 1

    if( IO_L ) write(IO_FID_LOG,*)
    if ( EXTIN_item_count == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** No external file specified.'
       return
    endif
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of external item : ', EXTIN_item_count

    ! allocate item container
    allocate( EXTIN_item(EXTIN_item_count) )

    ! read detail of external data from namelist
    rewind(IO_FID_CONF)
    do nid = 1, EXTIN_item_count

       ! set default
       step_limit            = EXTIN_step_limit
       basename              = ''
       varname               = ''
       step_fixed            = -1
       enable_periodic_year  = .false.
       enable_periodic_month = .false.
       enable_periodic_day   = .false.
       offset                = 0.0_RP
       defval                = UNDEF

       ! read namelist
       read(IO_FID_CONF,nml=EXTITEM)

       ! read from file
       call FileGetAllDatainfo( step_limit,               & ! [IN]
                                EXTIN_dim_limit,          & ! [IN]
                                basename,                 & ! [IN]
                                varname,                  & ! [IN]
                                PRC_myrank,               & ! [IN]
                                step_nmax,                & ! [OUT]
                                description,              & ! [OUT]
                                unit,                     & ! [OUT]
                                datatype,                 & ! [OUT]
                                dim_rank,                 & ! [OUT]
                                dim_name  (:),            & ! [OUT]
                                dim_size  (:),            & ! [OUT]
                                time_start(1:step_limit), & ! [OUT]
                                time_end  (1:step_limit), & ! [OUT]
                                time_units                ) ! [OUT]

       if ( step_nmax == 0 ) then
          write(*,*) 'xxx Data not found! basename,varname = ', trim(basename), ', ', trim(varname)
          call PRC_MPIstop
       endif

       do n=dim_rank+1, 3
          dim_size(n) = 1
       end do

       ! setup item
       EXTIN_item(nid)%basename    = basename
       EXTIN_item(nid)%varname     = varname
       EXTIN_item(nid)%dim_size(:) = dim_size(:)
       EXTIN_item(nid)%step_num    = step_nmax

       if ( enable_periodic_day ) then
          EXTIN_item(nid)%flag_periodic = I_periodic_day
       elseif( enable_periodic_month ) then
          EXTIN_item(nid)%flag_periodic = I_periodic_month
       elseif( enable_periodic_year ) then
          EXTIN_item(nid)%flag_periodic = I_periodic_year
       else
          EXTIN_item(nid)%flag_periodic = 0
       endif

       allocate( EXTIN_item(nid)%value( dim_size(1),dim_size(2),dim_size(3),2) )
       EXTIN_item(nid)%value(:,:,:,:) = defval
       EXTIN_item(nid)%offset         = offset

       allocate( EXTIN_item(nid)%time(step_nmax) )
       EXTIN_item(nid)%time(:) = 0.0_DP

       do t = 1, EXTIN_item(nid)%step_num
          cftime = 0.5_DP * ( time_start(t) + time_end(t) )

          EXTIN_item(nid)%time(t) = CALENDAR_CFunits2sec( cftime, time_units, TIME_OFFSET_year, TIME_STARTDAYSEC )
       enddo

       if ( EXTIN_item(nid)%step_num == 1 ) then
          step_fixed = 1
       endif



       if ( step_fixed > 0 ) then ! fixed time step mode

          EXTIN_item(nid)%fixed_step           = .true.
          EXTIN_item(nid)%data_steppos(I_prev) = step_fixed
          EXTIN_item(nid)%data_steppos(I_next) = step_fixed

       else

          EXTIN_item(nid)%fixed_step = .false.

          ! seek start position
          EXTIN_item(nid)%data_steppos(I_next) = 1
          do t = 1, EXTIN_item(nid)%step_num
             if ( EXTIN_item(nid)%time(t) > TIME_NOWDAYSEC ) exit
             EXTIN_item(nid)%data_steppos(I_next) = t + 1
          enddo

          EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%data_steppos(I_next) - 1

          if ( EXTIN_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

             if ( EXTIN_item(nid)%data_steppos(I_next) == 1 ) then ! between first-1 and first

                ! first-1 = last
                EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%step_num

             elseif( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then ! between last and last+1

                ! last+1 = first
                EXTIN_item(nid)%data_steppos(I_next) = 1

                ! update data time in periodic condition
                do t = 1, EXTIN_item(nid)%step_num
                   dataday     = 0
                   datasec     = EXTIN_item(nid)%time(t)
                   offset_year = 0
                   call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                   call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              offset_year  ) ! [IN]

                   if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
                      datadate(I_day)   = datadate(I_day)   + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
                      datadate(I_month) = datadate(I_month) + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
                      datadate(I_year)  = datadate(I_year)  + 1
                   endif

                   call CALENDAR_date2daysec( dataday,     & ! [OUT]
                                              datasec,     & ! [OUT]
                                              datadate(:), & ! [IN]
                                              datasubsec,  & ! [IN]
                                              offset_year  ) ! [IN]

                   EXTIN_item(nid)%time(t) = CALENDAR_combine_daysec( dataday, datasec )
                enddo

                if( IO_L ) write(IO_FID_LOG,*) '*** data time is updated.'
             endif

          else ! normal mode

             if (      EXTIN_item(nid)%data_steppos(I_next) == 1                          &
                  .OR. EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then
                write(*,*) 'xxx Current time is out of period of external data! ', EXTIN_item(nid)%varname
                call PRC_MPIstop
             endif

          endif

       endif

       !--- read first data
       if (       dim_size(1) >= 1 &
            .AND. dim_size(2) == 1 &
            .AND. dim_size(3) == 1 ) then ! 1D

          ! read prev
          call FileRead( EXTIN_item(nid)%value(:,1,1,I_prev),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_prev), &
                         PRC_myrank                            )
          ! read next
          call FileRead( EXTIN_item(nid)%value(:,1,1,I_next),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_next), &
                         PRC_myrank                            )

       elseif (       dim_size(1) >= 1 &
                .AND. dim_size(2) >  1 &
                .AND. dim_size(3) == 1 ) then ! 2D

          ! read prev
          call FileRead( EXTIN_item(nid)%value(:,:,1,I_prev),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_prev), &
                         PRC_myrank                            )
          ! read next
          call FileRead( EXTIN_item(nid)%value(:,:,1,I_next),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_next), &
                         PRC_myrank                            )

       elseif (       dim_size(1) >= 1 &
                .AND. dim_size(2) >  1 &
                .AND. dim_size(3) >  1 ) then ! 3D

          ! read prev
          call FileRead( EXTIN_item(nid)%value(:,:,:,I_prev),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_prev), &
                         PRC_myrank                            )
          ! read next
          call FileRead( EXTIN_item(nid)%value(:,:,:,I_next),  &
                         EXTIN_item(nid)%basename,             &
                         EXTIN_item(nid)%varname,              &
                         EXTIN_item(nid)%data_steppos(I_next), &
                         PRC_myrank                            )

       else
          write(*,*) 'xxx Unexpected dimsize: ', dim_size(:)
          call PRC_MPIstop
       endif

    enddo

    return
  end subroutine EXTIN_setup

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_1D( &
       var,          &
       varname,      &
       axistype,     &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:)       ! variable
    character(len=*), intent(in)  :: varname      ! item name
    character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: dim1_max, dim1_S, dim1_E, n1, nn1
    integer  :: n
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if( nid == 0 ) return

    if ( axistype == 'Z' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'X' ) then
       dim1_max = IMAXB
       dim1_S   = ISB
       dim1_E   = IEB
    elseif( axistype == 'Y' ) then
       dim1_max = JMAXB
       dim1_S   = JSB
       dim1_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( dim1_max /= EXTIN_item(nid)%dim_size(1) ) then
       write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
       write(*,*) 'xxx dim 1 (data,requested)    : ', EXTIN_item(nid)%dim_size(1), dim1_max
       call PRC_MPIstop
    endif

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var: ', trim(EXTIN_item(nid)%varname)

       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,1,1,I_next),  & ! [OUT]
                      EXTIN_item(nid)%basename,             & ! [IN]
                      EXTIN_item(nid)%varname,              & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next), & ! [IN]
                      PRC_myrank                            ) ! [IN]
    endif

    ! store data with weight
    do n1 = 1, dim1_max
       nn1 = n1 + dim1_S - 1

       var(nn1) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,1,1,I_prev) &
                + (        weight ) * EXTIN_item(nid)%value(n1,1,1,I_next)
    enddo

    error = .false.

    return
  end subroutine EXTIN_update_1D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_2D( &
       var,          &
       varname,      &
       axistype,     &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:)     ! variable
    character(len=*), intent(in)  :: varname      ! item name
    character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: dim1_max, dim1_S, dim1_E, n1, nn1
    integer  :: dim2_max, dim2_S, dim2_E, n2, nn2
    integer  :: n
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if( nid == 0 ) return

    if ( axistype == 'XY' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif( axistype == 'ZX' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if (      dim1_max /= EXTIN_item(nid)%dim_size(1) &
         .OR. dim2_max /= EXTIN_item(nid)%dim_size(2) ) then
       write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
       write(*,*) 'xxx dim 1 (data,requested)    : ', EXTIN_item(nid)%dim_size(1), dim1_max
       write(*,*) 'xxx dim 2 (data,requested)    : ', EXTIN_item(nid)%dim_size(2), dim2_max
       call PRC_MPIstop
    endif

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(EXTIN_item(nid)%varname)
       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,1,I_next),  & ! [OUT]
                      EXTIN_item(nid)%basename,             & ! [IN]
                      EXTIN_item(nid)%varname,              & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next), & ! [IN]
                      PRC_myrank                            ) ! [IN]
    endif

    ! store data with weight
    do n2 = 1, dim2_max
       nn2 = n2 + dim2_S - 1
    do n1 = 1, dim1_max
       nn1 = n1 + dim1_S - 1

       var(nn1,nn2) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,1,I_prev) &
                    + (        weight ) * EXTIN_item(nid)%value(n1,n2,1,I_next)
    enddo
    enddo

    error = .false.

    return
  end subroutine EXTIN_update_2D

  !-----------------------------------------------------------------------------
  !> Read data
  subroutine EXTIN_update_3D( &
       var,          &
       varname,      &
       axistype,     &
       current_time, &
       error         )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:)   ! variable
    character(len=*), intent(in)  :: varname      ! item name
    character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)
    real(DP),         intent(in)  :: current_time ! current time
    logical,          intent(out) :: error        ! error code

    integer  :: nid
    real(RP) :: weight
    logical  :: do_readfile

    integer  :: dim1_max, dim1_S, dim1_E, n1, nn1
    integer  :: dim2_max, dim2_S, dim2_E, n2, nn2
    integer  :: dim3_max, dim3_S, dim3_E, n3, nn3
    integer  :: n
    !---------------------------------------------------------------------------

    error = .true.

    ! searching the data ID
    nid = -1
    do n = 1, EXTIN_item_count
       if( varname == EXTIN_item(n)%varname ) nid = n
    enddo

    if( nid == 0 ) return

    if ( axistype == 'ZXY' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else if ( axistype == 'Land' ) then
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else if ( axistype == 'Urban' ) then
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if (      dim1_max /= EXTIN_item(nid)%dim_size(1) &
         .OR. dim2_max /= EXTIN_item(nid)%dim_size(2) &
         .OR. dim3_max /= EXTIN_item(nid)%dim_size(3) ) then
       write(*,*) 'xxx data length does not match! ', trim(axistype), ' item:', trim(varname)
       write(*,*) 'xxx dim 1 (data,requested)    : ', EXTIN_item(nid)%dim_size(1), dim1_max
       write(*,*) 'xxx dim 2 (data,requested)    : ', EXTIN_item(nid)%dim_size(2), dim2_max
       write(*,*) 'xxx dim 3 (data,requested)    : ', EXTIN_item(nid)%dim_size(3), dim3_max
       call PRC_MPIstop
    endif

    call EXTIN_time_advance( nid,          & ! [IN]
                             current_time, & ! [IN]
                             weight,       & ! [OUT]
                             do_readfile   ) ! [OUT]

    if ( do_readfile ) then
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(EXTIN_item(nid)%varname)
       ! next -> prev
       EXTIN_item(nid)%value(:,:,:,I_prev) = EXTIN_item(nid)%value(:,:,:,I_next)

       ! read next
       call FileRead( EXTIN_item(nid)%value(:,:,:,I_next),  & ! [OUT]
                      EXTIN_item(nid)%basename,             & ! [IN]
                      EXTIN_item(nid)%varname,              & ! [IN]
                      EXTIN_item(nid)%data_steppos(I_next), & ! [IN]
                      PRC_myrank                            ) ! [IN]
    endif

    ! store data with weight
    do n3 = 1, dim3_max
       nn3 = n3 + dim3_S - 1
    do n2 = 1, dim2_max
       nn2 = n2 + dim2_S - 1
    do n1 = 1, dim1_max
       nn1 = n1 + dim1_S - 1

       var(nn1,nn2,nn3) = ( 1.0_RP-weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_prev) &
                        + (        weight ) * EXTIN_item(nid)%value(n1,n2,n3,I_next)
    enddo
    enddo
    enddo

    error = .false.

    return
  end subroutine EXTIN_update_3D

  !-----------------------------------------------------------------------------
  !> Check time to read and calc weight
  subroutine EXTIN_time_advance( &
       nid,          &
       current_time, &
       weight,       &
       do_readfile   )
    use scale_calendar, only: &
       CALENDAR_adjust_daysec,  &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec, &
       I_year,                  &
       I_month,                 &
       I_day
    implicit none

    integer,  intent(in)  :: nid          ! item id
    real(DP), intent(in)  :: current_time ! current time
    real(RP), intent(out) :: weight       ! weight
    logical,  intent(out) :: do_readfile  ! read new data at this time?

    integer  :: datadate(6)   !< date
    real(DP) :: datasubsec    !< subsecond
    integer  :: dataday       !< absolute day
    real(DP) :: datasec       !< absolute second
    integer  :: offset_year   !< offset year

    real(DP) :: prev_time, next_time
    integer  :: t
    !---------------------------------------------------------------------------

    do_readfile = .false.

    if ( EXTIN_item(nid)%fixed_step ) then
       !--- no time-advance
    else
       ! time is passed?
       if ( current_time > EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) ) ) then

          do_readfile = .true.

          if( IO_L ) write(IO_FID_LOG,*) '*** Update external input : ', trim(EXTIN_item(nid)%varname)

          ! update step position
          EXTIN_item(nid)%data_steppos(I_prev) = EXTIN_item(nid)%data_steppos(I_prev) + 1
          EXTIN_item(nid)%data_steppos(I_next) = EXTIN_item(nid)%data_steppos(I_next) + 1

          if ( EXTIN_item(nid)%flag_periodic > 0 ) then ! periodic time step mode

             if ( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then

                ! last+1 = first
                EXTIN_item(nid)%data_steppos(I_next) = 1

                ! update data time in periodic condition
                do t = 1, EXTIN_item(nid)%step_num
                   dataday     = 0
                   datasec     = EXTIN_item(nid)%time(t)
                   offset_year = 0
                   call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

                   call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              offset_year  ) ! [IN]

                   if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
                      datadate(I_day)   = datadate(I_day)   + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
                      datadate(I_month) = datadate(I_month) + 1
                   elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
                      datadate(I_year)  = datadate(I_year)  + 1
                   endif

                   call CALENDAR_date2daysec( dataday,     & ! [IN]
                                              datasec,     & ! [IN]
                                              datadate(:), & ! [OUT]
                                              datasubsec,  & ! [OUT]
                                              offset_year  ) ! [IN]

                   EXTIN_item(nid)%time(t) = CALENDAR_combine_daysec( dataday, datasec )
                enddo
             endif

          else ! normal mode

             if ( EXTIN_item(nid)%data_steppos(I_next) == EXTIN_item(nid)%step_num+1 ) then
                write(*,*) 'xxx Current time is out of period of external data! '
                call PRC_MPIstop
             endif

          endif

       endif

    endif

    ! calc weight
    if ( EXTIN_item(nid)%fixed_step ) then

       weight = 0.0_RP

    elseif( EXTIN_item(nid)%data_steppos(I_next) == 1 ) then !<--- this case is only periodic one.

       dataday     = 0
       datasec     = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_prev) )
       offset_year = 0
       call CALENDAR_adjust_daysec( dataday, datasec ) ! [INOUT]

       call CALENDAR_daysec2date( datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  offset_year  ) ! [IN]

       if    ( EXTIN_item(nid)%flag_periodic == I_periodic_day   ) then
          datadate(I_day)   = datadate(I_day)   - 1
       elseif( EXTIN_item(nid)%flag_periodic == I_periodic_month ) then
          datadate(I_month) = datadate(I_month) - 1
       elseif( EXTIN_item(nid)%flag_periodic == I_periodic_year  ) then
          datadate(I_year)  = datadate(I_year)  - 1
       endif

       call CALENDAR_date2daysec( dataday,     & ! [IN]
                                  datasec,     & ! [IN]
                                  datadate(:), & ! [OUT]
                                  datasubsec,  & ! [OUT]
                                  offset_year  ) ! [IN]

       prev_time = CALENDAR_combine_daysec( dataday, datasec )

       next_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) )

       weight = ( current_time - prev_time ) &
              / ( next_time    - prev_time )

    else

       prev_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_prev) )
       next_time = EXTIN_item(nid)%time( EXTIN_item(nid)%data_steppos(I_next) )

       weight = ( current_time - prev_time ) &
              / ( next_time    - prev_time )

    endif

    return
  end subroutine EXTIN_time_advance

end module scale_external_input
!-------------------------------------------------------------------------------
