!-------------------------------------------------------------------------------
!> module MKTOPO
!!
!! @par Description
!!          subroutines for preparing topography data
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_mktopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_MPIstop
  use scale_grid, only: &
     CX => GRID_CX, &
     CY => GRID_CY
  use scale_topography, only: &
     TOPO_Zsfc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKTOPO_setup
  public :: MKTOPO

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, save      :: MKTOPO_TYPE = -1
  integer, public, parameter :: I_IGNORE    =  0

  integer, public, parameter :: I_BELLSHAPE =  1
  integer, public, parameter :: I_GTOPO     =  2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: BELL_setup
  private :: MKTOPO_bellshape
  !private :: MKTOPO_gtopo

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: z_sfc(:,:,:)

  real(RP), private, allocatable :: bell (:,:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKTOPO_setup
    implicit none

    character(len=H_SHORT) :: MKTOPO_initname = 'OFF'

    NAMELIST / PARAM_MKTOPO / &
       MKTOPO_initname

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MKTOPO]/Categ[MKTOPO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO,iostat=ierr)

     allocate( z_sfc(1,IA,JA) )
     allocate( bell (1,IA,JA) )

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKTOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKTOPO)

    select case(trim(MKTOPO_initname))
    case('OFF')
       MKTOPO_TYPE = I_IGNORE
    case('BELLSHAPE')
       MKTOPO_TYPE = I_BELLSHAPE
       call BELL_setup
    !case('GTOPO')
    !   MKTOPO_TYPE = I_GTOPO
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(MKTOPO_initname)
       call PRC_MPIstop
    endselect

    return
  end subroutine MKTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKTOPO
    use scale_const, only: &
       CONST_UNDEF8
    use scale_topography, only: &
       TOPO_write
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( MKTOPO_TYPE == I_IGNORE ) then
      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  MAKING BOUNDARY DATA ++++++'
    else
      do j = 1, JA
      do i = 1, IA
         z_sfc(1,i,j) = CONST_UNDEF8
      enddo
      enddo

      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING BOUNDARY DATA ++++++'

      select case(MKTOPO_TYPE)
      case(I_BELLSHAPE)
         call MKTOPO_bellshape
      case(I_GTOPO)
         !call MKTOPO_gtopo
      case default
         write(*,*) ' xxx Unsupported TYPE:', MKTOPO_TYPE
         call PRC_MPIstop
      endselect

      if( IO_L ) write(IO_FID_LOG,*) '++++++ END   MAKING BOUNDARY DATA ++++++'

       ! output topography file
       call TOPO_write
    endif

    return
  end subroutine MKTOPO

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine BELL_setup
    implicit none

    ! Bubble
    logical  :: BELL_eachnode = .false.  ! Arrange bubble at each node? [kg/kg]
    real(RP) :: BELL_CX       =  2.E3_RP ! center location [m]: x
    real(RP) :: BELL_CY       =  2.E3_RP ! center location [m]: y
    real(RP) :: BELL_RX       =  2.E3_RP ! bubble radius   [m]: x
    real(RP) :: BELL_RY       =  2.E3_RP ! bubble radius   [m]: y

    NAMELIST / PARAM_BELL / &
       BELL_eachnode, &
       BELL_CX,       &
       BELL_CY,       &
       BELL_RX,       &
       BELL_RY

    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: dist

    integer  :: ierr
    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BELL]/Categ[MKTOPO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_BELL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_BELL)

    if ( BELL_eachnode ) then
       CX_offset = CX(IS)
       CY_offset = CY(JS)
    else
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    do j = JS, JE
    do i = IS, IE

       ! make tracer bubble
       dist = ( (CX(i)-CX_offset-BELL_CX)/BELL_RX )**2 &
            + ( (CY(j)-CY_offset-BELL_CY)/BELL_RY )**2

       bell(1,i,j) = 1.0_RP / ( 1.0_RP + dist)

    enddo
    enddo

    return
  end subroutine BELL_setup

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  !-----------------------------------------------------------------------------
  subroutine MKTOPO_bellshape
    implicit none

    ! bell-shaped mountain parameter
    real(RP) :: MOUNTAIN_HEIGHT =  100.0_RP ! height of mountain [m]

    NAMELIST / PARAM_MKTOPO_BELLSHAPE / &
       MOUNTAIN_HEIGHT

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BELLSHAPE]/Categ[MKTOPO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_BELLSHAPE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKTOPO_BELLSHAPE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKTOPO_BELLSHAPE)

    ! make bell-shaped mountain
    do j = JS, JE
    do i = IS, IE
       TOPO_Zsfc(i,j) = MOUNTAIN_HEIGHT * bell(1,i,j)
    enddo
    enddo

    return
  end subroutine MKTOPO_bellshape

end module mod_mktopo
