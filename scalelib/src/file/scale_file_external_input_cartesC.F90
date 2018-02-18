!-------------------------------------------------------------------------------
!> module file / external_input_cartesC
!!
!! @par Description
!!          External file input module for the cartesian-C grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_file_external_input_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_ocean_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FILE_EXTERNAL_INPUT_CARTESC_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  private :: FILE_EXTERNAL_INPUT_CARTESC_get_dims1D
  private :: FILE_EXTERNAL_INPUT_CARTESC_get_dims2D
  private :: FILE_EXTERNAL_INPUT_CARTESC_get_dims3D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine FILE_EXTERNAL_INPUT_CARTESC_setup
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_setup, &
       FILE_EXTERNAL_INPUT_get_dims1D, &
       FILE_EXTERNAL_INPUT_get_dims2D, &
       FILE_EXTERNAL_INPUT_get_dims3D

    call FILE_EXTERNAL_INPUT_setup

    FILE_EXTERNAL_INPUT_get_dims1D => FILE_EXTERNAL_INPUT_CARTESC_get_dims1D
    FILE_EXTERNAL_INPUT_get_dims2D => FILE_EXTERNAL_INPUT_CARTESC_get_dims2D
    FILE_EXTERNAL_INPUT_get_dims3D => FILE_EXTERNAL_INPUT_CARTESC_get_dims3D

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> get_dims
  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims1D( &
       dim1_max, &
       dim1_S,   &
       dim1_E,   &
       varname,  &
       axistype  )
    use scale_process, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim1_E
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)

    select case ( axistype )
    case ( 'Z' )
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    case ( 'X' )
       dim1_max = IMAXB
       dim1_S   = ISB
       dim1_E   = IEB
    case ( 'Y' )
       dim1_max = JMAXB
       dim1_S   = JSB
       dim1_E   = JEB
    case ( 'OZ' )
       dim1_max = OKMAX
       dim1_S   = OKS
       dim1_E   = OKE
    case default
       write(*,*) 'xxx [FILE_EXTERNAL_INPUT_CARTESC_get_dims1D] unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims1D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims2D( &
       dim1_max,  &
       dim1_S,    &
       dim1_E,    &
       dim2_max,  &
       dim2_S,    &
       dim2_E,    &
       transpose, &
       varname,   &
       axistype   )
    use scale_process, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim1_E
    integer,          intent(out) :: dim2_max
    integer,          intent(out) :: dim2_S
    integer,          intent(out) :: dim2_E
    logical,          intent(out) :: transpose
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (XY/XZ/ZX)

    select case ( axistype )
    case ( 'XY' )
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       transpose = .false.
    case ( 'ZX' )
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       transpose = .false.
    case ( 'XZ' )
       dim1_max = IMAXB
       dim2_max = KMAX
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = KS
       dim2_E   = KE
       transpose = .true.
    case default
       write(*,*) 'xxx [FILE_EXTERNAL_INPUT_CARTESC_get_dims2D] unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims2D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims3D( &
       dim1_max,  &
       dim1_S,    &
       dim1_E,    &
       dim2_max,  &
       dim2_S,    &
       dim2_E,    &
       dim3_max,  &
       dim3_S,    &
       dim3_E,    &
       transpose, &
       varname,   &
       axistype   )
    use scale_process, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim1_E
    integer,          intent(out) :: dim2_max
    integer,          intent(out) :: dim2_S
    integer,          intent(out) :: dim2_E
    integer,          intent(out) :: dim3_max
    integer,          intent(out) :: dim3_S
    integer,          intent(out) :: dim3_E
    logical,          intent(out) :: transpose
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (ZXY/XYZ/LXY/XYL/UXY/XYU)

    select case ( axistype )
    case ( 'ZXY' )
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       transpose = .false.
    case ( 'XYZ' )
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim3_max = KMAX
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       dim3_S   = KS
       dim3_E   = KE
       transpose = .true.
    case ( 'OXY' ) ! OCEAN
       dim1_max = OKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = OKS
       dim1_E   = OKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       transpose = .false.
    case ( 'XYO' ) ! OCEAN
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim3_max = OKMAX
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       dim3_S   = OKS
       dim3_E   = OKE
       transpose = .true.
    case ( 'LXY' ) ! LAND
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       transpose = .false.
    case ( 'XYL' ) ! LAND
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim3_max = LKMAX
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       dim3_S   = LKS
       dim3_E   = LKE
       transpose = .true.
    case ( 'UXY' )
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       transpose = .false.
    case ( 'XYU' )
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim3_max = UKMAX
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       dim3_S   = UKS
       dim3_E   = UKE
       transpose = .true.
    case default
       write(*,*) 'xxx [FILE_EXTERNAL_INPUT_CARTESC_get_dims3D] unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims3D

end module scale_file_external_input_cartesC
