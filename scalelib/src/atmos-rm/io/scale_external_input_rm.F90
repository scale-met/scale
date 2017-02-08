!-------------------------------------------------------------------------------
!> module EXTERNAL INPUT
!!
!! @par Description
!!          External file input RM module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_external_input_rm
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: EXTIN_RM_get_dims_1D
  public :: EXTIN_RM_get_dims_2D
  public :: EXTIN_RM_get_dims_3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  !> get_dims
  subroutine EXTIN_RM_get_dims_1D( &
       dim1_max, &
       dim1_S,   &
       dim1_E,   &
       varname,  &
       axistype  )
    use scale_process, only: &
       PRC_MPIstop
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
    case default
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    end select

    return
  end subroutine EXTIN_RM_get_dims_1D

  subroutine EXTIN_RM_get_dims_2D( &
       dim1_max, &
       dim1_S,   &
       dim1_E,   &
       dim2_max, &
       dim2_S,   &
       dim2_E,   &
       varname,  &
       axistype, &
       transpose )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim1_E
    integer,          intent(out) :: dim2_max
    integer,          intent(out) :: dim2_S
    integer,          intent(out) :: dim2_E
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (XY/XZ/ZX)
    logical,          intent(in)  :: transpose

    integer :: tmp

    select case ( axistype )
    case ( 'XY' )
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    case ( 'ZX', 'XZ' )
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
    case default
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    end select

    if ( transpose ) then
       tmp = dim1_max
       dim1_max = dim2_max
       dim2_max = tmp
       tmp = dim1_S
       dim1_S = dim2_S
       dim2_S = tmp
       tmp = dim1_E
       dim1_E = dim2_E
       dim2_E = tmp
    end if

    return
  end subroutine EXTIN_RM_get_dims_2D

  subroutine EXTIN_RM_get_dims_3D( &
       dim1_max, &
       dim1_S,   &
       dim1_E,   &
       dim2_max, &
       dim2_S,   &
       dim2_E,   &
       dim3_max, &
       dim3_S,   &
       dim3_E,   &
       varname,  &
       axistype, &
       transpose )
    use scale_process, only: &
       PRC_MPIstop
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
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (ZXY/XYZ/Land/Urban)
    logical,          intent(in)  :: transpose

    integer :: tmp

    select case ( axistype )
    case ( 'ZXY', 'XYZ' )
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    case ( 'Land' )
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    case ( 'Urban' )
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    case default
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    end select

    if ( transpose ) then
       tmp = dim1_max
       dim1_max = dim2_max
       dim2_max = dim3_max
       dim3_max = tmp
       tmp = dim1_S
       dim1_S = dim2_S
       dim2_S = dim3_S
       dim3_S = tmp
       tmp = dim1_E
       dim1_E = dim2_E
       dim2_E = dim3_E
       dim3_E = tmp
    end if

    return
  end subroutine EXTIN_RM_get_dims_3D

end module scale_external_input_rm
