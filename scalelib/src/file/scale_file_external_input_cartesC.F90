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
#include "scalelib.h"
module scale_file_external_input_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_urban_grid_cartesC_index
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
       FILE_EXTERNAL_INPUT_get_dims3D, &
       FILE_EXTERNAL_INPUT_read_1D,    &
       FILE_EXTERNAL_INPUT_read_2D,    &
       FILE_EXTERNAL_INPUT_read_3D

    call FILE_EXTERNAL_INPUT_setup

    FILE_EXTERNAL_INPUT_get_dims1D => FILE_EXTERNAL_INPUT_CARTESC_get_dims1D
    FILE_EXTERNAL_INPUT_get_dims2D => FILE_EXTERNAL_INPUT_CARTESC_get_dims2D
    FILE_EXTERNAL_INPUT_get_dims3D => FILE_EXTERNAL_INPUT_CARTESC_get_dims3D

    FILE_EXTERNAL_INPUT_read_1D => FILE_EXTERNAL_INPUT_CARTESC_read_1D
    FILE_EXTERNAL_INPUT_read_2D => FILE_EXTERNAL_INPUT_CARTESC_read_2D
    FILE_EXTERNAL_INPUT_read_3D => FILE_EXTERNAL_INPUT_CARTESC_read_3D

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> get_dims
  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims1D( &
       dim1_size, &
       dim1_max,  &
       dim1_S,    &
       varname,   &
       axistype   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_size
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (Z/X/Y)

    select case ( axistype )
    case ( 'Z' )
       dim1_size = KA
       dim1_max  = KMAX
       dim1_S    = KS
    case ( 'X' )
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
    case ( 'Y' )
       dim1_size = JA
       dim1_max  = JMAXB
       dim1_S    = JSB
    case ( 'OZ' )
       dim1_size = OKA
       dim1_max  = OKMAX
       dim1_S    = OKS
    case default
       LOG_ERROR("FILE_EXTERNAL_INPUT_CARTESC_get_dims1D",*) 'unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims1D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims2D( &
       dim1_size, &
       dim1_max,  &
       dim1_S,    &
       dim2_size, &
       dim2_max,  &
       dim2_S,    &
       transpose, &
       varname,   &
       axistype   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_size
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim2_size
    integer,          intent(out) :: dim2_max
    integer,          intent(out) :: dim2_S
    logical,          intent(out) :: transpose
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (XY/XZ/ZX)

    select case ( axistype )
    case ( 'XY' )
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
       dim2_size = JA
       dim2_max  = JMAXB
       dim2_S    = JSB
       transpose = .false.
    case ( 'ZX' )
       dim1_size = KA
       dim1_max  = KMAX
       dim1_S    = KS
       dim2_size = IA
       dim2_max  = IMAXB
       dim2_S    = ISB
       transpose = .false.
    case ( 'XZ' )
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
       dim2_size = KA
       dim2_max  = KMAX
       dim2_S    = KS
       transpose = .true.
    case default
       LOG_ERROR("FILE_EXTERNAL_INPUT_CARTESC_get_dims2D",*) 'unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims2D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims3D( &
       dim1_size, &
       dim1_max,  &
       dim1_S,    &
       dim2_size, &
       dim2_max,  &
       dim2_S,    &
       dim3_size, &
       dim3_max,  &
       dim3_S,    &
       transpose, &
       varname,   &
       axistype   )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(out) :: dim1_size
    integer,          intent(out) :: dim1_max
    integer,          intent(out) :: dim1_S
    integer,          intent(out) :: dim2_size
    integer,          intent(out) :: dim2_max
    integer,          intent(out) :: dim2_S
    integer,          intent(out) :: dim3_size
    integer,          intent(out) :: dim3_max
    integer,          intent(out) :: dim3_S
    logical,          intent(out) :: transpose
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: axistype     ! axis type (ZXY/XYZ/LXY/XYL/UXY/XYU)

    select case ( axistype )
    case ( 'ZXY' )
       dim1_size = KA
       dim1_max  = KMAX
       dim1_S    = KS
       dim2_size = IA
       dim2_max  = IMAXB
       dim2_S    = ISB
       dim3_size = JA
       dim3_max  = JMAXB
       dim3_S    = JSB
       transpose = .false.
    case ( 'XYZ' )
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
       dim2_size = JA
       dim2_max  = JMAXB
       dim2_S    = JSB
       dim3_size = KA
       dim3_max  = KMAX
       dim3_S    = KS
       transpose = .true.
    case ( 'OXY' ) ! OCEAN
       dim1_size = OKA
       dim1_max  = OKMAX
       dim1_S    = OKS
       dim2_size = IA
       dim2_max  = IMAXB
       dim2_S    = ISB
       dim3_size = JA
       dim3_max  = JMAXB
       dim3_S    = JSB
       transpose = .false.
    case ( 'XYO' ) ! OCEAN
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
       dim2_size = JA
       dim2_max  = JMAXB
       dim2_S    = JSB
       dim3_size = OKA
       dim3_max  = OKMAX
       dim3_S    = OKS
       transpose = .true.
    case ( 'LXY' ) ! LAND
       dim1_size = LKA
       dim1_max  = LKMAX
       dim1_S    = LKS
       dim2_size = IA
       dim2_max  = IMAXB
       dim2_S    = ISB
       dim3_size = JA
       dim3_max  = JMAXB
       dim3_S    = JSB
       transpose = .false.
    case ( 'XYL' ) ! LAND
       dim1_size = IA
       dim1_max = IMAXB
       dim1_S   = ISB
       dim2_size = JA
       dim2_max  = JMAXB
       dim2_S    = JSB
       dim3_size = LKA
       dim3_max  = LKMAX
       dim3_S    = LKS
       transpose = .true.
    case ( 'UXY' )
       dim1_size = UKA
       dim1_max  = UKMAX
       dim1_S    = UKS
       dim2_size = IA
       dim2_max  = IMAXB
       dim2_S    = ISB
       dim3_size = JA
       dim3_max  = JMAXB
       dim3_S    = JSB
       transpose = .false.
    case ( 'XYU' )
       dim1_size = IA
       dim1_max  = IMAXB
       dim1_S    = ISB
       dim2_size = JA
       dim2_max  = JMAXB
       dim2_S    = JSB
       dim3_size = UKA
       dim3_max  = UKMAX
       dim3_S    = UKS
       transpose = .true.
    case default
       LOG_ERROR("FILE_EXTERNAL_INPUT_CARTESC_get_dims3D",*) 'unsupported axis type. Check! axistype:', trim(axistype), ', item:',trim(varname)
       call PRC_abort
    end select

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_get_dims3D


  subroutine FILE_EXTERNAL_INPUT_CARTESC_read_1D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step          )
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type
    real(RP),         intent(out) :: var(:)   !< value of the variable
    integer,          intent(in), optional :: step     !< step number

    call FILE_CARTESC_read( fid, varname, dim_type, & ! [IN]
                            var(:),                 & ! [OUT]
                            step = step             ) ! [IN]

    call FILE_CARTESC_flush( fid )

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_read_1D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_read_2D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step          )
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type
    real(RP),         intent(out) :: var(:,:) !< value of the variable
    integer,          intent(in), optional :: step     !< step number

    call FILE_CARTESC_read( fid, varname, dim_type, & ! [IN]
                            var(:,:),               & ! [OUT]
                            step = step             ) ! [IN]

    call FILE_CARTESC_flush( fid )

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_read_2D

  subroutine FILE_EXTERNAL_INPUT_CARTESC_read_3D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step          )
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    integer,          intent(in)  :: fid        !< file ID
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: dim_type   !< dimension type
    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    integer,          intent(in), optional :: step     !< step number

    call FILE_CARTESC_read( fid, varname, dim_type, & ! [IN]
                            var(:,:,:),             & ! [OUT]
                            step = step             ) ! [IN]

    call FILE_CARTESC_flush( fid )

    return
  end subroutine FILE_EXTERNAL_INPUT_CARTESC_read_3D

end module scale_file_external_input_cartesC
