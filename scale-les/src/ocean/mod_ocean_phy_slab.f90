!-------------------------------------------------------------------------------
!> module OCEAN / Physics Fixed-SST
!!
!! @par Description
!!          ocean physics module, fixed SST
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_phy_slab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_driver_setup
  public :: OCEAN_PHY_driver_first
  public :: OCEAN_PHY_driver_final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, save :: DZW    = 50.0_RP    !< water depth of slab ocean [m]
  logical,  private, save :: FLG_CR = .false.    !< is the fixed change rate used?
  real(RP), private, save :: CRATE  = 2.0E-5_RP  !< fixed change rate of water temperature [K/s]

  real(RP), private, save, allocatable :: WHFLX  (:,:)
  real(RP), private, save, allocatable :: PRECFLX(:,:)
  real(RP), private, save, allocatable :: QVFLX  (:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_ocean_vars, only: &
       OCEAN_TYPE_PHY
    implicit none

    real(RP) :: OCEAN_SLAB_DEPTH
    logical  :: OCEAN_SLAB_FLG_CR
    real(RP) :: OCEAN_SLAB_CRATE

    NAMELIST / PARAM_OCEAN_SLAB / &
       OCEAN_SLAB_DEPTH,  &
       OCEAN_SLAB_FLG_CR, &
       OCEAN_SLAB_CRATE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SLAB]/Categ[OCEAN]'

    allocate( WHFLX  (IA,JA) )
    allocate( PRECFLX(IA,JA) )
    allocate( QVFLX  (IA,JA) )

    OCEAN_SLAB_DEPTH  = DZW
    OCEAN_SLAB_FLG_CR = FLG_CR
    OCEAN_SLAB_CRATE  = CRATE

    if ( OCEAN_TYPE_PHY /= 'SLAB' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx OCEAN_TYPE_PHY is not SLAB. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_SLAB,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN_SLAB)

    DZW    = OCEAN_SLAB_DEPTH
    FLG_CR = OCEAN_SLAB_FLG_CR
    CRATE  = OCEAN_SLAB_CRATE

    return
  end subroutine OCEAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for ocean submodel
  subroutine OCEAN_PHY_driver_first
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_vars, only: &
       TW,                  &
       OCEAN_vars_fillhalo
    use mod_cpl_vars, only: &
       CPL_getCPL2Ocn
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean step: Slab'

    call CPL_getCPL2Ocn( WHFLX  (:,:), & ! [OUT]
                         PRECFLX(:,:), & ! [OUT]
                         QVFLX  (:,:)  ) ! [OUT]

    ! update water temperature
    if( FLG_CR ) then

      do j = JS, JE
      do i = IS, IE
        TW(i,j) = TW(i,j) + CRATE * dt
      end do
      end do

    else

      do j = JS, JE
      do i = IS, IE
        TW(i,j) = TW(i,j) - 2.0_RP * WHFLX(i,j) / ( DWATR * CL * DZW ) * dt
      end do
      end do

    end if

    call OCEAN_vars_fillhalo

    return
  end subroutine OCEAN_PHY_driver_first

  subroutine OCEAN_PHY_driver_final
    use mod_ocean_vars, only: &
       TW
    use mod_cpl_vars, only: &
       CPL_putOcn
    implicit none
    !---------------------------------------------------------------------------

    call CPL_putOcn( TW(:,:) ) ! [IN]

    return
  end subroutine OCEAN_PHY_driver_final

end module mod_ocean_phy_slab
