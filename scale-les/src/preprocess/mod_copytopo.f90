!-------------------------------------------------------------------------------
!> module Copy topography
!!
!! @par Description
!!          subroutines for preparing topography data (copy from parent domain)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_copytopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file_h
  use gtool_file, only: &
     FileRead
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COPYTOPO

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: CNVTOPO_TYPE = -1

  integer, public, parameter :: I_IGNORE     =  0
  integer, public, parameter :: I_GTOPO30    =  1
  integer, public, parameter :: I_DEM50M     =  2
  integer, public, parameter :: I_GMTED2010  =  3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: COPYTOPO_transgrid
  private :: COPYTOPO_setalpha
  private :: COPYTOPO_input_data
  private :: COPYTOPO_mix_data

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: handle = 1
  integer, private            :: itp_nh = 4

  character(len=H_LONG), private :: COPYTOPO_IN_BASENAME  = ''

  real(RP), private :: COPYTOPO_TRANSITION_DX = -1.0_RP  !< thickness of transition region [m]: x
  real(RP), private :: COPYTOPO_TRANSITION_DY = -1.0_RP  !< thickness of transition region [m]: y
  real(RP), private :: COPYTOPO_TRANSFACT     = -1.0_RP  !< stretch factor of transition region
  real(RP), private :: COPYTOPO_FRACX         =  1.0_RP  !< fraction of transition region (x) [0-1]
  real(RP), private :: COPYTOPO_FRACY         =  1.0_RP  !< fraction of transition region (y) [0-1]
  real(RP), private :: COPYTOPO_taux          =  1.0_RP  !< maximum value for mixing tau (x) [s]
  real(RP), private :: COPYTOPO_tauy          =  1.0_RP  !< maximum value for mixing tau (y) [s]

  logical,  private :: COPYTOPO_ENTIRE_REGION = .false.  !< copy parent topo over an entire region
  logical,  private :: COPYTOPO_LINEAR_H      = .true.   !< linear or non-linear profile of relax region
  real(RP), private :: COPYTOPO_EXP_H         = 2.0_RP   !< factor of non-linear profile of relax region

  real(RP), private, allocatable :: CTRX(:)              !< center buffer factor [0-1]: x
  real(RP), private, allocatable :: CTRY(:)              !< center buffer factor [0-1]: y
  real(RP), private, allocatable :: COPYTOPO_alpha(:,:)  !> damping coefficient  [0-1]
  real(RP), private, allocatable :: topo_pd(:,:)         !> topography of parent domain

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup and Main
  subroutine COPYTOPO( &
      topo_cd          ) ! [inout]
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       DX,         &
       DY,         &
       BUFFER_DX,  &
       BUFFER_DY,  &
       BUFFFACT
    use scale_grid_nest, only: &
       NEST_INTERP_LEVEL
    implicit none

    real(RP), intent(inout) :: topo_cd(:,:)

    NAMELIST / PARAM_COPYTOPO / &
       COPYTOPO_IN_BASENAME,    &
       COPYTOPO_TRANSITION_DX,  &
       COPYTOPO_TRANSITION_DY,  &
       COPYTOPO_TRANSFACT,      &
       COPYTOPO_FRACX,          &
       COPYTOPO_FRACY,          &
       COPYTOPO_taux,           &
       COPYTOPO_tauy,           &
       COPYTOPO_ENTIRE_REGION,  &
       COPYTOPO_LINEAR_H,       &
       COPYTOPO_EXP_H

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COPYTOPO]/Categ[COPYTOPO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COPYTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COPYTOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_COPYTOPO)

    if ( COPYTOPO_TRANSITION_DX < 0.0_RP ) then
       COPYTOPO_TRANSITION_DX = BUFFER_DX
    end if
    if ( COPYTOPO_TRANSITION_DY < 0.0_RP ) then
       COPYTOPO_TRANSITION_DY = BUFFER_DY
    end if
    if ( COPYTOPO_TRANSFACT < 0.0_RP ) then
       COPYTOPO_TRANSFACT = BUFFFACT
    end if

    allocate( CTRX(IA) )
    allocate( CTRY(JA) )
    allocate( COPYTOPO_alpha(IA,JA) )
    allocate( topo_pd       (IA,JA) )
    COPYTOPO_alpha(:,:) = 0.0_RP

    itp_nh = int( NEST_INTERP_LEVEL )

    ! copy topography from parent domain to transition region
    call COPYTOPO_transgrid

    call COPYTOPO_setalpha

    call COPYTOPO_input_data( topo_pd ) ! (out)

    call COPYTOPO_mix_data( topo_cd, &  ! (inout)
                            topo_pd  )  ! (in)

    return
  end subroutine COPYTOPO

  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Generate transition grid
  subroutine COPYTOPO_transgrid
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank
    use scale_les_process, only: &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       EPS    => CONST_EPS
    use scale_grid, only: &
       DX,                  &
       DY,                  &
       CXG   => GRID_CXG,   &
       FXG   => GRID_FXG,   &
       CYG   => GRID_CYG,   &
       FYG   => GRID_FYG,   &
       CBFXG => GRID_CBFXG, &
       CBFYG => GRID_CBFYG, &
       BUFFFACT
    implicit none

    real(RP), allocatable :: CTRXG(:) !< center buffer factor [0-1]: x, global
    real(RP), allocatable :: CTRYG(:) !< center buffer factor [0-1]: y, global

    real(RP), allocatable :: buffx(:),  buffy(:)
    real(RP)              :: bufftotx,  bufftoty
    real(RP), allocatable :: transx(:), transy(:)
    real(RP)              :: transtotx, transtoty

    integer :: IAG, JAG
    integer :: imain, ibuff, itrans
    integer :: jmain, jbuff, jtrans
    integer :: copy_is, copy_ie, copy_js, copy_je
    integer :: i, j, ii, jj
    !---------------------------------------------------------------------------

    ! array size (global domain)
    IAG = IHALO + IMAX*PRC_NUM_X + IHALO
    JAG = JHALO + JMAX*PRC_NUM_Y + JHALO

    allocate( buffx (0:IAG) )
    allocate( buffy (0:JAG) )
    allocate( transx(0:IAG) )
    allocate( transy(0:JAG) )
    allocate( CTRXG (  IAG) )
    allocate( CTRYG (  JAG) )

    ! X-direction
    ! calculate buffer grid size
    buffx(:)  = DX
    transx(0) = DX
    bufftotx  = 0.0_RP
    transtotx = 0.0_RP

    do i = IHALO+1, IAG
       if( abs(CBFXG(i) - 0.0_RP) < EPS ) exit
       buffx(i) = buffx(i-1) * BUFFFACT
       bufftotx = bufftotx + buffx(i)
    enddo
    ibuff = i - (IHALO+1)

    do i = 1, IAG
       if( transtotx >= COPYTOPO_TRANSITION_DX ) exit
       transx(i) = transx(i-1) * COPYTOPO_TRANSFACT
       transtotx = transtotx + transx(i)
    enddo
    itrans = i - 1
    imain  = IAG - 2*ibuff - 2*itrans - 2*IHALO

    if ( imain < 1 ) then
       write(*,*) 'xxx Not appropriate transition width for global domain(X).', COPYTOPO_TRANSITION_DX
       write(*,*) '    # of buffer region (one side)', ibuff
       write(*,*) '    # of transion region (one side)', itrans
       call PRC_MPIstop
    endif

    ! calc transition factor (global domaim)
    CTRXG(:) = 0.0_RP
    do i = 1, IHALO+ibuff
       CTRXG(i) = 1.0_RP
    enddo

    if ( itrans > 0 ) then
       copy_is = IHALO+ibuff+1
       copy_ie = IHALO+ibuff+itrans
       do i = copy_is, copy_ie
          CTRXG(i) = (transtotx+bufftotx+FXG(IHALO    )-CXG(i)) / transtotx
       enddo
       copy_is = IHALO+ibuff+itrans+imain+1
       copy_ie = IHALO+ibuff+itrans+imain+itrans+ibuff
       do i = copy_is, copy_ie
          CTRXG(i) = (transtotx+bufftotx-FXG(IAG-IHALO)+CXG(i)) / transtotx
       enddo
    endif

    copy_is = IHALO+ibuff+itrans+imain+itrans+ibuff+1
    copy_ie = IHALO+ibuff+itrans+imain+itrans+ibuff+IHALO
    do i = copy_is, copy_ie
       CTRXG(i) = 1.0_RP
    enddo
    CTRXG(:) = max( min( CTRXG(:), 1.0_RP ), 0.0_RP )

    ! Y-direction
    ! calculate buffer grid size
    buffy(:)  = DY
    transy(0) = DY
    bufftoty  = 0.0_RP
    transtoty = 0.0_RP

    do j = JHALO+1, JAG
       if( abs(CBFYG(j) - 0.0_RP) < EPS ) exit
       buffy(j) = buffy(j-1) * BUFFFACT
       bufftoty = bufftoty + buffy(j)
    enddo
    jbuff = j - (JHALO+1)

    do j = 1, JAG
       if( transtoty >= COPYTOPO_TRANSITION_DY ) exit
       transy(j) = transy(j-1) * COPYTOPO_TRANSFACT
       transtoty = transtoty + transy(j)
    enddo
    jtrans = j - 1
    jmain  = JAG - 2*jbuff - 2*jtrans - 2*JHALO

    if ( jmain < 1 ) then
       write(*,*) 'xxx Not appropriate transition width for global domain(Y).', COPYTOPO_TRANSITION_DY
       write(*,*) '    # of buffer region (one side)', jbuff
       write(*,*) '    # of transion region (one side)', jtrans
       call PRC_MPIstop
    endif

    ! calc transition factor (global domaim)
    CTRYG(:) = 0.0_RP
    do j = 1, JHALO+jbuff
       CTRYG(j) = 1.0_RP
    enddo

    if ( jtrans > 0 ) then
       copy_js = JHALO+jbuff+1
       copy_je = JHALO+jbuff+jtrans
       do j = copy_js, copy_je
          CTRYG(j) = (transtoty+bufftoty+FYG(JHALO    )-CYG(j)) / transtoty
       enddo
       copy_js = JHALO+jbuff+jtrans+jmain+1
       copy_je = JHALO+jbuff+jtrans+jmain+jtrans+jbuff
       do j = copy_js, copy_je
          CTRYG(j) = (transtoty+bufftoty-FYG(JAG-JHALO)+CYG(j)) / transtoty
       enddo
    endif

    copy_js = JHALO+jbuff+jtrans+jmain+jtrans+jbuff+1
    copy_je = JHALO+jbuff+jtrans+jmain+jtrans+jbuff+JHALO
    do j = copy_js, copy_je
       CTRYG(j) = 1.0_RP
    enddo
    CTRYG(:) = max( min( CTRYG(:), 1.0_RP ), 0.0_RP )

    ! horizontal coordinate (local domaim)
    do i = 1, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX
       CTRX(i) = CTRXG(ii)
    enddo
    do j = 1, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX
       CTRY(j) = CTRYG(jj)
    enddo

    deallocate( transx )
    deallocate( transy )

    return
  end subroutine COPYTOPO_transgrid


  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine COPYTOPO_setalpha
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait

    real(RP) :: coef_x, alpha_x1
    real(RP) :: coef_y, alpha_y1
    real(RP) :: ee1

    integer :: i, j
    !---------------------------------------------------------------------------

    ! check invalid fraction
    COPYTOPO_FRACX = max( min( COPYTOPO_FRACX, 1.0_RP ), EPS )
    COPYTOPO_FRACY = max( min( COPYTOPO_FRACY, 1.0_RP ), EPS )

    if ( COPYTOPO_taux <= 0.0_RP ) then ! invalid tau
       coef_x = 0.0_RP
    else
       coef_x = 1.0_RP / COPYTOPO_taux
    endif

    if ( COPYTOPO_tauy <= 0.0_RP ) then ! invalid tau
       coef_y = 0.0_RP
    else
       coef_y = 1.0_RP / COPYTOPO_tauy
    endif

    do j = 1, JA
    do i = 1, IA
       ee1 = CTRX(i)
       if ( ee1 <= 1.0_RP - COPYTOPO_FRACX ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + COPYTOPO_FRACX ) / COPYTOPO_FRACX
       endif

       if ( COPYTOPO_LINEAR_H ) then
          alpha_x1 = coef_x * ee1
       else
          alpha_x1 = coef_x * ee1 * exp( -(1.0_RP-ee1) * COPYTOPO_EXP_H )
       end if

       ee1 = CTRY(j)
       if ( ee1 <= 1.0_RP - COPYTOPO_FRACY ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + COPYTOPO_FRACY ) / COPYTOPO_FRACY
       endif

       if ( COPYTOPO_LINEAR_H ) then
          alpha_y1 = coef_y * ee1
       else
          alpha_y1 = coef_y * ee1 * exp( -(1.0_RP-ee1) * COPYTOPO_EXP_H )
       end if

       COPYTOPO_alpha(i,j) = max( alpha_x1, alpha_y1 )
    enddo
    enddo

    call COMM_vars8( COPYTOPO_alpha(:,:), 1 )
    call COMM_wait ( COPYTOPO_alpha(:,:), 1 )

    return
  end subroutine COPYTOPO_setalpha


  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine COPYTOPO_input_data( &
      topo_pd      ) ! (out)
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid_nest, only: &
       PARENT_IMAX,     &
       PARENT_JMAX,     &
       NEST_TILE_NUM_X, &
       NEST_TILE_NUM_Y, &
       NEST_TILE_ID,    &
       NEST_domain_shape
    use scale_interpolation_nest, only: &
       INTRPNEST_domain_compatibility,  &
       INTRPNEST_interp_fact_latlon,    &
       INTRPNEST_interp_2d
    use scale_grid_real, only: &
       LAT => REAL_LAT, &
       LON => REAL_LON
    implicit none

    real(RP), intent(out) :: topo_pd(:,:)

    real(RP)              :: dummy (1,1,1)
    real(RP), allocatable :: read2D(:,:)
    real(RP), allocatable :: lon_org (:,:)
    real(RP), allocatable :: lat_org (:,:)
    real(RP), allocatable :: topo_org(:,:)
    real(RP), allocatable :: hfact(:,:,:)
    integer,  allocatable :: igrd (:,:,:)
    integer,  allocatable :: jgrd (:,:,:)

    integer :: IALL, JALL  ! number of grids for whole domain
    integer :: PTI, PTJ    ! number of grids for a tile
    integer :: tilei, tilej
    integer :: rank
    integer :: i
    integer :: cxs, cxe, cys, cye  ! for child domain
    integer :: pxs, pxe, pys, pye  ! for parent domain
    !---------------------------------------------------------------------------

    PTI  = PARENT_IMAX(handle)
    PTJ  = PARENT_JMAX(handle)
    IALL = PTI * NEST_TILE_NUM_X
    JALL = PTJ * NEST_TILE_NUM_Y

    allocate( hfact   ( IA,   JA,  itp_nh ) )
    allocate( igrd    ( IA,   JA,  itp_nh ) )
    allocate( jgrd    ( IA,   JA,  itp_nh ) )
    allocate( lon_org ( IALL, JALL        ) )
    allocate( lat_org ( IALL, JALL        ) )
    allocate( topo_org( IALL, JALL        ) )

    do i = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(i)

       call NEST_domain_shape ( tilei, tilej,   & ! [out]
                                cxs,   cxe,     & ! [out]
                                cys,   cye,     & ! [out]
                                pxs,   pxe,     & ! [out]
                                pys,   pye,     & ! [out]
                                i               ) ! [in ]
       allocate( read2D  ( tilei,tilej ) )

       call FileRead( read2D(:,:), COPYTOPO_IN_BASENAME, "lon",  1, rank )
       lon_org (cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye) * D2R
       call FileRead( read2D(:,:), COPYTOPO_IN_BASENAME, "lat",  1, rank )
       lat_org (cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye) * D2R
       call FileRead( read2D(:,:), COPYTOPO_IN_BASENAME, "TOPO", 1, rank )
       topo_org(cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye)

       deallocate( read2D )
    end do

    call INTRPNEST_domain_compatibility( lon_org(:,:), lat_org(:,:), dummy(:,:,:), &
                                         LON(:,:),     LAT(:,:),     dummy(:,:,:), &
                                         skip_z=.true.                             )
    !call check_domain_compatibility( lon_org(:,:), lat_org(:,:), dummy(:,:,:), &
    !                                 LON(:,:),     LAT(:,:),     dummy(:,:,:), &
    !                                 skip_z=.true.                             )

    call INTRPNEST_interp_fact_latlon( hfact  (:,:,:),   & ! (out)
                                       igrd   (:,:,:),   & ! (out)
                                       jgrd   (:,:,:),   & ! (out)
                                       LAT    (:,:),     & ! (in)
                                       LON    (:,:),     & ! (in)
                                       IA,   JA,         & ! (in)
                                       lat_org(:,:),     & ! (in)
                                       lon_org(:,:),     & ! (in)
                                       IALL, JALL        ) ! (in)

    call INTRPNEST_interp_2d( topo_pd (:,:),    &
                              topo_org(:,:),    &
                              hfact   (:,:,:),  &
                              igrd    (:,:,:),  &
                              jgrd    (:,:,:),  &
                              IA, JA            )

    call COMM_vars8( topo_pd(:,:), 1 )
    call COMM_wait ( topo_pd(:,:), 1 )

    return
  end subroutine COPYTOPO_input_data

  !-----------------------------------------------------------------------------
  !> Mixing TOPO data using parent and child domains
  subroutine COPYTOPO_mix_data( &
      topo_cd,  &  ! (inout)
      topo_pd    ) ! (in)
    implicit none
    real(RP), intent(inout) :: topo_cd(:,:)  ! topography of child domain (mine)
    real(RP), intent(in)    :: topo_pd(:,:)  ! topography of parent domain

    real(RP) :: frac
    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( COPYTOPO_ENTIRE_REGION ) then
       topo_cd(:,:) = topo_pd(:,:)
    else
       do j = 1, JA
       do i = 1, IA
          frac = COPYTOPO_alpha(i,j)
          topo_cd(i,j) = topo_cd(i,j) * ( 1.0_RP - frac ) &
                       + topo_pd(i,j) * frac
       end do
       end do
    endif

    return
  end subroutine COPYTOPO_mix_data

end module mod_copytopo
