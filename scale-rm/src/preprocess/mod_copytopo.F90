!-------------------------------------------------------------------------------
!> module Copy topography
!!
!! @par Description
!!          subroutines for preparing topography data (copy from parent domain)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_copytopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
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

  character(len=H_LONG), private :: COPYTOPO_IN_BASENAME   = ''
  character(len=H_LONG), private :: COPYTOPO_IN_VARNAME    = 'topo'
  real(RP),              private :: COPYTOPO_TRANSITION_DX = -1.0_RP  !< thickness of transition region [m]: x
  real(RP),              private :: COPYTOPO_TRANSITION_DY = -1.0_RP  !< thickness of transition region [m]: y
  real(RP),              private :: COPYTOPO_FRACX         =  1.0_RP  !< fraction of transition region (x) (0-1)
  real(RP),              private :: COPYTOPO_FRACY         =  1.0_RP  !< fraction of transition region (y) (0-1)
  real(RP),              private :: COPYTOPO_taux          =  1.0_RP  !< maximum value for mixing tau (x) [s]
  real(RP),              private :: COPYTOPO_tauy          =  1.0_RP  !< maximum value for mixing tau (y) [s]
  logical,               private :: COPYTOPO_ENTIRE_REGION = .false.  !< copy parent topo over an entire region
  logical,               private :: COPYTOPO_LINEAR_H      = .true.   !< linear or non-linear profile of relax region
  real(RP),              private :: COPYTOPO_EXP_H         =  2.0_RP  !< factor of non-linear profile of relax region
  integer,               private :: COPYTOPO_FILTER_ORDER  = 2
  integer,               private :: COPYTOPO_FILTER_NITER  = 20

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup and Main
  subroutine COPYTOPO( &
       TOPO_child )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: TOPO_child(:,:) !< topography of child domain

    real(RP) :: CTRX       (IA)    !< transition factor   (0-1): x
    real(RP) :: CTRY       (JA)    !< transition factor   (0-1): y
    real(RP) :: alpha      (IA,JA) !< dumping coefficient (0-1)
    real(RP) :: TOPO_parent(IA,JA) !< topography of parent domain

    namelist / PARAM_COPYTOPO / &
       COPYTOPO_IN_BASENAME,    &
       COPYTOPO_IN_VARNAME,     &
       COPYTOPO_TRANSITION_DX,  &
       COPYTOPO_TRANSITION_DY,  &
       COPYTOPO_FRACX,          &
       COPYTOPO_FRACY,          &
       COPYTOPO_taux,           &
       COPYTOPO_tauy,           &
       COPYTOPO_ENTIRE_REGION,  &
       COPYTOPO_LINEAR_H,       &
       COPYTOPO_EXP_H,          &
       COPYTOPO_FILTER_ORDER,   &
       COPYTOPO_FILTER_NITER

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("COPYTOPO",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COPYTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("COPYTOPO",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("COPYTOPO",*) 'Not appropriate names in namelist PARAM_COPYTOPO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_COPYTOPO)

    ! copy topography from parent domain to transition region

    call COPYTOPO_transgrid( CTRX(:), CTRY(:) ) ! [OUT]

    call COPYTOPO_setalpha ( CTRX(:), CTRY(:), & ! [IN]
                             alpha(:,:)        ) ! [OUT]

    call COPYTOPO_input_data( TOPO_parent(:,:) ) ! [OUT]

    call COPYTOPO_mix_data( TOPO_parent(:,:), & ! [IN]
                            alpha      (:,:), & ! [IN]
                            TOPO_child (:,:)  ) ! [INOUT]

    return
  end subroutine COPYTOPO

  !-----------------------------------------------------------------------------
  !> Generate transition grid
  subroutine COPYTOPO_transgrid( &
       CTRX, &
       CTRY  )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    use scale_prc_cartesC, only: &
       PRC_2Drank
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_grid_cartesC, only: &
       CXG   => ATMOS_GRID_CARTESC_CXG,   &
       FXG   => ATMOS_GRID_CARTESC_FXG,   &
       CYG   => ATMOS_GRID_CARTESC_CYG,   &
       FYG   => ATMOS_GRID_CARTESC_FYG,   &
       CBFXG => ATMOS_GRID_CARTESC_CBFXG, &
       CBFYG => ATMOS_GRID_CARTESC_CBFYG
    implicit none

    real(RP), intent(out) :: CTRX(IA) !< transition factor (0-1): x, local
    real(RP), intent(out) :: CTRY(JA) !< transition factor (0-1): y, local

    real(RP) :: CTRXG (  IAG) !< transition factor (0-1): x, global
    real(RP) :: CTRYG (  JAG) !< transition factor (0-1): y, global

    real(RP) :: bufftotx, transtotx
    real(RP) :: bufftoty, transtoty
    integer  :: imain, ibuff, itrans
    integer  :: jmain, jbuff, jtrans
    integer  :: copy_is, copy_ie
    integer  :: copy_js, copy_je

    integer  :: i, j, ii, jj
    !---------------------------------------------------------------------------

    ! X-direction
    ! calculate buffer grid size

    do i = IHALO+1, IAG
       if( abs(CBFXG(i)) < EPS ) exit
    enddo
    ibuff = i - 1 - IHALO
    if ( ibuff == 0 ) then
       bufftotx = 0.0_RP
    else
       bufftotx = FXG(ibuff+IHALO) - FXG(IHALO)
    end if

    if ( COPYTOPO_TRANSITION_DX < 0.0_RP ) then
       COPYTOPO_TRANSITION_DX = bufftotx
    end if
    do i = ibuff+IHALO+1, IAG
       if( CXG(i) - bufftotx - FXG(IHALO) >= COPYTOPO_TRANSITION_DX ) exit
    enddo
    itrans = i - 1 - IHALO - ibuff
    if ( itrans == 0 ) then
       transtotx = 0.0_RP
    else
       transtotx = FXG(itrans+ibuff+IHALO) - FXG(IHALO) - bufftotx
    end if

    imain  = IAG - 2*ibuff - 2*itrans - 2*IHALO

    if ( imain < 1 ) then
       LOG_ERROR("COPYTOPO_transgrid",*) 'Not appropriate transition width for global domain(X).', COPYTOPO_TRANSITION_DX
       LOG_ERROR_CONT(*) '# of buffer   region (one side) = ', ibuff
       LOG_ERROR_CONT(*) '# of transion region (one side) = ', itrans
       call PRC_abort
    endif

    ! calc transition factor (global domaim)
    CTRXG(:) = 0.0_RP

    copy_is = 1
    copy_ie = IHALO+ibuff
    !$omp parallel do
    do i = copy_is, copy_ie
       CTRXG(i) = 1.0_RP
    enddo

    if ( itrans > 0 ) then
       copy_is = IHALO+ibuff+1
       copy_ie = IHALO+ibuff+itrans
       !$omp parallel do
       do i = copy_is, copy_ie
          CTRXG(i) = (transtotx+bufftotx+FXG(IHALO    )-CXG(i)) / transtotx
       enddo
       copy_is = IAG - IHALO - ibuff - itrans + 1
       copy_ie = IAG - IHALO - ibuff
       !$omp parallel do
       do i = copy_is, copy_ie
          CTRXG(i) = (transtotx+bufftotx-FXG(IAG-IHALO)+CXG(i)) / transtotx
       enddo
    endif

    copy_is = IAG - IHALO - ibuff + 1
    copy_ie = IAG
    !$omp parallel do
    do i = copy_is, copy_ie
       CTRXG(i) = 1.0_RP
    enddo

    CTRXG(:) = max( min( CTRXG(:), 1.0_RP ), 0.0_RP )

    ! Y-direction
    ! calculate buffer grid size
    do j = JHALO+1, JAG
       if( abs(CBFYG(j)) < EPS ) exit
    enddo
    jbuff = j - 1 - JHALO
    if ( jbuff == 0 ) then
       bufftoty = 0.0_RP
    else
       bufftoty = FYG(jbuff+JHALO) - FYG(JHALO)
    end if

    if ( COPYTOPO_TRANSITION_DY < 0.0_RP ) then
       COPYTOPO_TRANSITION_DY = bufftoty
    end if
    do j = jbuff+JHALO+1, JAG
       if( CYG(j) - bufftoty - FYG(JHALO) >= COPYTOPO_TRANSITION_DY ) exit
    enddo
    jtrans = j - 1 - JHALO - jbuff
    if ( jtrans == 0 ) then
       transtoty = 0.0_RP
    else
       transtoty = FYG(jtrans+jbuff+JHALO) - FYG(JHALO) - bufftoty
    end if

    jmain  = JAG - 2*jbuff - 2*jtrans - 2*JHALO

    if ( jmain < 1 ) then
       LOG_ERROR("COPYTOPO_transgrid",*) 'Not appropriate transition width for global domain(Y).', COPYTOPO_TRANSITION_DY
       LOG_ERROR_CONT(*) '# of buffer   region (one side)', jbuff
       LOG_ERROR_CONT(*) '# of transion region (one side)', jtrans
       call PRC_abort
    endif

    ! calc transition factor (global domaim)
    CTRYG(:) = 0.0_RP

    copy_js = 1
    copy_je = JHALO+jbuff
    !$omp parallel do
    do j = copy_js, copy_je
       CTRYG(j) = 1.0_RP
    enddo

    if ( jtrans > 0 ) then
       copy_js = JHALO+jbuff+1
       copy_je = JHALO+jbuff+jtrans
       !$omp parallel do
       do j = copy_js, copy_je
          CTRYG(j) = (transtoty+bufftoty+FYG(JHALO    )-CYG(j)) / transtoty
       enddo
       copy_js = JAG - JHALO - jbuff - jtrans + 1
       copy_je = JAG - JHALO - jbuff
       !$omp parallel do
       do j = copy_js, copy_je
          CTRYG(j) = (transtoty+bufftoty-FYG(JAG-JHALO)+CYG(j)) / transtoty
       enddo
    endif

    copy_js = JAG - JHALO - jbuff + 1
    copy_je = JAG
    !$omp parallel do
    do j = copy_js, copy_je
       CTRYG(j) = 1.0_RP
    enddo

    CTRYG(:) = max( min( CTRYG(:), 1.0_RP ), 0.0_RP )



    ! horizontal coordinate (local domaim)
    !$omp parallel do &
    !$omp private(ii)
    do i = 1, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX
       CTRX(i) = CTRXG(ii)
    enddo

    !$omp parallel do &
    !$omp private(jj)
    do j = 1, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX
       CTRY(j) = CTRYG(jj)
    enddo

    LOG_INFO("COPYTOPO_transgrid",*) 'COPYTOPO grid Information'
    LOG_INFO_CONT(*) 'size of buffer   region (x and y; one side) = ', bufftotx, bufftoty
    LOG_INFO_CONT(*) '#    of buffer   region (x and y; one side) = ', ibuff, jbuff
    LOG_INFO_CONT(*) 'size of transion region (x and y; one side) = ', transtotx, transtoty
    LOG_INFO_CONT(*) '#    of transion region (x and y; one side) = ', itrans, jtrans


    return
  end subroutine COPYTOPO_transgrid


  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine COPYTOPO_setalpha( &
       CTRX, &
       CTRY, &
       alpha )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(in)  :: CTRX (IA)    !< transition factor   (0-1): x
    real(RP), intent(in)  :: CTRY (JA)    !< transition factor   (0-1): y
    real(RP), intent(out) :: alpha(IA,JA) !< dumping coefficient (0-1)

    real(RP) :: coef_x, alpha_x1
    real(RP) :: coef_y, alpha_y1
    real(RP) :: ee1

    integer  :: i, j
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

    !$omp parallel do &
    !$omp private(ee1,alpha_x1,alpha_y1)
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
       endif

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
       endif

       alpha(i,j) = max( alpha_x1, alpha_y1 )
    enddo
    enddo

    call COMM_vars8( alpha(:,:), 1 )
    call COMM_wait ( alpha(:,:), 1 )

    return
  end subroutine COPYTOPO_setalpha


  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine COPYTOPO_input_data( &
       TOPO_parent )
    use scale_const, only: &
       EPS => CONST_EPS, &
       D2R => CONST_D2R
    use scale_file, only: &
       FILE_open, &
       FILE_read, &
       FILE_close
    use scale_interp, only: &
       INTERP_domain_compatibility, &
       INTERP_factor2d,             &
       INTERP_interp2d
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_real, only: &
       LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON => ATMOS_GRID_CARTESC_REAL_LON
    use scale_comm_cartesC_nest, only: &
       NEST_INTERP_LEVEL => COMM_CARTESC_NEST_INTERP_LEVEL, &
       NEST_domain_shape => COMM_CARTESC_NEST_domain_shape, &
       NEST_TILE_NUM_X   => COMM_CARTESC_NEST_TILE_NUM_X,   &
       NEST_TILE_NUM_Y   => COMM_CARTESC_NEST_TILE_NUM_Y,   &
       NEST_TILE_ID      => COMM_CARTESC_NEST_TILE_ID,      &
       PARENT_IMAX,       &
       PARENT_JMAX
    use scale_filter, only: &
       FILTER_hyperdiff
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    real(RP), intent(out) :: TOPO_parent(:,:)

    real(RP), allocatable :: LON_org (:,:)
    real(RP), allocatable :: LAT_org (:,:)
    real(RP), allocatable :: TOPO_org(:,:)
    real(RP), allocatable :: read2D  (:,:)

    real(RP) :: dummy(1,1)
    integer  :: idx_i(IA,JA,NEST_INTERP_LEVEL)
    integer  :: idx_j(IA,JA,NEST_INTERP_LEVEL)
    real(RP) :: hfact(IA,JA,NEST_INTERP_LEVEL)

    real(RP) :: TOPO_sign(IA,JA)

    integer :: IA_org, JA_org     ! number of grids for whole domain
    integer :: tilei, tilej
    integer :: cxs, cxe, cys, cye ! for child domain
    integer :: pxs, pxe, pys, pye ! for parent domain
    integer :: rank

    real(RP) :: ocean_flag

    integer :: fid
    integer :: n, i, j
    !---------------------------------------------------------------------------

    IA_org = PARENT_IMAX(handle) * NEST_TILE_NUM_X
    JA_org = PARENT_JMAX(handle) * NEST_TILE_NUM_Y

    allocate( LON_org (IA_org,JA_org) )
    allocate( LAT_org (IA_org,JA_org) )
    allocate( TOPO_org(IA_org,JA_org) )

    do n = 1, size( NEST_TILE_ID(:) )
       ! read data from split files
       rank = NEST_TILE_ID(n)

       call NEST_domain_shape( tilei, tilej, & ! [OUT]
                               cxs,   cxe,   & ! [OUT]
                               cys,   cye,   & ! [OUT]
                               pxs,   pxe,   & ! [OUT]
                               pys,   pye,   & ! [OUT]
                               n             ) ! [IN]

       allocate( read2D(tilei,tilej) )

       call FILE_open( COPYTOPO_IN_BASENAME,          & ! [IN]
                       fid,                           & ! [OUT]
                       aggregate=.false., rankid=rank ) ! [IN]

       call FILE_read( fid, "lon",  read2D(:,:) )
       LON_org (cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye) * D2R
       call FILE_read( fid, "lat",  read2D(:,:) )
       LAT_org (cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye) * D2R
       call FILE_read( fid, COPYTOPO_IN_VARNAME, read2D(:,:) )
       TOPO_org(cxs:cxe,cys:cye) = read2D(pxs:pxe,pys:pye)
       deallocate( read2D )

       call FILE_close( fid )

    enddo

    call INTERP_domain_compatibility( LON_org(:,:), LAT_org(:,:), dummy(:,:),             &
                                      LON(:,:),     LAT(:,:),     dummy(:,:), dummy(:,:), &
                                      skip_z=.true.                                       )

    call INTERP_factor2d( NEST_INTERP_LEVEL, & ! [IN]
                          IA_org, JA_org,    & ! [IN]
                          IA, JA,            & ! [IN]
                          LON_org(:,:),      & ! [IN]
                          LAT_org(:,:),      & ! [IN]
                          LON    (:,:),      & ! [IN]
                          LAT    (:,:),      & ! [IN]
                          idx_i  (:,:,:),    & ! [OUT]
                          idx_j  (:,:,:),    & ! [OUT]
                          hfact  (:,:,:)     ) ! [OUT]

    call INTERP_interp2d( NEST_INTERP_LEVEL,  & ! [IN]
                          IA_org, JA_org,     & ! [IN]
                          IA, JA,             & ! [IN]
                          idx_i      (:,:,:), & ! [IN]
                          idx_j      (:,:,:), & ! [IN]
                          hfact      (:,:,:), & ! [IN]
                          TOPO_org   (:,:),   & ! [IN]
                          TOPO_parent(:,:)    ) ! [OUT]

    if ( COPYTOPO_FILTER_NITER > 1 ) then
       !$omp parallel do &
       !$omp private(ocean_flag)
       do j = 1, JA
       do i = 1, IA
          ocean_flag = ( 0.5_RP + sign( 0.5_RP, LANDUSE_fact_ocean(i,j) - 1.0_RP + EPS ) ) & ! fact_ocean==1
                     * ( 0.5_RP + sign( 0.5_RP, EPS - abs(TOPO_parent(i,j)) ) )              ! |TOPO| > EPS
          TOPO_sign(i,j) = sign( 1.0_RP, TOPO_parent(i,j) ) * ( 1.0_RP - ocean_flag )
       end do
       end do

       call FILTER_hyperdiff( IA, ISB, IEB, JA, JSB, JEB, &
                              TOPO_parent(:,:), COPYTOPO_FILTER_ORDER, COPYTOPO_FILTER_NITER, &
                              limiter_sign = TOPO_sign(:,:) )
    end if


    call COMM_vars8( TOPO_parent(:,:), 1 )
    call COMM_wait ( TOPO_parent(:,:), 1 )


    return
  end subroutine COPYTOPO_input_data

  !-----------------------------------------------------------------------------
  !> Mixing TOPO data using parent and child domains
  subroutine COPYTOPO_mix_data( &
       TOPO_parent, &
       alpha,       &
       TOPO_child   )
    implicit none

    real(RP), intent(in)    :: TOPO_parent(:,:) ! topography of parent domain
    real(RP), intent(in)    :: alpha      (:,:) ! dumping coefficient (0-1)
    real(RP), intent(inout) :: TOPO_child (:,:) ! topography of child  domain

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( COPYTOPO_ENTIRE_REGION ) then
       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          TOPO_child(i,j) = TOPO_parent(i,j)
       enddo
       enddo
    else
       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          TOPO_child(i,j) = ( 1.0_RP-alpha(i,j) ) * TOPO_child (i,j) &
                          + (        alpha(i,j) ) * TOPO_parent(i,j)
       enddo
       enddo
    endif

    return
  end subroutine COPYTOPO_mix_data

end module mod_copytopo
