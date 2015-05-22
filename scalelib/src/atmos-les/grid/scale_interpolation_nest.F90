!-------------------------------------------------------------------------------
!> module INTERPOLATION (nesting system)
!!
!! @par Description
!!          INTERPOLATION module for nesting system
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-02-10 (R.Yoshida)  [new] rearranged sub-routines
!!
!<
!-------------------------------------------------------------------------------
module scale_interpolation_nest
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  use scale_index
  use scale_tracer
  use scale_const, only: &
     r_in_m => CONST_RADIUS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: INTRPNEST_setup
  public :: INTRPNEST_domain_compatibility
  public :: INTRPNEST_interp_fact_latlon
  public :: INTRPNEST_interp_fact_llz

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: INTRPNEST_search_nearest_block
  private :: INTRPNEST_search_horiz_1points
  private :: INTRPNEST_search_horiz_3points
  private :: INTRPNEST_search_horiz_4points
  private :: INTRPNEST_search_horiz_8points
  private :: INTRPNEST_search_horiz_12points
  private :: INTRPNEST_search_vert_offline
  private :: INTRPNEST_search_vert_online
  private :: INTRPNEST_interp_2d_1points
  private :: INTRPNEST_interp_3d_1points
  private :: INTRPNEST_interp_2d_3points
  private :: INTRPNEST_interp_3d_3points
  private :: INTRPNEST_interp_2d_4points
  private :: INTRPNEST_interp_3d_4points
  private :: INTRPNEST_interp_2d_8points
  private :: INTRPNEST_interp_3d_8points
  private :: INTRPNEST_interp_2d_12points
  private :: INTRPNEST_interp_3d_12points
  private :: INTRPNEST_haversine

  !-----------------------------------------------------------------------------
  abstract interface
     subroutine INTRPNEST_intfc_search_h(  &
          hfact,   & ! (out)
          igrd,    & ! (out)
          jgrd,    & ! (out)
          mylat,   & ! (in)
          mylon,   & ! (in)
          inlat,   & ! (in)
          inlon,   & ! (in)
          is,      & ! (in)
          ie,      & ! (in)
          js,      & ! (in)
          je       ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: hfact(:)      ! horizontal interp factor
       integer,  intent(out) :: igrd (:)      ! grid points of interp target
       integer,  intent(out) :: jgrd (:)      ! grid points of interp target
       real(RP), intent(in)  :: mylat         ! latitude data of mine
       real(RP), intent(in)  :: mylon         ! longitude data of mine
       real(RP), intent(in)  :: inlat(:,:)    ! latitude  data of you (input)
       real(RP), intent(in)  :: inlon(:,:)    ! longitude data of you (input)
       integer,  intent(in)  :: is            ! start index for x-direction
       integer,  intent(in)  :: ie            ! end   index for x-direction
       integer,  intent(in)  :: js            ! start index for y-direction
       integer,  intent(in)  :: je            ! end   index for y-direction
     end subroutine INTRPNEST_intfc_search_h
  end interface
  procedure(INTRPNEST_intfc_search_h), pointer :: INTRPNEST_search_horiz => NULL()
  private :: INTRPNEST_search_horiz

  !-----------------------------------------------------------------------------
  abstract interface
     subroutine INTRPNEST_intfc_search_v(  &
          vfact,   & ! (out)
          kgrd,    & ! (out)
          ncopy,   & ! (out)
          igrd,    & ! (in)
          jgrd,    & ! (in)
          myhgt,   & ! (in)
          inhgt,   & ! (in)
          iloc,    & ! (in)
          jloc,    & ! (in)
          ks,      & ! (in)
          ke,      & ! (in)
          inKA,    & ! (in)
          lndgrd   ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(inout) :: vfact(:,:,:,:,:)! vertical interp factor
       integer,  intent(inout) :: kgrd (:,:,:,:,:)! grid points of interp target
       integer,  intent(out) :: ncopy(:)          ! number of daughter's layers below the parent's lowest layer
       integer,  intent(in)  :: igrd(:)           ! grid points of interp target
       integer,  intent(in)  :: jgrd(:)           ! grid points of interp target
       real(RP), intent(in)  :: myhgt(:)          ! height data of mine
       real(RP), intent(in)  :: inhgt(:,:,:)      ! height data of you (input)
       integer,  intent(in)  :: iloc              ! locator index for x-direction
       integer,  intent(in)  :: jloc              ! locator index for y-direction
       integer,  intent(in)  :: ks                ! start index for z-direction
       integer,  intent(in)  :: ke                ! end   index for z-direction
       integer,  intent(in)  :: inKA              ! grid number of you (input)
       logical,  intent(in)  :: lndgrd            ! flag of land grid
     end subroutine INTRPNEST_intfc_search_v
  end interface
  procedure(INTRPNEST_intfc_search_v), pointer :: INTRPNEST_search_vert => NULL()
  private :: INTRPNEST_search_vert

  !-----------------------------------------------------------------------------
  abstract interface
     subroutine INTRPNEST_intfc_interp_2d(  &
          intp,    & ! (out)
          ref,     & ! (in)
          hfact,   & ! (in)
          igrd,    & ! (in)
          jgrd,    & ! (in)
          ia,      & ! (in)
          ja       ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: intp(:,:)      ! interpolated data
       real(RP), intent(in)  :: ref (:,:)      ! reference data
       real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
       integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
       integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
       integer,  intent(in)  :: ia             ! grid number of mine
       integer,  intent(in)  :: ja             ! grid number of mine
     end subroutine INTRPNEST_intfc_interp_2d
  end interface
  procedure(INTRPNEST_intfc_interp_2d), pointer :: INTRPNEST_interp_2d => NULL()
  public :: INTRPNEST_interp_2d

  !-----------------------------------------------------------------------------
  abstract interface
     subroutine INTRPNEST_intfc_interp_3d(  &
          intp,    & ! (out)
          ref,     & ! (in)
          hfact,   & ! (in)
          vfact,   & ! (in)
          kgrd,    & ! (in)
          igrd,    & ! (in)
          jgrd,    & ! (in)
          ia,      & ! (in)
          ja,      & ! (in)
          ks,      & ! (in)
          ke,      & ! (in)
          logwegt  ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: intp(:,:,:)       ! interpolated data
       real(RP), intent(in)  :: ref (:,:,:)       ! reference data
       real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
       real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
       integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
       integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
       integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
       integer,  intent(in)  :: ia                ! grid number of mine
       integer,  intent(in)  :: ja                ! grid number of mine
       integer,  intent(in)  :: ks                ! start grid number of mine
       integer,  intent(in)  :: ke                ! end grid number of mine
       logical,  intent(in), optional :: logwegt
     end subroutine INTRPNEST_intfc_interp_3d
  end interface
  procedure(INTRPNEST_intfc_interp_3d), pointer :: INTRPNEST_interp_3d => NULL()
  public :: INTRPNEST_interp_3d

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: large_number_1  = 9.999E+15_RP
  real(RP), private, parameter :: large_number_2  = 9.888E+15_RP
  real(RP), private, parameter :: large_number_3  = 9.777E+15_RP
  real(RP), private, parameter :: large_number_4  = 9.666E+15_RP
  real(RP), private, parameter :: large_number_5  = 9.555E+15_RP
  real(RP), private, parameter :: large_number_6  = 9.444E+15_RP
  real(RP), private, parameter :: large_number_7  = 9.333E+15_RP
  real(RP), private, parameter :: large_number_8  = 9.222E+15_RP
  real(RP), private, parameter :: large_number_9  = 9.111E+15_RP
  real(RP), private, parameter :: large_number_10 = 9.000E+15_RP
  real(RP), private, parameter :: large_number_11 = 8.999E+15_RP
  real(RP), private, parameter :: large_number_12 = 8.888E+15_RP

  integer,  private            :: divnum
  integer,  private            :: itp_nh

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTRPNEST_setup ( &
      interp_search_divnum,  &
      NEST_INTERP_LEVEL,     &
      OFFLINE     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer, intent(in) :: interp_search_divnum
    integer, intent(in) :: NEST_INTERP_LEVEL
    logical, intent(in) :: OFFLINE

    character(7) :: select_type
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NEST]/Categ[GRID INTERP]'

    divnum = interp_search_divnum

    select case ( NEST_INTERP_LEVEL )
    case ( 1 )
       INTRPNEST_search_horiz => INTRPNEST_search_horiz_1points
       INTRPNEST_interp_2d    => INTRPNEST_interp_2d_1points
       INTRPNEST_interp_3d    => INTRPNEST_interp_3d_1points
       itp_nh = 1

    case ( 3 )
       INTRPNEST_search_horiz => INTRPNEST_search_horiz_3points
       INTRPNEST_interp_2d    => INTRPNEST_interp_2d_3points
       INTRPNEST_interp_3d    => INTRPNEST_interp_3d_3points
       itp_nh = 3

    case ( 4 )
       INTRPNEST_search_horiz => INTRPNEST_search_horiz_4points
       INTRPNEST_interp_2d    => INTRPNEST_interp_2d_4points
       INTRPNEST_interp_3d    => INTRPNEST_interp_3d_4points
       itp_nh = 4

    case ( 8 )
       INTRPNEST_search_horiz => INTRPNEST_search_horiz_8points
       INTRPNEST_interp_2d    => INTRPNEST_interp_2d_8points
       INTRPNEST_interp_3d    => INTRPNEST_interp_3d_8points
       itp_nh = 8

    case ( 12 )
       INTRPNEST_search_horiz => INTRPNEST_search_horiz_12points
       INTRPNEST_interp_2d    => INTRPNEST_interp_2d_12points
       INTRPNEST_interp_3d    => INTRPNEST_interp_3d_12points
       itp_nh = 12

    case default
       write(*,*) 'xxx invarid NEST_INTERP_LEVEL (', NEST_INTERP_LEVEL, &
                  ') [setup: nest/interp]'
       call PRC_MPIstop
    end select

    if ( OFFLINE ) then
       select_type = "offline"
       INTRPNEST_search_vert => INTRPNEST_search_vert_offline
    else
       select_type = "online"
       INTRPNEST_search_vert => INTRPNEST_search_vert_online
    endif

    if( IO_L ) write(IO_FID_LOG,*) '+++ horizontal interpolation with ', &
                                   NEST_INTERP_LEVEL, " points."
    if( IO_L ) write(IO_FID_LOG,*) '+++ vertical interpolation for ', &
                                   trim(select_type)

    return
  end subroutine INTRPNEST_setup


  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon
  subroutine INTRPNEST_interp_fact_latlon( &
      hfact,      & ! (out)
      igrd,       & ! (out)
      jgrd,       & ! (out)
      mylat,      & ! (in)
      mylon,      & ! (in)
      myIA,       & ! (in)
      myJA,       & ! (in)
      inlat,      & ! (in)
      inlon,      & ! (in)
      inIA,       & ! (in)
      inJA        ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: hfact(:,:,:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:,:,:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:,:,:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat(:,:)         ! latitude data of mine
    real(RP), intent(in)  :: mylon(:,:)         ! longitude data of mine
    integer,  intent(in)  :: myIA               ! grid number of mine
    integer,  intent(in)  :: myJA               ! grid number of mine

    real(RP), intent(in)  :: inlat(:,:)         ! latitude data of you (input)
    real(RP), intent(in)  :: inlon(:,:)         ! longitude data of you (input)
    integer,  intent(in)  :: inIA               ! grid number of you (input)
    integer,  intent(in)  :: inJA               ! grid number of you (input)

    integer  :: i, j
    integer  :: is, ie
    integer  :: js, je
    !---------------------------------------------------------------------------

    hfact(:,:,:) = 0.0_RP

    do j = 1, myJA
    do i = 1, myIA
       ! nearest block search
       call INTRPNEST_search_nearest_block( is, ie, js, je,         &
                                            mylat(i,j), mylon(i,j), &
                                            inlat(:,:), inlon(:,:), &
                                            inIA,  inJA             )

       ! main search
       call INTRPNEST_search_horiz( hfact(i,j,:), &
                                    igrd(i,j,:),  &
                                    jgrd(i,j,:),  &
                                    mylat(i,j),   &
                                    mylon(i,j),   &
                                    inlat,        &
                                    inlon,        &
                                    is, ie,       &
                                    js, je        )
    enddo
    enddo

    return
  end subroutine INTRPNEST_interp_fact_latlon


  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon and Z-Height information
  subroutine INTRPNEST_interp_fact_llz( &
      hfact,      & ! (out)
      vfact,      & ! (out)
      kgrd,       & ! (out)
      igrd,       & ! (out)
      jgrd,       & ! (out)
      ncopy,      & ! (out)
      myhgt,      & ! (in)
      mylat,      & ! (in)
      mylon,      & ! (in)
      myKS,       & ! (in)
      myKE,       & ! (in)
      myIA,       & ! (in)
      myJA,       & ! (in)
      inhgt,      & ! (in)
      inlat,      & ! (in)
      inlon,      & ! (in)
      inKA,       & ! (in)
      inIA,       & ! (in)
      inJA,       & ! (in)
      landgrid    ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: hfact(:,:,:)       ! horizontal interp factor
    real(RP), intent(out) :: vfact(:,:,:,:,:)   ! vertical interp factor
    integer,  intent(out) :: kgrd (:,:,:,:,:)   ! grid points of interp target
    integer,  intent(out) :: igrd (:,:,:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:,:,:)       ! grid points of interp target
    integer,  intent(out) :: ncopy(:,:,:)       ! number of daughter's layers below parent lowest layer

    real(RP), intent(in)  :: myhgt(:,:,:)       ! height data of mine
    real(RP), intent(in)  :: mylat(:,:)         ! latitude data of mine
    real(RP), intent(in)  :: mylon(:,:)         ! longitude data of mine
    integer,  intent(in)  :: myKS               ! start grid number of mine
    integer,  intent(in)  :: myKE               ! end grid number of mine
    integer,  intent(in)  :: myIA               ! grid number of mine
    integer,  intent(in)  :: myJA               ! grid number of mine

    real(RP), intent(in)  :: inhgt(:,:,:)       ! height data of you (input)
    real(RP), intent(in)  :: inlat(:,:)         ! latitude data of you (input)
    real(RP), intent(in)  :: inlon(:,:)         ! longitude data of you (input)
    integer,  intent(in)  :: inKA               ! grid number of you (input)
    integer,  intent(in)  :: inIA               ! grid number of you (input)
    integer,  intent(in)  :: inJA               ! grid number of you (input)

    logical,  intent(in), optional :: landgrid

    integer  :: i, j, k
    integer  :: is, ie
    integer  :: js, je
    logical  :: lndgrd
    !---------------------------------------------------------------------------

    lndgrd = .false.
    if ( present(landgrid) ) then
    if ( landgrid ) then
       lndgrd = .true.
    endif
    endif

    hfact(:,:,:) = 0.0_RP
    vfact(:,:,:,:,:) = 0.0_RP
    ncopy(:,:,:) = 0

    do j = 1, myJA
    do i = 1, myIA

       ! nearest block search
       call INTRPNEST_search_nearest_block( is, ie, js, je,         &
                                            mylat(i,j), mylon(i,j), &
                                            inlat(:,:), inlon(:,:), &
                                            inIA,  inJA             )

       ! main search
       call INTRPNEST_search_horiz( hfact(i,j,:), &
                                    igrd(i,j,:),  &
                                    jgrd(i,j,:),  &
                                    mylat(i,j),   &
                                    mylon(i,j),   &
                                    inlat,        &
                                    inlon,        &
                                    is, ie,       &
                                    js, je        )


       call INTRPNEST_search_vert( vfact,         &
                                   kgrd,          &
                                   ncopy(i,j,:),  &
                                   igrd(i,j,:),   &
                                   jgrd(i,j,:),   &
                                   myhgt(:,i,j),  &
                                   inhgt,         &
                                   i,  j,         &
                                   myKS, myKE,    &
                                   inKA,          &
                                   lndgrd         )

    enddo
    enddo

    return
  end subroutine INTRPNEST_interp_fact_llz


  !-----------------------------------------------------------------------------
  ! search of nearest region for speed up of interpolation
  subroutine INTRPNEST_search_nearest_block( &
      is,      & ! (out)
      ie,      & ! (out)
      js,      & ! (out)
      je,      & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      inIA,    & ! (in)
      inJA     ) ! (in)
    implicit none

    integer,  intent(out)  :: is             ! start index for x-direction
    integer,  intent(out)  :: ie             ! end   index for x-direction
    integer,  intent(out)  :: js             ! start index for y-direction
    integer,  intent(out)  :: je             ! end   index for y-direction

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)
    integer,  intent(in)  :: inIA           ! grid number of you (input)
    integer,  intent(in)  :: inJA           ! grid number of you (input)

    real(RP) :: distance, dist
    integer  :: ii, jj
    integer  :: iinc, jinc
    integer  :: blk_i, blk_j
    !---------------------------------------------------------------------------

    iinc = max( (inIA + 1) / divnum, 1 )
    jinc = max( (inJA + 1) / divnum, 1 )
    dist = large_number_1

    jj = 1 + (jinc/2)
    do while (jj <= inJA)
       ii = 1 + (iinc/2)
       do while (ii <= inIA)
          distance = INTRPNEST_haversine( mylat,        mylon,       &
                                          inlat(ii,jj), inlon(ii,jj) )
          if( distance < dist )then
             dist = distance
             blk_i = ii
             blk_j = jj
          endif
          ii = ii + iinc
       enddo
       jj = jj + jinc
    enddo

    ! +- 3 is buffer for 12 points
    is = blk_i - (iinc/2) - 3
    if( is < 1 ) is = 1
    ie = blk_i + (iinc/2) + 3
    if( ie > inIA ) ie = inIA
    js = blk_j - (jinc/2) - 3
    if( js < 1 ) js = 1
    je = blk_j + (jinc/2) + 3
    if( je > inJA ) je = inJA

    return
  end subroutine INTRPNEST_search_nearest_block


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for one-points (nearest-neighbor)
  subroutine INTRPNEST_search_horiz_1points( &
      hfact,   & ! (out)
      igrd,    & ! (out)
      jgrd,    & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      is,      & ! (in)
      ie,      & ! (in)
      js,      & ! (in)
      je       ) ! (in)
    implicit none

    real(RP), intent(out) :: hfact(:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)

    integer,  intent(in)  :: is             ! start index for x-direction
    integer,  intent(in)  :: ie             ! end   index for x-direction
    integer,  intent(in)  :: js             ! start index for y-direction
    integer,  intent(in)  :: je             ! end   index for y-direction

    real(RP) :: distance
    real(RP) :: dist
    integer :: ii, jj
    !---------------------------------------------------------------------------

    dist    = large_number_1
    igrd(:) = -1
    jgrd(:) = -1

    do jj = js, je
    do ii = is, ie
       distance = INTRPNEST_haversine( mylat,mylon,inlat(ii,jj),inlon(ii,jj) )
       if ( distance <= dist ) then
          dist = distance;  igrd(1) = ii;       jgrd(1) = jj
       endif
    enddo
    enddo

    hfact(1) = 1.0_RP

    return
  end subroutine INTRPNEST_search_horiz_1points


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for three-points
  subroutine INTRPNEST_search_horiz_3points( &
      hfact,   & ! (out)
      igrd,    & ! (out)
      jgrd,    & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      is,      & ! (in)
      ie,      & ! (in)
      js,      & ! (in)
      je       ) ! (in)
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    real(RP), intent(out) :: hfact(:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)

    integer,  intent(in)  :: is             ! start index for x-direction
    integer,  intent(in)  :: ie             ! end   index for x-direction
    integer,  intent(in)  :: js             ! start index for y-direction
    integer,  intent(in)  :: je             ! end   index for y-direction

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(3)
    integer :: ii, jj
    !---------------------------------------------------------------------------

    dist(1) = large_number_3
    dist(2) = large_number_2
    dist(3) = large_number_1
    igrd(:) = -1
    jgrd(:) = -1

    do jj = js, je
    do ii = is, ie
       distance = INTRPNEST_haversine( mylat,mylon,inlat(ii,jj),inlon(ii,jj) )
       if ( distance <= dist(1) ) then
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = dist(1);   igrd(2) = igrd(1);  jgrd(2) = jgrd(1)
          dist(1) = distance;  igrd(1) = ii;       jgrd(1) = jj
       elseif ( dist(1) < distance .and. distance <= dist(2) ) then
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = distance;  igrd(2) = ii;       jgrd(2) = jj
       elseif ( dist(2) < distance .and. distance <= dist(3) ) then
          dist(3) = distance;  igrd(3) = ii;       jgrd(3) = jj
       endif
    enddo
    enddo

    if ( abs(dist(1)) < eps ) then
       hfact(1) = 1.0_RP
       hfact(2) = 0.0_RP
       hfact(3) = 0.0_RP
    else
       denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) + (1.0_RP/dist(3)) )
       hfact(1) = ( 1.0_RP/dist(1) ) * denom
       hfact(2) = ( 1.0_RP/dist(2) ) * denom
       hfact(3) = ( 1.0_RP/dist(3) ) * denom
    endif

    return
  end subroutine INTRPNEST_search_horiz_3points


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for four-points
  subroutine INTRPNEST_search_horiz_4points( &
      hfact,   & ! (out)
      igrd,    & ! (out)
      jgrd,    & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      is,      & ! (in)
      ie,      & ! (in)
      js,      & ! (in)
      je       ) ! (in)
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    real(RP), intent(out) :: hfact(:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)

    integer,  intent(in)  :: is             ! start index for x-direction
    integer,  intent(in)  :: ie             ! end   index for x-direction
    integer,  intent(in)  :: js             ! start index for y-direction
    integer,  intent(in)  :: je             ! end   index for y-direction

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(4)
    integer :: ii, jj
    !---------------------------------------------------------------------------

    dist(1) = large_number_4
    dist(2) = large_number_3
    dist(3) = large_number_2
    dist(4) = large_number_1
    igrd(:) = -1
    jgrd(:) = -1

    do jj = js, je
    do ii = is, ie
       distance = INTRPNEST_haversine( mylat,mylon,inlat(ii,jj),inlon(ii,jj) )
       if ( distance <= dist(1) ) then
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = dist(1);   igrd(2) = igrd(1);  jgrd(2) = jgrd(1)
          dist(1) = distance;  igrd(1) = ii;       jgrd(1) = jj
       elseif ( dist(1) < distance .and. distance <= dist(2) ) then
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = distance;  igrd(2) = ii;       jgrd(2) = jj
       elseif ( dist(2) < distance .and. distance <= dist(3) ) then
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = distance;  igrd(3) = ii;       jgrd(3) = jj
       elseif ( dist(3) < distance .and. distance <= dist(4) ) then
          dist(4) = distance;  igrd(4) = ii;       jgrd(4) = jj
       endif
    enddo
    enddo

    if ( abs(dist(1)) < eps ) then
       hfact(1) = 1.0_RP
       hfact(2) = 0.0_RP
       hfact(3) = 0.0_RP
       hfact(4) = 0.0_RP
    else
       denom = 1.0_RP / (  (1.0_RP/dist(1)) + (1.0_RP/dist(2)) &
                         + (1.0_RP/dist(3)) + (1.0_RP/dist(4)) )
       hfact(1) = ( 1.0_RP/dist(1) ) * denom
       hfact(2) = ( 1.0_RP/dist(2) ) * denom
       hfact(3) = ( 1.0_RP/dist(3) ) * denom
       hfact(4) = ( 1.0_RP/dist(4) ) * denom
    endif

    return
  end subroutine INTRPNEST_search_horiz_4points


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for eight-points
  subroutine INTRPNEST_search_horiz_8points( &
      hfact,   & ! (out)
      igrd,    & ! (out)
      jgrd,    & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      is,      & ! (in)
      ie,      & ! (in)
      js,      & ! (in)
      je       ) ! (in)
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    real(RP), intent(out) :: hfact(:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)

    integer,  intent(in)  :: is             ! start index for x-direction
    integer,  intent(in)  :: ie             ! end   index for x-direction
    integer,  intent(in)  :: js             ! start index for y-direction
    integer,  intent(in)  :: je             ! end   index for y-direction

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(8)
    integer :: ii, jj
    !---------------------------------------------------------------------------

    dist(1) = large_number_8
    dist(2) = large_number_7
    dist(3) = large_number_6
    dist(4) = large_number_5
    dist(5) = large_number_4
    dist(6) = large_number_3
    dist(7) = large_number_2
    dist(8) = large_number_1
    igrd(:) = -1
    jgrd(:) = -1

    do jj = js, je
    do ii = is, ie
       distance = INTRPNEST_haversine( mylat,mylon,inlat(ii,jj),inlon(ii,jj) )
       if ( distance <= dist(1) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = dist(5);   igrd(6) = igrd(5);  jgrd(6) = jgrd(5)
          dist(5) = dist(4);   igrd(5) = igrd(4);  jgrd(5) = jgrd(4)
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = dist(1);   igrd(2) = igrd(1);  jgrd(2) = jgrd(1)
          dist(1) = distance;  igrd(1) = ii;       jgrd(1) = jj
       elseif ( dist(1) < distance .and. distance <= dist(2) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = dist(5);   igrd(6) = igrd(5);  jgrd(6) = jgrd(5)
          dist(5) = dist(4);   igrd(5) = igrd(4);  jgrd(5) = jgrd(4)
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = dist(2);   igrd(3) = igrd(2);  jgrd(3) = jgrd(2)
          dist(2) = distance;  igrd(2) = ii;       jgrd(2) = jj
       elseif ( dist(2) < distance .and. distance <= dist(3) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = dist(5);   igrd(6) = igrd(5);  jgrd(6) = jgrd(5)
          dist(5) = dist(4);   igrd(5) = igrd(4);  jgrd(5) = jgrd(4)
          dist(4) = dist(3);   igrd(4) = igrd(3);  jgrd(4) = jgrd(3)
          dist(3) = distance;  igrd(3) = ii;       jgrd(3) = jj
       elseif ( dist(3) < distance .and. distance <= dist(4) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = dist(5);   igrd(6) = igrd(5);  jgrd(6) = jgrd(5)
          dist(5) = dist(4);   igrd(5) = igrd(4);  jgrd(5) = jgrd(4)
          dist(4) = distance;  igrd(4) = ii;       jgrd(4) = jj
       elseif ( dist(4) < distance .and. distance <= dist(5) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = dist(5);   igrd(6) = igrd(5);  jgrd(6) = jgrd(5)
          dist(5) = distance;  igrd(5) = ii;       jgrd(5) = jj
       elseif ( dist(5) < distance .and. distance <= dist(6) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = dist(6);   igrd(7) = igrd(6);  jgrd(7) = jgrd(6)
          dist(6) = distance;  igrd(6) = ii;       jgrd(6) = jj
       elseif ( dist(6) < distance .and. distance <= dist(7) ) then
          dist(8) = dist(7);   igrd(8) = igrd(7);  jgrd(8) = jgrd(7)
          dist(7) = distance;  igrd(7) = ii;       jgrd(7) = jj
       elseif ( dist(7) < distance .and. distance <= dist(8) ) then
          dist(8) = distance;  igrd(8) = ii;       jgrd(8) = jj
       endif
    enddo
    enddo

    if ( abs(dist(1)) < eps ) then
       hfact(:) = 0.0_RP
       hfact(1) = 1.0_RP
    else
       denom = 1.0_RP / (   (1.0_RP/dist(1)) + (1.0_RP/dist(2)) &
                          + (1.0_RP/dist(3)) + (1.0_RP/dist(4)) &
                          + (1.0_RP/dist(5)) + (1.0_RP/dist(6)) &
                          + (1.0_RP/dist(7)) + (1.0_RP/dist(8)) )
       hfact(1) = ( 1.0_RP/dist(1) ) * denom
       hfact(2) = ( 1.0_RP/dist(2) ) * denom
       hfact(3) = ( 1.0_RP/dist(3) ) * denom
       hfact(4) = ( 1.0_RP/dist(4) ) * denom
       hfact(5) = ( 1.0_RP/dist(5) ) * denom
       hfact(6) = ( 1.0_RP/dist(6) ) * denom
       hfact(7) = ( 1.0_RP/dist(7) ) * denom
       hfact(8) = ( 1.0_RP/dist(8) ) * denom
    endif

    return
  end subroutine INTRPNEST_search_horiz_8points


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for twelve-points
  subroutine INTRPNEST_search_horiz_12points( &
      hfact,   & ! (out)
      igrd,    & ! (out)
      jgrd,    & ! (out)
      mylat,   & ! (in)
      mylon,   & ! (in)
      inlat,   & ! (in)
      inlon,   & ! (in)
      is,      & ! (in)
      ie,      & ! (in)
      js,      & ! (in)
      je       ) ! (in)
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    real(RP), intent(out) :: hfact(:)       ! horizontal interp factor
    integer,  intent(out) :: igrd (:)       ! grid points of interp target
    integer,  intent(out) :: jgrd (:)       ! grid points of interp target

    real(RP), intent(in)  :: mylat          ! latitude data of mine
    real(RP), intent(in)  :: mylon          ! longitude data of mine
    real(RP), intent(in)  :: inlat(:,:)     ! latitude  data of you (input)
    real(RP), intent(in)  :: inlon(:,:)     ! longitude data of you (input)

    integer,  intent(in)  :: is             ! start index for x-direction
    integer,  intent(in)  :: ie             ! end   index for x-direction
    integer,  intent(in)  :: js             ! start index for y-direction
    integer,  intent(in)  :: je             ! end   index for y-direction

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(12)
    integer :: ii, jj
    !---------------------------------------------------------------------------

    dist(1 ) = large_number_12
    dist(2 ) = large_number_11
    dist(3 ) = large_number_10
    dist(4 ) = large_number_9
    dist(5 ) = large_number_8
    dist(6 ) = large_number_7
    dist(7 ) = large_number_6
    dist(8 ) = large_number_5
    dist(9 ) = large_number_4
    dist(10) = large_number_3
    dist(11) = large_number_2
    dist(12) = large_number_1
    igrd(:) = -1
    jgrd(:) = -1

    do jj = js, je
    do ii = is, ie
       distance = INTRPNEST_haversine( mylat,mylon,inlat(ii,jj),inlon(ii,jj) )
       if ( distance <= dist(1) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = dist(5 );   igrd(6 ) = igrd(5 );  jgrd(6 ) = jgrd(5 )
          dist(5 ) = dist(4 );   igrd(5 ) = igrd(4 );  jgrd(5 ) = jgrd(4 )
          dist(4 ) = dist(3 );   igrd(4 ) = igrd(3 );  jgrd(4 ) = jgrd(3 )
          dist(3 ) = dist(2 );   igrd(3 ) = igrd(2 );  jgrd(3 ) = jgrd(2 )
          dist(2 ) = dist(1 );   igrd(2 ) = igrd(1 );  jgrd(2 ) = jgrd(1 )
          dist(1 ) = distance;   igrd(1 ) = ii;        jgrd(1 ) = jj
       elseif ( dist(1) < distance .and. distance <= dist(2) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = dist(5 );   igrd(6 ) = igrd(5 );  jgrd(6 ) = jgrd(5 )
          dist(5 ) = dist(4 );   igrd(5 ) = igrd(4 );  jgrd(5 ) = jgrd(4 )
          dist(4 ) = dist(3 );   igrd(4 ) = igrd(3 );  jgrd(4 ) = jgrd(3 )
          dist(3 ) = dist(2 );   igrd(3 ) = igrd(2 );  jgrd(3 ) = jgrd(2 )
          dist(2 ) = distance;   igrd(2 ) = ii;        jgrd(2 ) = jj
       elseif ( dist(2) < distance .and. distance <= dist(3) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = dist(5 );   igrd(6 ) = igrd(5 );  jgrd(6 ) = jgrd(5 )
          dist(5 ) = dist(4 );   igrd(5 ) = igrd(4 );  jgrd(5 ) = jgrd(4 )
          dist(4 ) = dist(3 );   igrd(4 ) = igrd(3 );  jgrd(4 ) = jgrd(3 )
          dist(3 ) = distance;   igrd(3 ) = ii;        jgrd(3 ) = jj
       elseif ( dist(3) < distance .and. distance <= dist(4) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = dist(5 );   igrd(6 ) = igrd(5 );  jgrd(6 ) = jgrd(5 )
          dist(5 ) = dist(4 );   igrd(5 ) = igrd(4 );  jgrd(5 ) = jgrd(4 )
          dist(4 ) = distance;   igrd(4 ) = ii;        jgrd(4 ) = jj
       elseif ( dist(4) < distance .and. distance <= dist(5) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = dist(5 );   igrd(6 ) = igrd(5 );  jgrd(6 ) = jgrd(5 )
          dist(5 ) = distance;   igrd(5 ) = ii;        jgrd(5 ) = jj
       elseif ( dist(5) < distance .and. distance <= dist(6) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = dist(6 );   igrd(7 ) = igrd(6 );  jgrd(7 ) = jgrd(6 )
          dist(6 ) = distance;   igrd(6 ) = ii;        jgrd(6 ) = jj
       elseif ( dist(6) < distance .and. distance <= dist(7) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = dist(7 );   igrd(8 ) = igrd(7 );  jgrd(8 ) = jgrd(7 )
          dist(7 ) = distance;   igrd(7 ) = ii;        jgrd(7 ) = jj
       elseif ( dist(7) < distance .and. distance <= dist(8) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = dist(8 );   igrd(9 ) = igrd(8 );  jgrd(9 ) = jgrd(8 )
          dist(8 ) = distance;   igrd(8 ) = ii;        jgrd(8 ) = jj
       elseif ( dist(8) < distance .and. distance <= dist(9) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = dist(9 );   igrd(10) = igrd(9 );  jgrd(10) = jgrd(9 )
          dist(9 ) = distance;   igrd(9 ) = ii;        jgrd(9 ) = jj
       elseif ( dist(9) < distance .and. distance <= dist(10) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = dist(10);   igrd(11) = igrd(10);  jgrd(11) = jgrd(10)
          dist(10) = distance;   igrd(10) = ii;        jgrd(10) = jj
       elseif ( dist(10) < distance .and. distance <= dist(11) ) then
          dist(12) = dist(11);   igrd(12) = igrd(11);  jgrd(12) = jgrd(11)
          dist(11) = distance;   igrd(11) = ii;        jgrd(11) = jj
       elseif ( dist(11) < distance .and. distance <= dist(12) ) then
          dist(12) = distance;   igrd(12) = ii;        jgrd(12) = jj
       endif
    enddo
    enddo

    if ( abs(dist(1)) < eps ) then
       hfact(:) = 0.0_RP
       hfact(1) = 1.0_RP
    else
       denom = 1.0_RP / (  (1.0_RP/dist(1 )) + (1.0_RP/dist(2 )) &
                         + (1.0_RP/dist(3 )) + (1.0_RP/dist(4 )) &
                         + (1.0_RP/dist(5 )) + (1.0_RP/dist(6 )) &
                         + (1.0_RP/dist(7 )) + (1.0_RP/dist(8 )) &
                         + (1.0_RP/dist(9 )) + (1.0_RP/dist(10)) &
                         + (1.0_RP/dist(11)) + (1.0_RP/dist(12)) )
       hfact(1 ) = ( 1.0_RP/dist(1 ) ) * denom
       hfact(2 ) = ( 1.0_RP/dist(2 ) ) * denom
       hfact(3 ) = ( 1.0_RP/dist(3 ) ) * denom
       hfact(4 ) = ( 1.0_RP/dist(4 ) ) * denom
       hfact(5 ) = ( 1.0_RP/dist(5 ) ) * denom
       hfact(6 ) = ( 1.0_RP/dist(6 ) ) * denom
       hfact(7 ) = ( 1.0_RP/dist(7 ) ) * denom
       hfact(8 ) = ( 1.0_RP/dist(8 ) ) * denom
       hfact(9 ) = ( 1.0_RP/dist(9 ) ) * denom
       hfact(10) = ( 1.0_RP/dist(10) ) * denom
       hfact(11) = ( 1.0_RP/dist(11) ) * denom
       hfact(12) = ( 1.0_RP/dist(12) ) * denom
    endif

    return
  end subroutine INTRPNEST_search_horiz_12points


  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points (online)
  subroutine INTRPNEST_search_vert_online( &
      vfact,   & ! (inout)
      kgrd,    & ! (inout)
      ncopy,   & ! (out)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      myhgt,   & ! (in)
      inhgt,   & ! (in)
      iloc,    & ! (in)
      jloc,    & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      inKA,    & ! (in)
      lndgrd   ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    real(RP), intent(inout) :: vfact(:,:,:,:,:) ! vertical interp factor
    integer,  intent(inout) :: kgrd (:,:,:,:,:) ! grid points of interp target
    integer,  intent(out) :: ncopy(:)           ! number of daughter's layers below inKS

    integer,  intent(in)  :: igrd(:)            ! grid points of interp target
    integer,  intent(in)  :: jgrd(:)            ! grid points of interp target
    real(RP), intent(in)  :: myhgt(:)           ! height data of mine
    real(RP), intent(in)  :: inhgt(:,:,:)       ! height data of you (input)
    integer,  intent(in)  :: iloc               ! locator index for x-direction
    integer,  intent(in)  :: jloc               ! locator index for y-direction
    integer,  intent(in)  :: ks                 ! start index for z-direction
    integer,  intent(in)  :: ke                 ! end   index for z-direction
    integer,  intent(in)  :: inKA               ! grid number of you (input)
    logical,  intent(in)  :: lndgrd             ! flag of land grid

    real(RP) :: distance
    real(RP) :: dist(2)
    logical  :: dflag                           ! flag: data found or not
    integer  :: ii, jj, idx
    integer  :: k, kk
    integer  :: inKS, inKE
    logical  :: copy
    !---------------------------------------------------------------------------

    if ( lndgrd ) then
       write(*,*) 'xxx internal error [interporation: nest/interp]'
       write(*,*) '    land grid is not araviable in online'
       call PRC_MPIstop
    endif

    inKS = 1 + KHALO
    inKE = inKA - KHALO

    do idx = 1, itp_nh
       ii = igrd(idx)
       jj = jgrd(idx)

       do k = ks, ke
          dist(1) = large_number_2
          dist(2) = large_number_1
          copy    = .false.

          kgrd(k,iloc,jloc,idx,:) = -1

          dflag   = .false.

          if( myhgt(k) < inhgt(inKS,ii,jj) ) then
             copy       = .true.
             ncopy(idx) = ncopy(idx) + 1

             kgrd(k,iloc,jloc,idx,:)  = inKS

             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP

             dflag = .true.
          else
             do kk = inKS, inKE
                dist(1) = myhgt(k) - inhgt(kk  ,ii,jj)
                dist(2) = myhgt(k) - inhgt(kk+1,ii,jj)

                if( dist(1) >= 0.0_RP .and. dist(2) < 0.0_RP ) then
                   kgrd(k,iloc,jloc,idx,1) = kk
                   kgrd(k,iloc,jloc,idx,2) = kk+1

                   vfact(k,iloc,jloc,idx,1) = abs(dist(2)) / ( abs(dist(1)) + abs(dist(2)) )
                   vfact(k,iloc,jloc,idx,2) = abs(dist(1)) / ( abs(dist(1)) + abs(dist(2)) )

                   dflag = .true.

                   exit ! loop end
                endif
             enddo
          endif

          if( .not. dflag ) then
             write(*,*) 'xxx internal error [INTRPNEST_search_vert_online]'
             write(*,*) 'xxx data for interpolation was not found.'
             write(*,*) 'xxx iloc=',iloc,' jloc=',jloc,' k=',k,' idx=',idx
             call PRC_MPIstop
          endif
       enddo

    enddo

    return
  end subroutine INTRPNEST_search_vert_online

  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points (offline)
  subroutine INTRPNEST_search_vert_offline( &
      vfact,   & ! (inout)
      kgrd,    & ! (inout)
      ncopy,   & ! (out)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      myhgt,   & ! (in)
      inhgt,   & ! (in)
      iloc,    & ! (in)
      jloc,    & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      inKA,    & ! (in)
      lndgrd   ) ! (in)
    use scale_const, only: &
       eps => CONST_EPS
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: vfact(:,:,:,:,:) ! vertical interp factor
    integer,  intent(inout) :: kgrd (:,:,:,:,:) ! grid points of interp target
    integer,  intent(out) :: ncopy(:)           ! number of daughter's layers below inKS

    integer,  intent(in)  :: igrd(:)            ! grid points of interp target
    integer,  intent(in)  :: jgrd(:)            ! grid points of interp target
    real(RP), intent(in)  :: myhgt(:)           ! height data of mine
    real(RP), intent(in)  :: inhgt(:,:,:)       ! height data of you (input)
    integer,  intent(in)  :: iloc               ! locator index for x-direction
    integer,  intent(in)  :: jloc               ! locator index for y-direction
    integer,  intent(in)  :: ks                 ! start index for z-direction
    integer,  intent(in)  :: ke                 ! end   index for z-direction
    integer,  intent(in)  :: inKA               ! grid number of you (input)
    logical,  intent(in)  :: lndgrd             ! flag of land grid

    real(RP) :: distance
    real(RP) :: denom
    real(RP) :: dist(2)
    logical  :: dflag                           ! flag: data found or not
    integer  :: ii, jj
    integer  :: k, kk, kks, kke
    integer  :: idx
    !---------------------------------------------------------------------------

    if ( lndgrd ) then
       kks = 1;     kke = LKMAX
    else
       kks = ks;  kke = ke
    endif

    ncopy = 0  ! dummy
    do idx = 1, itp_nh
       ii = igrd(idx)
       jj = jgrd(idx)

       do k = kks, kke
          dist(1) = large_number_2
          dist(2) = large_number_1
          kgrd(k,iloc,jloc,idx,:) = -1
          dflag = .false.

          if( myhgt(k) < inhgt(1,ii,jj) )then ! copy
             kgrd(k,iloc,jloc,idx,:)  = 1
             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP
             dflag = .true.
          else if( abs(inhgt(inKA,ii,jj)-myhgt(k))<eps )then
             kgrd(k,iloc,jloc,idx,:)  = inKA
             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP
             dflag = .true.
          else if( inhgt(inKA,ii,jj) < myhgt(k) )then
             write(*,*) 'xxx internal error [INTRPNEST_search_vert_offline]'
             write(*,*) 'xxx data level is beyond parent data'
             write(*,*) 'in',ii,jj,inKA,inhgt(inKA,ii,jj),'my',k,myhgt(k)
             call PRC_MPIstop
          else

             do kk = 1, inKA-1
                if( (inhgt(kk,ii,jj)<=myhgt(k)).and.(myhgt(k)<inhgt(kk+1,ii,jj)) )then
                   kgrd(k,iloc,jloc,idx,1) = kk
                   kgrd(k,iloc,jloc,idx,2) = kk+1
                   dist(1) = abs( myhgt(k) - inhgt(kk,ii,jj)   )
                   dist(2) = abs( myhgt(k) - inhgt(kk+1,ii,jj) )
                   dflag = .true.
                   if ( abs(dist(1))<eps )then
                      vfact(k,iloc,jloc,idx,1) = 1.0_RP
                      vfact(k,iloc,jloc,idx,2) = 0.0_RP
                   else
                      denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
                      vfact(k,iloc,jloc,idx,1) = ( 1.0_RP/dist(1) ) * denom
                      vfact(k,iloc,jloc,idx,2) = ( 1.0_RP/dist(2) ) * denom
                   endif
                   exit
                endif
             enddo
          endif

          if( .not. dflag )then
             write(*,*) 'xxx internal error [INTRPNEST_search_vert_offline]'
             write(*,*) 'xxx data for interpolation was not found.'
             write(*,*) 'xxx iloc=',iloc,' jloc=',jloc,' k=',k,' idx=',idx
             call PRC_MPIstop
          endif
       enddo

    enddo

    return
  end subroutine INTRPNEST_search_vert_offline

  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 2D data (nearest-neighbor)
  subroutine INTRPNEST_interp_2d_1points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja       ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: ia             ! grid number of mine
    integer,  intent(in)  :: ja             ! grid number of mine

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, ja
    do i = 1, ia
       intp(i,j) = ref(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1)
    end do
    end do

    return
  end subroutine INTRPNEST_interp_2d_1points


  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 3D data (nearest-neighbor)
  subroutine INTRPNEST_interp_3d_1points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      vfact,   & ! (in)
      kgrd,    & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja,      & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      logwegt  ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:,:)       ! interpolated data

    real(RP), intent(in)  :: ref (:,:,:)       ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
    integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine
    integer,  intent(in)  :: ks                ! start grid number of mine
    integer,  intent(in)  :: ke                ! end grid number of mine

    logical,  intent(in), optional :: logwegt

    integer :: i, j, k
    logical :: logarithmic
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(logwegt) ) then
    if ( logwegt ) then
       logarithmic = .true.
    endif
    endif

       ! linear interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,1)              &
                      + ref(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,2)
       end do
       end do
       end do

    ! logarithmic weighting (for pres, dens)
    if ( logarithmic ) then
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( intp(k,i,j) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_1points


  !-----------------------------------------------------------------------------
  ! interpolation using three-points for 2D data
  subroutine INTRPNEST_interp_2d_3points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja       ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: ia             ! grid number of mine
    integer,  intent(in)  :: ja             ! grid number of mine

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, ja
    do i = 1, ia
       intp(i,j) = ref(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1)  &
                 + ref(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2)  &
                 + ref(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)
    end do
    end do

    return
  end subroutine INTRPNEST_interp_2d_3points


  !-----------------------------------------------------------------------------
  ! interpolation using three-points for 3D data
  subroutine INTRPNEST_interp_3d_3points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      vfact,   & ! (in)
      kgrd,    & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja,      & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      logwegt  ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:,:)       ! interpolated data

    real(RP), intent(in)  :: ref (:,:,:)       ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
    integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine
    integer,  intent(in)  :: ks                ! start grid number of mine
    integer,  intent(in)  :: ke                ! end grid number of mine

    logical,  intent(in), optional :: logwegt

    integer :: i, j, k
    logical :: logarithmic
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(logwegt) ) then
    if ( logwegt ) then
       logarithmic = .true.
    endif
    endif

       ! linear interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,1)              &
                      + ref(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,1)              &
                      + ref(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,1)              &
                      + ref(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,2)              &
                      + ref(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,2)              &
                      + ref(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,2)
       end do
       end do
       end do

    ! logarithmic weighting (for pres, dens)
    if ( logarithmic ) then
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( intp(k,i,j) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_3points


  !-----------------------------------------------------------------------------
  ! interpolation using four-points for 2D data
  subroutine INTRPNEST_interp_2d_4points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja       ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: ia             ! grid number of mine
    integer,  intent(in)  :: ja             ! grid number of mine

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, ja
    do i = 1, ia
       intp(i,j) = ref(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1)  &
                 + ref(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2)  &
                 + ref(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)  &
                 + ref(igrd(i,j,4),jgrd(i,j,4)) * hfact(i,j,4)
    end do
    end do

    return
  end subroutine INTRPNEST_interp_2d_4points


  !-----------------------------------------------------------------------------
  ! interpolation using four-points for 3D data
  subroutine INTRPNEST_interp_3d_4points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      vfact,   & ! (in)
      kgrd,    & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja,      & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      logwegt  ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:,:)       ! interpolated data

    real(RP), intent(in)  :: ref (:,:,:)       ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
    integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine
    integer,  intent(in)  :: ks                ! start grid number of mine
    integer,  intent(in)  :: ke                ! end grid number of mine

    logical,  intent(in), optional :: logwegt

    integer :: i, j, k
    logical :: logarithmic
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(logwegt) ) then
    if ( logwegt ) then
       logarithmic = .true.
    endif
    endif


       ! linear interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,1)              &
                      + ref(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,1)              &
                      + ref(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,1)              &
                      + ref(kgrd(k,i,j,4,1),igrd(i,j,4),jgrd(i,j,4)) &
                      * hfact(i,j,4) * vfact(k,i,j,4,1)              &
                      + ref(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,2)              &
                      + ref(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,2)              &
                      + ref(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,2)              &
                      + ref(kgrd(k,i,j,4,2),igrd(i,j,4),jgrd(i,j,4)) &
                      * hfact(i,j,4) * vfact(k,i,j,4,2)
       end do
       end do
       end do

    ! logarithmic weighting (for pres, dens)
    if ( logarithmic ) then
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( intp(k,i,j) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_4points


  !-----------------------------------------------------------------------------
  ! interpolation using eight-points for 2D data
  subroutine INTRPNEST_interp_2d_8points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja       ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: ia             ! grid number of mine
    integer,  intent(in)  :: ja             ! grid number of mine

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, ja
    do i = 1, ia
       intp(i,j) = ref(igrd(i,j,1),jgrd(i,j,1)) * hfact(i,j,1)  &
                 + ref(igrd(i,j,2),jgrd(i,j,2)) * hfact(i,j,2)  &
                 + ref(igrd(i,j,3),jgrd(i,j,3)) * hfact(i,j,3)  &
                 + ref(igrd(i,j,4),jgrd(i,j,4)) * hfact(i,j,4)  &
                 + ref(igrd(i,j,5),jgrd(i,j,5)) * hfact(i,j,5)  &
                 + ref(igrd(i,j,6),jgrd(i,j,6)) * hfact(i,j,6)  &
                 + ref(igrd(i,j,7),jgrd(i,j,7)) * hfact(i,j,7)  &
                 + ref(igrd(i,j,8),jgrd(i,j,8)) * hfact(i,j,8)
    end do
    end do

    return
  end subroutine INTRPNEST_interp_2d_8points


  !-----------------------------------------------------------------------------
  ! interpolation using eight-points for 3D data
  subroutine INTRPNEST_interp_3d_8points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      vfact,   & ! (in)
      kgrd,    & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja,      & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      logwegt  ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:,:)       ! interpolated data

    real(RP), intent(in)  :: ref (:,:,:)       ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
    integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine
    integer,  intent(in)  :: ks                ! start grid number of mine
    integer,  intent(in)  :: ke                ! end grid number of mine

    logical,  intent(in), optional :: logwegt

    integer :: i, j, k
    logical :: logarithmic
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(logwegt) ) then
    if ( logwegt ) then
       logarithmic = .true.
    endif
    endif

       ! linear interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,1)              &
                      + ref(kgrd(k,i,j,2,1),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,1)              &
                      + ref(kgrd(k,i,j,3,1),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,1)              &
                      + ref(kgrd(k,i,j,4,1),igrd(i,j,4),jgrd(i,j,4)) &
                      * hfact(i,j,4) * vfact(k,i,j,4,1)              &
                      + ref(kgrd(k,i,j,5,1),igrd(i,j,5),jgrd(i,j,5)) &
                      * hfact(i,j,5) * vfact(k,i,j,5,1)              &
                      + ref(kgrd(k,i,j,6,1),igrd(i,j,6),jgrd(i,j,6)) &
                      * hfact(i,j,6) * vfact(k,i,j,6,1)              &
                      + ref(kgrd(k,i,j,7,1),igrd(i,j,7),jgrd(i,j,7)) &
                      * hfact(i,j,7) * vfact(k,i,j,7,1)              &
                      + ref(kgrd(k,i,j,8,1),igrd(i,j,8),jgrd(i,j,8)) &
                      * hfact(i,j,8) * vfact(k,i,j,8,1)              &
                      + ref(kgrd(k,i,j,1,2),igrd(i,j,1),jgrd(i,j,1)) &
                      * hfact(i,j,1) * vfact(k,i,j,1,2)              &
                      + ref(kgrd(k,i,j,2,2),igrd(i,j,2),jgrd(i,j,2)) &
                      * hfact(i,j,2) * vfact(k,i,j,2,2)              &
                      + ref(kgrd(k,i,j,3,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,3) * vfact(k,i,j,3,2)              &
                      + ref(kgrd(k,i,j,4,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,4) * vfact(k,i,j,4,2)              &
                      + ref(kgrd(k,i,j,5,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,5) * vfact(k,i,j,5,2)              &
                      + ref(kgrd(k,i,j,6,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,6) * vfact(k,i,j,6,2)              &
                      + ref(kgrd(k,i,j,7,2),igrd(i,j,3),jgrd(i,j,3)) &
                      * hfact(i,j,7) * vfact(k,i,j,7,2)              &
                      + ref(kgrd(k,i,j,8,2),igrd(i,j,8),jgrd(i,j,8)) &
                      * hfact(i,j,8) * vfact(k,i,j,8,2)
       end do
       end do
       end do

    ! logarithmic weighting (for pres, dens)
    if ( logarithmic ) then
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( intp(k,i,j) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_8points


  !-----------------------------------------------------------------------------
  ! interpolation using twelve-points for 2D data
  subroutine INTRPNEST_interp_2d_12points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja       ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)   ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)   ! grid points of interp target
    integer,  intent(in)  :: ia             ! grid number of mine
    integer,  intent(in)  :: ja             ! grid number of mine

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, ja
    do i = 1, ia
       intp(i,j) = ref(igrd(i,j,1), jgrd(i,j,1))  * hfact(i,j,1)  &
                 + ref(igrd(i,j,2), jgrd(i,j,2))  * hfact(i,j,2)  &
                 + ref(igrd(i,j,3), jgrd(i,j,3))  * hfact(i,j,3)  &
                 + ref(igrd(i,j,4), jgrd(i,j,4))  * hfact(i,j,4)  &
                 + ref(igrd(i,j,5), jgrd(i,j,5))  * hfact(i,j,5)  &
                 + ref(igrd(i,j,6), jgrd(i,j,6))  * hfact(i,j,6)  &
                 + ref(igrd(i,j,7), jgrd(i,j,7))  * hfact(i,j,7)  &
                 + ref(igrd(i,j,8), jgrd(i,j,8))  * hfact(i,j,8)  &
                 + ref(igrd(i,j,9), jgrd(i,j,9))  * hfact(i,j,9)  &
                 + ref(igrd(i,j,10),jgrd(i,j,10)) * hfact(i,j,10) &
                 + ref(igrd(i,j,11),jgrd(i,j,11)) * hfact(i,j,11) &
                 + ref(igrd(i,j,12),jgrd(i,j,12)) * hfact(i,j,12)
    end do
    end do

    return
  end subroutine INTRPNEST_interp_2d_12points


  !-----------------------------------------------------------------------------
  ! interpolation using twelve-points for 3D data
  subroutine INTRPNEST_interp_3d_12points( &
      intp,    & ! (out)
      ref,     & ! (in)
      hfact,   & ! (in)
      vfact,   & ! (in)
      kgrd,    & ! (in)
      igrd,    & ! (in)
      jgrd,    & ! (in)
      ia,      & ! (in)
      ja,      & ! (in)
      ks,      & ! (in)
      ke,      & ! (in)
      logwegt  ) ! (in)
    implicit none

    real(RP), intent(out) :: intp(:,:,:)      ! interpolated data

    real(RP), intent(in)  :: ref (:,:,:)      ! reference data
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    real(RP), intent(in)  :: vfact(:,:,:,:,:)  ! vertical interp factor
    integer,  intent(in)  :: kgrd (:,:,:,:,:)  ! grid points of interp target
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine
    integer,  intent(in)  :: ks                ! start grid number of mine
    integer,  intent(in)  :: ke                ! end grid number of mine

    logical,  intent(in), optional :: logwegt

    integer :: i, j, k
    logical :: logarithmic
    !---------------------------------------------------------------------------

    logarithmic = .false.
    if ( present(logwegt) ) then
    if ( logwegt ) then
       logarithmic = .true.
    endif
    endif

       ! linear interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = ref(kgrd(k,i,j,1, 1),igrd(i,j,1 ),jgrd(i,j,1 )) &
                      * hfact(i,j,1 ) * vfact(k,i,j,1, 1)               &
                      + ref(kgrd(k,i,j,2, 1),igrd(i,j,2 ),jgrd(i,j,2 )) &
                      * hfact(i,j,2 ) * vfact(k,i,j,2, 1)               &
                      + ref(kgrd(k,i,j,3, 1),igrd(i,j,3 ),jgrd(i,j,3 )) &
                      * hfact(i,j,3 ) * vfact(k,i,j,3, 1)               &
                      + ref(kgrd(k,i,j,4, 1),igrd(i,j,4 ),jgrd(i,j,4 )) &
                      * hfact(i,j,4 ) * vfact(k,i,j,4, 1)               &
                      + ref(kgrd(k,i,j,5, 1),igrd(i,j,5 ),jgrd(i,j,5 )) &
                      * hfact(i,j,5 ) * vfact(k,i,j,5, 1)               &
                      + ref(kgrd(k,i,j,6, 1),igrd(i,j,6 ),jgrd(i,j,6 )) &
                      * hfact(i,j,6 ) * vfact(k,i,j,6, 1)               &
                      + ref(kgrd(k,i,j,7, 1),igrd(i,j,7 ),jgrd(i,j,7 )) &
                      * hfact(i,j,7 ) * vfact(k,i,j,7, 1)               &
                      + ref(kgrd(k,i,j,8, 1),igrd(i,j,8 ),jgrd(i,j,8 )) &
                      * hfact(i,j,8 ) * vfact(k,i,j,8, 1)               &
                      + ref(kgrd(k,i,j,9, 1),igrd(i,j,9 ),jgrd(i,j,9 )) &
                      * hfact(i,j,9 ) * vfact(k,i,j,9, 1)               &
                      + ref(kgrd(k,i,j,10,1),igrd(i,j,10),jgrd(i,j,10)) &
                      * hfact(i,j,10) * vfact(k,i,j,10,1)               &
                      + ref(kgrd(k,i,j,11,1),igrd(i,j,11),jgrd(i,j,11)) &
                      * hfact(i,j,11) * vfact(k,i,j,11,1)               &
                      + ref(kgrd(k,i,j,12,1),igrd(i,j,12),jgrd(i,j,12)) &
                      * hfact(i,j,12) * vfact(k,i,j,12,1)               &
                      + ref(kgrd(k,i,j,1, 2),igrd(i,j,1 ),jgrd(i,j,1 )) &
                      * hfact(i,j,1 ) * vfact(k,i,j,1, 2)               &
                      + ref(kgrd(k,i,j,2, 2),igrd(i,j,2 ),jgrd(i,j,2 )) &
                      * hfact(i,j,2 ) * vfact(k,i,j,2, 2)               &
                      + ref(kgrd(k,i,j,3, 2),igrd(i,j,3 ),jgrd(i,j,3 )) &
                      * hfact(i,j,3 ) * vfact(k,i,j,3, 2)               &
                      + ref(kgrd(k,i,j,4, 2),igrd(i,j,4 ),jgrd(i,j,4 )) &
                      * hfact(i,j,4 ) * vfact(k,i,j,4, 2)               &
                      + ref(kgrd(k,i,j,5, 2),igrd(i,j,5 ),jgrd(i,j,5 )) &
                      * hfact(i,j,5 ) * vfact(k,i,j,5, 2)               &
                      + ref(kgrd(k,i,j,6, 2),igrd(i,j,6 ),jgrd(i,j,6 )) &
                      * hfact(i,j,6 ) * vfact(k,i,j,6, 2)               &
                      + ref(kgrd(k,i,j,7, 2),igrd(i,j,7 ),jgrd(i,j,7 )) &
                      * hfact(i,j,7 ) * vfact(k,i,j,7, 2)               &
                      + ref(kgrd(k,i,j,8, 2),igrd(i,j,8 ),jgrd(i,j,8 )) &
                      * hfact(i,j,8 ) * vfact(k,i,j,8, 2)               &
                      + ref(kgrd(k,i,j,9, 2),igrd(i,j,9 ),jgrd(i,j,9 )) &
                      * hfact(i,j,9 ) * vfact(k,i,j,9, 2)               &
                      + ref(kgrd(k,i,j,10,2),igrd(i,j,10),jgrd(i,j,10)) &
                      * hfact(i,j,10) * vfact(k,i,j,10,2)               &
                      + ref(kgrd(k,i,j,11,2),igrd(i,j,11),jgrd(i,j,11)) &
                      * hfact(i,j,11) * vfact(k,i,j,11,2)               &
                      + ref(kgrd(k,i,j,12,2),igrd(i,j,12),jgrd(i,j,12)) &
                      * hfact(i,j,12) * vfact(k,i,j,12,2)
       end do
       end do
       end do

    ! logarithmic weighting (for pres, dens)
    if ( logarithmic ) then
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( intp(k,i,j) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_12points


  !-----------------------------------------------------------------------------
  subroutine INTRPNEST_domain_compatibility( &
      lon_org,     & ! (in)
      lat_org,     & ! (in)
      lev_org,     & ! (in)
      lon_loc,     & ! (in)
      lat_loc,     & ! (in)
      lev_loc,     & ! (in)
      skip_x,      & ! (in)
      skip_y,      & ! (in)
      skip_z       ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none
    real(RP), intent(in) :: lon_org(:,:)
    real(RP), intent(in) :: lat_org(:,:)
    real(RP), intent(in) :: lev_org(:,:,:)
    real(RP), intent(in) :: lon_loc(:,:)
    real(RP), intent(in) :: lat_loc(:,:)
    real(RP), intent(in) :: lev_loc(:,:,:)
    logical,  intent(in), optional :: skip_x
    logical,  intent(in), optional :: skip_y
    logical,  intent(in), optional :: skip_z

    real(RP) :: max_ref, min_ref
    real(RP) :: max_loc, min_loc

    logical :: do_xdirec
    logical :: do_ydirec
    logical :: do_zdirec
    !---------------------------------------------------------------------------

    do_xdirec = .true.
    if ( present(skip_x) ) then
    if ( skip_x ) then
       do_xdirec = .false.
    endif
    endif

    do_ydirec = .true.
    if ( present(skip_y) ) then
    if ( skip_y ) then
       do_ydirec = .false.
    endif
    endif

    do_zdirec = .true.
    if ( present(skip_z) ) then
    if ( skip_z ) then
       do_zdirec = .false.
    endif
    endif

    if ( do_xdirec ) then
       max_ref = maxval( lon_org(:,:) / D2R )
       min_ref = minval( lon_org(:,:) / D2R )
       max_loc = maxval( lon_loc(:,:) / D2R )
       min_loc = minval( lon_loc(:,:) / D2R )

       if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- LONGITUDINAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    if ( do_ydirec ) then
       max_ref = maxval( lat_org(:,:) / D2R )
       min_ref = minval( lat_org(:,:) / D2R )
       max_loc = maxval( lat_loc(:,:) / D2R )
       min_loc = minval( lat_loc(:,:) / D2R )

       if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- LATITUDINAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    if ( do_zdirec ) then
       max_ref = maxval( lev_org(:,:,:) )
       !max_loc = maxval( lev_loc(KS-1:KE+1,:,:) ) ! HALO + 1
       max_loc = maxval( lev_loc(:,:,:) ) ! HALO + 1
       !min_ref = minval( lev_org(:,:,:) )
       !min_loc = minval( lev_loc(3:KA,:,:) ) ! HALO + 1

       if ( max_ref < max_loc ) then
       !if ( max_ref < max_loc .or. min_ref > min_loc ) then
          write(*,*) 'xxx ERROR: REQUESTED DOMAIN IS TOO MUCH BROAD'
          write(*,*) 'xxx -- VERTICAL direction over the limit'
          write(*,*) 'xxx -- reference max: ', max_ref
          !write(*,*) 'xxx -- reference min: ', min_ref
          write(*,*) 'xxx --     local max: ', max_loc
          !write(*,*) 'xxx --     local min: ', min_loc
          call PRC_MPIstop
       endif
    endif

    return
  end subroutine INTRPNEST_domain_compatibility

  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function INTRPNEST_haversine( &
      la0,       &
      lo0,       &
      la,        &
      lo )       &
      result( d )
    implicit none
    real(RP), intent(in) :: la0, lo0, la, lo   ! la,la0: Lat, lo,lo0: Lon; [rad]
    real(RP) :: d, dlon, dlat, work1, work2
    !---------------------------------------------------------------------------

    ! output unit : [m]
    dlon = lo0 - lo
    dlat = la0 - la
    work1 = (sin(dlat/2.0_RP))**2.0_RP + &
            cos(la0) * cos(la) * (sin(dlon/2.0_RP))**2.0_RP
    work2 = 2.0_RP * asin(min( 1.0_RP, sqrt(work1) ))
    d = r_in_m * work2

  end function INTRPNEST_haversine

end module scale_interpolation_nest
