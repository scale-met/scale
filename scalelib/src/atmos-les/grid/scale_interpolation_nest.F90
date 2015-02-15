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
  private :: INTRPNEST_search_horiz_3points
  private :: INTRPNEST_search_horiz_4points
  private :: INTRPNEST_search_vert_offline
  private :: INTRPNEST_search_vert_online
  private :: INTRPNEST_interp_2d_3points
  private :: INTRPNEST_interp_3d_3points
  private :: INTRPNEST_interp_2d_4points
  private :: INTRPNEST_interp_3d_4points
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
  public :: INTRPNEST_search_horiz

  !-----------------------------------------------------------------------------
  abstract interface
     subroutine INTRPNEST_intfc_search_v(  &
          vfact,   & ! (out)
          kgrd,    & ! (out)
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

       real(RP), intent(out) :: vfact(:,:,:,:,:)  ! vertical interp factor
       integer,  intent(out) :: kgrd (:,:,:,:,:)  ! grid points of interp target
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
  public :: INTRPNEST_search_vert

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
  real(RP), private, parameter :: large_number_1 = 9.999E+15_RP
  real(RP), private, parameter :: large_number_2 = 8.888E+15_RP
  real(RP), private, parameter :: large_number_3 = 7.777E+15_RP
  real(RP), private, parameter :: large_number_4 = 6.666E+15_RP

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

    real(RP) :: distance
    real(RP) :: dist
    integer  :: i, j, ii, jj
    integer  :: idx
    integer  :: is, ie
    integer  :: js, je
    integer  :: iinc, jinc
    integer  :: blk_i, blk_j
    !---------------------------------------------------------------------------

    hfact(:,:,:) = 0.0_RP

    do j = 1, myJA
    do i = 1, myIA
       ! nearest block search
       iinc = max( (inIA + 1) / divnum, 1 )
       jinc = max( (inJA + 1) / divnum, 1 )
       dist = large_number_1
       jj = 1 + (jinc/2)
       do while (jj <= inJA)
          ii = 1 + (iinc/2)
          do while (ii <= inIA)
             distance = INTRPNEST_haversine( mylat(i,j),   mylon(i,j),  &
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
       is = blk_i - (iinc/2) - 1
       if( is < 1 ) is = 1
       ie   = blk_i + (iinc/2) + 1
       if( ie  > inIA ) ie   = inIA
       js = blk_j - (jinc/2) - 1
       if( js < 1 ) js = 1
       je   = blk_j + (jinc/2) + 1
       if( je  > inJA ) je   = inJA

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

    real(RP) :: distance
    real(RP) :: dist
    integer  :: i, j, k, ii, jj, kk
    integer  :: idx
    integer  :: is, ie
    integer  :: js, je
    integer  :: iinc, jinc
    integer  :: blk_i, blk_j
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

    do j = 1, myJA
    do i = 1, myIA

       ! nearest block search
       iinc = max( (inIA + 1) / divnum, 1 )
       jinc = max( (inJA + 1) / divnum, 1 )
       dist = large_number_1
       jj = 1 + (jinc/2)
       do while (jj <= inJA)
          ii = 1 + (iinc/2)
          do while (ii <= inIA)
             distance = INTRPNEST_haversine( mylat(i,j),   mylon(i,j),  &
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
       is = blk_i - (iinc/2) - 1
       if( is < 1 ) is = 1
       ie   = blk_i + (iinc/2) + 1
       if( ie  > inIA ) ie   = inIA
       js = blk_j - (jinc/2) - 1
       if( js < 1 ) js = 1
       je   = blk_j + (jinc/2) + 1
       if( je  > inJA ) je   = inJA

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
  ! vertical search of interpolation points for two-points (online)
  subroutine INTRPNEST_search_vert_online( &
      vfact,   & ! (out)
      kgrd,    & ! (out)
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

    real(RP), intent(out) :: vfact(:,:,:,:,:)   ! vertical interp factor
    integer,  intent(out) :: kgrd (:,:,:,:,:)   ! grid points of interp target

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
    integer  :: ii, jj, idx
    integer  :: k, kk
    integer  :: KS_local
    integer  :: ncopy
    logical  :: copy
    !---------------------------------------------------------------------------

    if ( lndgrd ) then
       write(*,*) 'xxx internal error [interporation: nest/interp]'
       write(*,*) '    land grid is not araviable in online'
       call PRC_MPIstop
    endif

    KS_local = 1 + KHALO
    do idx = 1, itp_nh
       ii = igrd(idx)
       jj = jgrd(idx)
       ncopy = 0

       do k = ks, ke
          dist(1) = large_number_2
          dist(2) = large_number_1
          kgrd(k,iloc,jloc,idx,:) = -1
          copy = .false.

          do kk = 1+KHALO, inKA-KHALO
             distance = abs( myhgt(k) - inhgt(kk,ii,jj) )
             if ( distance <= dist(1) ) then
                dist(2) = dist(1);     kgrd(k,iloc,jloc,idx,2) = kgrd(k,iloc,jloc,idx,1)
                dist(1) = distance;    kgrd(k,iloc,jloc,idx,1) = kk
             elseif ( dist(1) < distance .and. distance <= dist(2) ) then
                dist(2) = distance;    kgrd(k,iloc,jloc,idx,2) = kk
             endif
          enddo

          if( inhgt(KS_local,ii,jj) > myhgt(k) ) then
             copy = .true.
          endif

          if( copy ) then
             kgrd(k,iloc,jloc,idx,1)  = KS_local
             kgrd(k,iloc,jloc,idx,2)  = KS_local + 1 ! not used
             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP
             ncopy = ncopy + 1
          elseif( (.NOT. copy) .and. abs(dist(1)) < eps ) then
             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP
          elseif( .NOT. copy ) then
             denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
             vfact(k,iloc,jloc,idx,1) = ( 1.0_RP/dist(1) ) * denom
             vfact(k,iloc,jloc,idx,2) = ( 1.0_RP/dist(2) ) * denom
          else
             write(*,*) 'xxx internal error [interporation: nest/interp]'
             call PRC_MPIstop
          endif
       enddo

       ! not to output message: tentative
       !if( ncopy > 1 )then ! copy is allowed only one time.
       !   write(*,*) 'xxx ERROR: times of copying is exceeded allowed times'
       !   write(*,*) 'xxx domain number: ', ONLINE_DOMAIN_NUM
       !   write(*,*) 'xxx copy times: ', ncopy
       !   call PRC_MPIstop
       !endif
    enddo

    return
  end subroutine INTRPNEST_search_vert_online


  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points (offline)
  subroutine INTRPNEST_search_vert_offline( &
      vfact,   & ! (out)
      kgrd,    & ! (out)
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
    implicit none

    real(RP), intent(out) :: vfact(:,:,:,:,:)   ! vertical interp factor
    integer,  intent(out) :: kgrd (:,:,:,:,:)   ! grid points of interp target

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
    integer  :: ii, jj
    integer  :: k, kk, kks, kke
    integer  :: idx
    !---------------------------------------------------------------------------

    if ( lndgrd ) then
       kks = 1;     kke = LKMAX
    else
       kks = ks-1;  kke = ke+1
    endif

    do idx = 1, itp_nh
       ii = igrd(idx)
       jj = jgrd(idx)

       do k = kks, kke
          dist(1) = large_number_2
          dist(2) = large_number_1
          kgrd(k,iloc,jloc,idx,:) = -1

          do kk = 1, inKA
             distance = abs( myhgt(k) - inhgt(kk,ii,jj) )
             if ( distance <= dist(1) ) then
                dist(2) = dist(1);     kgrd(k,iloc,jloc,idx,2) = kgrd(k,iloc,jloc,idx,1)
                dist(1) = distance;    kgrd(k,iloc,jloc,idx,1) = kk
             elseif ( dist(1) < distance .and. distance <= dist(2) ) then
                dist(2) = distance;    kgrd(k,iloc,jloc,idx,2) = kk
             endif
          enddo

          if ( abs(dist(1)) < eps ) then
             vfact(k,iloc,jloc,idx,1) = 1.0_RP
             vfact(k,iloc,jloc,idx,2) = 0.0_RP
          else
             denom = 1.0_RP / ( (1.0_RP/dist(1)) + (1.0_RP/dist(2)) )
             vfact(k,iloc,jloc,idx,1) = ( 1.0_RP/dist(1) ) * denom
             vfact(k,iloc,jloc,idx,2) = ( 1.0_RP/dist(2) ) * denom
          endif
       enddo
    enddo

    return
  end subroutine INTRPNEST_search_vert_offline


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
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine

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

    if ( .not. logarithmic ) then
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
    else
       ! logarithmic weighted interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
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
                           * hfact(i,j,3) * vfact(k,i,j,3,2) )
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
    real(RP), intent(in)  :: hfact(:,:,:)      ! horizontal interp factor
    integer,  intent(in)  :: igrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: jgrd (:,:,:)      ! grid points of interp target
    integer,  intent(in)  :: ia                ! grid number of mine
    integer,  intent(in)  :: ja                ! grid number of mine

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

    if ( .not. logarithmic ) then
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
    else
       ! logarithmic weighted interpolation
       do j = 1, ja
       do i = 1, ia
       do k = ks, ke
          intp(k,i,j) = exp( ref(kgrd(k,i,j,1,1),igrd(i,j,1),jgrd(i,j,1)) &
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
                           * hfact(i,j,4) * vfact(k,i,j,4,2) )
       end do
       end do
       end do
    endif

    return
  end subroutine INTRPNEST_interp_3d_4points


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
