!-------------------------------------------------------------------------------
!> Module mkgrd
!!
!! @par Description
!!          Making horizontal grid systems based on the icosahedral grid configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_mkgrd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_grd, only: &
     GRD_XDIR,  &
     GRD_YDIR,  &
     GRD_ZDIR,  &
     GRD_x,     &
     GRD_x_pl,  &
     GRD_xt,    &
     GRD_xt_pl, &
     GRD_s,     &
     GRD_s_pl,  &
     GRD_st,    &
     GRD_st_pl
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKGRD_setup
  public :: MKGRD_standard
  public :: MKGRD_spring
  public :: MKGRD_prerotate
  public :: MKGRD_stretch
  public :: MKGRD_shrink
  public :: MKGRD_rotate
  public :: MKGRD_gravcenter
  public :: MKGRD_diagnosis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_LONG),  public :: MKGRD_IN_BASENAME  = ''
  character(len=H_LONG),  public :: MKGRD_OUT_BASENAME = ''
  character(len=H_SHORT), public :: MKGRD_IN_io_mode   = 'ADVANCED'
  character(len=H_SHORT), public :: MKGRD_OUT_io_mode  = 'ADVANCED'

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: MKGRD_DOSPRING         = .true.
  logical,  private :: MKGRD_DOPREROTATE      = .false.
  logical,  private :: MKGRD_DOSTRETCH        = .false.
  logical,  private :: MKGRD_DOSHRINK         = .false.
  logical,  private :: MKGRD_DOROTATE         = .false.

  real(RP), private :: MKGRD_spring_beta      = 1.15_RP ! parameter beta for spring dynamics
  real(RP), private :: MKGRD_prerotation_tilt =  0.0_RP ! [deg]
  real(RP), private :: MKGRD_stretch_alpha    = 1.00_RP ! parameter alpha for stretch
  integer,  private :: MKGRD_shrink_level     = 0       ! shrink level (only for 1-diamond experiment)
  real(RP), private :: MKGRD_rotation_lon     =  0.0_RP ! [deg]
  real(RP), private :: MKGRD_rotation_lat     = 90.0_RP ! [deg]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKGRD_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF  => CONST_UNDEF
    use mod_adm, only: &
       ADM_nxyz,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_TI,      &
       ADM_TJ
    implicit none

    namelist / PARAM_MKGRD / &
      MKGRD_DOSPRING,         &
      MKGRD_DOPREROTATE,      &
      MKGRD_DOSTRETCH,        &
      MKGRD_DOSHRINK,         &
      MKGRD_DOROTATE,         &
      MKGRD_IN_BASENAME,      &
      MKGRD_IN_io_mode,       &
      MKGRD_OUT_BASENAME,     &
      MKGRD_OUT_io_mode,      &
      MKGRD_spring_beta,      &
      MKGRD_prerotation_tilt, &
      MKGRD_stretch_alpha,    &
      MKGRD_shrink_level,     &
      MKGRD_rotation_lon,     &
      MKGRD_rotation_lat

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Program[mkgrd]/Category[prep]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKGRD,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** PARAM_MKGRD is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist PARAM_MKGRD. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist PARAM_MKGRD. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKGRD)

#ifndef _FIXEDINDEX_
    allocate( GRD_x    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz) )
    allocate( GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz) )
    allocate( GRD_xt   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz) )
    allocate( GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz) )

    allocate( GRD_s    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              2) )
    allocate( GRD_s_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              2) )
    allocate( GRD_st   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,2) )
    allocate( GRD_st_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              2) )
#endif
    GRD_x    (:,:,:,:)   = UNDEF
    GRD_x_pl (:,:,:,:)   = UNDEF
    GRD_xt   (:,:,:,:,:) = UNDEF
    GRD_xt_pl(:,:,:,:)   = UNDEF

    GRD_s    (:,:,:,:)   = UNDEF
    GRD_s_pl (:,:,:,:)   = UNDEF
    GRD_st   (:,:,:,:,:) = UNDEF
    GRD_st_pl(:,:,:,:)   = UNDEF

    return
  end subroutine MKGRD_setup

  !-----------------------------------------------------------------------------
  !> Make standard grid system
  subroutine MKGRD_standard
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_rlevel,  &
       ADM_glevel,  &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_gslf_pl, &
       ADM_NPL,     &
       ADM_SPL
    use scale_const, only: &
       PI => CONST_PI
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP), allocatable :: r0(:,:,:)
    real(RP), allocatable :: r1(:,:,:)
    real(RP), allocatable :: g0(:,:,:)
    real(RP), allocatable :: g1(:,:,:)

    real(RP) :: alpha2, phi

    integer  :: rgnid, dmd
    real(RP) :: rdmd

    integer  :: rgn_all_1d, rgn_all
    integer  :: rgnid_dmd, ir, jr

    integer  :: nmax, nmax_prev, rl, gl
    integer  :: i, j, ij, k, l
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Make standard grid system'

    k = ADM_KNONE

    alpha2 = 2.0_RP * PI / 5.0_RP
    phi    = asin( cos(alpha2) / (1.0_RP-cos(alpha2) ) )

    rgn_all_1d = 2**ADM_rlevel
    rgn_all    = rgn_all_1d * rgn_all_1d

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       nmax = 2
       allocate( r0(nmax,nmax,3) )
       allocate( r1(nmax,nmax,3) )

       dmd = (rgnid-1) / rgn_all + 1

       if ( dmd <= 5 ) then ! northern hemisphere
          rdmd = real(dmd-1,kind=RP)

          r0(1,1,GRD_XDIR) = cos( phi) * cos(alpha2*rdmd)
          r0(1,1,GRD_YDIR) = cos( phi) * sin(alpha2*rdmd)
          r0(1,1,GRD_ZDIR) = sin( phi)

          r0(2,1,GRD_XDIR) = cos(-phi) * cos(alpha2*(rdmd+0.5_RP))
          r0(2,1,GRD_YDIR) = cos(-phi) * sin(alpha2*(rdmd+0.5_RP))
          r0(2,1,GRD_ZDIR) = sin(-phi)

          r0(1,2,GRD_XDIR) =  0.0_RP
          r0(1,2,GRD_YDIR) =  0.0_RP
          r0(1,2,GRD_ZDIR) =  1.0_RP

          r0(2,2,GRD_XDIR) = cos( phi) * cos(alpha2*(rdmd+1.0_RP))
          r0(2,2,GRD_YDIR) = cos( phi) * sin(alpha2*(rdmd+1.0_RP))
          r0(2,2,GRD_ZDIR) = sin( phi)
       else ! southern hemisphere
          rdmd = real(dmd-6,kind=RP)

          r0(1,1,GRD_XDIR) = cos(-phi) * cos(-alpha2*(rdmd+0.5_RP))
          r0(1,1,GRD_YDIR) = cos(-phi) * sin(-alpha2*(rdmd+0.5_RP))
          r0(1,1,GRD_ZDIR) = sin(-phi)

          r0(2,1,GRD_XDIR) =  0.0_RP
          r0(2,1,GRD_YDIR) =  0.0_RP
          r0(2,1,GRD_ZDIR) = -1.0_RP

          r0(1,2,GRD_XDIR) = cos( phi) * cos(-alpha2*rdmd)
          r0(1,2,GRD_YDIR) = cos( phi) * sin(-alpha2*rdmd)
          r0(1,2,GRD_ZDIR) = sin( phi)

          r0(2,2,GRD_XDIR) = cos(-phi) * cos(-alpha2*(rdmd-0.5_RP))
          r0(2,2,GRD_YDIR) = cos(-phi) * sin(-alpha2*(rdmd-0.5_RP))
          r0(2,2,GRD_ZDIR) = sin(-phi)
       endif

       do rl = 1, ADM_rlevel
          nmax_prev = nmax
          nmax = 2 * (nmax-1) + 1

          deallocate( r1 )
          allocate( r1(nmax,nmax,3) )

          call decomposition( nmax_prev, & ! [IN]
                              r0(:,:,:), & ! [IN]
                              nmax,      & ! [IN]
                              r1(:,:,:)  ) ! [OUT]

          deallocate( r0 )
          allocate( r0(nmax,nmax,3) )

          r0(:,:,:) = r1(:,:,:)
       enddo

       nmax = 2
       allocate( g0(nmax,nmax,3) )
       allocate( g1(nmax,nmax,3) )

       rgnid_dmd = mod(rgnid-1,rgn_all) + 1
       ir        = mod(rgnid_dmd-1,rgn_all_1d) + 1
       jr        = (rgnid_dmd-ir) / rgn_all_1d + 1

       g0(1,1,:) = r0(ir  ,jr  ,:)
       g0(2,1,:) = r0(ir+1,jr  ,:)
       g0(1,2,:) = r0(ir  ,jr+1,:)
       g0(2,2,:) = r0(ir+1,jr+1,:)

       do gl = ADM_rlevel+1, ADM_glevel
          nmax_prev = nmax
          nmax = 2 * (nmax-1) + 1

          deallocate( g1 )
          allocate( g1(nmax,nmax,3) )

          call decomposition( nmax_prev, & ! [IN]
                              g0(:,:,:), & ! [IN]
                              nmax,      & ! [IN]
                              g1(:,:,:)  ) ! [OUT]

          deallocate( g0 )
          allocate( g0(nmax,nmax,3) )

          g0(:,:,:) = g1(:,:,:)
       enddo

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          GRD_x(ij,k,l,:) = g0(i-1,j-1,:)
       enddo
       enddo

       deallocate( r0 )
       deallocate( r1 )
       deallocate( g0 )
       deallocate( g1 )
    enddo

    ij = ADM_gslf_pl

    GRD_x_pl(ij,k,ADM_NPL,GRD_XDIR) =  0.0_RP
    GRD_x_pl(ij,k,ADM_NPL,GRD_YDIR) =  0.0_RP
    GRD_x_pl(ij,k,ADM_NPL,GRD_ZDIR) =  1.0_RP

    GRD_x_pl(ij,k,ADM_SPL,GRD_XDIR) =  0.0_RP
    GRD_x_pl(ij,k,ADM_SPL,GRD_YDIR) =  0.0_RP
    GRD_x_pl(ij,k,ADM_SPL,GRD_ZDIR) = -1.0_RP

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_standard

  !-----------------------------------------------------------------------------
  !> Apply spring dynamics
  subroutine MKGRD_spring
    use scale_const, only: &
       PI => CONST_PI
    use scale_vector, only: &
       VECTR_cross, &
       VECTR_dot,   &
       VECTR_abs,   &
       VECTR_angle
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_KNONE,    &
       ADM_have_sgp, &
       ADM_glevel,   &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_gmin,     &
       ADM_gmax
    use mod_comm, only: &
       COMM_data_transfer
    use mod_gm_statistics, only: &
       GTL_max, &
       GTL_min
    implicit none

    integer,  parameter :: var_vindex = 8
    integer,  parameter :: I_Rx   = 1
    integer,  parameter :: I_Ry   = 2
    integer,  parameter :: I_Rz   = 3
    integer,  parameter :: I_Wx   = 4
    integer,  parameter :: I_Wy   = 5
    integer,  parameter :: I_Wz   = 6
    integer,  parameter :: I_Fsum = 7
    integer,  parameter :: I_Ek   = 8

    real(RP) :: var   ( ADM_gall,   ADM_KNONE,ADM_lall,   var_vindex)
    real(RP) :: var_pl( ADM_gall_pl,ADM_KNONE,ADM_lall_pl,var_vindex)

    real(RP), parameter :: dump_coef = 1.0_RP   !> friction coefficent in spring dynamics
    real(RP), parameter :: dt        = 2.E-2_RP !> delta t for solution of spring dynamics
    real(RP), parameter :: criteria  = 1.E-4_RP !> criteria of convergence
    real(RP)            :: lambda, dbar

    real(RP)            :: P(ADM_nxyz,0:6,ADM_gall)
    real(RP)            :: F(ADM_nxyz,1:6,ADM_gall)
    real(RP), parameter :: o(3) = 0.0_RP
    real(RP)            :: fixed_point(3)
    real(RP)            :: P0Pm(3), P0PmP0(3), Fsum(3), R0(3), W0(3)
    real(RP)            :: length, distance, E

    integer,  parameter :: itelim = 10000000
    integer             :: ite
    real(RP)            :: Fsum_max, Ek_max

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, k0, l, m
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSPRING ) return

    k0 = ADM_KNONE

    lambda = 2.0_RP*PI / ( 10.0_RP*2.0_RP**(ADM_glevel-1) )
    dbar   = MKGRD_spring_beta * lambda

    if( IO_L ) write(IO_FID_LOG,*) '*** Apply grid modification with spring dynamics'
    if( IO_L ) write(IO_FID_LOG,*) '*** spring factor beta  = ', MKGRD_spring_beta
    if( IO_L ) write(IO_FID_LOG,*) '*** length lambda       = ', lambda
    if( IO_L ) write(IO_FID_LOG,*) '*** delta t             = ', dt
    if( IO_L ) write(IO_FID_LOG,*) '*** conversion criteria = ', criteria
    if( IO_L ) write(IO_FID_LOG,*) '*** dumping coefficient = ', dump_coef
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(3(A16))') 'itelation', 'max. Kinetic E', 'max. forcing'

    var   (:,:,:,:) = 0.0_RP
    var_pl(:,:,:,:) = 0.0_RP

    var   (:,:,:,I_Rx:I_Rz) = GRD_x   (:,:,:,GRD_XDIR:GRD_ZDIR)
    var_pl(:,:,:,I_Rx:I_Rz) = GRD_x_pl(:,:,:,GRD_XDIR:GRD_ZDIR)

    !--- Solving spring dynamics
    do ite = 1, itelim

       do l = 1, ADM_lall
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij     = suf(i  ,j  )
             ip1j   = suf(i+1,j  )
             ip1jp1 = suf(i+1,j+1)
             ijp1   = suf(i  ,j+1)
             im1j   = suf(i-1,j  )
             im1jm1 = suf(i-1,j-1)
             ijm1   = suf(i  ,j-1)

             P(GRD_XDIR,0,ij) = var(ij    ,k0,l,I_Rx)
             P(GRD_XDIR,1,ij) = var(ip1j  ,k0,l,I_Rx)
             P(GRD_XDIR,2,ij) = var(ip1jp1,k0,l,I_Rx)
             P(GRD_XDIR,3,ij) = var(ijp1  ,k0,l,I_Rx)
             P(GRD_XDIR,4,ij) = var(im1j  ,k0,l,I_Rx)
             P(GRD_XDIR,5,ij) = var(im1jm1,k0,l,I_Rx)
             P(GRD_XDIR,6,ij) = var(ijm1  ,k0,l,I_Rx)

             P(GRD_YDIR,0,ij) = var(ij    ,k0,l,I_Ry)
             P(GRD_YDIR,1,ij) = var(ip1j  ,k0,l,I_Ry)
             P(GRD_YDIR,2,ij) = var(ip1jp1,k0,l,I_Ry)
             P(GRD_YDIR,3,ij) = var(ijp1  ,k0,l,I_Ry)
             P(GRD_YDIR,4,ij) = var(im1j  ,k0,l,I_Ry)
             P(GRD_YDIR,5,ij) = var(im1jm1,k0,l,I_Ry)
             P(GRD_YDIR,6,ij) = var(ijm1  ,k0,l,I_Ry)

             P(GRD_ZDIR,0,ij) = var(ij    ,k0,l,I_Rz)
             P(GRD_ZDIR,1,ij) = var(ip1j  ,k0,l,I_Rz)
             P(GRD_ZDIR,2,ij) = var(ip1jp1,k0,l,I_Rz)
             P(GRD_ZDIR,3,ij) = var(ijp1  ,k0,l,I_Rz)
             P(GRD_ZDIR,4,ij) = var(im1j  ,k0,l,I_Rz)
             P(GRD_ZDIR,5,ij) = var(im1jm1,k0,l,I_Rz)
             P(GRD_ZDIR,6,ij) = var(ijm1  ,k0,l,I_Rz)
          enddo
          enddo

          if ( ADM_have_sgp(l) ) then ! pentagon
             P(:,6,suf(ADM_gmin,ADM_gmin)) = P(:,1,suf(ADM_gmin,ADM_gmin))
          endif

          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             do m = 1, 6
                call VECTR_cross( P0Pm  (:), o(:), P(:,0,ij), o(:), P(:,m,ij) ) ! P0 X Pm
                call VECTR_cross( P0PmP0(:), o(:), P0Pm(:),   o(:), P(:,0,ij) ) ! ( P0 X Pm ) X P0
                call VECTR_abs  ( length, P0PmP0(:) )

                call VECTR_angle( distance, P(:,0,ij), o(:), P(:,m,ij) )

                F(:,m,ij) = ( distance - dbar ) * P0PmP0(:) / length
             enddo
          enddo
          enddo

          if ( ADM_have_sgp(l) ) then ! pentagon
             F(:,6,suf(ADM_gmin,ADM_gmin)) = 0.0_RP

             ! save value of fixed point
             fixed_point(:) = var(suf(ADM_gmin,ADM_gmin),k0,l,I_Rx:I_Rz)
          endif

          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             R0(:) = var(ij,k0,l,I_Rx:I_Rz)
             W0(:) = var(ij,k0,l,I_Wx:I_Wz)

             Fsum(:) = F(:,1,ij) + F(:,2,ij) + F(:,3,ij) + F(:,4,ij) + F(:,5,ij) + F(:,6,ij)

             ! update R0
             R0(:) = R0(:) + W0(:) * dt

             ! normalize
             call VECTR_abs( length, R0(:) )
             R0(:) = R0(:) / length

             ! update W0
             W0(:) = W0(:) + ( Fsum(:) - dump_coef * W0(:) ) * dt

             ! horizontalize
             call VECTR_dot( E, o(:), R0(:), o(:), W0(:) )
             W0(:) = W0(:) - E * R0(:)

             var(ij,k0,l,I_Rx:I_Rz) = R0(:)
             var(ij,k0,l,I_Wx:I_Wz) = W0(:)

             ! check dw0/dt
             call VECTR_abs( length, Fsum(:) )
             var(ij,k0,l,I_Fsum) = length / lambda

             ! kinetic energy
             call VECTR_dot( E, o(:), W0(:), o(:), W0(:) )
             var(ij,k0,l,I_Ek) = 0.5_RP * E
          enddo
          enddo

          ! restore value of fixed point
          if ( ADM_have_sgp(l) ) then ! pentagon
             var(suf(ADM_gmin,ADM_gmin),k0,l,:)         = 0.0_RP
             var(suf(ADM_gmin,ADM_gmin),k0,l,I_Rx:I_Rz) = fixed_point(:)
          endif

       enddo ! l loop

       call COMM_data_transfer( var(:,:,:,:), var_pl(:,:,:,:) )

       Fsum_max = GTL_max( var(:,:,:,I_Fsum), var_pl(:,:,:,I_Fsum), 1, 1, 1 )
       Ek_max   = GTL_max( var(:,:,:,I_Ek),   var_pl(:,:,:,I_Ek)  , 1, 1, 1 )

       if( IO_L ) write(IO_FID_LOG,'(I16,4(E16.8))') ite, Ek_max, Fsum_max

       if( Fsum_max < criteria ) exit

    enddo ! itelation loop

    GRD_x   (:,:,:,GRD_XDIR:GRD_ZDIR) = var   (:,:,:,I_Rx:I_Rz)
    GRD_x_pl(:,:,:,GRD_XDIR:GRD_ZDIR) = var_pl(:,:,:,I_Rx:I_Rz)

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_spring

  !-----------------------------------------------------------------------------
  !> Apply rotation before stretching, for 1-diamond grid system
  subroutine MKGRD_prerotate
    use scale_const, only: &
       PI => CONST_PI
    use scale_vector, only: &
       VECTR_rotation, &
       I_Yaxis,        &
       I_Zaxis
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP) :: g(3)
    real(RP) :: angle_y, angle_z, angle_tilt
    real(RP) :: alpha2

    real(RP) :: d2r
    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOPREROTATE ) return

    k = ADM_KNONE

    d2r        = PI / 180.0_RP
    alpha2     = 2.0_RP * PI / 5.0_RP
    angle_z    = alpha2 / 2.0_RP
    angle_y    = 0.25_RP*PI * ( 3.0_RP - sqrt(3.0_RP) )
    angle_tilt = MKGRD_prerotation_tilt * d2r

    if( IO_L ) write(IO_FID_LOG,*) '*** Apply pre-rotation'
    if( IO_L ) write(IO_FID_LOG,*) '*** Diamond tilting factor = ', MKGRD_prerotation_tilt
    if( IO_L ) write(IO_FID_LOG,*) '*** angle_z   (deg) = ', angle_z    / d2r
    if( IO_L ) write(IO_FID_LOG,*) '*** angle_y   (deg) = ', angle_y    / d2r
    if( IO_L ) write(IO_FID_LOG,*) '*** angle_tilt(deg) = ', angle_tilt / d2r

    do l = 1, ADM_lall
       do ij = 1, ADM_gall
          g(:) = GRD_x(ij,k,l,:)

          ! align lowermost vertex of diamond to x-z coordinate plane
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_z, & ! [IN]
                               I_Zaxis  ) ! [IN]
          ! rotate around y-axis, for fitting the center of diamond to north pole
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_y, & ! [IN]
                               I_Yaxis  ) ! [IN]
          ! rotate the diamond around z-axis
          call VECTR_rotation( g(:),       & ! [INOUT]
                               angle_tilt, & ! [IN]
                               I_Zaxis     ) ! [IN]

          GRD_x(ij,k,l,:) = g(:)
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          g(:) = GRD_x_pl(ij,k,l,:)

          ! align lowermost vertex of diamond to x-z coordinate plane
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_z, & ! [IN]
                               I_Zaxis  ) ! [IN]
          ! rotate around y-axis, for fitting the center of diamond to north pole
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_y, & ! [IN]
                               I_Yaxis  ) ! [IN]
          ! rotate the diamond around z-axis
          call VECTR_rotation( g(:),       & ! [INOUT]
                               angle_tilt, & ! [IN]
                               I_Zaxis     ) ! [IN]

          GRD_x_pl(ij,k,l,:) = g(:)
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_prerotate

  !-----------------------------------------------------------------------------
  !> Apply stretching to grid system
  subroutine MKGRD_stretch
    use scale_const, only: &
       PI => CONST_PI
    use scale_vector, only: &
       VECTR_xyz2latlon, &
       VECTR_latlon2xyz
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP) :: lat, lon, lat_trans

    real(RP), parameter :: criteria = 1.E-10_RP

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSTRETCH ) return

    if( IO_L ) write(IO_FID_LOG,*) '*** Apply stretch'
    if( IO_L ) write(IO_FID_LOG,*) '*** Stretch factor = ', MKGRD_stretch_alpha

    k = ADM_KNONE

    do l = 1, ADM_lall
       do ij = 1, ADM_gall

          call VECTR_xyz2latlon( GRD_x(ij,k,l,GRD_XDIR), & ! [IN]
                                 GRD_x(ij,k,l,GRD_YDIR), & ! [IN]
                                 GRD_x(ij,k,l,GRD_ZDIR), & ! [IN]
                                 lat,                    & ! [OUT]
                                 lon                     ) ! [OUT]

          if ( 0.5_RP*PI-abs(lat) > criteria ) then
             lat_trans = asin( ( MKGRD_stretch_alpha*(1.0_RP+sin(lat)) / (1.0_RP-sin(lat)) - 1.0_RP ) &
                             / ( MKGRD_stretch_alpha*(1.0_RP+sin(lat)) / (1.0_RP-sin(lat)) + 1.0_RP ) )
          else
             lat_trans = lat
          endif

          call VECTR_latlon2xyz( lat_trans,              & ! [IN]
                                 lon,                    & ! [IN]
                                 GRD_x(ij,k,l,GRD_XDIR), & ! [OUT]
                                 GRD_x(ij,k,l,GRD_YDIR), & ! [OUT]
                                 GRD_x(ij,k,l,GRD_ZDIR), & ! [OUT]
                                 1.0_RP                  ) ! [IN]
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl

          call VECTR_xyz2latlon( GRD_x_pl(ij,k,l,GRD_XDIR), & ! [IN]
                                 GRD_x_pl(ij,k,l,GRD_YDIR), & ! [IN]
                                 GRD_x_pl(ij,k,l,GRD_ZDIR), & ! [IN]
                                 lat,                       & ! [OUT]
                                 lon                        ) ! [OUT]

          if ( 0.5_RP*PI-abs(lat) > criteria ) then
             lat_trans = asin( ( MKGRD_stretch_alpha*(1.0_RP+sin(lat)) / (1.0_RP-sin(lat)) - 1.0_RP ) &
                             / ( MKGRD_stretch_alpha*(1.0_RP+sin(lat)) / (1.0_RP-sin(lat)) + 1.0_RP ) )
          else
             lat_trans = lat
          endif

          call VECTR_latlon2xyz( lat_trans,                 & ! [IN]
                                 lon,                       & ! [IN]
                                 GRD_x_pl(ij,k,l,GRD_XDIR), & ! [OUT]
                                 GRD_x_pl(ij,k,l,GRD_YDIR), & ! [OUT]
                                 GRD_x_pl(ij,k,l,GRD_ZDIR), & ! [OUT]
                                 1.0_RP                     ) ! [IN]
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_stretch

  !-----------------------------------------------------------------------------
  !> Apply shrinkng to grid system
  subroutine MKGRD_shrink
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP) :: o(3), g(3), len

    integer  :: ij, k, l, ite
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOSHRINK ) return

    if( IO_L ) write(IO_FID_LOG,*) '*** Apply shrink'
    if( IO_L ) write(IO_FID_LOG,*) '*** Shrink level = ', MKGRD_shrink_level

    k = ADM_KNONE

    o(GRD_XDIR) = 0.0_RP
    o(GRD_YDIR) = 0.0_RP

    do ite = 1, MKGRD_shrink_level
       do l  = 1, ADM_lall
       do ij = 1, ADM_gall
          o(GRD_ZDIR) = sign(1.0_RP,GRD_x(ij,k,l,GRD_ZDIR))

          g(GRD_XDIR) = GRD_x(ij,k,l,GRD_XDIR) + o(GRD_XDIR)
          g(GRD_YDIR) = GRD_x(ij,k,l,GRD_YDIR) + o(GRD_YDIR)
          g(GRD_ZDIR) = GRD_x(ij,k,l,GRD_ZDIR) + o(GRD_ZDIR)

          len = ( g(GRD_XDIR)*g(GRD_XDIR) &
                + g(GRD_YDIR)*g(GRD_YDIR) &
                + g(GRD_ZDIR)*g(GRD_ZDIR) )

          GRD_x(ij,k,l,GRD_XDIR) = g(GRD_XDIR) / len
          GRD_x(ij,k,l,GRD_YDIR) = g(GRD_YDIR) / len
          GRD_x(ij,k,l,GRD_ZDIR) = g(GRD_ZDIR) / len
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then
    do ite = 1, MKGRD_shrink_level-1
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          o(GRD_ZDIR) = sign(1.0_RP,GRD_x_pl(ij,k,l,GRD_ZDIR))

          g(GRD_XDIR) = GRD_x_pl(ij,k,l,GRD_XDIR) + o(GRD_XDIR)
          g(GRD_YDIR) = GRD_x_pl(ij,k,l,GRD_YDIR) + o(GRD_YDIR)
          g(GRD_ZDIR) = GRD_x_pl(ij,k,l,GRD_ZDIR) + o(GRD_ZDIR)

          len = ( g(GRD_XDIR)*g(GRD_XDIR) &
                + g(GRD_YDIR)*g(GRD_YDIR) &
                + g(GRD_ZDIR)*g(GRD_ZDIR) )

          GRD_x_pl(ij,k,l,GRD_XDIR) = g(GRD_XDIR) / len
          GRD_x_pl(ij,k,l,GRD_YDIR) = g(GRD_YDIR) / len
          GRD_x_pl(ij,k,l,GRD_ZDIR) = g(GRD_ZDIR) / len
       enddo
       enddo
    enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_shrink

  !-----------------------------------------------------------------------------
  !> Apply rotation to grid system
  subroutine MKGRD_rotate
    use scale_const, only: &
       PI => CONST_PI
    use scale_vector, only: &
       VECTR_rotation, &
       I_Yaxis,        &
       I_Zaxis
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl
    use mod_comm, only: &
       COMM_data_transfer
    implicit none

    real(RP) :: g(3)
    real(RP) :: angle_y, angle_z

    real(RP) :: d2r
    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    if( .NOT. MKGRD_DOROTATE ) return

    if( IO_L ) write(IO_FID_LOG,*) '*** Apply rotation'
    if( IO_L ) write(IO_FID_LOG,*) '*** North pole -> Longitude(deg) = ', MKGRD_rotation_lon
    if( IO_L ) write(IO_FID_LOG,*) '*** North pole -> Latitude (deg) = ', MKGRD_rotation_lat

    k = ADM_KNONE

    d2r = PI / 180.0_RP
    angle_y = ( MKGRD_rotation_lat - 90.0_RP ) * d2r
    angle_z = - MKGRD_rotation_lon * d2r

    do l = 1, ADM_lall
       do ij = 1, ADM_gall
          g(:) = GRD_x(ij,k,l,:)

          ! rotate around y-axis
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_y, & ! [IN]
                               I_Yaxis  ) ! [IN]
          ! rotate around z-axis
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_z, & ! [IN]
                               I_Zaxis  ) ! [IN]

          GRD_x(ij,k,l,:) = g(:)
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l  = 1, ADM_lall_pl
       do ij = 1, ADM_gall_pl
          g(:) = GRD_x_pl(ij,k,l,:)

          ! rotate around y-axis
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_y, & ! [IN]
                               I_Yaxis  ) ! [IN]
          ! rotate around z-axis
          call VECTR_rotation( g(:),    & ! [INOUT]
                               angle_z, & ! [IN]
                               I_Zaxis  ) ! [IN]

          GRD_x_pl(ij,k,l,:) = g(:)
       enddo
       enddo
    endif

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_rotate

  !-----------------------------------------------------------------------------
  !> Arrange gravitational center
  subroutine MKGRD_gravcenter
    use mod_comm, only: &
       COMM_data_transfer
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Calc gravitational center'

    if( IO_L ) write(IO_FID_LOG,*) '*** center -> vertex'
    call MKGRD_center2vertex

    if( IO_L ) write(IO_FID_LOG,*) '*** vertex -> center'
    call MKGRD_vertex2center

    call COMM_data_transfer( GRD_x(:,:,:,:), GRD_x_pl(:,:,:,:) )

    return
  end subroutine MKGRD_gravcenter

  !-----------------------------------------------------------------------------
  !> Diagnose grid property
  subroutine MKGRD_diagnosis
    use scale_const, only: &
       PI     => CONST_PI,     &
       RADIUS => CONST_RADIUS
    use scale_vector, only: &
       VECTR_cross, &
       VECTR_dot,   &
       VECTR_abs
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_KNONE,    &
       ADM_glevel,   &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_gmax,     &
       ADM_gmin
    use mod_gmtr, only: &
       GMTR_p_AREA, &
       GMTR_p,      &
       GMTR_p_pl
    use mod_gm_statistics, only: &
       GTL_global_sum_srf, &
       GTL_max,            &
       GTL_min
    implicit none

    real(RP) :: angle    (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: angle_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: length   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: length_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: sqarea   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: sqarea_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: dummy    (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP) :: dummy_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)

    real(RP) :: len(6), ang(6)
    real(RP) :: p(ADM_nxyz,0:7)
    real(RP) :: nvlenC, nvlenS, nv(3)

    real(RP) :: nlen, len_tot
    real(RP) :: l_mean, area, temp
    real(RP) :: sqarea_avg, sqarea_max, sqarea_min
    real(RP) :: angle_max,  length_max, length_avg

    real(RP) :: global_area
    integer  :: global_grid

    integer  :: i, j, ij, k, l, m
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Diagnose grid property'

    k = ADM_KNONE

    angle    (:,:,:) = 0.0_RP
    angle_pl (:,:,:) = 0.0_RP
    length   (:,:,:) = 0.0_RP
    length_pl(:,:,:) = 0.0_RP

    nlen    = 0.0_RP
    len_tot = 0.0_RP

    do l = 1, ADM_lall
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ij = suf(i,j)

       if (       ADM_have_sgp(l) &
            .AND. i == ADM_gmin   &
            .AND. j == ADM_gmin   ) then ! Pentagon

          p(:,0) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
          p(:,1) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)
          p(:,2) = GRD_xt(suf(i,  j  ),k,l,ADM_TJ,:)
          p(:,3) = GRD_xt(suf(i-1,j  ),k,l,ADM_TI,:)
          p(:,4) = GRD_xt(suf(i-1,j-1),k,l,ADM_TJ,:)
          p(:,5) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
          p(:,6) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)

          len(:) = 0.0_RP
          ang(:) = 0.0_RP
          do m = 1, 5
             ! vector length of Pm->Pm-1, Pm->Pm+1
             call VECTR_dot( len(m), p(:,m), p(:,m-1), p(:,m), p(:,m-1) )
             len(m) = sqrt( len(m) )

             ! angle of Pm-1->Pm->Pm+1
             call VECTR_dot( nvlenC, p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
             call VECTR_cross( nv(:), p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
             call VECTR_abs( nvlenS, nv(:) )

             ang(m) = atan2( nvlenS, nvlenC )
          enddo

          ! maximum/minimum ratio of angle between the cell vertexes
          angle(ij,k,l) = maxval( ang(1:5) ) / minval( ang(1:5) ) - 1.0_RP

          ! l_mean: side length of regular pentagon =sqrt(area/1.7204774005)
          area   = GMTR_p(ij,k,l,GMTR_p_AREA)
          l_mean = sqrt( 4.0_RP / sqrt( 25.0_RP + 10.0_RP*sqrt(5.0_RP)) * area )

          temp = 0.0_RP
          do m = 1, 5
             nlen    = nlen + 1.0_RP
             len_tot = len_tot + len(m)

             temp = temp + (len(m)-l_mean) * (len(m)-l_mean)
          enddo
          ! distortion of side length from l_mean
          length(ij,k,l) = sqrt( temp/5.0_RP ) / l_mean

       else ! Hexagon

          p(:,0) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
          p(:,1) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)
          p(:,2) = GRD_xt(suf(i,  j  ),k,l,ADM_TJ,:)
          p(:,3) = GRD_xt(suf(i-1,j  ),k,l,ADM_TI,:)
          p(:,4) = GRD_xt(suf(i-1,j-1),k,l,ADM_TJ,:)
          p(:,5) = GRD_xt(suf(i-1,j-1),k,l,ADM_TI,:)
          p(:,6) = GRD_xt(suf(i,  j-1),k,l,ADM_TJ,:)
          p(:,7) = GRD_xt(suf(i,  j  ),k,l,ADM_TI,:)

          len(:) = 0.0_RP
          ang(:) = 0.0_RP
          do m = 1, 6
             ! vector length of Pm->Pm-1, Pm->Pm+1
             call VECTR_dot( len(m), p(:,m), p(:,m-1), p(:,m), p(:,m-1) )
             len(m) = sqrt( len(m) )

             ! angle of Pm-1->Pm->Pm+1
             call VECTR_dot( nvlenC, p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
             call VECTR_cross( nv(:), p(:,m), p(:,m-1), p(:,m), p(:,m+1) )
             call VECTR_abs( nvlenS, nv(:) )

             ang(m) = atan2( nvlenS, nvlenC )
          enddo

          ! maximum/minimum ratio of angle between the cell vertexes
          angle(ij,k,l) = maxval( ang(:) ) / minval( ang(:) ) - 1.0_RP

          ! l_mean: side length of equilateral triangle
          area   = GMTR_p(ij,k,l,GMTR_p_AREA)
          l_mean = sqrt( 4.0_RP / sqrt(3.0_RP) / 6.0_RP * area )

          temp = 0.0_RP
          do m = 1, 6
             nlen = nlen + 1.0_RP
             len_tot = len_tot + len(m)

             temp = temp + (len(m)-l_mean)*(len(m)-l_mean)
          enddo
          ! distortion of side length from l_mean
          length(ij,k,l) = sqrt( temp/6.0_RP ) / l_mean

       endif
    enddo
    enddo
    enddo

    dummy    (:,:,:) = 1.0_RP
    dummy_pl (:,:,:) = 1.0_RP
    global_area = GTL_global_sum_srf( dummy(:,:,:), dummy_pl(:,:,:) )
    global_grid = 10*4**ADM_glevel + 2
    sqarea_avg = sqrt( global_area / real(global_grid,kind=RP) )

    sqarea   (:,:,:) = sqrt( GMTR_p   (:,:,:,GMTR_p_AREA) )
    sqarea_pl(:,:,:) = sqrt( GMTR_p_pl(:,:,:,GMTR_p_AREA) )
    sqarea_max = GTL_max ( sqarea(:,:,:), sqarea_pl(:,:,:), 1, 1, 1 )
    sqarea_min = GTL_min ( sqarea(:,:,:), sqarea_pl(:,:,:), 1, 1, 1 )

    length_avg = len_tot / nlen
    length_max = GTL_max( length(:,:,:), length_pl(:,:,:), 1, 1, 1 )
    angle_max  = GTL_max( angle (:,:,:), angle_pl (:,:,:), 1, 1, 1 )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '------ Diagnosis result ---'
    if( IO_L ) write(IO_FID_LOG,*) '--- ideal  global surface area  = ', 4.0_RP*PI*RADIUS*RADIUS*1.E-6_RP,' [km2]'
    if( IO_L ) write(IO_FID_LOG,*) '--- actual global surface area  = ', global_area*1.E-6_RP,' [km2]'
    if( IO_L ) write(IO_FID_LOG,*) '--- global total number of grid = ', global_grid
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '--- average grid interval       = ', sqarea_avg * 1.E-3_RP,' [km]'
    if( IO_L ) write(IO_FID_LOG,*) '--- max grid interval           = ', sqarea_max * 1.E-3_RP,' [km]'
    if( IO_L ) write(IO_FID_LOG,*) '--- min grid interval           = ', sqarea_min * 1.E-3_RP,' [km]'
    if( IO_L ) write(IO_FID_LOG,*) '--- ratio max/min grid interval = ', sqarea_max / sqarea_min
    if( IO_L ) write(IO_FID_LOG,*) '--- average length of arc(side) = ', length_avg * 1.E-3_RP,' [km]'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '--- max length distortion       = ', length_max * 1.D-3,' [km]'
    if( IO_L ) write(IO_FID_LOG,*) '--- max angle distortion        = ', angle_max*180.0_RP/PI,' [deg]'

    return
  end subroutine MKGRD_diagnosis

  !---------------------------------------------------------------------------==
  subroutine decomposition( &
      n0, &
      g0, &
      n1, &
      g1  )
    implicit none

    integer,  intent(in)  :: n0
    real(RP), intent(in)  :: g0(n0,n0,3)
    integer,  intent(in)  :: n1
    real(RP), intent(out) :: g1(n1,n1,3)

    real(RP) :: r
    integer  :: i, j, inew, jnew
    !---------------------------------------------------------------------------

    do i = 1, n0
    do j = 1, n0
       inew = 2 * i - 1
       jnew = 2 * j - 1

       g1(inew,jnew,:) = g0(i,j,:)

       if ( i < n0 ) then
          g1(inew+1,jnew  ,:) = g0(i+1,j  ,:) + g0(i,j,:)
       endif
       if ( j < n0 ) then
          g1(inew  ,jnew+1,:) = g0(i  ,j+1,:) + g0(i,j,:)
       endif
       if ( i < n0 .AND. j < n0 ) then
          g1(inew+1,jnew+1,:) = g0(i+1,j+1,:) + g0(i,j,:)
       endif
    enddo
    enddo

    do i = 1, n1
    do j = 1, n1
       r = sqrt( g1(i,j,1)*g1(i,j,1) &
               + g1(i,j,2)*g1(i,j,2) &
               + g1(i,j,3)*g1(i,j,3) )

       g1(i,j,1) = g1(i,j,1) / r
       g1(i,j,2) = g1(i,j,2) / r
       g1(i,j,3) = g1(i,j,3) / r
    enddo
    enddo

    return
  end subroutine decomposition

  !-----------------------------------------------------------------------------
  !> gnomonic projection
  subroutine MISC_latlon2gnom( &
      x,          &
      y,          &
      lat,        &
      lon,        &
      lat_center, &
      lon_center  )
    implicit none

    real(RP), intent(out) :: x          !> gnomonic, x
    real(RP), intent(out) :: y          !> gnomonic, y
    real(RP), intent(in)  :: lat        !> spheric, latitude
    real(RP), intent(in)  :: lon        !> spheric, longitude
    real(RP), intent(in)  :: lat_center !> projection center, latitude
    real(RP), intent(in)  :: lon_center !> projection center, longitude

    real(RP) :: cosc
    !---------------------------------------------------------------------------

    cosc = sin(lat_center) * sin(lat) &
         + cos(lat_center) * cos(lat) * cos(lon-lon_center)

    x = ( cos(lat) * sin(lon-lon_center) ) / cosc
    y = ( cos(lat_center) * sin(lat)                       &
        - sin(lat_center) * cos(lat) * cos(lon-lon_center) ) / cosc

  end subroutine MISC_latlon2gnom

  !-----------------------------------------------------------------------------
  !> gnomonic projection (inverse)
  subroutine MISC_gnom2latlon( &
      lat,        &
      lon,        &
      x,          &
      y,          &
      lat_center, &
      lon_center  )
    implicit none

    real(RP), intent(out) :: lat        !> spheric, latitude
    real(RP), intent(out) :: lon        !> spheric, longitude
    real(RP), intent(in)  :: x          !> gnomonic, x
    real(RP), intent(in)  :: y          !> gnomonic, y
    real(RP), intent(in)  :: lat_center !> projection center, latitude
    real(RP), intent(in)  :: lon_center !> projection center, longitude

    real(RP) :: rho, c
    !---------------------------------------------------------------------------

    rho = sqrt( x*x + y*y )

    if ( rho == 0.0_RP ) then ! singular point
       lat = lat_center
       lon = lon_center
       return
    endif

    c = atan( rho )

    lon = lon_center + atan2( x*sin(c), ( rho*cos(lat_center)*cos(c) - y*sin(lat_center)*sin(c) ) )
    lat = asin( cos(c)*sin(lat_center) + y*sin(c)*cos(lat_center) / rho )

    return
  end subroutine MISC_gnom2latlon

  !-----------------------------------------------------------------------------
  !> Make center grid -> vertex grid
  subroutine MKGRD_center2vertex
    use mod_adm, only: &
       ADM_nxyz,     &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_KNONE,    &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use scale_vector, only: &
       VECTR_cross, &
       VECTR_dot,   &
       VECTR_abs
    implicit none

    real(RP) :: wk   (ADM_nxyz,4,ADM_gall,ADM_TI:ADM_TJ)
    real(RP) :: wk_pl(ADM_nxyz,4)

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP) :: r(3), gc(3)
    real(RP) :: r_lenS, r_lenC, gc_len

    integer  :: ij
    integer  :: ip1j, ip1jp1, ijp1

    integer  :: i, j, k0, l, d, v, n, t, m
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij     = suf(i  ,j  )
          ip1j   = suf(i+1,j  )
          ip1jp1 = suf(i+1,j+1)
          ijp1   = suf(i  ,j+1)

          do d = 1, ADM_nxyz
             wk(d,1,ij,ADM_TI) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,ADM_TI) = GRD_x(ip1j  ,k0,l,d)
             wk(d,3,ij,ADM_TI) = GRD_x(ip1jp1,k0,l,d)
             wk(d,4,ij,ADM_TI) = GRD_x(ij    ,k0,l,d)

             wk(d,1,ij,ADM_TJ) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,ADM_TJ) = GRD_x(ip1jp1,k0,l,d)
             wk(d,3,ij,ADM_TJ) = GRD_x(ijp1  ,k0,l,d)
             wk(d,4,ij,ADM_TJ) = GRD_x(ij    ,k0,l,d)
          enddo
       enddo
       enddo

       !--- treat unused triangle
       wk(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TI) = wk(:,:,suf(ADM_gmax,ADM_gmin-1),ADM_TJ)
       wk(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TJ) = wk(:,:,suf(ADM_gmin-1,ADM_gmax),ADM_TI)

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,:,suf(ADM_gmin-1,ADM_gmin-1),ADM_TI) = wk(:,:,suf(ADM_gmin,ADM_gmin-1),ADM_TJ)
       endif

       do t = ADM_TI, ADM_TJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          gc(:) = 0.0_RP
          do m = 1, 3
             call VECTR_dot  ( r_lenC, o(:), wk(:,m,ij,t), o(:), wk(:,m+1,ij,t) )
             call VECTR_cross( r(:),   o(:), wk(:,m,ij,t), o(:), wk(:,m+1,ij,t) )
             call VECTR_abs  ( r_lenS, r(:) )

             r(:) = r(:) / r_lenS * atan2( r_lenS, r_lenC )

             gc(:) = gc(:) + r(:)
          enddo

          call VECTR_abs( gc_len, gc(:) )

          GRD_xt(ij,k0,l,t,:) = gc(:) / gc_len
       enddo
       enddo
       enddo

    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1,ADM_lall_pl
       do v = ADM_gmin_pl, ADM_gmax_pl
          ij   = v
          ijp1 = v + 1
          if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

          do d = 1, ADM_nxyz
             wk_pl(:,1) = GRD_x_pl(n   ,k0,l,:)
             wk_pl(:,2) = GRD_x_pl(ij  ,k0,l,:)
             wk_pl(:,3) = GRD_x_pl(ijp1,k0,l,:)
             wk_pl(:,4) = GRD_x_pl(n   ,k0,l,:)
          enddo

          gc(:) = 0.0_RP
          do m = 1, 3
             call VECTR_dot  ( r_lenC, o(:), wk_pl(:,m), o(:), wk_pl(:,m+1) )
             call VECTR_cross( r(:),   o(:), wk_pl(:,m), o(:), wk_pl(:,m+1) )
             call VECTR_abs  ( r_lenS, r(:) )

             r(:) = r(:) / r_lenS * atan2( r_lenS, r_lenC )

             gc(:) = gc(:) + r(:)
          enddo

          call VECTR_abs( gc_len, gc(:) )

          GRD_xt_pl(v,k0,l,:) = -gc(:) / gc_len
       enddo
       enddo
    endif

    return
  end subroutine MKGRD_center2vertex

  !-----------------------------------------------------------------------------
  !> Make vertex grid -> center grid
  subroutine MKGRD_vertex2center
    use scale_const, only: &
       EPS => CONST_EPS
    use mod_adm, only : &
       ADM_nxyz,     &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_KNONE,    &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_vlink,    &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl
    use scale_vector, only: &
       VECTR_cross, &
       VECTR_dot,   &
       VECTR_abs
    implicit none

    real(RP) :: wk   (ADM_nxyz,7,ADM_gall)
    real(RP) :: wk_pl(ADM_nxyz,ADM_vlink+1)

    real(RP), parameter :: o(3) = 0.0_RP
    real(RP) :: r(3), gc(3)
    real(RP) :: r_lenS, r_lenC, gc_len
    real(RP) :: zerosw

    integer  :: ij
    integer  :: im1j, im1jm1, ijm1

    integer  :: i, j, k0, l, d, v, n, m
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = suf(i  ,j  )
          im1j   = suf(i-1,j  )
          im1jm1 = suf(i-1,j-1)
          ijm1   = suf(i  ,j-1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(ijm1  ,k0,l,ADM_TJ,d)
             wk(d,2,ij) = GRD_xt(ij    ,k0,l,ADM_TI,d)
             wk(d,3,ij) = GRD_xt(ij    ,k0,l,ADM_TJ,d)
             wk(d,4,ij) = GRD_xt(im1j  ,k0,l,ADM_TI,d)
             wk(d,5,ij) = GRD_xt(im1jm1,k0,l,ADM_TJ,d)
             wk(d,6,ij) = GRD_xt(im1jm1,k0,l,ADM_TI,d)
             wk(d,7,ij) = GRD_xt(ijm1  ,k0,l,ADM_TJ,d)
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,6,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
          wk(:,7,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
       endif

       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          gc(:) = 0.0_RP
          do m = 1, 6
             call VECTR_dot  ( r_lenC, o(:), wk(:,m,ij), o(:), wk(:,m+1,ij) )
             call VECTR_cross( r(:),   o(:), wk(:,m,ij), o(:), wk(:,m+1,ij) )
             call VECTR_abs  ( r_lenS, r(:) )

             zerosw = 0.5_RP - sign(0.5_RP,abs(r_lenS)-EPS)
             r(:) = r(:) * ( 1.0_RP - zerosw ) / ( r_lenS + zerosw ) * atan2( r_lenS, r_lenC )

             gc(:) = gc(:) + r(:)
          enddo

          call VECTR_abs( gc_len, gc(:) )

          GRD_x(ij,k0,l,:) = gc(:) / gc_len
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1,ADM_lall_pl
          do d = 1, ADM_nxyz
             do v = 1, ADM_vlink ! (ICO=5)
                wk_pl(d,v) = GRD_xt_pl(v+1,k0,l,d)
             enddo
             wk_pl(d,ADM_vlink+1) = wk_pl(d,1)
          enddo

          gc(:) = 0.0_RP
          do v = 1, ADM_vlink ! (ICO=5)
             call VECTR_dot  ( r_lenC, o(:), wk_pl(:,v), o(:), wk_pl(:,v+1) )
             call VECTR_cross( r(:),   o(:), wk_pl(:,v), o(:), wk_pl(:,v+1) )
             call VECTR_abs  ( r_lenS, r(:) )

             r(:) = r(:) / r_lenS * atan2( r_lenS, r_lenC )

             gc(:) = gc(:) + r(:)
          enddo

          call VECTR_abs( gc_len, gc(:) )

          GRD_x_pl(n,k0,l,:) = -gc(:) / gc_len
       enddo
    endif

    return
  end subroutine MKGRD_vertex2center

  !-----------------------------------------------------------------------------
  !> suffix calculation
  !> @return suf
  function suf(i,j) result(suffix)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: suffix
    integer :: i, j
    !---------------------------------------------------------------------------

    suffix = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_mkgrd
