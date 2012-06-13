!-------------------------------------------------------------------------------
!
!+  NICAM Double moment Water 6 scheme
!
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for the Spectral Bin Model
  !
  !       
  !++ Current Corresponding Author : Y.Sato, K. Suzuki
  ! 
  !++ History: NDW6
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !        0.00   12/06/04 Y.Sato, import from version 4.1 of original code
  !               
  !      -----------------------------------------------------------------------
  !      Reference:  -- Journals 
  !                   Suzuki et al.(2006): SOLA,vol.2,pp.116-119(original code)
  !                   Suzuki et al.(2010): J. Atmos. Sci.,vol.67,pp.1126-1141
  !                   Sato et al. (2009): J. Geophys. Res.114,doi
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_const, only: &
     pi => CONST_PI, &
     CONST_CPdry, &
     CONST_DWATR, &
     CONST_GRAV, &
     CONST_Rvap, &
     CONST_Rdry, &
     CONST_LH0, &
     CONST_LHS0, &
     CONST_EMELT, &
     CONST_TEM00, &
     CONST_TMELT, &
     CONST_PSAT0, &
     CONST_PRE00, &
     CONST_CPovCV, &
     CONST_RovCP
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  !++ Private parameters
  !
  !--- Indeces
  integer, parameter :: il = 1

  !--- bin information of hydrometeors
  real(8) :: xctr( nbin )         !--- log( m ) value of bin center 
  real(8) :: xbnd( nbin+1 )       !--- log( m ) value of bin boundary
  real(8) :: radc( nbin )         !--- radius of hydrometeor at bin center [m]
  real(8) :: dxmic                !--- d( log(m) ) of hydrometeor bin
  real(8) :: cctr( 7,nbin )       !--- capacitance of hydrometeor at bin center
  real(8) :: cbnd( 7,nbin+1 )     !--- capacitance of hydrometeor at bin boundary
  real(8) :: ck( 7,7,nbin,nbin )  !-- collection kernel
  real(8) :: vt( 7,nbin )         !--- terminal velocity of hydrometeor [m/s]
  real(8) :: br( 7,nbin )         !--- bulk density of hydrometeor [kg/m^3] 
  !--- bin information of aerosol (not supported)
  real(8) :: xactr( nccn )        !--- log( ma ) value of bin center
  real(8) :: xabnd( nccn+1 )      !--- log( ma ) value of bin boundary
  real(8) :: dxaer                !--- d( log(ma) ) of aerosol bin
  real(8) :: xasta
  real(8) :: xaend
  real(8) :: rada( nccn )
  !--- constant for bin
  real(8), parameter :: cldmin = 1.0D-10, eps = 1.0D-30, OneovThird=1.D0/3.D0
  !--- constant for aerosol
  real(8) :: rhoa   = 2.25D+03 ! density of aerosol ( NaCl )
  real(8) :: emaer  = 58.D0    ! molecular weight of aerosol ( salt )
  real(8) :: emwtr  = 18.D0    ! molecular weight of water
  real(8) :: rasta  = 1.D-08   ! minimum radius of aerosol (m)
  real(8) :: raend  = 1.D-06   ! maximum radius of aerosol (m)
  real(8) :: r0a    = 1.D-07   ! average radius of aerosol (m)
  logical :: flg_regeneration=.false.  ! flag regeneration of aerosol
  logical :: flg_nucl=.false.  ! flag regeneration of aerosol
  logical :: flg_sf_aero =.false.  ! flag surface flux of aerosol
  real(8) :: marate( nccn )    ! mass rate of each aerosol bin to total aerosol mass
  integer, private, save      :: K10_1, K10_2        ! scaling factor for 10m value (momentum)
  real(8), private, save      :: R10M1, R10M2        ! scaling factor for 10m value (momentum)
  real(8), private, save      :: R10H1, R10H2        ! scaling factor for 10m value (heat)
  real(8), private, save      :: R10E1, R10E2        ! scaling factor for 10m value (tracer)

  character(11),parameter :: fname_micpara="micpara.dat"
  integer(4) :: fid_micpara
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_stdio, only: &
      IO_FID_CONF
    use mod_process, only: &
      PRC_MPIstop
    use mod_grid, only: &
      CDZ => GRID_CDZ, &
      CZ  => GRID_CZ,  &
      FZ  => GRID_FZ
    implicit none
    !---------------------------------------------------------------------------

    real(8) :: ATMOS_PHY_MP_RHOA  !--- density of aerosol
    real(8) :: ATMOS_PHY_MP_EMAER !--- moleculer weight of aerosol
    real(8) :: ATMOS_PHY_MP_RAMIN !--- minimum radius of aerosol (um)
    real(8) :: ATMOS_PHY_MP_RAMAX !--- maximum radius of aerosol (um)
    real(8) :: ATMOS_PHY_MP_R0A   !--- maximum radius of aerosol (um)
    logical :: ATMOS_PHY_MP_FLAG_REGENE  !--- flag of regeneration
    logical :: ATMOS_PHY_MP_FLAG_NUCLEAT !--- flag of regeneration
    logical :: ATMOS_PHY_MP_FLAG_SFAERO  !--- flag of surface flux of aeorol

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       ATMOS_PHY_MP_RHOA,  &
       ATMOS_PHY_MP_EMAER, &
       ATMOS_PHY_MP_RAMIN, &
       ATMOS_PHY_MP_RAMAX, &
       ATMOS_PHY_MP_FLAG_REGENE,  &
       ATMOS_PHY_MP_FLAG_NUCLEAT, &
       ATMOS_PHY_MP_FLAG_SFAERO,  &
       ATMOS_PHY_MP_R0A

    integer :: nnspc, nnbin
    integer :: nn, mm, mmyu, nnyu
    integer :: myu, nyu, i, j, k, n, ierr

    ATMOS_PHY_MP_RHOA = rhoa
    ATMOS_PHY_MP_EMAER = emaer
    ATMOS_PHY_MP_RAMIN = rasta
    ATMOS_PHY_MP_RAMAX = raend
    ATMOS_PHY_MP_R0A = r0a
    ATMOS_PHY_MP_FLAG_REGENE = flg_regeneration
    ATMOS_PHY_MP_FLAG_NUCLEAT = flg_nucl
    ATMOS_PHY_MP_FLAG_SFAERO = flg_sf_aero

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for SBM (warm cloud)'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP, Check!'
     call PRC_MPIstop
    end if
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    rhoa = ATMOS_PHY_MP_RHOA
    emaer = ATMOS_PHY_MP_EMAER
    rasta = ATMOS_PHY_MP_RAMIN
    raend = ATMOS_PHY_MP_RAMAX
    r0a   = ATMOS_PHY_MP_R0A
    flg_regeneration = ATMOS_PHY_MP_FLAG_REGENE
    flg_nucl = ATMOS_PHY_MP_FLAG_NUCLEAT
    flg_sf_aero = ATMOS_PHY_MP_FLAG_SFAERO

    call fio_register_file(n,fname_micpara)
    fid_micpara = n
    !--- cloud
    open ( fid_micpara, file = fname_micpara, form = 'formatted', status = 'old' )

    read( fid_micpara,* ) nnspc, nnbin

    ! grid parameter
    if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
    do n = 1, nbin
      read( fid_micpara,* ) nn, xctr( n ), radc( n )
      if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)')  "Radius of ", n, "th cloud bin= ", radc( n ) , "[m]"
    end do
    do n = 1, nbin+1
      read( fid_micpara,* ) nn, xbnd( n )
    end do
    read( fid_micpara,* ) dxmic
    if( IO_L ) write(IO_FID_LOG,*)  '*** Width of Cloud SDF= ', dxmic

    ! capacity
    do myu = 1, 7
     do n = 1, nbin
      read( fid_micpara,* ) mmyu, nn, cctr( myu,n )
     end do
     do n = 1, nbin+1
      read( fid_micpara,* ) mmyu, nn, cbnd( myu,n )
     end do
    end do

    ! collection kernel 
    do myu = 1, 7
     do nyu = 1, 7
      do i = 1, nbin
       do j = 1, nbin
        read( fid_micpara,* ) mmyu, nnyu, mm, nn, ck( myu,nyu,i,j )
       end do
      end do
     end do
    end do

    ! terminal velocity
    do myu = 1, 7
     do n = 1, nbin
      read( fid_micpara,* ) mmyu, nn, vt( myu,n )
     end do
    end do

    ! bulk density
    do myu = 1, 7
     do n = 1, nbin
      read( fid_micpara,* ) mmyu, nn, br( myu,n )
     end do
    end do

    close ( fid_micpara )

    !--- aerosol ( CCN ) (not supported)
    xasta = log( rhoa*4.D0/3.D0*pi * ( rasta )**3 )
    xaend = log( rhoa*4.D0/3.D0*pi * ( raend )**3 )

    dxaer = ( xaend-xasta )/nccn

    do n = 1, nccn+1
     xabnd( n ) = xasta + dxaer*( n-1 )
    end do
    do n = 1, nccn
     xactr( n ) = ( xabnd( n )+xabnd( n+1 ) )*0.5D0
     rada( n )  = ( exp( xactr( n ) )*3.D0/4.D0/pi/rhoa )**( OneovThird )
    end do

    if( flg_sf_aero ) then
     if ( CZ(KS) >= 10.D0 ) then
          R10M1 = 10.D0 / CZ(KS) * 0.5D0 ! scale with height
          R10M2 = 10.D0 / CZ(KS) * 0.5D0 ! scale with height
          R10H1 = 1.D0 * 0.5D0
          R10H2 = 1.D0 * 0.5D0
          R10E1 = 1.D0 * 0.5D0
          R10E2 = 1.D0 * 0.5D0
          K10_1 = KS
          K10_2 = KS
     else
       k = 1
       do while ( CZ(k) < 10.D0 )
          k = k + 1
          K10_1 = k
          K10_2 = k + 1
          R10M1 = ( CZ(k+1) - 10.D0 ) / CDZ(k)
          R10M2 = ( 10.D0   - CZ(k) ) / CDZ(k)
          R10H1 = ( CZ(k+1) - 10.D0 ) / CDZ(k)
          R10H2 = ( 10.D0   - CZ(k) ) / CDZ(k)
          R10E1 = ( CZ(k+1) - 10.D0 ) / CDZ(k)
          R10E2 = ( 10.D0   - CZ(k) ) / CDZ(k)
       enddo
     endif
    end if

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_process, only: &
       PRC_myrank
    use mod_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP, &
       TIME_NOWSEC
    use mod_grid, only : &
       GRID_CZ,  &
       GRID_FZ,  &
       GRID_CDZ, &
       GRID_FDZ
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       ATMOS_vars_total,   &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    implicit none

    real(8) :: dz (KA)
    real(8) :: dzh(KA)
    real(8) :: dt
    real(8) :: ct

    real(8) :: rho      (IA,JA,KA)
    real(8) :: temp_mp  (IA,JA,KA)
    real(8) :: pres_mp  (IA,JA,KA)
    real(8) :: gdgc     (IA,JA,KA,nbin) !-- SDF of hydrometeors [kg/m^3/unit ln(r)]
    real(8) :: qv_mp    (IA,JA,KA)      !-- Qv [kg/kg]
    real(8) :: wfall( KA )  
    real(8) :: gdga     (IA,JA,KA,nccn) !-- SDF of aerosol (not supported)
    
    real(8) :: ssliq, ssice, sum1, mixrate, rtotal, sum2
    integer :: n, myu, k, i, j, ij, iq, iv
    logical, save :: ofirst_sdfa = .true.

    real(8) :: VELX(IA,JA)
    real(8) :: VELY(IA,JA)
    real(8) :: SFLX_AERO(IA,JA,nccn)
    real(8) :: Uabs, bparam
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("MICROPHYSICS")
#endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(SBM-liquid only)'

    if( flg_sf_aero ) then
     do j = JS-2, JE+2
     do i = IS-2, IE+1
       VELX(i,j) = MOMX(K10_1,i,j) / ( DENS(K10_1,i+1,j)+DENS(K10_1,i,j) ) * R10M1 &
                 + MOMX(K10_2,i,j) / ( DENS(K10_2,i+1,j)+DENS(K10_2,i,j) ) * R10M2
     enddo
     enddo

     do j = JS-2, JE+1
     do i = IS-2, IE+2
       VELY(i,j) = MOMY(K10_1,i,j) / ( DENS(K10_1,i,j+1)+DENS(K10_1,i,j) ) * R10M1 &
                 + MOMY(K10_2,i,j) / ( DENS(K10_2,i,j+1)+DENS(K10_2,i,j) ) * R10M2
     enddo
     enddo
    end if

    call TIME_rapstart('MPX ijkconvert')
    dz (:) = GRID_CDZ(:)
    dzh(1) = GRID_FDZ(1)
    dzh(2:KA) = GRID_FDZ(1:KA-1)

    dt = TIME_DTSEC_ATMOS_PHY_MP
    ct = TIME_NOWSEC

    gdgc(:,:,:,:) = 0.d0
    gdga(:,:,:,:) = 0.d0
    pres_mp(:,:,:) = 0.d0
    temp_mp(:,:,:) = 0.d0
    qv_mp(:,:,:) = 0.d0
    do j = JS-1, JE
    do i = IS-1, IE
       do k = KS, KE
          mixrate = 0.d0
          do n = 1, QA
           mixrate = mixrate + QTRC(k,i,j,n)
          end do
          mixrate = 1.d0 - mixrate
          rtotal = CONST_Rdry*mixrate + CONST_Rvap*QTRC(k,i,j,I_QV)
          rho    (i,j,k) = DENS(k,i,j)*mixrate
          pres_mp(i,j,k) = CONST_PRE00 * &
             ( RHOT(k,i,j) * rtotal / CONST_PRE00 )**CONST_CPovCV
          temp_mp(i,j,k) = &
                  pres_mp(i,j,k)/( DENS(k,i,j)*rtotal )
          qv_mp(i,j,k)   = QTRC(k,i,j,I_QV)
          do n = 1, nbin
           gdgc( i,j,k,n ) = &
               QTRC(k,i,j,n+1)*DENS(k,i,j)*mixrate/dxmic
          end do
          do n = 1, nccn
           gdga( i,j,k,n ) = &
               QTRC(k,i,j,n+nbin+1)*DENS(k,i,j)*mixrate/dxaer
          end do
          !--- store initial SDF of aerosol 
          if( ofirst_sdfa ) then
           sum2 = 0.d0
           do n = 1, nccn
             marate( n ) = gdga(i,j,k,n)/exp( xactr( n ) )
             sum2 = sum2 + gdga(i,j,k,n)/exp( xactr( n ) )
           end do
           if( sum2 /= 0.d0 ) then
            marate( 1:nccn ) = marate( 1:nccn )/sum2
            ofirst_sdfa = .false.
           end if
          end if
       enddo
       do k = 1, KS-1
          mixrate = 0.d0
          do n = 1, QA
           mixrate = mixrate + QTRC(KS,i,j,n)
          end do
          mixrate = 1.d0 - mixrate
          rtotal = CONST_Rdry*mixrate + CONST_Rvap*QTRC(KS,i,j,I_QV)
          rho    (i,j,k) = DENS(KS,i,j)*mixrate
          pres_mp(i,j,k) = CONST_PRE00 * &
             ( RHOT(KS,i,j) * rtotal / CONST_PRE00 )**CONST_CPovCV
          temp_mp(i,j,k) = &
                  pres_mp(i,j,k)/( DENS(KS,i,j)*rtotal )
          qv_mp(i,j,k)   = QTRC(KS,i,j,I_QV)
          do n = 1, nbin
           gdgc( i,j,k,n ) = &
                QTRC(KS,i,j,n+1)*DENS(KS,i,j)*mixrate/dxmic
          end do
          do n = 1, nccn
           gdga( i,j,k,n ) = &
                QTRC(KS,i,j,n+nbin+1)*DENS(KS,i,j)*mixrate/dxaer
          end do
       enddo
       do k = KE+1, KA
          mixrate = 0.d0
          do n = 1, QA
           mixrate = mixrate + QTRC(KE,i,j,n)
          end do
          mixrate = 1.d0 - mixrate
          rtotal = CONST_Rdry*mixrate + CONST_Rvap*QTRC(KE,i,j,I_QV)
          rho   (i,j,k) = DENS(KE,i,j)*mixrate
          pres_mp(i,j,k) = CONST_PRE00 * &
             ( RHOT(KE,i,j) * rtotal / CONST_PRE00 )**CONST_CPovCV
          temp_mp(i,j,k) = &
                  pres_mp(i,j,k)/( DENS(KE,i,j)*rtotal )
          qv_mp(i,j,k) = QTRC(KE,i,j,I_QV)
          do n = 1, nbin
           gdgc( i,j,k,n ) = &
              QTRC(KE,i,j,n+1)*DENS(KE,i,j)*mixrate/dxmic
          end do
          do n = 1, nccn
           gdga( i,j,k,n ) = &
              QTRC(KE,i,j,n+nbin+1)*DENS(KE,i,j)*mixrate/dxaer
          end do
       enddo
    enddo
    enddo

    call TIME_rapend  ('MPX ijkconvert')

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE

      call getsups &
        (  qv_mp(i,j,k), temp_mp(i,j,k), pres_mp(i,j,k), &
           ssliq, ssice )
      sum1 = 0.d0
      do n = 1, nbin
        sum1 = sum1 + gdgc(i,j,k,n)*dxmic
      end do

      if( ssliq > 0.d0 .or. ssice > 0.d0 .or. sum1 > cldmin ) then
       call mp_hbinw_evolve                    &
            ( pres_mp(i,j,k), rho(i,j,k), dt,  &  !--- in
              gdgc(i,j,k,1:nbin),              &  !--- inout
              gdga(i,j,k,1:nccn),              &  !--- inout  for aerosol tracer
              qv_mp(i,j,k), temp_mp(i,j,k)     )  !--- inout
      end if
    
    end do
    end do
    end do

    !--- gravitational falling
     do j = JS, JE
      do i = IS, IE
       sum1 = 0.d0
       do k = KS, KE
        do n = 1, nbin
         sum1 = sum1 + gdgc(i,j,k,n)*dxmic
        end do
       end do
       if( sum1 > cldmin ) then
        do n = 1, nbin
         wfall( 1:KA ) = -vt( il,n )
         call  advec_1d             &
              ( dz, dzh,            & !--- in
                wfall,              & !--- in
                dt,                 & !--- in
                gdgc( i,j,1:KA,n )  ) !--- inout
         do k = KS, KE
          if( gdgc(i,j,k,n) < cldmin ) then
           gdgc(i,j,k,n) = max( gdgc(i,j,k,n),0.d0 ) 
          end if
         end do
        end do
       end if
      end do
     end do

    !--- SURFACE FLUX by Monahan et al. (1986)
    if( flg_sf_aero ) then
     do j = JS-1, JE
     do i = IS-1, IE
       Uabs = sqrt(  ( ( VELX(i,j) + VELX(i-1,j  ) ) * 0.5d0 )**2 &
                   + ( ( VELY(i,j) + VELY(i  ,j-1) ) * 0.5d0 )**2 )
       do n = 1, nccn 
        if( rada( n ) <= 2.0d-5 .and. rada( n ) >= 3.0d-7 ) then
         bparam = ( 0.38D0 - log( rada( n ) ) )/0.65D0
         SFLX_AERO(i,j,n) = 1.373d0 * Uabs**( 3.41d0 ) * rada( n )**( -3.D0 ) &
                          * ( 1.D0 + 0.057D0 * rada( n )**( 1.05D0 ) ) &
                          * 10.d0**( 1.19D0 * exp( -bparam*bparam ) ) 
         ! convert from [#/m^2/um/s] -> [kg/m^3/unit log (m)]
         SFLX_AERO(i,j,n) = SFLX_AERO(i,j,n) / rho(i,j,KS) &
                          / GRID_CDZ(KS) * rada( n ) / 3.D0 * dt * exp( xactr( n ) )
         gdga(i,j,KS,n) = gdga(i,j,KS,n)+SFLX_AERO(i,j,n)/dxaer
        end if     
       end do     
     end do
     end do    
    end if

    call TIME_rapstart('MPX ijkconvert')
    do j = JS, JE
     do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv_mp(i,j,k)  
          do n = 1, nbin
            QTRC(k,i,j,n+1)= &
                gdgc(i,j,k,n)/rho(i,j,k)*dxmic
            if( QTRC(k,i,j,n+1) <= eps ) then 
              QTRC(k,i,j,n+1) = 0.d0
            end if
          end do
          do n = 1, nccn
            QTRC(k,i,j,n+1+nbin)= &
                gdga(i,j,k,n)/rho(i,j,k)*dxaer
          end do
          mixrate = 0.d0
          do n = 1, QA
           mixrate = mixrate + QTRC(k,i,j,n)
          end do
          mixrate = 1.d0 - mixrate
          DENS(k,i,j) = rho(i,j,k)/mixrate
          RHOT(k,i,j) = temp_mp(i,j,k)* &
            ( CONST_PRE00/pres_mp(i,j,k) )**( CONST_RovCP )*DENS(k,i,j)
       enddo
     enddo
    enddo

    call TIME_rapend  ('MPX ijkconvert')

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( MOMZ(:,:,:), 2 )
    call COMM_vars8( MOMX(:,:,:), 3 )
    call COMM_vars8( MOMY(:,:,:), 4 )
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( DENS(:,:,:), 1 )
    call COMM_wait ( MOMZ(:,:,:), 2 )
    call COMM_wait ( MOMX(:,:,:), 3 )
    call COMM_wait ( MOMY(:,:,:), 4 )
    call COMM_wait ( RHOT(:,:,:), 5 )

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("MICROPHYSICS")
#endif

    call ATMOS_vars_total

    return
  end subroutine ATMOS_PHY_MP
  !-----------------------------------------------------------------------------
  subroutine getsups        &
       ( qvap, temp, pres,  & !--- in
         ssliq, ssice )       !--- out
  !
  real(8), intent(in) :: qvap  !  specific humidity [ kg/kg ]
  real(8), intent(in) :: temp  !  temperature [ K ]
  real(8), intent(in) :: pres  !  pressure [ Pa ]
  !
  real(8), intent(out) :: ssliq
  real(8), intent(out) :: ssice
  !
  real(8) :: epsl, rr, evap, esatl, esati
  !
  epsl = CONST_Rdry/CONST_Rvap
  !
  rr = qvap / ( 1.D0-qvap )
  evap = rr*pres/( epsl+rr )

  esatl = fesatl( temp )
  esati = fesati( temp )

  ssliq = evap/esatl - 1.D0
  ssice = evap/esati - 1.D0

  return

  end subroutine getsups
  !-----------------------------------------------------------------------------
  subroutine mp_hbinw_evolve        &
      ( pres, dens,                 & !--- in
        dtime,                      & !--- in
        gc,                         & !--- inout
        ga,                         & !--- inout
        qvap, temp                  ) !--- inout

  real(8), intent(in) :: pres   !  pressure
  real(8), intent(in) :: dens   !  density
  real(8), intent(in) :: dtime  !  time interval

  real(8), intent(inout) :: gc( nbin )
  real(8), intent(inout) :: ga( nccn )  !--- aerosol SDF (not supported)
  real(8), intent(inout) :: qvap  !  specific humidity
  real(8), intent(inout) :: temp  !  temperature
  integer(4) :: n
  !
  !
  !--- nucleat
  call nucleat                  &
         ( dens, pres, dtime,   & !--- in
           gc, ga, qvap, temp )   !--- inout

  !--- condensation / evaporation
  call cndevpsbl                &
         ( dtime,               & !--- in
           dens, pres,          & !--- in
           gc, ga, qvap, temp   ) !--- inout

  !--- collision-coagulation
  call  collmain                &
          ( dtime, temp,        & !--- in
            gc                  ) !--- inout

  return

  end subroutine mp_hbinw_evolve 
  !-----------------------------------------------------------------------------
  subroutine nucleat             &
      ( dens, pres, dtime,       & !--- in
        gc, ga, qvap, temp  )      !--- inout
  !
  !  liquid nucleation from aerosol particle
  !
  !
  real(8), intent(in) :: dens   !  density  [ kg/m3 ]
  real(8), intent(in) :: pres   !  pressure [ Pa ]
  real(8), intent(in) :: dtime
  !
  real(8), intent(inout) :: gc( nbin )  !  SDF ( hydrometeors )
  real(8), intent(inout) :: ga( nccn )  !  SDF ( aerosol ) : mass
  real(8), intent(inout) :: qvap  !  specific humidity [ kg/kg ]
  real(8), intent(inout) :: temp  !  temperature [ K ]
  !
  real(8), parameter :: rhow  = CONST_DWATR
  real(8), parameter :: qlevp = CONST_LH0
  real(8), parameter :: cp    = CONST_CPdry
  real(8), parameter :: rvap  = CONST_Rvap
  !--- local
  real(8) :: gan( nccn )  !  SDF ( aerosol ) : number
  real(8) :: ssliq, ssice, delcld
  real(8) :: sumold, sumnew, acoef, bcoef, xcrit, rcrit
  real(8) :: ractr, rcld, xcld, part, vdmp, dmp
  integer :: n, nc, ncrit
  integer, allocatable, save :: ncld( : )
  logical, save :: ofirst = .true.
  !
  !--- use for aerosol coupled model
  real(8), parameter :: sigma = 7.5D-02  ! water surface tension [ N/m2 ]
  real(8), parameter :: vhfct = 2.D0     ! van't hoff factor
  !
  !--- supersaturation
  call  getsups               &
          ( qvap, temp, pres, & !--- in
            ssliq, ssice  )     !--- out

  if ( ssliq <= 0.D0 ) return
  !--- use for aerosol coupled model
  !--- mass -> number
  do n = 1, nccn
    gan( n ) = ga( n )/exp( xactr( n ) )
  end do

  acoef = 2.D0*sigma/rvap/rhow/temp
  bcoef = vhfct* rhoa/rhow * emwtr/emaer

  !--- relationship of bin number
  if ( ofirst ) then
    allocate ( ncld( 1:nccn ) )
    do n = 1, nccn
      ractr = ( exp( xactr( n ) )*3.D0/4.D0/pi/rhoa )**( OneovThird )
      rcld  = sqrt( 3.D0*bcoef*ractr**3 / acoef )
      xcld  = log( rhow * 4.D0*pi/3.D0 * rcld**3 )
     if( flg_nucl ) then
      ncld( n ) = 1
     else
      ncld( n ) = int( ( xcld-xctr( 1 ) )/dxmic ) + 1
      ncld( n ) = min( max( ncld( n ),1 ),nbin )
     end if
    end do
    ofirst = .false.
  end if

  !--- nucleation
  do n = nccn, 1, -1

    call  getsups                &
            ( qvap, temp, pres,  & !--- in
              ssliq, ssice )       !--- out

    if ( ssliq <= 0.D0 ) exit
    !--- use for aerosol coupled model
    acoef = 2.D0*sigma/rvap/rhow/temp
    rcrit = acoef/3.D0 * ( 4.D0/bcoef )**( OneovThird ) / ssliq**( 2.D0/3.D0 )
    xcrit = log( rhoa * 4.D0*pi/3.D0 * rcrit**3 )
    ncrit = int( ( xcrit-xabnd( 1 ) )/dxaer ) + 1

    if ( n == ncrit ) then
      part = ( xabnd( ncrit+1 )-xcrit )/dxaer
    else if ( n > ncrit ) then
      part = 1.D0
    else
      exit
    end if

    nc = ncld( n )
    dmp = part*gan( n )*dxaer*exp( xctr( nc ) )
    dmp = min( dmp,qvap*dens )
    gc( nc ) = gc( nc ) + dmp/dxmic
    gan( n ) = gan( n ) - dmp/dxaer/exp( xctr( nc ) )
    gan( n ) = max( gan( n ), 0.D0 )
    qvap = qvap - dmp/dens
    qvap = max( qvap,0.D0 )
    temp = temp + dmp/dens*qlevp/cp
  end do

  !--- number -> mass
  do n = 1, nccn
    ga( n ) = gan( n )*exp( xactr( n ) )
  end do
  !
  return
  !
  end subroutine nucleat
  !-----------------------------------------------------------------------------
  subroutine  cndevpsbl     &
      ( dtime,              & !--- in
        dens, pres,         & !--- in
        gc, ga, qvap, temp  ) !--- inout
  !
  real(8), intent(in) :: dtime
  real(8), intent(in) :: dens   !  atmospheric density [ kg/m3 ]
  real(8), intent(in) :: pres   !  atmospheric pressure [ Pa ]
  real(8), intent(inout) :: gc( nbin )  ! Size Distribution Function
  real(8), intent(inout) :: ga( nccn )  !  SDF ( aerosol ) : mass
  real(8), intent(inout) :: qvap    !  specific humidity [ kg/kg ]
  real(8), intent(inout) :: temp    !  temperature [ K ]
  !
  !--- local variables
  integer :: iflg( il ), n, iliq, iice
  real(8) :: csum( il )
  real(8) :: regene_gcn
  !
  !
  iflg( : ) = 0
  csum( : ) = 0.D0
  regene_gcn = 0.d0
  do n = 1, nbin
    csum( il ) = csum( il )+gc( n )*dxmic
  end do

  if ( csum( il ) > cldmin ) iflg( il ) = 1

  iliq = iflg( il )

  if ( iliq == 1 ) then
      call  liqphase            &
              ( dtime, iliq,    & !--- in
                dens, pres,     & !--- in
                gc, qvap, temp, & !--- inout
                regene_gcn      ) !--- out
     !--- regeneration of aerosol
      if( flg_regeneration ) then
       call faero( regene_gcn,  & !--- in
                   ga           ) !--- inout
      end if
  end if
  !
  end subroutine cndevpsbl
  !-----------------------------------------------------------------------------
  subroutine liqphase   &
      ( dtime, iflg,    & !--- in
        dens, pres,     & !--- in
        gc, qvap, temp, & !--- inout
        regene_gcn      ) !--- out
  !
  real(8), intent(in) :: dtime
  integer, intent(in) :: iflg
  real(8), intent(in) :: dens, pres
  real(8), intent(inout) :: gc( nbin ), qvap, temp
  real(8), intent(out) :: regene_gcn
  !
  real(8), parameter :: rhow  = CONST_DWATR
  real(8), parameter :: qlevp = CONST_LH0
  real(8), parameter :: rvap  = CONST_Rvap
  real(8), parameter :: cp    = CONST_CPdry
  !
  !--- local variables
  integer :: n, nloop, ncount
  real(8) :: gclold, gclnew, gtliq, umax, dtcnd
  real(8) :: sumliq, cefliq, a, sliqtnd, cndmss, ssliq, ssice
  real(8) :: gcn( nbin ), gdu( nbin+1 ), gcnold( nbin )
  real(8), parameter :: cflfct = 0.5D0
  real(8) :: old_sum_gcn, new_sum_gcn
  !
  gclold = 0.D0
  do n = 1, nbin
    gclold = gclold + gc( n )
  end do
  gclold = gclold * dxmic
  !
  !------- mass -> number
  gcn( 1:nbin ) = gc( 1:nbin ) / exp( xctr( 1:nbin ) )

  !
  !------- CFL condition
  call  getsups                &
          ( qvap, temp, pres,  & !--- in
            ssliq, ssice )       !--- out  

  gtliq = gliq( pres,temp )
  umax = cbnd( il,1 )/exp( xbnd( 1 ) )*gtliq*abs( ssliq )
  dtcnd = cflfct*dxmic/umax
  nloop = int( dtime/dtcnd ) + 1
  dtcnd = dtime / nloop
  !
  ncount = 0
  regene_gcn = 0.d0
  !------- loop
  1000 continue

  !----- matrix for supersaturation tendency
  call  getsups                &
          ( qvap, temp, pres,  & !--- in
            ssliq, ssice )       !--- out

  gtliq = gliq( pres,temp )
  sumliq = 0.D0
  old_sum_gcn = 0.D0
  do n = 1, nbin
    sumliq = sumliq + gcn( n )*cctr( il,n )
  end do
  sumliq = sumliq * dxmic
  cefliq = ( ssliq+1.D0 )*( 1.D0/qvap + qlevp*qlevp/cp/rvap/temp/temp )
  a = - cefliq*sumliq*gtliq/dens
  !
  !----- supersaturation tendency
  if ( abs( a*dtcnd ) >= 0.1D0 ) then
    sliqtnd = ssliq*( exp( a*dtcnd )-1.D0 )/( a*dtcnd )
  else
    sliqtnd = ssliq
  end if
  !
  !----- change of SDF
  gdu( 1:nbin+1 ) = cbnd( il,1:nbin+1 )/exp( xbnd( 1:nbin+1 ) )*gtliq*sliqtnd
  gcnold( : ) = gcn( : )
  call  advection              &
          ( dtcnd,             & !--- in
            gdu( 1:nbin+1 ),   & !--- in
            gcn( 1:nbin ),     & !--- inout
            regene_gcn         ) !--- inout
  !
  !----- new mass
  gclnew = 0.D0
  new_sum_gcn = 0.d0
  do n = 1, nbin
    gclnew = gclnew + gcn( n )*exp( xctr( n ) )
    old_sum_gcn = old_sum_gcn + gcnold( n )*dxmic
    new_sum_gcn = new_sum_gcn + gcn( n )*dxmic
  end do
!  if( sliqtnd < 0.d0 .and. gcnold( 1 ) >= cldmin ) then
!   regene_gcn = regene_gcn + ( old_sum_gcn-new_sum_gcn )
!   regene_gcn = max( regene_gcn,0.d0 )
!  end if
     
  gclnew = gclnew*dxmic
  !
  !----- change of humidity and temperature
  cndmss = gclnew - gclold
  qvap = qvap - cndmss/dens
  temp = temp + cndmss/dens*qlevp/cp
  !
  gclold = gclnew
  !
  !----- continue/end
  ncount = ncount + 1
  if ( ncount < nloop ) go to 1000
  !
  !------- number -> mass
  gc( 1:nbin ) = gcn( 1:nbin )*exp( xctr( 1:nbin ) )
  ! 
  end subroutine liqphase
  !-------------------------------------------------------------------------------
  subroutine advection  &
       ( dtime,         & !--- in
         gdu,           & !--- in
         gdq, regene    ) !--- inout
  !
  real(8), intent(in) :: dtime, gdu( 1:nbin+1 )
  real(8), intent(inout) :: gdq( 1:nbin ), regene
  !
  !--- local variables
  real(8) :: delx 
  real(8) :: qadv( -1:nbin+2 ), uadv( 0:nbin+2 )
  real(8) :: flq( 1:nbin+1 )
  integer, parameter :: ldeg = 2
  real(8), parameter :: epsl = 1.D-10
  real(8) :: acoef( 0:nbin+1,0:ldeg ), sum
  real(8) :: crn( 0:nbin+2 )
  real(8) :: aip( 0:nbin+1 ), aim( 0:nbin+1 ), ai( 0:nbin+1 )
  real(8) :: cmins, cplus
  integer :: i, n
  !
  !
  delx = dxmic
  do i = 1, nbin
    qadv( i ) = gdq( i )
  end do
  qadv( -1 )     = 0.D0
  qadv( 0  )     = 0.D0
  qadv( nbin+1 ) = 0.D0
  qadv( nbin+2 ) = 0.D0

  do i = 1, nbin+1
    uadv( i ) = gdu( i )
  end do
  uadv( 0  ) = 0.D0
  uadv( nbin+2 ) = 0.D0

  !--- flux
  do i = 0, nbin+1
    acoef( i,0 ) = - ( qadv( i+1 )-26.D0*qadv( i )+qadv( i-1 ) ) / 24.D0
    acoef( i,1 ) = ( qadv( i+1 )-qadv( i-1 ) ) * 0.5D0
    acoef( i,2 ) = ( qadv( i+1 )-2.D0*qadv( i )+qadv( i-1 ) ) * 0.5D0
  end do

  crn( 0:nbin+2 ) = uadv( 0:nbin+2 )*dtime/delx

  do i = 0, nbin+1
    cplus = ( crn( i+1 ) + abs( crn( i+1 ) ) ) * 0.5D0
    sum = 0.D0
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.D0**( n+1 )  &
                *( 1.D0-( 1.D0-2.D0*cplus )**( n+1 ) )
    end do
    aip( i ) = max( sum,0.D0 )
  end do

  do i = 0, nbin+1
    cmins = - ( crn( i ) - abs( crn( i ) ) ) * 0.5D0
    sum = 0.D0
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.D0**( n+1 ) * (-1)**n &
                *( 1.D0-( 1.D0-2.D0*cmins )**( n+1 ) )
    end do
    aim( i ) = max( sum,0.D0 )
  end do

  do i = 0, nbin+1
    sum = 0.D0
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.D0**( n+1 ) * ( (-1)**n+1 )
    end do
    ai( i ) = max( sum,aip( i )+aim( i )+epsl )
  end do

  do i = 1, nbin+1
    flq( i ) = ( aip( i-1 )/ai( i-1 )*qadv( i-1 )  &
                -aim( i   )/ai( i   )*qadv( i   ) )*delx/dtime
  end do

  if( gdu( 1 ) < 0.d0 ) then
   regene = regene+( -flq( 1 )*dtime/delx )
  end if

  do i = 1, nbin
    gdq( i ) = gdq( i ) - ( flq( i+1 )-flq( i ) )*dtime/delx
  end do

  end subroutine advection
  !-------------------------------------------------------------------------------
  function  gliq( pres,temp )
  !
  use mod_const, only : rair  => CONST_Rdry,  &
                        rvap  => CONST_Rvap,  &
                        qlevp => CONST_LH0
  !
  real(8), intent(in) :: pres, temp
  real(8) :: gliq
  !
  real(8) :: emu, dens, cefd, cefk, f1, f2
  real(8), parameter :: fct = 1.4D3
  !
  emu = fmyu( temp )
  dens = pres/rair/temp
  cefd = emu/dens
  cefk = fct*emu

  f1 = rvap*temp/fesatl( temp )/cefd
  f2 = qlevp/cefk/temp*( qlevp/rvap/temp - 1.D0 )

  gliq = 4.D0*pi/( f1+f2 )

  end function gliq
  !-------------------------------------------------------------------------------
  function fmyu( temp )
  !
  !
  real(8), intent(in) :: temp
  real(8) :: fmyu
  real(8), parameter :: tmlt = CONST_TMELT
  !
  real(8), parameter :: a = 1.72D-05, b = 3.93D2, c = 1.2D02
  !
  fmyu = a*( b/( temp+c ) )*( temp/tmlt )**1.5D0
  !
  end function fmyu
  !-------------------------------------------------------------------------------
  function fesatl( temp )
  !
  real(8), intent(in) :: temp
  real(8) :: fesatl
  !
  real(8), parameter :: qlevp = CONST_LH0
  real(8), parameter :: rvap  = CONST_Rvap
  real(8), parameter :: esat0 = CONST_PSAT0
  real(8), parameter :: temp0 = CONST_TEM00
  !
  fesatl = esat0*exp( qlevp/rvap*( 1.D0/temp0 - 1.D0/temp ) )
  !
  return

  end function fesatl
  !-----------------------------------------------------------------------
  function fesati( temp )
  !
  real(8), intent(in) :: temp
  real(8) :: fesati
  !
  real(8), parameter :: qlsbl = CONST_LHS0
  real(8), parameter :: rvap  = CONST_Rvap
  real(8), parameter :: esat0 = CONST_PSAT0
  real(8), parameter :: temp0 = CONST_TEM00
  !
  fesati = esat0*exp( qlsbl/rvap*( 1.D0/temp0 - 1.D0/temp ) )
  !
  return

  end function fesati
  !-------------------------------------------------------------------------------
 subroutine collmain &
      ( dtime, temp, & !--- in
        gc           ) !--- inout

  !--- Y.Sato added on 2012/3/2
!  use mod_randomset, only : mbinp, mspcp,      &
!                            wgtbin, wgtspc,    &
!                            blrg, bsml,        &
!                            slrg, ssml,        &
!                            mbinp, msetp,      &
!                            ntimep, rndm_flgp, &
!                            i_grid, k_grid
  !
  real(8), intent(in) :: dtime
  real(8), intent(in) :: temp
  real(8), intent(inout) :: gc( nbin )
  !
  !--- local variables
  integer :: iflg( 1 ), n, isml, ilrg, irsl
  real(8) :: csum( 1 )
  real(8), parameter :: tcrit = 271.15D0
  !--- Y.Sato added on 2012/3/2
  integer :: det, s
!  real(8) :: snums( mspcp ), snuml( mspcp )
  !
  !--- judgement of particle existence
    iflg( : ) = 0
    csum( : ) = 0.D0
    do n = 1, nbin
      csum( il ) = csum( il ) + gc( n )*dxmic
    end do
    if ( csum( il ) > cldmin ) iflg( il ) = 1
  !
  !--- interaction
!   if( rndm_flgp <= 1 ) then

     if ( iflg( il ) == 1 ) then
!           if ( rndm_flgp == 1 ) then  !--- stochastic method
!            !--- Y.Sato added on 2012/3/2
!            call r_collcoag                    &
!                    ( nspc, nbin,              & !--- in
!                      xctr, dxmic, dtime, ck,  & !--- in
!                      isml, ilrg, irsl,        & !--- in
!                      i_grid, k_grid, wgtbin,  & !--- in
!                      gc   )                     !--- inout
!           else                       !--- default method
        call  collcoag        &
                ( dtime,      & !--- in
                  il, il, il, & !--- in
                  gc          ) !--- inout
!           end if
     end if

!   else if( rndm_flgp == 2 ) then
!    !--- Y.Sato added on 2012/3/2 ( not supported )
!    !--- use random number for both bin and spc
!    det = mod( ntimep*i_grid*k_grid, msetp ) + 1
!    snums( 1:mspcp ) = ssml( det, 1:mspcp )
!    snuml( 1:mspcp ) = slrg( det, 1:mspcp )
!    do s = 1, mspcp
!     isml = snums( s )
!      if ( iflg( isml ) == 1 ) then
!        ilrg = snuml( s )
!         if ( iflg( ilrg ) == 1 ) then
!           call r_collcoag                    &
!             ( dtime,                         & !--- in
!               isml, ilrg, irsl,              & !--- in
!               i_grid, k_grid, wgtbin*wgtspc, & !--- in
!               gc   )                           !--- inout 
!         end if
!      end if
!    end do
!
!   end if
!  !
  !
  return
  !
  end subroutine collmain
  !-------------------------------------------------------------------------------
  subroutine  collcoag            &
       ( dtime,                   & !--- in
         isml, ilrg, irsl,        & !--- in
         gc   )                     !--- inout
  !
  real(8), intent(in) :: dtime
  integer, intent(in) :: isml, ilrg, irsl
  real(8), intent(inout) :: gc( nbin )
  !
  !--- local variables
  integer :: i, j, k, l
  real(8) :: xi, xj, xnew, dmpi, dmpj, frci, frcj
  real(8) :: gprime, gprimk, wgt, crn, sum, flux
  integer, parameter :: ldeg = 2
  real(8) :: acoef( 0:ldeg )
  real(8), parameter :: dmpmin = 1.D-01
  !
  small : do i = 1, nbin-1
    if ( gc( i ) <= cldmin ) cycle small
  large : do j = i+1, nbin
    if ( gc( j ) <= cldmin ) cycle large

    xi = exp( xctr( i ) )
    xj = exp( xctr( j ) )
    xnew = log( xi+xj )
    k = int( ( xnew-xctr( 1 ) )/dxmic ) + 1
    k = max( max( k,j ),i )
    if ( k >= nbin ) cycle small

    dmpi = ck( isml,ilrg,i,j )*gc( j )/xj*dxmic*dtime
    dmpj = ck( isml,ilrg,i,j )*gc( i )/xi*dxmic*dtime

    if ( dmpi <= dmpmin ) then
      frci = gc( i )*dmpi
    else
      frci = gc( i )*( 1.D0-exp( -dmpi ) )
    end if

    if ( dmpj <= dmpmin ) then
      frcj = gc( j )*dmpj
    else
      frcj = gc( j )*( 1.D0-exp( -dmpj ) )
    end if

    gc( i ) = max( gc( i )-frci, 0.D0 )
    gc( j ) = max( gc( j )-frcj, 0.D0 )

    gprime = frci+frcj
    if ( gprime <= 0.D0 ) cycle large
    gprimk = gc( k ) + gprime
    wgt = gprime / gprimk
    crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )

    acoef( 0 ) = -( gc( k+1 )-26.D0*gprimk+gc( k-1 ) )/24.D0
    acoef( 1 ) = ( gc( k+1 )-gc( k-1 ) ) *0.5D0
    acoef( 2 ) = ( gc( k+1 )-2.D0*gprimk+gc( k-1 ) ) *0.5D0

    sum = 0.D0
    do l = 0, ldeg
      sum = sum + acoef( l )/( l+1 )/2.D0**( l+1 )   &
                *( 1.D0-( 1.D0-2.D0*crn )**( l+1 ) )
    end do

    flux = wgt*sum
    flux = min( max( flux,0.D0 ),gprime )

    gc( k ) = gprimk - flux
    gc( k+1 ) = gc( k+1 ) + flux

  end do large
  end do small
  !  
  return
  !
  end subroutine collcoag
  !-------------------------------------------------------------------------------
  subroutine  advec_1d     &
      ( delxa, delxb,      & !--- in
        gdu,               & !--- in
        dtime,             & !--- in
        gdq              )   !--- inout

  real(8), intent(in) :: delxa( KMAX ), delxb( KMAX )
  real(8), intent(in) :: gdu( KMAX )
  real(8), intent(in) :: dtime
  real(8), intent(inout) :: gdq( KMAX )
!  real(8), intent(in) :: delxa( KA ), delxb( KA )
!  real(8), intent(in) :: gdu( KA )
!  real(8), intent(in) :: dtime
!  real(8), intent(inout) :: gdq( KA )

  !--- local
  integer :: i
  real(8) :: fq( KS:KE+1 ), tmp(KMAX)
  real(8) :: dqr, dql, dq, qstar
  !
  !
  !--- reset of fluxes
  fq( : ) = 0.D0
  tmp(:) = gdq(:)

  !--- tracer flux
  do i = KS, KE+1
    if ( gdu( i ) > 0.D0 ) then
      dqr = ( gdq( i   )-gdq( i-1 ) )/delxb( i )
      dql = ( gdq( i-1 )-gdq( i-2 ) )/delxb( i-1 )
      if ( dqr*dql > 0.D0 ) then
        dq = 2.D0 / ( 1.D0/dqr+1.D0/dql )
      else
        dq = 0.D0
      end if
      qstar = gdq( i-1 )+( delxa( i-1 )-gdu( i )*dtime )*dq*0.5D0
    else
      dqr = ( gdq( i+1 )-gdq( i   ) )/delxb( i+1 )
      dql = ( gdq( i   )-gdq( i-1 ) )/delxb( i )
      if ( dqr*dql > 0.D0 ) then
        dq = 2.D0 / ( 1.D0/dqr+1.D0/dql )
      else
        dq = 0.D0
      end if
      qstar = gdq( i )-( delxa( i )+gdu( i )*dtime )*dq*0.5D0
    end if
    fq( i ) = qstar*gdu( i )
  end do

  !--- change of concentration by flux convergence
  do i = KS, KE
    gdq( i ) = gdq( i ) - dtime/delxa( i )*( fq( i+1 )-fq( i ) )
  end do
  !
  !
  return

  end subroutine  advec_1d
  !-----------------------------------------------------------------------------
  subroutine faero( f0,ga )

  real(8), intent(in) ::  f0
  real(8), intent(inout) :: ga( nccn )
  real(8) :: gaero( nccn ), f1, radmax, radmin
  real(8), parameter :: alpha = 3.D0
  integer :: n

!  radmin = ( exp( xactr( 1 ) )*3.D0/4.D0/pi/rhoa )**( OneovThird )
!  radmax = ( exp( xactr( nccn ) )*3.D0/4.D0/pi/rhoa )**( OneovThird )
!  f1 = 2.d0*f0/r0a/r0a/r0a* &
!        ( radmax*radmax*radmin*radmin/( radmax*radmax-radmin*radmin ) )
  do n = 1, nccn
!   gaero( n ) = f1*( rada( n )/r0a )**( -alpha )*exp( xactr( n ) )/dxaer
   gaero( n ) = f0*marate( n )*exp( xactr( n ) )/dxaer
   ga( n ) = ga( n )+gaero( n )
  end do

  return

  end subroutine faero
  !-------------------------------------------------------------------------------
  !
  ! + Y. Sato added for stochastic method 
  ! + Reference Sato et al. (2009) JGR, doi:10.1029/2008JD011247
  !
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
end module mod_atmos_phy_mp
!-------------------------------------------------------------------------------
