!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coalescence subroutines for the SDM
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-12 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_coalescence
  use scale_precision

  implicit none
  private
  public :: sdm_coales

contains
  subroutine sdm_coales(sdm_colkrnl,sdm_colbrwn,               &
                        sdm_aslset,sdm_aslrho,                 &
                        sdm_dtcol,            &
                        pres_scale, t_scale,                   &
                        zph_crs,                &
                        ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl, &
                        sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,&
                        sort_id,sort_key,sort_freq,sort_tag,   &
                        sd_rng,sd_rand,                        &
                        sort_tag0,fsort_id,icp,sd_perm,c_rate  )
    use gadg_algorithm, only: &
         gadg_count_sort
    use rng_uniform_mt, only: &
         c_rng_uniform_mt, gen_rand_array
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_const, only:  &
         t0 => CONST_TEM00, &    ! Gas Constant of vapor [J/K/kg]
         rw => CONST_Rvap,  &    ! Gas Constant of vapor [J/K/kg]
         rd => CONST_Rdry,  &    ! Gas Constant of dry air [J/K/kg]
         cp => CONST_CPdry, &
         p0 => CONST_PRE00       ! Reference Pressure [Pa]
    use m_sdm_common, only: &
         VALID2INVALID,INVALID,knum_sdm, &
         rho_amsul,rho_nacl,ONE_PI,m2micro,r0col,ratcol,ecoll,micro2m,dxiv_sdm,dyiv_sdm,F_THRD,O_THRD,rrst,boltz,mass_air
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort, sdm_getperm
    !  Input variables
    integer, intent(in) :: sdm_colkrnl   ! Kernel type for coalescence process
    integer, intent(in) :: sdm_colbrwn   ! Control flag of Brownian Coagulation and Scavenging process
    integer, intent(in) :: sdm_aslset    ! Control flag to set species and way of chemical material as water-soluble aerosol
    real(RP),intent(in) :: sdm_aslrho(20)! User specified density of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sdm_dtcol   ! tims step of {stochastic coalescence} process
    real(RP), intent(in) :: pres_scale(KA,IA,JA)  ! Pressure
    real(RP), intent(in) :: t_scale(KA,IA,JA)    ! Temperature
    real(RP), intent(in) :: zph_crs(KA,IA,JA) ! z physical coordinate
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer, intent(in) :: sd_numasl ! Number of kind of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    ! Input and output variables
    type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
    real(RP),intent(inout) :: sd_rand(1:sd_num) ! random numbers
    integer, intent(inout) :: sort_id(1:sd_num) ! super-droplets sorted by SD-grids
    integer, intent(inout) :: sort_key(1:sd_num) ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
    integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2) ! accumulated number of super-droplets in each SD-grid
    integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Output variables
    integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
    integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    integer, intent(out) :: icp(1:sd_num) ! index of coalescence pair
    integer, intent(out) :: sd_perm(1:sd_num) ! random permutations
    real(RP), intent(out) :: c_rate(1:sd_num) ! coalescence probability
    ! Internal shared variables
    real(RP) :: sd_aslrho(1:22) ! Density of chemical material contained as water-soluble aerosol in super droplets
    integer(RP) :: sd_ncol ! how many times coalescence occurs
    integer :: freq_max ! get the maximum number of super-droplets in each grid
    integer :: hfreq_max ! hfreq_max / 2
    integer :: ipremium  ! premium coef. for coalescence
    ! Work variables
    real(RP) :: dmask(1:22) ! mask for vactorization
    real(RP) :: sd_asl1(1:sd_numasl)
    real(RP) :: sd_asl2(1:sd_numasl) ! aerosol mass of super-droplets with large/small multiplicity
    real(RP) :: sd_cc1  ! slip correction of super-droplets
    real(RP) :: sd_cc2  ! with large/small multiplicity
    real(RP) :: sd_dia1 ! diameter of super-droplets
    real(RP) :: sd_dia2 ! with large/small multiplicity
    real(RP) :: sd_lmd1 ! mean free path of super-droplets
    real(RP) :: sd_lmd2 ! with large/small multiplicity
    real(RP) :: sd_m1   ! mass of super-droplets
    real(RP) :: sd_m2   ! with large/small multiplicity
    real(RP) :: sd_r1   ! radius of super-droplets
    real(RP) :: sd_r2   ! with large/small multiplicity
    real(RP) :: sd_rk1  ! index[k/real] of super-droplets with large multiplicity
    real(RP) :: sd_rw1  ! radius of water parts in super-droplets
    real(RP) :: sd_rw2  ! with large/small multiplicity
    real(RP) :: sd_tmasl1
    real(RP) :: sd_tmasl2 ! total mass of aerosol part in super droplets with large/small multiplicity
    real(RP) :: sd_tvasl1 ! total volume of aerosol part in super
    real(RP) :: sd_tvasl2 ! droplets with large/small multiplicity
    real(RP) :: sd_v1   ! volume of super-droplets
    real(RP) :: sd_v2   ! with large/small multiplicity
    real(RP) :: sd_c1   ! temporary
    real(RP) :: sd_c2
    real(RP) :: sd_d1
    real(RP) :: sd_d2
    real(RP) :: sd_g1
    real(RP) :: sd_g2
    real(RP) :: lmd_crs ! air mean free path
    real(RP) :: p_crs   ! pressure
    real(RP) :: pt_crs  ! potential temperarure
    real(RP) :: t_crs   ! temperarure
    real(RP) :: vis_crs ! dynamic viscosity
    real(RP) :: sumdia  ! sum of variables of a pair of droplets
    real(RP) :: sumd
    real(RP) :: sumc
    real(RP) :: sumg
    real(RP) :: sumr
    real(RP) :: k12     ! brownian coagulation coefficient
    real(RP) :: dvz     ! difference in terminal velocity of a pair of super-droplets
    real(RP) :: dtmp    ! temporary
    real(RP) :: frac    ! fraction parts
    real(RP) :: ivvol(KA,IA,JA)   ! inverse of a grid volume
    real(RP) :: tdeg    ! temperature in degree
    real(RP) :: rq      ! radius ratio of a pair of super-droplets
    real(RP) :: ek      ! temporary
    real(RP) :: p       ! temporary
    real(RP) :: q       ! temporary

    integer(DP) :: sd_nmax  ! maximum multiplicity
    integer(DP) :: sd_n1    ! multiplicity of super-droplets with large multiplicity
    integer(DP) :: sd_n2    ! multiplicity of super-droplets  with small multiplicity

    integer, allocatable :: fsort_tag(:) ! buffer for sorting
    integer, allocatable :: fsort_freq(:) ! buffer for sorting

    integer :: idx_nasl(1:22)  ! index for vactorization


    integer :: gnum          ! grid number

    integer :: in1, in2, in3, in4 ! index
    integer :: irr, iqq           ! index of coef.
    integer :: i, j, k, m, n, s   ! index
    integer :: t, tc, tp          ! index
    integer :: ix, jy

    integer :: sort_tag0m
    integer :: sort_freqm
    integer :: icptc, icptp
    !--------------------------------------------------------------------

    ! Section specification for fapp profiler
    call fapp_start("sdm_coales",1,1)

    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    ! ni_sdm=IE-IS+1, nj_sdm=JE-JS+1, knum_sdm=floor(rkumax)+1)-KS+1
    gnum = ni_sdm * nj_sdm * knum_sdm

    freq_max = 1

    !### aerosol type ###!
    ! why we need sd_aslrho? Isn't it same to sdm_aslrho?? Check later
    do n=1,22
       sd_aslrho(n) = 1.0_RP
    end do

    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then

          !### numasl=1 @ init : NaCl ###!

          sd_aslrho(1) = real(rho_nacl,kind=RP)

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslrho(1) = real(rho_amsul,kind=RP)
          sd_aslrho(2) = real(rho_nacl,kind=RP)

       end if
    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)
       sd_aslrho(2) = real(rho_nacl,kind=RP)

       !         do n=1,20
       !            call getrname( id_sdm_aslrho + (n-1), sd_aslrho(n+2) )
       !         end do
       do n=1,20
          sd_aslrho(n+2) = sdm_aslrho(n)
       end do

    end if

    ! Sorting super-droplets.
    
    call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
         sort_id,sort_key,sort_freq,sort_tag,'valid')

    ! Initialize
    
    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do

    do n=1,gnum+2
       sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
    end do

    do n=1,sd_num
       c_rate(n) = 0.0_RP
    end do

    ! Get the maximum number of super-droplets in each grid
    
    do m=1,gnum
       freq_max = max( freq_max, sort_freq(m) )
    end do

    hfreq_max = int(freq_max/2)

    ! Sorting the grids by the number of super-droplets in each grid

    allocate( fsort_tag(0:freq_max+1) )
    allocate( fsort_freq(0:freq_max)  )

    call gadg_count_sort( sort_freq(1:gnum), 0, freq_max, &
         fsort_freq, fsort_tag, fsort_id )

    fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

    ! Get random number using random number generator
    
    call gen_rand_array( sd_rng, sd_rand )

    ! Get random permutation layout of super-droplets in each grid

    call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
         &                 sort_tag0,fsort_tag,fsort_id,sd_rand,sd_perm)

    ! Get random number using random number generator

    call gen_rand_array( sd_rng, sd_rand )

    ! Select a pair of super-droples in random permutation layout

!!$      do n=1,hfreq_max
!!$
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle
       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          in1 = 2 * n
          in2 = in1 - 1

          !### select pair of super-droplets in each grid ###!

!!$            in3 = sd_perm( sort_tag0(m) + in1 )
!!$            in4 = sd_perm( sort_tag0(m) + in2 )
          in3 = sd_perm( sort_tag0m + in1 )
          in4 = sd_perm( sort_tag0m + in2 )

          !### set the random index ###!
!!$            tc = sort_tag0(m) + n
          tc = sort_tag0m + n
!!$            tp = tc + int(sort_freq(m)/2)
          tp = tc + int(sort_freqm/2)

!!$            icp(tc) = sort_id( sort_tag0(m) + in3 )
!!$            icp(tp) = sort_id( sort_tag0(m) + in4 )
          icp(tc) = sort_id( sort_tag0m + in3 )
          icp(tp) = sort_id( sort_tag0m + in4 )

       end do

    end do

    ! Get effective collision probability for "Gravitational Settling" 
    
    if( sdm_colkrnl==0 ) then

       !### Golovin's kernel (m^3/s?) ###!
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = sd_r(icptc)
             sd_r2 = sd_r(icptp)

             c_rate(tc) = 1500.0_RP * 1.333333_RP * ONE_PI              &
                  * ( sd_r1*sd_r1*sd_r1 + sd_r2*sd_r2*sd_r2 )

          end do

       end do

    else if( sdm_colkrnl==1 ) then

       !### Long's kernel (m^2)  ###!

!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = max( sd_r(icptc), sd_r(icptp) )  !! large
             sd_r2 = min( sd_r(icptc), sd_r(icptp) )  !! small

             if( sd_r1 <= 5.E-5_RP ) then
                c_rate(tc) = 4.5E+8_RP * ( sd_r1*sd_r1 )                 &
                     * ( 1.0_RP - 3.E-6_RP/(max(3.01E-6_RP,sd_r1)) )
             else
                c_rate(tc) = 1.d0
             end if

             sumr = sd_r1 + sd_r2

             c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

          end do

       end do

    else if( sdm_colkrnl==2 ) then

       !### Hall's kernel (m^2) ###!

!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + int(sort_freq(m)/2)

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sd_r1 = max( sd_r(icptc), sd_r(icptp) )  !! large
             sd_r2 = min( sd_r(icptc), sd_r(icptp) )  !! small

             rq    = sd_r2 / sd_r1
             sd_r1 = sd_r1 * m2micro    !! [m] => [micro-m]
             sd_r2 = sd_r2 * m2micro    !! [m] => [micro-m]

             !! Get index of the array {r0,rat}.

             if( sd_r1 <= r0col(1) ) then
                irr = 1
             else if( sd_r1 <= r0col(2) ) then
                irr = 2
             else if( sd_r1 <= r0col(3) ) then
                irr = 3
             else if( sd_r1 <= r0col(4) ) then
                irr = 4
             else if( sd_r1 <= r0col(5) ) then
                irr = 5
             else if( sd_r1 <= r0col(6) ) then
                irr = 6
             else if( sd_r1 <= r0col(7) ) then
                irr = 7
             else if( sd_r1 <= r0col(8) ) then
                irr = 8
             else if( sd_r1 <= r0col(9) ) then
                irr = 9
             else if( sd_r1 <= r0col(10) ) then
                irr = 10
             else if( sd_r1 <= r0col(11) ) then
                irr = 11
             else if( sd_r1 <= r0col(12) ) then
                irr = 12
             else if( sd_r1 <= r0col(13) ) then
                irr = 13
             else if( sd_r1 <= r0col(14) ) then
                irr = 14
             else if( sd_r1 <= r0col(15) ) then
                irr = 15
             else
                irr = 16
             end if

             if( rq <= ratcol(2) ) then
                iqq = 2
             else if( rq <= ratcol(3) ) then
                iqq = 3
             else if( rq <= ratcol(4) ) then
                iqq = 4
             else if( rq <= ratcol(5) ) then
                iqq = 5
             else if( rq <= ratcol(6) ) then
                iqq = 6
             else if( rq <= ratcol(7) ) then
                iqq = 7
             else if( rq <= ratcol(8) ) then
                iqq = 8
             else if( rq <= ratcol(9) ) then
                iqq = 9
             else if( rq <= ratcol(10) ) then
                iqq = 10
             else if( rq <= ratcol(11) ) then
                iqq = 11
             else if( rq <= ratcol(12) ) then
                iqq = 12
             else if( rq <= ratcol(13) ) then
                iqq = 13
             else if( rq <= ratcol(14) ) then
                iqq = 14
             else if( rq <= ratcol(15) ) then
                iqq = 15
             else if( rq <= ratcol(16) ) then
                iqq = 16
             else if( rq <= ratcol(17) ) then
                iqq = 17
             else if( rq <= ratcol(18) ) then
                iqq = 18
             else if( rq <= ratcol(19) ) then
                iqq = 19
             else if( rq <= ratcol(20) ) then
                iqq = 20
             else
                iqq = 21
             end if
             !! Get c_rate

             if( irr>=16 ) then

                q  = (rq-ratcol(iqq-1)) / (ratcol(iqq)-ratcol(iqq-1))
                ek = (1.0_RP-q)*ecoll(15,iqq-1) + q*ecoll(15,iqq)

                c_rate(tc) = min( ek, 1.0_RP )

             else if( irr>=2 .and. irr<16 ) then

                p = (sd_r1-r0col(irr-1))/(r0col(irr)-r0col(irr-1))
                q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                c_rate(tc) = (1.0_RP-p)*(1.0_RP-q)*ecoll(irr-1,iqq-1)     &
                     + p*(1.0_RP-q)*ecoll(irr,iqq-1)              &
                     + q*(1.0_RP-p)*ecoll(irr-1,iqq)              &
                     + p*q*ecoll(irr,iqq)

             else

                q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                c_rate(tc) = (1.0_RP-q)*ecoll(1,iqq-1) + q*ecoll(1,iqq)

             end if

             sd_r1 = sd_r1 * micro2m    !! [micro-m] => [m]
             sd_r2 = sd_r2 * micro2m    !! [micro-m] => [m]

             !! Get c_rate

             sumr = sd_r1 + sd_r2

             c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

          end do

       end do

    else if( sdm_colkrnl==3 ) then

       !### no coalescence effeciency hydrodynamic kernel (m^2) ###!
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2

!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             sumr = sd_r(icptc) + sd_r(icptp)

             c_rate(tc) = ONE_PI * (sumr*sumr)

          end do

       end do

    end if

    ! -----

    if( sdm_colkrnl==0 ) then

       !### Golovin's kernel [-] ###!
!!$         do n=1,hfreq_max
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2
             
             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             !! Get location of Super-Droplets
             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do
       end do

    else

       !### Long's kernel, Hall's kernel,       ###!
       !### no col_effi hydrodynamic kernel [-] ###!
       
!!$         do n=1,hfreq_max
!!$
!!$            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$               m = fsort_id(t)
!!$
!!$               tc = sort_tag0(m) + n
!!$               tp = tc + sort_freq(m)/2
       
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)

          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             dvz = abs( sd_vz(icptc) - sd_vz(icptp) )

             !! Get location of Super-Droplets

             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP)        &
                  * ivvol(k,i,j) * dvz

          end do

       end do

    end if

    ! Get effective collision for "Brownian Coagulation and Scavenging
    ! (Seinfeld & Pandis,2006)" 
    ! This is mechanisim is effective for droplets less than micrometer-size ( below 1um )
    if( sdm_colbrwn>0 ) then
!!$        do n=1,hfreq_max
!!$
!!$          do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)
          
          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             !### information of super-droplets ###

             !! radius of water parts in droplets

             sd_rw1 = sd_r(icptc)    !! [m]
             sd_rw2 = sd_r(icptp)

             !! mass and volume of aerosol parts in droplets

             sd_tmasl1 = 0.0_RP
             sd_tmasl2 = 0.0_RP
             sd_tvasl1 = 0.0_RP
             sd_tvasl2 = 0.0_RP

             do k=1,22

                s = idx_nasl(k)

                sd_tmasl1 = sd_tmasl1 + sd_asl(icptc,s) * dmask(k)
                sd_tmasl2 = sd_tmasl2 + sd_asl(icptp,s) * dmask(k)

                sd_tvasl1 = sd_tvasl1                                    &
                     + sd_asl(icptc,s)/sd_aslrho(s) * dmask(k)
                sd_tvasl2 = sd_tvasl2                                    &
                     + sd_asl(icptp,s)/sd_aslrho(s) * dmask(k)

             end do

             sd_tmasl1 = sd_tmasl1  * 1.E-3_RP    !! [g]=>[kg]
             sd_tmasl2 = sd_tmasl2  * 1.E-3_RP

             !! diameter and mass and volume of droplets

             dtmp = ONE_PI * F_THRD

             sd_v1 = sd_tvasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1)
             sd_v2 = sd_tvasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2)

             sd_m1 = sd_tmasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1) * rw
             sd_m2 = sd_tmasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2) * rw

             sd_dia1 = (6.0_RP*sd_v1/ONE_PI)**O_THRD
             sd_dia2 = (6.0_RP*sd_v2/ONE_PI)**O_THRD

             !### location of super-droplets ###

             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             p_crs = pres_scale(k,i,j) !! [Pa]
             t_crs = t_scale(k,i,j)    !! [K]

             !### dynamic viscosity [Pa*s]  ###!
             !### (Pruppacher & Klett,1997) ###!
             tdeg = t_crs - t0     !! [K] => [degC]
             if( tdeg>=0.0_RP ) then
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg                        &
                     -1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             !### air mean free path [m] ###!
             dtmp = dsqrt(8.0_RP*mass_air*1.E-3_RP/(ONE_PI*rrst*t_crs))
             lmd_crs = (2.0_RP*vis_crs)/(p_crs*dtmp)

             !### slip correction of droplets [-]  ###!
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia1/lmd_crs)
             sd_cc1 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia1)
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia2/lmd_crs)
             sd_cc2 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia2)

             !### diffusion term [m*m/s] ###!
             dtmp = (boltz*t_crs)/(3.0_RP*ONE_PI*vis_crs)
             sd_d1 = dtmp * (sd_cc1/sd_dia1)
             sd_d2 = dtmp * (sd_cc2/sd_dia2)

             !### velocity term [m/s] ###!
             dtmp = (8.0_RP*boltz*t_crs)/ONE_PI
             sd_c1 = dsqrt(dtmp/sd_m1)
             sd_c2 = dsqrt(dtmp/sd_m2)

             !### mean free path of droplets [m] ###!
             dtmp = 8.0_RP/ONE_PI
             sd_lmd1 = dtmp * (sd_d1/sd_c1)
             sd_lmd2 = dtmp * (sd_d2/sd_c2)

             !### length term [m] ###!
             dtmp = (sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)&
                  - (sd_dia1*sd_dia1+sd_lmd1*sd_lmd1)**1.50_RP
             sd_g1 = dtmp/(3.0_RP*sd_dia1*sd_lmd1) - sd_dia1
             dtmp = (sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)&
                  - (sd_dia2*sd_dia2+sd_lmd2*sd_lmd2)**1.50_RP
             sd_g2 = dtmp/(3.0_RP*sd_dia2*sd_lmd2) - sd_dia2

             !### Brownian Coagulation Coefficient K12 [m3/s] ###!
             sumdia = sd_dia1 + sd_dia2
             sumd   = sd_d1   + sd_d2
             sumc   = dsqrt( sd_c1*sd_c1 + sd_c2*sd_c2 )
             sumg   = dsqrt( sd_g1*sd_g1 + sd_g2*sd_g2 )

             dtmp = sumdia/(sumdia+2.0_RP*sumg) + (8.0_RP*sumd)/(sumdia*sumc)
             k12 = 2.0_RP*ONE_PI * sumdia*sumd/dtmp

             !### add effective collision [-] ###!
             c_rate(tc) = c_rate(tc)                                     &
                  + k12 * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do

       end do

    end if

    ! Get total effective collision of droplets
    
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          sd_nmax  = max( sd_n(icptc), sd_n(icptp) )
          ! maximum multiplicity
          ipremium = sort_freqm - 1 + iand(sort_freqm,1)
          ! IAND(sort_freq(i),1) => even:0, odd:1

          c_rate(tc) = c_rate(tc) * real( sd_nmax*ipremium, kind=RP )
       end do
    end do

    ! Stochastic coalescence process.
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          !### set coalescence count ###!
          sd_ncol = int( c_rate(tc), kind=RP )
          frac = c_rate(tc) - real( sd_ncol, kind=RP )

          ! judge coalescence by random number and fractional part

          if( sd_rand(tc) < frac ) then
             sd_ncol =  sd_ncol + 1
          end if

          if( sd_ncol<=0 ) cycle  !! no coalesecense

          !### coalescence procudure ###!

          if( sd_n(icptc) > sd_n(icptp) ) then

             sd_n1  = sd_n( icptc )
             sd_r1  = sd_r( icptc )
             sd_rk1 = sd_rk( icptc )
             sd_m1  = sd_r1 * sd_r1 * sd_r1

             sd_n2  = sd_n( icptp )
             sd_r2  = sd_r( icptp )
             sd_m2  = sd_r2 * sd_r2 * sd_r2

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptc,s )
                sd_asl2(s) = sd_asl( icptp,s )
             end do

          else

             sd_n1  = sd_n( icptp )
             sd_r1  = sd_r( icptp )
             sd_rk1 = sd_rk( icptp )
             sd_m1  = sd_r1 * sd_r1 * sd_r1

             sd_n2  = sd_n( icptc )
             sd_r2  = sd_r( icptc )
             sd_m2  = sd_r2 * sd_r2 * sd_r2
             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptp,s )
                sd_asl2(s) = sd_asl( icptc,s )
             end do

          end if

          sd_ncol = min( sd_ncol, int(sd_n1/sd_n2,kind=RP) )

          if( sd_n1 > sd_n2*sd_ncol ) then

             sd_n1 = sd_n1 - sd_n2*sd_ncol
             sd_r2 = exp( O_THRD                                      &
                  * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
             !ORG           sd_r2 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD

             do k=1,22
                s = idx_nasl(k)
                dtmp = sd_asl1(s) * real(sd_ncol,kind=RP)
                sd_asl2(s) = sd_asl2(s) + dmask(k) * dtmp
             end do

          else

             !! coalescent SDs with same multiplicity
             !!  - do not change the order of equations for
             !!  - vectorization

             sd_n1 = int( sd_n2/2, kind=RP )
             sd_n2 = sd_n2 - sd_n1

             sd_r1 = exp( O_THRD                                      &
                  * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
             !ORG           sd_r1 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD
             sd_r2 = sd_r1
             do k=1,22

                s = idx_nasl(k)
                dtmp = sd_asl1(s)*real(sd_ncol,kind=RP) + sd_asl2(s)
                sd_asl1(s) = sd_asl1(s) + dmask(k)*(dtmp-sd_asl1(s))
                sd_asl2(s) = sd_asl1(s)

             end do

             !! invalid by collisions between SDs with
             !! sd_n1=sd_n2*sd_ncol and sd_n2=1

             if( sd_n1==0 ) then
                sd_rk1 = INVALID
             end if

          end if

          !! This never happens
!!$            !! check muliplicity
!!$
!!$            if( sd_n1>(2.0_RP**63._RP) .or. sd_n2>(2.0_RP**63._RP) ) then
!!$               iexced = -1
!!$               cycle
!!$            end if

          if( sd_n(icptc) > sd_n(icptp) ) then
             sd_n( icptc )  = sd_n1
             sd_r( icptc )  = sd_r1
             sd_rk( icptc ) = sd_rk1

             sd_n( icptp )  = sd_n2
             sd_r( icptp )  = sd_r2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptc,s ) = sd_asl1(s)
                sd_asl( icptp,s ) = sd_asl2(s)
             end do

          else

             sd_n( icptp )  = sd_n1
             sd_r( icptp )  = sd_r1
             sd_rk( icptp ) = sd_rk1

             sd_n( icptc )  = sd_n2
             sd_r( icptc )  = sd_r2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptp,s ) = sd_asl1(s)
                sd_asl( icptc,s ) = sd_asl2(s)
             end do

          end if

       end do

    end do

    ! Deallocate
    deallocate( fsort_tag  )
    deallocate( fsort_freq )

    ! Section specification for fapp profiler
    call fapp_stop("sdm_coales",1,1)

    return
  end subroutine sdm_coales
end module m_sdm_coalescence
