!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Utility subroutines dealing super-droplet IDs
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
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_idutil
  use scale_precision

  implicit none
  private
  public :: sdm_sort,sdm_getperm

contains
  subroutine sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,    &
                         sort_tag0,fsort_tag,fsort_id,            &
                         sd_rand,sd_perm)

    ! Input variables
    integer, intent(in) :: freq_max ! maximum number of SD in each sd-grid
    integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
    integer, intent(in) :: sd_num   ! Number of super-droplets
    integer, intent(in) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! = sort_tag(n) - 1
    ! sort_tag(m) : accumulated number of
    !               super-droplets in each
    !               SDM-grid
    integer, intent(in) :: fsort_tag(0:freq_max+1)
    ! accumulated number with respect to
    ! the number contaiend super-droplets
    integer, intent(in) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! super-droplets sorted by the number
    ! contaiend super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
    ! Output variables
    integer, intent(out) :: sd_perm(1:sd_num)    ! random permutations
    ! Work variables
    integer :: m, n, t, ss, tt         ! index
    !---------------------------------------------------------------

    ! Initialize for calculating random permutations
    
    do t=fsort_tag(2),fsort_tag(freq_max+1)-1

       m = fsort_id(t)

       tt = sort_tag0(m) + 1
       sd_perm(tt) = 1

    end do

    ! Calculate random permutations

    do n=2,freq_max

       do t=fsort_tag(n),fsort_tag(freq_max+1)-1

          m = fsort_id(t)

          tt = sort_tag0(m) + n
          !ORG        ss = sort_tag0(m) + ( int(sd_rand(tt)*n) + 1 )
          ss = sort_tag0(m) + ( int(sd_rand(tt)*real(n,kind=RP)) + 1 )

          !### swap data ###!

          sd_perm(tt) = sd_perm(ss)
          sd_perm(ss) = n

       end do

    end do

    return
  end subroutine sdm_getperm
  !----------------------------------------------------------------------------
  subroutine sdm_sort(ni_sdm,nj_sdm,nk_sdm,                       &
                      sd_num,sd_n,sd_ri,sd_rj,sd_rk,                &
                      sort_id,sort_key,sort_freq,sort_tag,jdgtype)
    use scale_grid_index, only: &
         IS,JS,KS
    use m_sdm_common, only: &
         VALID2INVALID,knum_sdm
    use gadg_algorithm, only: &
         gadg_count_sort
    ! Input variables
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    real(RP), intent(in) :: sd_ri(1:sd_num)  ! index[i/real] of super-droplets
    real(RP), intent(in) :: sd_rj(1:sd_num)  ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)  ! index[k/real] of super-droplets
    character(len=5), intent(in) :: jdgtype   ! flag for sorting
    ! Output variables
    integer, intent(out) :: sort_id(1:sd_num)   ! id that super-droplets sorted by sd-grids
    integer, intent(out) :: sort_key(1:sd_num)  ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each sd-grid
    integer, intent(out) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)     ! accumulated number of super-droplets in each sd-grid
    integer :: adr                ! grid id
    integer :: max_key            ! total grid number
    integer :: i, j, k, n         ! index
    integer :: ix, jy
    !--------------------------------------------------------------------

    ! Initialize

    max_key = ni_sdm * nj_sdm * knum_sdm + 1

    ! Sorting [step-1] -- initialize IKEY
    if( jdgtype .eq. 'valid' ) then

       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID ) then

             ! i \in {0,...,ni_sdm-1}, j \in {0,...,nj_sdm-1}, k \in {0,...,knum_sdm-1}
             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do

    else if( jdgtype .eq. 'multi' ) then

       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID .and. sd_n(n)>1 ) then

             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do

    end if

    
    ! Sorting [step-2] -- counting sort
    call gadg_count_sort( sort_key, 1, max_key,                       &
         sort_freq, sort_tag, sort_id )

    sort_tag(max_key+1) = sort_tag(max_key) + sort_freq(max_key)

    return
  end subroutine sdm_sort
end module m_sdm_idutil
