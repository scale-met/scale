!-------------------------------------------------------------------------------
!> module atmosphere / physics / chemistry / RN222
!!
!! @par Description
!!          General component for rn222 tracer
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2017-02-28 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_ch_rn222
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CH_rn222_setup
  public :: ATMOS_PHY_CH_rn222_tendency

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: QA_CH = 1

  character(len=H_SHORT), public, target :: ATMOS_PHY_CH_rn222_NAME(QA_CH)
  character(len=H_MID)  , public, target :: ATMOS_PHY_CH_rn222_DESC(QA_CH)
  character(len=H_SHORT), public, target :: ATMOS_PHY_CH_rn222_UNIT(QA_CH)

  data ATMOS_PHY_CH_rn222_NAME / &
                 'RN222'  /

  data ATMOS_PHY_CH_rn222_DESC / &
                 'Mixing Ratio of Rn222'   /

  data ATMOS_PHY_CH_rn222_UNIT / &
                 'Bq/kg'  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  public, parameter :: I_ch_rn222 = 1
  integer,  public            :: QS_CH
  integer,  public            :: QE_CH

  real(RP),          private :: Rn222_decay_ratio                       ! Decay constant [/s]

  character(len=64), private :: Rn222_emission_type        = 'CONST'    ! Emission type
  real(RP),          private :: Rn222_const_emission_land  = 20.8E-3_RP ! Surface flux from land  [Bq/m2/s]
  real(RP),          private :: Rn222_const_emission_ocean = 0.14E-3_RP ! Surface flux from ocean [Bq/m2/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_rn222_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    real(RP) :: Rn222_half_life = 3.30048E+5_RP

    namelist / PARAM_ATMOS_PHY_CH_RN222 / &
       Rn222_half_life,           &
       Rn222_emission_type,       &
       Rn222_const_emission_land, &
       Rn222_const_emission_ocean

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Chemistry] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** rn222 process'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CH_RN222,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CH_RN222. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_CH_RN222)

    Rn222_decay_ratio = log(2.0_RP) / Rn222_half_life

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(A)')           ' *** Characteristics of Rn222'
    if( IO_L ) write(IO_FID_LOG,'(A,E16.6)')     ' *** Half life   [s]      : ', rn222_half_life
    if( IO_L ) write(IO_FID_LOG,'(A,E16.6)')     ' *** Decay ratio [1/s]    : ', Rn222_decay_ratio
    if( IO_L ) write(IO_FID_LOG,*)               ' *** Type of emission     : ', trim(RN222_emission_type)
    if ( Rn222_emission_type == 'CONST' ) then
       if( IO_L ) write(IO_FID_LOG,'(A,ES16.6)') ' *** From land  [Bq/m2/s] : ', Rn222_const_emission_land
       if( IO_L ) write(IO_FID_LOG,'(A,ES16.6)') ' *** From ocean [Bq/m2/s] : ', Rn222_const_emission_ocean
    else
       write(*,*) 'xxx Not supported type of Rn222 emission! Stop.'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_CH_rn222_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_CH_rn222_tendency( &
       KA, KS, KE,                 &
       IA, IS, IE,                 &
       JA, JS, JE,                 &
       QQA,                        &
       DENS,                       &
       QTRC,                       &
       ATMOS_GRID_CARTESC_REAL_FZ, &
       LANDUSE_fact_land,          &
       RHOQ_t                      )
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_tracer, only: &
       QA
    implicit none

    integer,  intent(in)    :: KA
    integer,  intent(in)    :: KS
    integer,  intent(in)    :: KE
    integer,  intent(in)    :: IA
    integer,  intent(in)    :: IS
    integer,  intent(in)    :: IE
    integer,  intent(in)    :: JA
    integer,  intent(in)    :: JS
    integer,  intent(in)    :: JE
    integer,  intent(in)    :: QQA
    real(RP), intent(in)    :: DENS  (KA,IA,JA)
    real(RP), intent(in)    :: QTRC  (KA,IA,JA,1)
    real(RP), intent(in)    :: ATMOS_GRID_CARTESC_REAL_FZ (0:KA,IA,JA)
    real(RP), intent(in)    :: LANDUSE_fact_land (IA,JA)
    real(RP), intent(inout) :: RHOQ_t(KA,IA,JA,QQA)

    real(RP) :: emission(IA,JA)

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Chemistry(Rn222)'

    !iq = QS_CH - 1 + I_ch_rn222

    !--- Decay based on half life

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,I_ch_rn222) = RHOQ_t(k,i,j,I_ch_rn222) &
                                - DENS(k,i,j) * QTRC(k,i,j,1) * Rn222_decay_ratio ! [Bq/m3/s]
    enddo
    enddo
    enddo

    !--- Surface emission

    if ( RN222_emission_type == "CONST" ) then

       do j  = JS, JE
       do i  = IS, IE
          emission(i,j) = ( 1.0_RP-LANDUSE_fact_land(i,j) ) * Rn222_const_emission_ocean &
                        + (        LANDUSE_fact_land(i,j) ) * Rn222_const_emission_land
       enddo
       enddo

    endif

    do j  = JS, JE
    do i  = IS, IE
       RHOQ_t(KS,i,j,I_ch_rn222) = RHOQ_t(KS,i,j,I_ch_rn222) + emission(i,j) / ( ATMOS_GRID_CARTESC_REAL_FZ(KS,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(KS-1,i,j) )
    enddo
    enddo

    call FILE_HISTORY_in( emission(:,:), 'EMIT_RN222', 'Emission Rn222', 'Bq/m2/s' )

    return
  end subroutine ATMOS_PHY_CH_rn222_tendency

end module scale_atmos_phy_ch_rn222
