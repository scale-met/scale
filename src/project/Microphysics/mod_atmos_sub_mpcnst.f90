module mod_atmos_cnst
  implicit none
  public
  !--------------------------------------------------------------------------------------
  ! physical constant parameters imported from NICAM
  !--------------------------------------------------------------------------------------
  real(8), parameter :: CNST_RAIR  = 287.04D0    ! gas const. of dry air
  real(8), parameter :: CNST_RVAP  = 461.50D0    ! gas const. of vapor
  real(8), parameter :: CNST_EPSV  = CNST_RAIR / CNST_RVAP
  real(8), parameter :: CNST_EPSVT = 1.0D0/CNST_EPSV - 1.0D0
  real(8), parameter :: CNST_CP    = 1004.6D0    ! Specific heat of air   (consant pressure)
  real(8), parameter :: CNST_CPV   = 1850.0D0    ! specific heat of vapor (consant pressure)
  real(8), parameter :: CNST_CV    = CNST_CP  - CNST_RAIR ! specific heat of air   (consant volume)
  real(8), parameter :: CNST_CVV   = CNST_CPV - CNST_RVAP ! specific heat of vapor (consant volume)
  real(8), parameter :: CNST_CL    = 4218.0D0    ! specific heat of water 
  real(8), parameter :: CNST_CI    = 2006.0D0
  real(8), parameter :: CNST_KAPPA = CNST_RAIR / CNST_CP
  real(8), parameter :: CNST_TEM00 = 273.15d0    ! 0 degree
  real(8), parameter :: CNST_EMELT = 3.40D+5
  real(8), parameter :: CNST_TMELT = 273.15D0
  !
  real(8), parameter :: CNST_PRES0 = 101325.0D0  ! surface pressure
  real(8), parameter :: CNST_PRE00 = 100000.0D0  ! standard pressure
  real(8), parameter :: CNST_LH0   = 2.5008D6    ! latent heat of vaporizaion at 0 degre  
  !                                              ! latent heat of vaporizaion at 0 K
  real(8), parameter :: CNST_LH00  = CNST_LH0  - ( CNST_CPV - CNST_CL ) * CNST_TEM00
  real(8), parameter :: CNST_LHS0  = 2.8342D6
  real(8), parameter :: CNST_LHS00 = CNST_LHS0 - ( CNST_CPV - CNST_CI ) * CNST_TEM00
  real(8), parameter :: CNST_LHF0  = CNST_LHS0 - CNST_LH0
  real(8), parameter :: CNST_LHF00 = CNST_LHF0 - ( CNST_CL  - CNST_CI ) * CNST_TEM00
  real(8), parameter :: CNST_PSAT0 = 610.7D0     ! saturate pressure of water vapor at 0C
  real(8), parameter :: CNST_DWATR = 1000.d0     ! density of water
  real(8), parameter :: CNST_DICE  =  916.8D0
  real(8), parameter :: CNST_EGRAV = 9.80616D0
  real(8), parameter :: CNST_PI    = 3.14159265358979323846264d0  
  real(8), parameter :: CNST_ERADIUS  = 6.37122D+6 
  real(8), parameter :: CNST_UNDEF    = -99.9E+33
  real(4), parameter :: CNST_UNDEF4   = -99.9E+33
  real(8), parameter :: CNST_VMISS    = 0.0D0
  !--------------------------------------------------------------------------------------
  ! experimental settings
  !--------------------------------------------------------------------------------------
  integer, parameter :: I_QV       = 1
  integer, parameter :: I_QC       = 2
  integer, parameter :: I_QR       = 3
  integer, parameter :: I_QI       = 4
  integer, parameter :: I_QS       = 5
  integer, parameter :: I_QG       = 6
  integer, parameter :: I_NC       = 7
  integer, parameter :: I_NR       = 8
  integer, parameter :: I_NI       = 9
  integer, parameter :: I_NS       = 10
  integer, parameter :: I_NG       = 11
  integer, parameter :: TRC_VMAX   = I_NG ! total number of tracers
  integer, parameter :: HYDRO_MAX  = 6    ! total number of mixing ratio of water 
  integer, parameter :: NQW_STR    = I_QV ! 
  integer, parameter :: NQW_END    = I_QG ! 
  integer, parameter :: NNW_STR    = I_NC ! 
  integer, parameter :: NNW_END    = I_NG ! 
  integer, parameter :: NQS_STR    = I_QI ! solid start
  integer, parameter :: NQS_END    = I_QG ! solid end
  integer, parameter :: NQL_STR    = I_QC ! liquid start
  integer, parameter :: NQL_END    = I_QR ! liquid end
  character(len=32)  :: WLABEL(TRC_VMAX)=(/&
       "VAPOR", "CLOUD", "RAIN", "ICE", "SNOW", "GRAUPEL", &
       "CLOUD_NUM", "RAIN_NUM", "ICE_NUM", "SNOW_NUM", "GRAUPEL_NUM" /)
  !-------------------------------------------------------------------------------
  ! specific constants for each hydrometeor
  !-------------------------------------------------------------------------------
  !                                      I_QV     I_QC    I_QR    I_QI    I_QS    I_QG
  real(8), parameter :: CVW(HYDRO_MAX)=(/CNST_CVV,CNST_CL,CNST_CL,CNST_CI,CNST_CI,CNST_CI/)
  real(8), parameter :: CPW(HYDRO_MAX)=(/CNST_CPV,CNST_CL,CNST_CL,CNST_CI,CNST_CI,CNST_CI/)
  real(8), parameter :: LHV = CNST_LH00
  real(8), parameter :: LHS = CNST_LHS00
  real(8), parameter :: LHF = CNST_LHF00
  !
end module mod_atmos_cnst
