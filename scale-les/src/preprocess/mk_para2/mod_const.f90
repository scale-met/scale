!-------------------------------------------------------------------
!
!+  Module : definition of constants for general use
!
!-------------------------------------------------------------------
!
!--- History
!
!    2005/09/17  K.Suzuki  created for F90 version
!
!-------------------------------------------------------------------
module mod_const

implicit none
private

real(8), public, save :: CNST_RHOW  = 1.0D+03     !  density of water
real(8), public, save :: CNST_RHOA  = 2.25D+03    !  density of aerosol
real(8), public, save :: CNST_PI    = 3.141592D0  !  circle constant
real(8), public, save :: CNST_CP    = 1004.D0     !  specific heat
real(8), public, save :: CNST_GRAV  = 9.8D0       !  gravity
real(8), public, save :: CNST_RVAP  = 461.D0      !  gas constant for water vapor
real(8), public, save :: CNST_RAIR  = 287.04D0    !  gas constant for dry air
real(8), public, save :: CNST_QLEVP = 2.50D+06    !  latent heat of evaporation
real(8), public, save :: CNST_QLSBL = 2.83D+06    !  latent heat of sublimation
real(8), public, save :: CNST_QLMLT = 3.3D+05     !  latent heat of melting
real(8), public, save :: CNST_TMP0  = 273.15D0    !  freezing temperature
real(8), public, save :: CNST_TMLT  = 273.15D0    !  melting temperature
real(8), public, save :: CNST_ESAT0 = 611.D0      !  saturation vapor pressure at freezing level
real(8), public, save :: CNST_UNDEF = -999.D0     !  undefined value

end module mod_const
