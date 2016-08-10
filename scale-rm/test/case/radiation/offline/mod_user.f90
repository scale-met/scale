!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Kusaka et al. (2000,BLM)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  !++ Private parameters & variables
  !
  logical, private :: USER_do   = .false. !< do user step?
  integer, private :: ITIME     =  1      !< record number of input data

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_admin, only: &
       ATMOS_do,          &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       ITIME

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    if( IO_L ) write(IO_FID_LOG,*)

    ! atmosphric model set to off
    !  ATMOS_do         = .false.
       ATMOS_sw_dyn     = .false.
       ATMOS_sw_phy_mp  = .false.
       ATMOS_sw_phy_ae  = .false.
       ATMOS_sw_phy_ch  = .false.
    !  ATMOS_sw_phy_rd  = .false.
       ATMOS_sw_phy_sf  = .false.
       ATMOS_sw_phy_tb  = .false.
       ATMOS_sw_phy_cp  = .false.

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
   use scale_const, only: &
       TEM00 => CONST_TEM00, &
       Rvap  => CONST_Rvap        !< gas constant (water vapor) [J/kg/K]
    use scale_grid_real, only: &
       REAL_lon
    use scale_time, only: &
       dt_RD => TIME_DTSEC_ATMOS_PHY_RD, &
       TIME_NOWDATE
    use scale_landuse, only: &
       LANDUSE_fact_ocean, &
       LANDUSE_fact_land,  &
       LANDUSE_fact_urban
    use mod_atmos_vars, only: &
       TEMP,              &
       PRES,              &
       DENS,              &
       RHOT,              &
       QTRC,              &
       RHOT_t => RHOT_tp
   use scale_atmos_phy_rd_common, only: &
       I_SW, &
       I_LW
    use mod_admin_time, only: &
       TIME_DOATMOS_PHY_RD !< execute physics(radiation)
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,  &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo
    implicit none

    real(RP) :: WORK_3d(KA,IA,JA)
    real(RP) :: WORK_2d(IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do .AND. TIME_DOATMOS_PHY_RD ) then

       call read_rad_inputdata('input/DENS.grd',DENS)
       call read_rad_inputdata('input/RHOT.grd',RHOT)
       call read_rad_inputdata('input/T.grd'   ,TEMP)
       call read_rad_inputdata('input/PRES.grd',PRES)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QV.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QV) = WORK_3d(:,:,:)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QC.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QC) = WORK_3d(:,:,:)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QR.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QR) = WORK_3d(:,:,:)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QI.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QI) = WORK_3d(:,:,:)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QS.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QS) = WORK_3d(:,:,:)
       WORK_3d = 0.0_RP
       call read_rad_inputdata('input/QG.grd'  ,WORK_3d)
       QTRC(:,:,:,I_QG) = WORK_3d(:,:,:)

       call read_rad_inputdata_2d('input/SFC_TEMP.grd',SFC_TEMP)
       call read_rad_inputdata_2d('input/SFC_ALB_LW.grd',WORK_2d)
       SFC_albedo(:,:,I_LW) = WORK_2d(:,:)
       call read_rad_inputdata_2d('input/SFC_ALB_SW.grd',WORK_2d)
       SFC_albedo(:,:,I_SW) = WORK_2d(:,:)

       print *,"DENS",DENS(KA/2,IA/2,JA/2)
       print *,"RHOT",RHOT(KA/2,IA/2,JA/2)
       print *,"T   ",TEMP(KA/2,IA/2,JA/2)
       print *,"PRES",PRES(KA/2,IA/2,JA/2)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QV)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QC)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QR)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QI)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QS)
       print *,"QTRC",QTRC(KA/2,IA/2,JA/2,I_QG)
       print *,"SFC_TEMP",SFC_TEMP(IA/2,JA/2)
       print *,"ALB_LW",SFC_albedo(IA/2,JA/2,I_LW)
       print *,"ALB_SW",SFC_albedo(IA/2,JA/2,I_SW)

    endif

    return
  end subroutine USER_step

!----------------------------------------------------------
  subroutine read_rad_inputdata(filename,var)

    implicit none
    character(*),intent(in)  :: filename
    real(RP),intent(out)     :: var (KA,IA,JA)
    real(SP)                 :: work(IMAX,JMAX,KMAX)
    integer                  :: k, i, j, irecl

     irecl=(IE-IS+1)*(JE-JS+1)*(KE-KS+1)*4
     open(50,file=trim(filename),status='old',access='direct', &
             form='unformatted',recl=irecl)
     read(50,rec=ITIME) WORK
     close(50)

     if(RP==SP)then  ! single
       do k=KS,KE
         do j=JS,JE
           do i=IS,IE
            var(k,i,j) = work(i-IHALO,j-JHALO,k-KHALO)
           enddo
         enddo
       enddo
     else if(RP==DP)then  ! single
       do k=KS,KE
         do j=JS,JE
           do i=IS,IE
            var(k,i,j) = dble(work(i-IHALO,j-JHALO,k-KHALO))
           enddo
         enddo
       enddo
     endif

       do k=KS,KE
         do j=JS,JE
           do i=1,IS-1
            var(k,i,j) = var(k,IS,j)
           enddo
           do i=IE+1,IA
            var(k,i,j) = var(k,IE,j)
           enddo
         enddo
         do j=1,JS-1
           var(k,:,j) = var(k,:,JS)
         enddo
         do j=JE+1,JA
           var(k,:,j) = var(k,:,JE)
         enddo
       enddo
       do k=1,KS-1
         var(k,:,:) = var(KS,:,:)
       enddo
       do k=KE+1,KA
         var(k,:,:) = var(KE,:,:)
       enddo

    return
  end subroutine read_rad_inputdata
!----------------------------------------------------------
  subroutine read_rad_inputdata_2d(filename,var)

    implicit none
    character(*),intent(in)  :: filename
    real(RP),intent(out)     :: var (IA,JA)
    real(SP)                 :: work(IMAX,JMAX)
    integer                  :: i, j, irecl

     irecl=(IE-IS+1)*(JE-JS+1)*4
     open(50,file=trim(filename),status='old',access='direct', &
             form='unformatted',recl=irecl)
     read(50,rec=ITIME) WORK
     close(50)

     if(RP==SP)then  ! single
       do j=JS,JE
       do i=IS,IE
        var(i,j) = work(i-IHALO,j-JHALO)
       enddo
       enddo
     else if(RP==DP)then  ! single
       do j=JS,JE
       do i=IS,IE
        var(i,j) = dble(work(i-IHALO,j-JHALO))
       enddo
       enddo
     endif

      do j=JS,JE
         do i=1,IS-1
            var(i,j) = var(IS,j)
         enddo
         do i=IE+1,IA
            var(i,j) = var(IE,j)
         enddo
      enddo
      do j=1,JS-1
           var(:,j) = var(:,JS)
      enddo
      do j=JE+1,JA
           var(:,j) = var(:,JE)
      enddo

    return
  end subroutine read_rad_inputdata_2d

end module mod_user
