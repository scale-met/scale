!-------------------------------------------------------------------------------
!> module ADMIN VERSIONCHECK
!!
!! @par Description
!!         Check obsolete namelistes
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module mod_admin_versioncheck
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ADMIN_versioncheck

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ADMIN_versioncheck
    use scale_prc, only: &
       PRC_abort
    implicit none

    logical :: dummy

    namelist / PARAM_TIME /                  dummy
    namelist / PARAM_PRC /                   dummy
    namelist / PARAM_COMM /                  dummy
    namelist / PARAM_NEST /                  dummy

    namelist / PARAM_FILEIO /                dummy
    namelist / PARAM_HISTORY /               dummy
    namelist / PARAM_HIST /                  dummy
    namelist / HISTITEM /                    dummy
    namelist / EXTITEM /                     dummy
    namelist / MONITITEM /                   dummy

    namelist / PARAM_INDEX /                 dummy
    namelist / PARAM_GRID /                  dummy
    namelist / PARAM_LAND_INDEX /            dummy
    namelist / PARAM_LAND_GRID /             dummy
    namelist / PARAM_URBAN_INDEX /           dummy
    namelist / PARAM_URBAN_GRID /            dummy

    namelist / PARAM_TRACER /                dummy
    namelist / PARAM_TRACER_KAJINO13 /       dummy

    namelist / PARAM_MAPPROJ /               dummy
    namelist / PARAM_GTRANS /                dummy

    namelist / PARAM_ATMOS_PHY_MP_BIN2BULK / dummy
    namelist / PARAM_BIN /                   dummy
    namelist / NM_MP_SN14_COLLECTION /       dummy
    namelist / NM_MP_SN14_CONDENSATION /     dummy
    namelist / NM_MP_SN14_INIT /             dummy
    namelist / NM_MP_SN14_NUCLEATION /       dummy
    namelist / NM_MP_SN14_PARTICLES /        dummy
    namelist / PARAM_ATMOS_PHY_SF /          dummy
    namelist / PARAM_ATMOS_PHY_SF_BULKCOEF / dummy
    namelist / PARAM_ATMOS_PHY_TB_HYBRID /   dummy
    namelist / PARAM_ATMOS_PHY_TB_MYNN /     dummy

    namelist / PARAM_LAND_PHY_SLAB /         dummy
    namelist / PARAM_LAND_SFC_SLAB /         dummy
    namelist / PARAM_LAND_SFC_THIN_SLAB /    dummy
    namelist / PARAM_LAND_SFC_THICK_SLAB /   dummy

    namelist / PARAM_OCEAN_PHY_SLAB /        dummy
    namelist / PARAM_OCEAN_PHY_FILE /        dummy
    namelist / PARAM_ROUGHNESS /             dummy
    namelist / PARAM_ROUGHNESS_MILLER92 /    dummy
    namelist / PARAM_ROUGHNESS_MOON07 /      dummy

    namelist / PARAM_URBAN_PHY_SLC /         dummy

    namelist / PARAM_SBMAERO /               dummy
    namelist / PARAM_MKINIT_INTERPORATION /  dummy

    logical :: stop4error
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ADMIN_versioncheck",*) 'Check version'

    stop4error = .false.

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PRC,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_PRC is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_PRC_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_COMM is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_COMM_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_NEST,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_NEST is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_COMM_CARTESC_NEST?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_FILEIO,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_FILEIO is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_FILE_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HISTORY,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_HISTORY is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_FILE_HISTORY?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HIST,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_HIST is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_FILE_HISTORY_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=HISTITEM,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist HISTITEM is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use HISTORY_ITEM?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=EXTITEM,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist EXTITEM is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use EXTERNAL_ITEM?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=MONITITEM,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist MONITITEM is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use MONITOR_ITEM?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INDEX,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_INDEX is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_GRID_CARTESC_INDEX?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GRID,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_GRID is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_GRID_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_INDEX,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_INDEX is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_GRID_CARTESC_INDEX?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_GRID,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_GRID is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_GRID_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_INDEX,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_URBAN_INDEX is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_URBAN_GRID_CARTESC_INDEX?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_GRID,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_URBAN_GRID is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_URBAN_GRID_CARTESC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TRACER,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_TRACER is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TRACER_KAJINO13,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_TRACER_KAJINO13 is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MAPPROJ,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_MAPPROJ is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_MAPPROJECTION?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GTRANS,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_GTRANS is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_GRID_CARTESC_METRIC?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_BIN2BULK,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ATMOS_PHY_MP_BIN2BULK is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BIN,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_BIN is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SUZUKI10_bin?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_MP_SN14_COLLECTION,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist NM_MP_SN14_COLLECTION is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SN14_collection?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_MP_SN14_CONDENSATION,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist NM_MP_SN14_CONDENSATION is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SN14_condensation?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_MP_SN14_INIT,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist NM_MP_SN14_INIT is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SN14_init?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_MP_SN14_NUCLEATION,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist NM_MP_SN14_NUCLEATION is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SN14_nucleation?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NM_MP_SN14_PARTICLES,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist NM_MP_SN14_PARTICLES is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_MP_SN14_particles?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ATMOS_PHY_SF is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_ATMOS_PHY_SF_BULK?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_BULKCOEF,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ATMOS_PHY_SF_BULKCOEF is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_HYBRID,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ATMOS_PHY_TB_HYBRID is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_MYNN,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ATMOS_PHY_TB_MYNN is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_SLAB,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_PHY_SLAB is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_DYN_BUCKET?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_SLAB,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_SFC_SLAB is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_SFC_SKIN?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_THIN_SLAB,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_SFC_THIN_SLAB is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_SFC_SKIN?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_THICK_SLAB,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_LAND_SFC_THICK_SLAB is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_LAND_SFC_FIXED_TEMP?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_SLAB,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_OCEAN_PHY_SLAB is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_OCEAN_DYN_SLAB?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_FILE,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_OCEAN_PHY_FILE is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_OCEAN_DYN_SLAB?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ROUGHNESS is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_OCEAN_ROUGHNESS?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS_MILLER92,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ROUGHNESS_MILLER92 is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_OCEAN_PHY_ROUGHNESS_MILLER92?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS_MOON07,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_ROUGHNESS_MOON07 is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_OCEAN_PHY_ROUGHNESS_MOON07?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_PHY_SLC,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_URBAN_PHY_SLC is found.'
       LOG_WARN_CONT(*)                 '=> You expected to use PARAM_URBAN_DYN_KUSAKA01?'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SBMAERO,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_SBMAERO is found.'
       stop4error = .true.
    endif

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_INTERPORATION,iostat=ierr)
    if ( ierr >= 0 ) then !--- exists
       LOG_WARN("ADMIN_versioncheck",*) 'Obsolete namelist PARAM_MKINIT_INTERPORATION is found.'
       stop4error = .true.
    endif

    if ( stop4error ) then !--- exists
       LOG_ERROR("ADMIN_versioncheck",*) 'Obsolete namelist found. Check the log file.'
       call PRC_abort
    else
       LOG_INFO("ADMIN_versioncheck",*) 'Obsolete namelists were not found. OK.'
    endif

    return
  end subroutine ADMIN_versioncheck

end module mod_admin_versioncheck
