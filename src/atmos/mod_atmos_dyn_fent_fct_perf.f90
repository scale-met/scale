module mod_atmos_dyn_fent_fct_perf
  include 'inc_index.h'
  ! Follow the Folding@Home convention. See
  ! http://folding.stanford.edu/English/FAQ-flops and
  ! http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
  integer, parameter :: EXP_OPS = 61
  
contains

  integer*8 function calc_ops_2d(is_off, ie_off, js_off, &
       je_off)
    implicit none
    integer, intent(in) :: is_off, ie_off, js_off, je_off
    calc_ops_2d = ((IE+ie_off) - (IS+is_off) + 1) * ((JE+je_off) - (JS+js_off) + 1)
    return
  end function calc_ops_2d

  subroutine update_ops_2d(is_off, ie_off, js_off, &
       je_off,fp_ops, fp_total_ops, &
       ld_ops, ld_total_ops, st_ops, st_total_ops)
    implicit none
    integer, intent(in) :: is_off, ie_off, js_off, je_off, &
         fp_ops, ld_ops, st_ops
    integer*8, intent(inout) :: fp_total_ops, ld_total_ops, st_total_ops
    integer*8 num_ops

    num_ops = calc_ops_2d(is_off, ie_off, js_off, je_off)
    fp_total_ops = fp_total_ops + num_ops * fp_ops
    ld_total_ops = ld_total_ops + num_ops * ld_ops
    st_total_ops = st_total_ops + num_ops * st_ops        
  end subroutine update_ops_2d

  integer*8 function calc_ops_3d(is_off, ie_off, js_off, &
       je_off, ks_off, ke_off)
    implicit none
    integer, intent(in) :: is_off, ie_off, js_off, je_off, &
         ks_off, ke_off
    calc_ops_3d = ((IE+ie_off) - (IS+is_off) + 1) * ((JE+je_off) - (JS+js_off) + 1) &
         * ((KE+ke_off) - (KS+ks_off) + 1)
    return
  end function calc_ops_3d

  subroutine update_ops_3d(is_off, ie_off, js_off, &
       je_off, ks_off, ke_off, fp_ops, fp_total_ops, &
       ld_ops, ld_total_ops, st_ops, st_total_ops)
    implicit none
    integer, intent(in) :: is_off, ie_off, js_off, je_off, &
         ks_off, ke_off, fp_ops, ld_ops, st_ops
    integer*8, intent(inout) :: fp_total_ops, ld_total_ops, st_total_ops
    integer*8 num_ops

    num_ops = calc_ops_3d(is_off, ie_off, js_off, je_off, &
         ks_off, ke_off)
    fp_total_ops = fp_total_ops + num_ops * fp_ops
    ld_total_ops = ld_total_ops + num_ops * ld_ops
    st_total_ops = st_total_ops + num_ops * st_ops        
  end subroutine update_ops_3d
  
  subroutine rk_ops(fp_ops, ld_ops, st_ops)
    integer*8, intent(out) :: fp_ops, ld_ops, st_ops
    
    fp_ops = 0
    ld_ops = 0
    st_ops = 0

    ! momentum -> velocity
    call update_ops_3d(0, 1, 0, 1, 0, -1, 3, fp_ops, 3, ld_ops, 1, st_ops)
    call update_ops_3d(-1, 1, 0, 1, 0, 0, 3, fp_ops, 3, ld_ops, 1, st_ops)
    call update_ops_3d(0, 1, -1, 1, 0, 0, 3, fp_ops, 3, ld_ops, 1, st_ops)

    ! pressure, pott. temp.
    call update_ops_3d(-2, 2, -2, 2, 0, 0, EXP_OPS + 4, fp_ops, 4, ld_ops, &
         2, st_ops)

    ! continuity equation
    ! at (x, y, interface)
    call update_ops_3d(0, 0, 0, 0, 0, -1, 2, fp_ops, 2, ld_ops, 1, st_ops)
    !!!!!!!!!!!!!!!!!!!
    ! at (u, y, layer)
    call update_ops_3d(-1, 0, 0, 0, 0, 0, 2, fp_ops, 2, ld_ops, 1, st_ops)
    ! at (x, v, layer)
    call update_ops_3d(0, 0, -1, 0, 0, 0, 2, fp_ops, 2, ld_ops, 1, st_ops)

    ! update density
    call update_ops_3d(0, 0, 0, 0, 0, 0, 10, fp_ops, 10, ld_ops, 1, st_ops)

    ! momentum equation (z)
    ! at (x, y, layer)
    call update_ops_3d(0, 0, 0, 0, 0, 0, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! at (u, y, interface)
    call update_ops_3d(-1, 0, 0, 0, 0, -1, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! at (x, v, interface)
    call update_ops_3d(0, 0, -1, 0, 0, -1, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! update momentum (z)
    call update_ops_3d(0, 0, 0, 0, 0, -1, 18, fp_ops, 16, ld_ops, 1, st_ops)

    ! momentum equation (x)
    ! at (u, y, interface)
    call update_ops_3d(0, 0, 0, 0, 1, -2, 10, fp_ops, 7, ld_ops, 1, st_ops)

    call update_ops_2d(0, 0, 0, 0, 12, fp_ops, 9, ld_ops, 4, st_ops)

    ! at (x, y, layer)
    call update_ops_3d(0, 1, 0, 0, 0, 0, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! at (u, v, layer)
    call update_ops_3d(0, 0, -1, 0, 0, 0, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! update momentum (x)
    call update_ops_3d(0, 0, 0, 0, 0, 0, 21, fp_ops, 20, ld_ops, 1, st_ops)

    ! momentum equation (y)
    ! at (x, v, interface)
    call update_ops_3d(0, 0, 0, 0, 1, -2, 10, fp_ops, 7, ld_ops, 1, st_ops)

    call update_ops_2d(0, 0, 0, 0, 12, fp_ops, 10, ld_ops, 4, st_ops)

    ! at (u, v, layer)
    call update_ops_3d(-1, 0, 0, 0, 0, 0, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! at (x, y, layer)
    call update_ops_3d(0, 0, 0, 1, 0, 0, 10, fp_ops, 7, ld_ops, 1, st_ops)

    ! update momentum (y)
    call update_ops_3d(0, 0, 0, 0, 0, 0, 21, fp_ops, 20, ld_ops, 1, st_ops)

    ! Thermodynamic Equation

    ! at (x, y, interface)
    call update_ops_3d(0, 0, 0, 0, 1, -2, 9, fp_ops, 6, ld_ops, 1, st_ops)

    call update_ops_2d(0, 0, 0, 0, 10, fp_ops, 8, ld_ops, 4, st_ops)

    ! at (u, y, layer)
    call update_ops_3d(-1, 0, 0, 0, 0, 0, 9, fp_ops, 6, ld_ops, 1, st_ops)

    ! at (x, v, layer)
    call update_ops_3d(0, 0, -1, 0, 0, 0,  9, fp_ops, 6, ld_ops, 1, st_ops)

    ! update rho*theta
    call update_ops_3d(0, 0, 0, 0, 0, 0, 11, fp_ops, 11, ld_ops, 1, st_ops)

    !write (*,*) "#FP ops: ", fp_ops
    !write (*,*) "#LD ops: ", ld_ops
    !write (*,*) "#ST ops: ", st_ops        
    
  end subroutine rk_ops

  subroutine rk_min_ops(ld_ops, st_ops)
    integer*8, intent(out) :: ld_ops, st_ops
    
    ld_ops = 0
    st_ops = 0

    ! MOMZ
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, -2, 1)
    ! DENS
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! VELZ
    st_ops = st_ops + calc_ops_3d(0, 1, 0, 1, 0, -1)
    ld_ops = ld_ops + calc_ops_3d(0, 1, 0, 1, -1, 0)
    ! MOMX
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! VELX
    st_ops = st_ops + calc_ops_3d(-1, 1, 0, 1, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(-1, 1, 0, 1, 0, 0)
    ! MOMY
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! VELY
    st_ops = st_ops + calc_ops_3d(0, 1, -1, 1, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 1, -1, 1, 0, 0)
    ! RHOT
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! RTOT
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! PRES
    st_ops = st_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 1, 0, 1, 0, 0)
    ! POTT
    st_ops = st_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(-2, 2, -2, 2, 0, 0)
    ! mflx_hi z
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ! mflx_hi x
    st_ops = st_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ! mflx_hi y
    st_ops = st_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ! DENS0
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! DENS_RK
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! qflx_hi z
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ! qflx_hi x
    st_ops = st_ops + calc_ops_3d(-1, 0, 0, 0, 0, -1)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, -1)
    st_ops = st_ops + calc_ops_3d(0, 1, 0, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 1, 0, 0, 0, 0)
    st_ops = st_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    st_ops = st_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ! qflx_hi y
    st_ops = st_ops + calc_ops_3d(0, 0, -1, 0, 0, -1)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, -1)    
    st_ops = st_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 1, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 1, 0, 0)
    st_ops = st_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ! num_diff(,,,I_DENS,ZDIR)
    ld_ops = ld_ops + calc_ops_3d(0, -1, 0, 0, 0, 0)
    ! num_diff(,,,I_DENS,XDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ! num_diff(,,,I_DENS,YDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, -1, 0)
    ! num_diff(,,,I_MOMZ,ZDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! num_diff(,,,I_MOMZ,XDIR)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, -1)
    ! num_diff(,,,I_MOMZ,YDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, -1)
    ! num_diff(,,,I_MOMX,ZDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! num_diff(,,,I_MOMX,XDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 1, 0, 0, 0, 0)
    ! num_diff(,,,I_MOMX,YDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ! num_diff(,,,I_MOMY,ZDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! num_diff(,,,I_MOMY,XDIR)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ! num_diff(,,,I_MOMY,YDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 1, 0, 0)
    ! num_diff(,,,I_RHOT,ZDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! num_diff(,,,I_RHOT,XDIR)
    ld_ops = ld_ops + calc_ops_3d(-1, 0, 0, 0, 0, 0)
    ! num_diff(,,,I_RHOT,YDIR)
    ld_ops = ld_ops + calc_ops_3d(0, 0, -1, 0, 0, 0)
    ! RCDZ
    ld_ops = ld_ops + KE-KS+1
    ! RCDX
    ld_ops = ld_ops + IE-IS+1
    ! RCDY
    ld_ops = ld_ops + JE-JS+1
    ! RFDZ
    ld_ops = ld_ops + KE-KS
    ! RFDX
    ld_ops = ld_ops + IE-IS+1
    ! RFDY
    ld_ops = ld_ops + JE-JS+1
    ! ray_damp(,,,I_MOMZ)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! ray_damp(,,,I_MOMX)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)    
    ! ray_damp(,,,I_MOMY)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)    
    ! ray_damp(,,,I_RHOT)
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! CORIOLI
    ld_ops = ld_ops + (IE+1-IS+1)*(JE+1-JS+1)
    ! MOMZ0
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! MOMZ_RK
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, -1)
    ! MOMX0
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! MOMX_RK
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! MOMY0
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! MOMY_RK
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! RHOT0
    ld_ops = ld_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
    ! RHOT_RK
    st_ops = st_ops + calc_ops_3d(0, 0, 0, 0, 0, 0)
  end subroutine rk_min_ops


end module mod_atmos_dyn_fent_fct_perf

