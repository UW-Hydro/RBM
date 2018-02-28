SUBROUTINE flow_subroutine (res_no, nyear, nd)

    use Block_Hydro
    use Block_Network
    use Block_Reservoir

    implicit none

    !real :: flow_in_epi_x, flow_in_hyp_x, flow_epi_hyp_x, flow_out_hyp_x,flow_out_epi_x
    real :: ratio_sp, ratio_pen
    real :: Q1, Q2
    !real :: res_vol_delta_x, vol_change_hyp_x, vol_change_epi_x
    integer :: res_no, nnd, nyear, nd
 
 
    !*************************************************************************
    !      read inflow from vic/rvic simulations
    !*************************************************************************
    if (.not.adjust_timestep) then
        res_storage_pre = res_storage(res_no,2)
        res_storage_post = res_storage(res_no,1)
    end if
    !res_storage_pre = res_capacity_mcm(res_no) * (10**6) * 0.95 ! TESTINFLOW
    !res_storage_post = res_capacity_mcm(res_no) * (10**6) * 0.95 ! TESTINFLOW

    Q1 = Q_res_inflow(res_no) * ftsec_to_msec * dt_res
    ! converts ft^3/sec to m^3/sec, multiplys by seconds per time step

    if ( density_in(res_no) .le. density_hypo(res_no) ) then
        flow_in_hyp_x = Q1 * 0
        flow_in_epi_x = Q1 * 1   ! all flow goes to epil.
    else
        flow_in_hyp_x = Q1 * 1
        flow_in_epi_x = Q1 * 0   ! all flow goes to hypo
    end if
    !
    !     calculate the reservoir storage
    !
    Q2 = Q_res_outflow(res_no) * ftsec_to_msec * dt_res
    !
    !     Initialization of reservoir storage
    !
    if (initial_storage(res_no)) then
        water_withdrawal(res_no) = 0
        res_depth_meter(res_no) =  res_storage(res_no,1) / surface_area(res_no)
        depth_e(res_no) = res_depth_meter(res_no) * depth_e_frac
        depth_h(res_no) = res_depth_meter(res_no) * depth_h_frac
        volume_e_x(res_no) = surface_area(res_no) * depth_e(res_no)
        volume_h_x(res_no) = surface_area(res_no) * depth_h(res_no)
        res_storage_pre = res_storage(res_no,1)
        initial_storage(res_no)=.FALSE.
    end if
    !
    !     calculate water withdrawal based on inflow/outflow and storage change
    !
    if (nsub.eq.1) water_withdrawal(res_no) = Q1 - Q2 - (res_storage_post - res_storage_pre)/numsub
    flow_out_hyp_x = Q2 ! * ftsec_to_msec * dt_res
    flow_out_epi_x = 0
    flow_epi_hyp_x = flow_in_epi_x
    !
    !     if res_storage_post > res_storage_pre, when reservoir is storing water
    !
    if (Q1 - Q2 .gt. 0 .and. flow_in_epi_x .gt. 0 &
        .and. volume_e_x(res_no) .lt. res_capacity_mcm(res_no)*(10**6)*0.3) then
        flow_epi_hyp_x = flow_out_hyp_x
    end if
    !
    ! based on inflow and outflow
    !
    vol_change_epi_x = flow_in_epi_x - flow_out_epi_x - flow_epi_hyp_x
    vol_change_hyp_x = flow_in_hyp_x + flow_epi_hyp_x - flow_out_hyp_x
    !
    !     Check whether hypolimnion volume is smaller than minimum hypolimnion volume
    !
    !volume_h_min(res_no) = 0 !res_capacity_mcm(res_no) * (10**6) * 0.1
    volume_h_min(res_no) = res_capacity_mcm(res_no) * (10**6) * 0.05
    if ((volume_h_x(res_no) + vol_change_hyp_x) .lt. volume_h_min(res_no)) then
        vol_change_epi_x =  vol_change_epi_x + vol_change_hyp_x + & 
                           (volume_h_x(res_no) - volume_h_min(res_no))
        vol_change_hyp_x = - (volume_h_x(res_no) - volume_h_min(res_no))
    end if
    
    volume_e_min(res_no) = res_capacity_mcm(res_no) * (10**6) * 0.05
    !if (res_no .eq. 19) write(*,*) volume_h_min(res_no), volume_e_min(res_no) 
    if ((volume_e_x(res_no) + vol_change_epi_x) .lt. volume_e_min(res_no)) then
        vol_change_hyp_x = vol_change_epi_x + vol_change_hyp_x + &
                           (volume_e_x(res_no) - volume_e_min(res_no))
        vol_change_epi_x = - (volume_e_x(res_no) - volume_e_min(res_no))
    end if 
    !----- update epilimnion and hypolimnion volume  -------
    volume_e_x(res_no) = volume_e_x(res_no) + vol_change_epi_x
    volume_h_x(res_no) = volume_h_x(res_no) + vol_change_hyp_x
    outflow_x = flow_out_epi_x + flow_out_hyp_x
    ! ------------------------- calculate dV/dt terms ---------------------------
    dV_dt_epi(res_no) = vol_change_epi_x * T_epil(res_no)
    dV_dt_hyp(res_no) = vol_change_hyp_x * T_hypo(res_no)

    depth_e(res_no) = volume_e_x(res_no) / surface_area(res_no)
    depth_h(res_no) = volume_h_x(res_no) / surface_area(res_no)
    !if(res_no .eq. 11) write(*,*) vol_change_epi_x, &
    !        vol_change_hyp_x, depth_e(res_no), depth_h(res_no) 
END SUBROUTINE flow_subroutine
