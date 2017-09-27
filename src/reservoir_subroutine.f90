SUBROUTINE reservoir_subroutine(res_no,q_surf)
    use Block_Network
    use Block_Reservoir
    !
    implicit none
    !
    real :: dayx, q_surf, log_K_z, n_stability, density_dif
    real :: dif_epi_x, dif_hyp_x, advec_in_epix, advec_in_hypx
    real :: advec_epi_hyp, energy_x
    real :: epix,hypox,volume_tot
    real :: water_withdrawal_epi, water_withdrawal_hyp

    integer :: nd,  res_no

    ! ---------------- turnover loop driven only by T_epil and T_hyp ----------
    density_dif = density_epil(res_no) - density_hypo(res_no)
    if(density_dif .gt. -0.00001 .and. density_dif .lt. 0.00001) density_dif= 0.00001
    n_stability = (-1) * (gravity/density_epil(res_no)) * &
        ((density_dif)/((depth_e(res_no) + depth_h(res_no))/2))
    if(n_stability < 0.000001) n_stability = 0.000001
    log_K_z = log10(n_stability) * (-1)  - 5.699 ! high scenario - 2  w/ adjusted intercept based on empirical equation in Quay et al. 1980, Fig 11
    K_z(res_no) = 10**log_K_z
    ! if (K_z(res_no) > 2) K_z(res_no) = 2
    ! ONLY for no stratification run:
    write(57, *) res_no,K_z(res_no)



    ! -------------------- calculate temperature terms  -------------------------
    dif_epi_x  = K_z(res_no) * surface_area(res_no) *  (T_hypo(res_no) - T_epil(res_no)) / volume_e_x(res_no)
    dif_hyp_x  = K_z(res_no) * surface_area(res_no) *  (T_epil(res_no) - T_hypo(res_no)) / volume_h_x(res_no)

    ! --------------------- calculate advection terms ---------------------------
    if(flow_in_epi_x .lt. volume_e_x(res_no)) then
        advec_in_epix  = flow_in_epi_x * (T_res_inflow(res_no) - T_epil(res_no)) /volume_e_x(res_no)
    else
        advec_in_epix  =  volume_e_x(res_no) * (T_res_inflow(res_no) - T_epil(res_no))/volume_e_x(res_no)
    end if

    if(flow_in_hyp_x .lt. volume_h_x(res_no)) then
        advec_in_hypx = flow_in_hyp_x * (T_res_inflow(res_no) - T_hypo(res_no)) /volume_h_x(res_no)
    else
        advec_in_hypx = volume_h_x(res_no) * (T_res_inflow(res_no) - T_hypo(res_no))/volume_h_x(res_no)
    end if

    advec_epi_hyp = flow_epi_hyp_x *  (T_epil(res_no) - T_hypo(res_no)) /volume_h_x(res_no)

    if(res_no .eq. 28) then
    write(*,*) 'advec_in_epi',advec_in_epix, 'flow_in_epi',flow_in_epi_x, 'T_res_inflow', &
            T_res_inflow(res_no),'T_epil)', T_epil(res_no), &
            'vol',volume_e_x(res_no)
    end if
    ! ------------------- calculate change in temperature  ---------------------

    ! ---------------- epilimnion -----------

    ! ------------ calculate total energy ----------
    energy_x  = (q_surf * dt_comp ) / (depth_e(res_no) * density * heat_c_kcal ) ! kcal/sec*m2 to C/day

    temp_change_ep(res_no) = advec_in_epix + dif_epi_x + energy_x

    !----- update epilimnion volume for next time step -------
    T_epil(res_no) = T_epil(res_no) + temp_change_ep(res_no)
    if(T_epil(res_no) .lt. 0) T_epil(res_no) =  0.5

    ! ------------------ hypolimnion ----------------

    ! ------------ calculate total energy ----------
    temp_change_hyp(res_no) = advec_in_hypx + advec_epi_hyp + dif_hyp_x !  - advec_out_hypx - dV_dt_hyp(res_no)

    !----- update epilimnion volume for next time step -------
    T_hypo(res_no) = T_hypo(res_no) +  temp_change_hyp(res_no)
    if(T_hypo(res_no) .lt. 0) T_hypo(res_no) = 0.5
    !---------- calculate combined (hypo. and epil.) temperature of outflow -----
    epix = T_epil(res_no)*(flow_out_epi_x/outflow_x)  ! portion of temperature from epilim.
    hypox= T_hypo(res_no)*(flow_out_hyp_x/outflow_x)  ! portion of temperature from hypol.
    temp_out(res_no) = epix + hypox   ! average outflow temperature
    volume_tot = volume_e_x(res_no)  + volume_h_x(res_no)
    !
    !     Update reservoir storage based on water withdrawal
    !
    water_withdrawal_epi = (volume_e_x(res_no)/volume_tot) * water_withdrawal(res_no)
    water_withdrawal_hyp = (volume_h_x(res_no)/volume_tot) * water_withdrawal(res_no)

    volume_e_x(res_no) = volume_e_x(res_no) - water_withdrawal_epi
    volume_h_x(res_no) = volume_h_x(res_no) - water_withdrawal_hyp
    volume_tot = volume_e_x(res_no)  + volume_h_x(res_no)

    T_res(res_no) = (T_epil(res_no) * (volume_e_x(res_no)/volume_tot)) + &
        (T_hypo(res_no)*(volume_h_x(res_no)/volume_tot) ) ! weighted averge temp


    ! ---------------- write data for reservoir_file ----------
    diffusion_tot(res_no) = dif_hyp_x
    advec_hyp_tot(res_no) = advec_in_hypx
    advec_epi_tot(res_no) = advec_in_epix
    qsurf_tot(res_no) = energy_x


end subroutine reservoir_subroutine
