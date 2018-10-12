SUBROUTINE reservoir_subroutine(nresx, nd, q_surf,time, nd_year, nyear)
    use Block_Reservoir
    use Block_Flow
    use Block_Network

    implicit none

    real :: dayx, q_surf, log_K_z, n_stability, density_dif
    real(8):: time
    integer :: nd,  nresx, nyear, nd_year
    ! ---------------- turnover loop driven only by T_epil and T_hyp ----------
    density_dif = density_epil(nresx) - density_hypo(nresx)
    if(density_dif .gt. -0.00001 .and. density_dif .lt. 0.00001) density_dif= 0.00001
    n_stability = (-1) * (gravity/density_epil(nresx)) * & 
        ((density_dif)/((depth_e(nresx) + depth_h(nresx))/2))
    if(n_stability < 0.000001) n_stability = 0.000001
    log_K_z = log10(n_stability) * (-0.65)  - 3.1 ! high scenario - 2  w/ adjusted intercept based on empirical equation in Quay et al. 1980, Fig 11
    K_z(nresx) = 10**log_K_z
    if (K_z(nresx) > 2) K_z(nresx) = 2
    write(57, *) nresx,K_z(nresx)

    ! -------------------- calculate temperature terms  -------------------------
    dif_epi_x  = K_z(nresx) * surface_area(nresx) *  (T_hypo(nresx) - T_epil(nresx)) / volume_e_x(nresx)
    dif_hyp_x  = K_z(nresx) * surface_area(nresx) *  (T_epil(nresx) - T_hypo(nresx)) / volume_h_x(nresx)

    ! --------------------- calculate advection terms --------------------------- 
    if(flow_in_epi_x .lt. volume_e_x(nresx)) then
        advec_in_epix  = flow_in_epi_x * (T_res_in(nresx) - T_epil(nresx)) /volume_e_x(nresx)
    else
        advec_in_epix  =  volume_e_x(nresx) * (T_res_in(nresx) - T_epil(nresx))/volume_e_x(nresx)
    end if

    if(flow_in_hyp_x .lt. volume_h_x(nresx)) then
        advec_in_hypx = flow_in_hyp_x * (T_res_in(nresx) - T_hypo(nresx)) /volume_h_x(nresx)
    else
        advec_in_hypx = volume_h_x(nresx) * (T_res_in(nresx) - T_hypo(nresx))/volume_h_x(nresx)
    end if

    advec_epi_hyp = flow_epi_hyp_x *  (T_epil(nresx) - T_hypo(nresx)) /volume_h_x(nresx)
    ! ------------------- calculate change in temperature  ---------------------

    ! ---------------- epilimnion -----------

    ! ------------ calculate total energy ----------

    energy_x  = (q_surf * dt_comp ) / (depth_e(nresx) * density * heat_c_kcal ) ! kcal/sec*m2 to C/day
    temp_change_ep(nresx) = advec_in_epix + dif_epi_x + energy_x

    !----- update epilimnion volume for next time step -------

    T_epil(nresx) = T_epil(nresx) + temp_change_ep(nresx)
    if(T_epil(nresx) .lt. 0) T_epil(nresx) =  0.5 

    ! ------------------ hypolimnion ----------------

    ! ------------ calculate total energy ----------
    temp_change_hyp(nresx) = advec_in_hypx + advec_epi_hyp + dif_hyp_x !  - advec_out_hypx - dV_dt_hyp(nresx)

    !----- update epilimnion volume for next time step -------
    T_hypo(nresx) = T_hypo(nresx) +  temp_change_hyp(nresx)
    if(T_hypo(nresx) .lt. 0) T_hypo(nresx) = 0.5 
    !---------- calculate combined (hypo. and epil.) temperature of outflow -----
    epix = T_epil(nresx)*(flow_out_epi_x/outflow_x)  ! portion of temperature from epilim. 
    hypox= T_hypo(nresx)*(flow_out_hyp_x/outflow_x)  ! portion of temperature from hypol.
    temp_out(nresx) = epix + hypox   ! average outflow temperature
    volume_tot = volume_e_x(nresx)  + volume_h_x(nresx)
    T_res(nresx) = (T_epil(nresx) * (volume_e_x(nresx)/volume_tot)) + &
        (T_hypo(nresx)*(volume_h_x(nresx)/volume_tot) ) ! weighted averge temp
    ! ---------------- write data for reservoir_file ----------
    diffusion_tot(nresx) = dif_hyp_x
    advec_hyp_tot(nresx) = advec_in_hypx
    advec_epi_tot(nresx) = advec_in_epix
    qsurf_tot(nresx) = energy_x
    ! non-essential - only to print out specific calculated variables
    if(nresx.eq.3) then
 
        write(48,*) time,T_res_in(nresx),T_epil(nresx), advec_in_epix, energy_x, dif_epi_x, flow_in_epi_x, volume_e_x(nresx) &
            ,T_hypo(nresx), advec_in_hypx, advec_epi_hyp, dif_hyp_x, flow_in_hyp_x,  volume_h_x(nresx) &
            ,  flow_epi_hyp_x,  temp_out(nresx)
    end if

    if(nresx.eq.21) then

        write(89,*) time,T_res_in(nresx),T_epil(nresx), advec_in_epix, energy_x,dif_epi_x, flow_in_epi_x, volume_e_x(nresx) &
            ,T_hypo(nresx), advec_in_hypx, advec_epi_hyp,dif_hyp_x, flow_in_hyp_x,  volume_h_x(nresx) &
            ,  flow_epi_hyp_x,  temp_out(nresx)


    end if


end subroutine reservoir_subroutine
