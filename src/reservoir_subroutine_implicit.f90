SUBROUTINE reservoir_subroutine_implicit(res_no,q_surf)
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
    real :: coeff_e1, coeff_e2, coeff_h1, coeff_h2, const_1, const_2
    real :: temp_epi, temp_hyp

    integer :: nd,  res_no

    ! ---------------- turnover loop driven only by T_epil and T_hyp ----------
    density_dif = density_epil(res_no) - density_hypo(res_no)
    if(density_dif .gt. -0.00001 .and. density_dif .lt. 0.00001) density_dif= 0.00001
    n_stability = (-1) * (gravity/density_epil(res_no)) * &
        ((density_dif)/((depth_e(res_no) + depth_h(res_no))/2))
    if(n_stability < 0.000001) n_stability = 0.000001
    log_K_z = log10(n_stability) * (-1)  - 5.699 ! high scenario - 2  w/ adjusted intercept based on empirical equation in Quay et al. 1980, Fig 11
    K_z(res_no) = 10**log_K_z
    if (res_no .eq. 19) write(*,*) &
        'K_z',log_K_z, K_z(res_no), depth_e(res_no), depth_h(res_no)
    ! if (K_z(res_no) > 2) K_z(res_no) = 2
    ! ONLY for no stratification run:
    write(57, *) res_no,K_z(res_no)
    !
    !     Implicit method to solve the energy balance equation
    !
    coeff_e1 = volume_e_x(res_no) + flow_in_epi_x + K_z(res_no) * surface_area(res_no)
    coeff_h1 = - K_z(res_no) * surface_area(res_no)
    const_1  = volume_e_x(res_no)*T_epil(res_no) + flow_in_epi_x*T_res_inflow(res_no) + &
                (q_surf * dt_comp * surface_area(res_no)) / (density * heat_c_kcal)
    coeff_e2 = -(flow_epi_hyp_x + K_z(res_no) * surface_area(res_no))
    coeff_h2 = volume_h_x(res_no) + flow_in_hyp_x + flow_epi_hyp_x  &
                + K_z(res_no) * surface_area(res_no)
    const_2  = volume_h_x(res_no)*T_hypo(res_no) + flow_in_hyp_x*T_res_inflow(res_no)

    temp_epi = (const_1 * coeff_h2 - const_2 * coeff_h1) / &
               (coeff_e1 * coeff_h2 - coeff_e2 * coeff_h1)
    temp_hyp = (const_1 * coeff_e2 - const_2 * coeff_e1) / &
               (coeff_e2 * coeff_h1 - coeff_e1 * coeff_h2)
    if (res_no .eq. 19) write(*,*) &
        coeff_e1, volume_e_x(res_no), flow_in_epi_x, K_z(res_no), surface_area(res_no)
    !
    !     If T of hypolimnion is higher than T of epilimnion, K_z is too large.
    !     Thus, we adjust the temperature to make T_epil equals to T_hypo
    !
    if (temp_hyp .gt. temp_epi) then
        temp_epi = (volume_h_x(res_no)*T_hypo(res_no) + volume_e_x(res_no)*T_epil(res_no) + &
                   (flow_in_epi_x + flow_in_hyp_x) * T_res_inflow(res_no) + &
                   (q_surf * dt_comp * surface_area(res_no)) / (density * heat_c_kcal)) / &
                   (volume_h_x(res_no) + volume_e_x(res_no) + flow_in_epi_x + flow_in_hyp_x)
        temp_hyp = temp_epi
    end if

    T_epil(res_no) = temp_epi
    T_hypo(res_no) = temp_hyp
    !
    !     Update reservoir storage based on water withdrawal
    !
    volume_tot = volume_e_x(res_no)  + volume_h_x(res_no)
    water_withdrawal_epi = (volume_e_x(res_no)/volume_tot) * water_withdrawal(res_no)
    water_withdrawal_hyp = (volume_h_x(res_no)/volume_tot) * water_withdrawal(res_no)
    !
    !     Check whether hypolimnion volume is smaller than minimum hypolimnion volume
    !
    if ((volume_h_x(res_no) + water_withdrawal_hyp) .lt. volume_h_min(res_no)) then
        water_withdrawal_epi = water_withdrawal_epi + water_withdrawal_hyp - &
                                (volume_h_x(res_no) - volume_h_min(res_no))
        water_withdrawal_hyp = volume_h_x(res_no) - volume_h_min(res_no)
    end if
    !
    volume_e_x(res_no) = volume_e_x(res_no) - water_withdrawal_epi
    volume_h_x(res_no) = volume_h_x(res_no) - water_withdrawal_hyp
    volume_tot = volume_e_x(res_no)  + volume_h_x(res_no)
    !
    T_res(res_no) = (T_epil(res_no) * (volume_e_x(res_no)/volume_tot)) + &
        (T_hypo(res_no)*(volume_h_x(res_no)/volume_tot) ) ! weighted averge temp
end subroutine reservoir_subroutine_implicit
