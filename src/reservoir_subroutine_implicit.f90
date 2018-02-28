SUBROUTINE reservoir_subroutine_implicit(res_no,q_surf,nd,tair)
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
    real :: temp_epi_pre, temp_hyp_pre

    integer :: nd,  res_no, tair

    logical :: mixed

    mixed=.FALSE.
    temp_epi_pre = T_epil(res_no)
    temp_hyp_pre = T_hypo(res_no)
    ! ---------------- turnover loop driven only by T_epil and T_hyp ----------
    density_dif = density_epil(res_no) - density_hypo(res_no)
    !if(density_dif .gt. -0.00001 .and. density_dif .lt. 0.00001) density_dif= 0.00001
    !if(density_dif .gt. -0.00001) density_dif= 0.00001
    n_stability = (-1) * (gravity/density_epil(res_no)) * &
        ((density_dif)/((depth_e(res_no) + depth_h(res_no))/2))
    if(n_stability .le. 0) then
        n_stability = 1e-10
    end if
    log_K_z = log10(n_stability) * (-1)  - 5.699 ! high scenario - 2  w/ adjusted intercept based on empirical equation in Quay et al. 1980, Fig 11
    K_z(res_no) = 10**log_K_z
    ! ONLY for no stratification run:
    write(57, *) res_no,K_z(res_no)
    !
    !     Estimate the numerical error in reservoir temperature 
    !
!    if (temp_epi_pre-temp_hyp_pre>0.01) then
!        m11 = -(flow_in_epi_x+K_z(res_no)*surface_area(res_no))/volume_e_x(res_no)/dt_res
!        m12 = K_z(res_no)*surface_area(res_no)/volume_e_x(res_no)/dt_res
!        m21 = (flow_epi_hyp_x+K_z(res_no)*surface_area(res_no))/volume_h_x(res_no)/dt_res
!        m22 = -(flow_in_hyp_x+flow_epi_hyp_x+K_z(res_no)*surface_area(res_no))/volume_h_x(res_no)/dt_res
!        const1 = ((q_surf*dt_res*surface_area(res_no))/(density*heat_c_kcal) + &
!                 flow_in_epi_x*T_res_inflow(res_no))/volume_e_x(res_no)/dt_res
!        const2 = flow_in_hyp_x*T_res_inflow(res_no)/volume_h_x(res_no)/dt_res
!        error_e = -0.5*((m11**2+m12*m21)*temp_epi_pre + &
!                        (m11*m12+m12*m22)*temp_hyp_pre + &
!                         m11*const1 + m12*const2)*dt_res**2
!        error_h = -0.5*((m21*m11+m22*m21)*temp_epi_pre + &
!                        (m21*m12+m22**2)*temp_hyp_pre + &
!                         m21*const1 + m22*const2)*dt_res**2
!    else
!        m11 = -(flow_in_epi_x+flow_in_hyp_x)/(volume_e_x(res_no)+volume_h_x(res_no))/dt_res
!        const1 = ((flow_in_epi_x+flow_in_hyp_x)*T_res_inflow(res_no) + &
!                 (q_surf*dt_res*surface_area(res_no))/(density*heat_c_kcal))&
!                 /(volume_e_x(res_no)+volume_h_x(res_no))/dt_res
!        error_e = -0.5*(m11**2*temp_epi_pre+m11*const1)
!        error_h = error_e
!    end if
!    if (res_no.eq.66) write(*,*) res_no,nd,dt_res,error_e,error_h
!    if (error_e.gt.error_threshold.or.error_h.gt.error_threshold) then
!        exceed_error_bound=.TRUE.
!        go to 500
!    end if
    !
    !     Implicit method to solve the energy balance equation
    !
    coeff_e1 = volume_e_x(res_no) + flow_in_epi_x + K_z(res_no) * surface_area(res_no)
    coeff_h1 = - K_z(res_no) * surface_area(res_no)
    const_1  = volume_e_x(res_no)*temp_epi_pre + flow_in_epi_x*T_res_inflow(res_no) + &
                (q_surf * dt_res * surface_area(res_no)) / (density * heat_c_kcal)
    coeff_e2 = -(flow_epi_hyp_x + K_z(res_no) * surface_area(res_no))
    coeff_h2 = volume_h_x(res_no) + flow_in_hyp_x + flow_epi_hyp_x  &
                + K_z(res_no) * surface_area(res_no)
    const_2  = volume_h_x(res_no)*temp_hyp_pre + flow_in_hyp_x*T_res_inflow(res_no)

    temp_epi = (const_1 * coeff_h2 - const_2 * coeff_h1) / &
               (coeff_e1 * coeff_h2 - coeff_e2 * coeff_h1)
    temp_hyp = (const_1 * coeff_e2 - const_2 * coeff_e1) / &
               (coeff_e2 * coeff_h1 - coeff_e1 * coeff_h2)
    if (temp_epi - temp_hyp .gt. 0.2 .and. .NOT. res_stratif_start(res_no)) then
        res_stratif_start(res_no) = .TRUE.
    end if
    
    if (temp_epi - temp_hyp .lt. 1 .and. nd.gt.180 .and. .NOT. res_turnover(res_no)) then
        res_turnover(res_no) = .TRUE.
    end if
    !
    !     If T of hypolimnion is higher than T of epilimnion, K_z is too large.
    !     Thus, we adjust the temperature to make T_epil equals to T_hypo
    !
    !if (.NOT. res_stratif_start(res_no)) then
    if (temp_epi - temp_hyp .le. 0 .or. res_turnover(res_no)) then
        temp_epi = (volume_h_x(res_no)*T_hypo(res_no) + volume_e_x(res_no)*T_epil(res_no) + &
                   (flow_in_epi_x + flow_in_hyp_x) * T_res_inflow(res_no) + &
                   (q_surf * dt_res * surface_area(res_no)) / (density * heat_c_kcal)) / &
                   (volume_h_x(res_no) + volume_e_x(res_no) + flow_in_epi_x + flow_in_hyp_x)
        temp_hyp = temp_epi
    end if
    !if (res_no.eq.66) write(99,'(8f30.6)') &
    !    volume_h_x(res_no),volume_e_x(res_no),flow_in_epi_x,flow_in_hyp_x, &
    !    volume_h_x(res_no)*T_hypo(res_no), volume_e_x(res_no)*T_epil(res_no),&
    !    (flow_in_epi_x + flow_in_hyp_x) * T_res_inflow(res_no), &
    !    q_surf * dt_res * surface_area(res_no)
    T_epil(res_no) = temp_epi
    T_hypo(res_no) = temp_hyp
    !if (res_no .eq. 3) write(49,*) flow_in_hyp_x+flow_in_epi_x, T_res_inflow(res_no)
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
    depth_e(res_no) = volume_e_x(res_no) / surface_area(res_no)
    depth_h(res_no) = volume_h_x(res_no) / surface_area(res_no)
    
    T_res(res_no) = (T_epil(res_no) * (volume_e_x(res_no)/volume_tot)) + &
        (T_hypo(res_no)*(volume_h_x(res_no)/volume_tot) ) ! weighted averge temp

    if(T_epil(res_no).lt.0.5) T_epil(res_no) = 0.5
    if(T_hypo(res_no).lt.0.5) T_hypo(res_no) = 0.5
    !if(res_no.eq.1) write(*,*) & 
    !    'volume',nsub,volume_e_x(res_no),volume_h_x(res_no),water_withdrawal(res_no),&
    !    depth_e(res_no),depth_h(res_no)
    !
    !if(res_no .eq. 19) write(*,*) nd, T_epil(res_no), T_hypo(res_no),&
    !    volume_e_x(res_no)/surface_area(res_no), &
    !    volume_h_x(res_no)/surface_area(res_no)
end subroutine reservoir_subroutine_implicit
