SUBROUTINE Error_estimate(nd,res_no)
    use Block_Network
    use Block_Reservoir
    !
    implicit none
    !
    real :: dayx, q_surf, log_K_z, n_stability, density_dif
    real :: temp_epi_pre, temp_hyp_pre
    integer :: nd,res_no

    temp_epi_pre = T_epil(res_no)
    temp_hyp_pre = T_hypo(res_no)
    ! ---------------- turnover loop driven only by T_epil and T_hyp
    ! ----------
    density_dif = density_epil(res_no) - density_hypo(res_no)
    n_stability = (-1) * (gravity/density_epil(res_no)) * &
        ((density_dif)/((depth_e(res_no) + depth_h(res_no))/2))
    if(n_stability .le. 0) then
        n_stability = 1e-10
    end if
    log_K_z = log10(n_stability) * (-1)  - 5.699 ! high scenario - 2  w/
    K_z(res_no) = 10**log_K_z
    ! ONLY for no stratification run:
    write(57, *) res_no,K_z(res_no)
    !
    !     Estimate the numerical error in reservoir temperature 
    !
    if (temp_epi_pre-temp_hyp_pre>0.01) then
        m11 = -(flow_in_epi_x+K_z(res_no)*surface_area(res_no))/volume_e_x(res_no)/dt_res
        m12 = K_z(res_no)*surface_area(res_no)/volume_e_x(res_no)/dt_res
        m21 = (flow_epi_hyp_x+K_z(res_no)*surface_area(res_no))/volume_h_x(res_no)/dt_res
        m22 = -(flow_in_hyp_x+flow_epi_hyp_x+K_z(res_no)*surface_area(res_no))/volume_h_x(res_no)/dt_res
        const1 = ((q_surf*dt_res*surface_area(res_no))/(density*heat_c_kcal) + &
                 flow_in_epi_x*T_res_inflow(res_no))/volume_e_x(res_no)/dt_res
        const2 = flow_in_hyp_x*T_res_inflow(res_no)/volume_h_x(res_no)/dt_res
        error_e = -0.5*((m11**2+m12*m21)*temp_epi_pre + &
                        (m11*m12+m12*m22)*temp_hyp_pre + &
                         m11*const1 + m12*const2)*dt_res**2
        error_h = -0.5*((m21*m11+m22*m21)*temp_epi_pre + &
                        (m21*m12+m22**2)*temp_hyp_pre + &
                         m21*const1 + m22*const2)*dt_res**2
    else
        m11 = -(flow_in_epi_x+flow_in_hyp_x)/(volume_e_x(res_no)+volume_h_x(res_no))/dt_res
        const1 = ((flow_in_epi_x+flow_in_hyp_x)*T_res_inflow(res_no) + &
                 (q_surf*dt_res*surface_area(res_no))/(density*heat_c_kcal))&
                 /(volume_e_x(res_no)+volume_h_x(res_no))/dt_res
        error_e = -0.5*(m11**2*temp_epi_pre+m11*const1)*dt_res**2
        error_h = error_e
    end if
    if (abs(error_e) > error_threshold.or. abs(error_h)>error_threshold) then
        exceed_error_bound=.TRUE.
    end if
end subroutine Error_estimate 
