SUBROUTINE reservoir_single_subroutine(res_no,q_surf,nd)
    use Block_Network
    use Block_Reservoir
    !
    implicit none
    !
    real :: q_surf, volume_tot
    real :: temp_epi, temp_hyp
    integer :: nd, res_no

temp_epi = (volume_h_x(res_no)*T_hypo(res_no) + volume_e_x(res_no)*T_epil(res_no) + &
           (flow_in_epi_x + flow_in_hyp_x) * T_res_inflow(res_no) + &
           (q_surf * dt_comp * surface_area(res_no)) / (density * heat_c_kcal)) / &
           (volume_h_x(res_no) + volume_e_x(res_no) + flow_in_epi_x + flow_in_hyp_x)
temp_hyp = temp_epi

T_epil(res_no) = temp_epi
T_hypo(res_no) = temp_hyp

volume_tot = volume_e_x(res_no)  + volume_h_x(res_no) - water_withdrawal(res_no)

volume_e_x(res_no) = volume_tot / 2.0
volume_h_x(res_no) = volume_tot / 2.0

end subroutine reservoir_single_subroutine

