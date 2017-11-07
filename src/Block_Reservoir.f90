!
! Module for reservoir characteristic and parameters
!
module Block_Reservoir
    !
    ! Logical
    !
    logical :: reservoir ! whether there is reservoir in the system
    logical, dimension(:),   allocatable  :: res_pres ! whether reservoir exists at specific grid cell
    logical, dimension(:),   allocatable  :: res_start !
    logical, dimension(:),   allocatable  :: initial_storage !
    logical, dimension(:),   allocatable  :: res_trib_calc !
    logical, dimension(:),   allocatable  :: res_stratif_start !
    logical, dimension(:),   allocatable  :: res_turnover !
    !
    !
    real, parameter :: depth_e_frac=0.4, depth_h_frac=0.6
    real, parameter :: ftsec_to_msec = 0.0283168, gravity = 9.8
    real, parameter :: density = 1000, heat_c = 4180
    real, parameter :: heat_c_kcal = 1
    real :: flow_in_epi_x, flow_in_hyp_x
    real :: flow_out_epi_x,flow_out_hyp_x
    real :: flow_epi_hyp_x, outflow_x
    real :: res_vol_delta_x, vol_change_hyp_x, vol_change_epi_x
    real :: res_storage_pre, res_storage_post
    real, dimension(:),   allocatable  :: K_z
    real, dimension(:),   allocatable  :: depth_e, depth_h
    real, dimension(:),   allocatable  :: density_epil, density_hypo
    real, dimension(:),   allocatable  :: density_in
    real, dimension(:),   allocatable  :: volume_e_x
    real, dimension(:),   allocatable  :: volume_h_x
    real, dimension(:),   allocatable  :: dV_dt_epi
    real, dimension(:),   allocatable  :: dV_dt_hyp
    real, dimension(:),   allocatable  :: dam_grid_lat
    real, dimension(:),   allocatable  :: dam_grid_lon
    real, dimension(:),   allocatable  :: res_depth_meter
    real, dimension(:),   allocatable  :: res_width_meter
    real, dimension(:),   allocatable  :: res_area_km2
    real, dimension(:),   allocatable  :: surface_area
    real, dimension(:),   allocatable  :: Q_res_inflow
    real, dimension(:),   allocatable  :: Q_res_outflow
    real, dimension(:),   allocatable  :: T_res_inflow
    real, dimension(:),   allocatable  :: T_res_in
    real, dimension(:),   allocatable  :: Q_res_in
    real, dimension(:),   allocatable  :: water_withdrawal
    real, dimension(:),   allocatable  :: temp_change_ep, temp_change_hyp
    real, dimension(:),   allocatable  :: T_epil
    real, dimension(:),   allocatable  :: T_hypo
    real, dimension(:),   allocatable  :: temp_out
    real, dimension(:),   allocatable  :: T_res
    real, dimension(:),   allocatable  :: diffusion_tot
    real, dimension(:),   allocatable  :: advec_hyp_tot
    real, dimension(:),   allocatable  :: advec_epi_tot
    real, dimension(:),   allocatable  :: qsurf_tot
    real, dimension(:),   allocatable  :: res_capacity_mcm
    real, dimension(:),   allocatable  :: volume_h_min
    real, dimension(:),   allocatable  :: volume_e_min
    real, dimension(:,:),   allocatable  :: res_storage
    !
    !
    integer :: nres
    integer, parameter :: km_to_m = 1000
    integer, dimension(:),   allocatable  :: res_num
    integer, dimension(:),   allocatable  :: dam_number
    integer, dimension(:),   allocatable  :: res_start_node
    integer, dimension(:),   allocatable  :: res_end_node
    integer, dimension(:),   allocatable  :: ns_res_start
    integer, dimension(:),   allocatable  :: ns_res_end


end module Block_Reservoir
