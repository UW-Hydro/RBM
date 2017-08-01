Module Block_Reservoir



        implicit none

        ! --------------------------------------------------------------------
        !       reservoir to RBM variables
        !----------------------------------------------------------------------

        !------------------------read in reservoir
        !info--------------------------
        real, dimension (:), allocatable:: dam_lat, dam_lon, res_grid_lat,res_grid_lon
        real, dimension (:), allocatable:: res_depth_feet, res_width_feet,res_length_feet
        integer, dimension (:), allocatable:: dam_number, start_operating_year
        real, dimension (:), allocatable:: res_top_vol, res_bot_vol,res_max_flow, res_min_flow
        integer, dimension (:), allocatable:: res_start_node, res_end_node
       !  integer, dimension (:,:), allocatable:: nodes_x   !for each reach,
       !  what are all the nodes
        real , dimension(:), allocatable  :: rmile_node
        integer, dimension(:,:), allocatable  :: res_num, nodes_x, x_dist_res
        logical, dimension(:,:), allocatable  :: res_pres, res_upstream
        logical, dimension(:), allocatable :: flag_turnover
        logical :: reservoir, res_upstreamx ! the first is TRUE or FALSE in fifth line of _Network file whether reserovirs are present 
        integer, dimension (:), allocatable :: xres, resx
        integer :: xres2, nres, nm_start, ncell0res, resx2
        real, dimension(:), allocatable :: dx_res, dt_res,reservoir_storage,reservoir_storage_prev
        real ::  advec_tot,q_surf_tot
        real, dimension (:),allocatable::diffusion_tot,advec_hyp_tot,advec_epi_tot, qsurf_tot

        !----------------------------------------------------------------------
        !       variables from reservoir model
        !----------------------------------------------------------------------

        ! ----------------------- constants --------------------------
        real, parameter ::  prcnt_flow_epil = 1, prcnt_flow_hypo = 0, gravity = 9.8
        real, parameter :: density = 1000, heat_c = 4180, ftsec_to_msec =0.0283168
        real, parameter :: J_to_kcal = 0.00023885, kPa_to_mb = 10, ft_to_m =0.3048
        real, parameter :: acrefeet_to_m3 = 1233.48
        real, parameter :: heat_c_kcal = 1  ! heat capacity in kcal/kg*C
        real, dimension (:), allocatable :: K_z    !diffusion coefficient (m^2/sec)

        ! -------------------- temperature and meterological variables ------
        real, dimension (:), allocatable :: temp_change_ep, temp_change_hyp,temp_out, temp_out_i ! energy
        real, dimension (:), allocatable :: T_epil,T_hypo, stream_T_in
        real, dimension (:), allocatable :: density_epil, density_hypo,density_in

        ! -------------------- reservoinr information --------------
        real, dimension (:), allocatable ::  surface_area, depth_total, depth_e,depth_h
        real, dimension (:), allocatable ::  delta_vol_e_x, delta_vol_h_x
        real, dimension (:), allocatable :: delta_vol_e_T_x, delta_vol_h_T_x,dV_dt_epi,dV_dt_hyp
        real, parameter :: depth_e_frac=0.4, depth_h_frac=0.6
        real, dimension (:), allocatable :: Q_tot, Q_pen, Q_spill
        real, dimension (:), allocatable :: depth_e_inital, volume_e_initial,depth_h_inital, volume_h_initial
        real, dimension (:), allocatable :: volume_e_x,volume_h_x, T_res,T_res_in, T_trib_tot, Q_res_in
        logical, dimension (:), allocatable :: res_run, res_start, trib_res !logical to only get start of reservoir and model entire reservoir once each loop
        real :: outflow_x, volume_tot, T_res_in_x, Q_trib_tot_x, T_trib_in_x
        ! -------------------- energy terms -----------
        real, dimension (:), allocatable :: area
        ! real  :: flow_in_hyp_x, flow_in_epi_x, flow_out_epi_x, flow_out_hyp_x
        ! real ::  flow_epi_hyp_x
        real  :: epix, hypox, dif_epi_x, dif_hyp_x, energy_x
        real :: advec_in_epix, advec_out_epix, advec_in_hypx, advec_out_hypx
        real :: advec_epi_hyp,tntrp_x

        ! --------------- path and directories of input and output files -------
        character(len=300 ) :: reservoir_file
        character(len=300 ) :: releases_file

end module Block_Reservoir
