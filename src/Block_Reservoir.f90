module Block_Reservoir

        implicit none

        !------------------------read in reservoir info--------------------------
        real, dimension (:), allocatable:: dam_lat, dam_lon, res_grid_lat, res_grid_lon
        real, dimension (:), allocatable:: res_depth_feet, res_width_feet, res_length_feet, temp_out
        integer, dimension (:), allocatable:: dam_number, start_operating_year
        integer, dimension (:), allocatable:: res_top_vol, res_bot_vol, res_max_flow, res_min_flow
        integer, dimension (:), allocatable:: res_start_node, res_end_node
       !  integer, dimension (:,:), allocatable:: nodes_x   !for each reach, what are all the nodes
        real , dimension(:), allocatable  :: rmile_node
        integer, dimension(:,:), allocatable  :: res_num, nodes_x, x_dist_res
        logical, dimension(:,:), allocatable  :: res_pres, res_upstream
        logical :: reservoir ! the TRUE or FALSE in fifth line of _Network file whether reserovirs are present 
        integer, dimension (:), allocatable :: xres
        integer :: xres2
        real, dimension(:), allocatable :: dx_res, dt_res
end module Block_Reservoir
