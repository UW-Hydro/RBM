!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
    integer, dimension(2000):: no_dt,nstrt_elm
    integer :: nsource,  sourcex, sourcex2  ! number of thermal point-source inputs
    real, dimension(2000)   :: dt_part,x_part
    !
    real, dimension(:),   allocatable  :: depth
    real, dimension(:),   allocatable  :: width
    real, dimension(:),   allocatable  :: u
    real, dimension(:),   allocatable  :: dt
    real, dimension(:),   allocatable  :: dx
    real, dimension(:),   allocatable  :: Q_in
    real, dimension(:),   allocatable  :: Q_trib
    real, dimension(:),   allocatable  :: Q_out
    real, dimension(:),   allocatable  :: Q_diff
    real, dimension(:,:), allocatable  :: Q_nps
    real, dimension(:,:), allocatable  :: temp_trib
    real, dimension(:,:), allocatable  :: temp_nps,thermal
    real, dimension(:,:), allocatable  :: x_dist
    real, dimension(:,:,:), allocatable :: temp
    real, dimension(:), allocatable :: Q_trib_tot
    real, dimension(:), allocatable :: T_trib, T_head
    real, dimension(:), allocatable :: flow_source, source_num_cell

    ! Leopold coefficients (constant throughout basin)
    ! <a_d> and <b_d>: a and b coefficients in flow depth estimation: z=aQ^b, where Q is flow discharge [cfs], z is flow depth [ft]
    ! <a_w> and <b_w>: a and b coefficients in channel width estimation: w=aQ^b, where Q is flow discharge [cfs], w is channel width [ft]
    real, parameter :: a_d  =  0.408, b_d= 0.392, a_w = 4.346, b_w = 0.520

    integer, dimension(:), allocatable :: source_cell
    logical, dimension(:), allocatable :: source_cell_tf
    logical :: source
    character(len=300 ) :: source_list
end module Block_Hydro
