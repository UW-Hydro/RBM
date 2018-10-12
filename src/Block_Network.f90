Module Block_Network
    !
    ! Module with stream topology variables
    !
    integer, dimension(:), allocatable  ::no_celm,no_cells,ndelta,no_tribs
    integer, dimension(:), allocatable  ::head_cell
    !
    integer, dimension(:,:), allocatable::segment_cell,trib, cell_segment
    !
    !
    ! Integer variables 
    !
    integer:: flow_cells,heat_cells
    integer:: ndays,nreach,ntrb,nwpd
    integer,parameter::ns_max=1000
    integer:: start_year,end_year
    integer:: n_default=2
    !
    ! Real variables
    !
    real:: delta_n,dt_comp
    real,dimension(4):: ta,xa
end module Block_Network
