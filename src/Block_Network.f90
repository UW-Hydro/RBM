Module Block_Network
!
! Module with stream topology variables
!
    integer, dimension(:), allocatable  ::no_celm,no_cells,no_tribs
    integer, dimension(:), allocatable  ::head_cell
!
    integer, dimension(:,:), allocatable::conflnce,reach_cell,segment_cell,trib
!
    integer, dimension(:,:,:), allocatable::nseg_out
!
! Integer variables 
!
    integer:: flow_cells,heat_cells
    integer:: ndays,nreach,ntrb,nwpd
    integer,parameter::ns_max=3000
    integer,parameter::nseg_out_num=2
    integer:: start_year,start_month,start_day
    integer:: end_year,end_month,end_day
    integer:: numsub !number of subdaily timestep
    integer:: nsub
!
! Real variables
!
    real   :: delta_n,n_default=2
    real   :: dt_comp
    real   :: dt_res
    real, dimension(:), allocatable  :: ndelta
!
!      Logical variables
!



end module Block_Network
