Module Block_Network
!
! Module with stream topology variables
!
    integer, dimension(:), allocatable  ::no_celm,no_cells,no_tribs
    integer, dimension(:), allocatable  ::head_cell
!
    integer, dimension(:,:), allocatable::conflnce,reach_cell,segment_cell,trib
!
!
! Integer variables 
!
    integer:: flow_cells,heat_cells
    integer:: ndays,nreach,ntrb,nwpd
    integer,parameter::ns_max=200
    integer:: start_year,start_month,start_day
    integer:: end_year,end_month,end_day
!
! Real variables
!
    real   :: delta_n,n_default=2
    real   :: dt_comp
    real, dimension(:), allocatable  :: ndelta
!
!      Logical variables
!
    logical :: source  ! presence of thermal plant



end module Block_Network
