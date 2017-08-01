Subroutine BEGIN(param_file,spatial_file)
!
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_Reservoir
!
implicit none
!    
    character (len=8) :: end_date,start_date     
    character (len=8) :: lat
    character (len=10):: long
    character (len=200):: param_file,source_file,spatial_file
    character :: dam_name
!
    integer:: Julian
    integer:: head_name,trib_cell
    integer:: jul_start,main_stem,nyear1,nyear2,nc,ncell,nseg
    integer:: ns_max_test,node,ncol,nrow,nr,cum_sgmnt
    integer:: nreservoir, node_res
!
    logical:: first_cell, res_presx
!
    real :: nndlta
    real :: rmile0,rmile1,xwpd
!
    real,parameter   :: miles_to_ft=5280.
!
!   Mohseni parameters, if used
!
!
!
!     Card Group I
!
read(90,*) start_date,end_date
read(start_date,'(i4,2i2)') start_year,start_month,start_day
read(end_date,  '(i4,2i2)') end_year,end_month,end_day
nyear1=start_year
nyear2=end_year
write(*,'(2(2x,i4,2i2))')  &
 start_year,start_month,start_day,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
!
jul_start = Julian(start_year,start_month,start_day)
!
read(90,*) nreach, flow_cells, heat_cells, source, nsource, reservoir,  nres
!
! Allocate dynamic arrays
!
 allocate(ndelta(heat_cells))
 allocate(mu(nreach))
 allocate(alphamu(nreach))
 allocate(beta(nreach))
 allocate(gmma(nreach))
 allocate (smooth_param(nreach))
 allocate(dx(heat_cells))
 allocate(no_celm(nreach))
 no_celm=0
 allocate(no_cells(nreach))
 no_cells=0
 allocate(no_tribs(heat_cells))
 no_tribs=0
 allocate(trib(heat_cells,10))
 trib=0
 allocate (conflnce(heat_cells,10))
 conflnce=0
 allocate(reach_cell(nreach,ns_max)) 
 allocate(head_cell(nreach))
 allocate(segment_cell(nreach,ns_max))
 allocate(x_dist(nreach,0:ns_max))
! Allocate reservoir info
 allocate(dam_lat(nres))
 allocate(dam_lon(nres))
 allocate(res_grid_lat(nres))
 allocate(res_grid_lon(nres))
 allocate(res_depth_feet(nres))
 allocate(res_width_feet(nres))
 allocate(res_length_feet(nres))
 allocate(dam_number(nres))
 allocate(start_operating_year(nres))
 allocate(res_top_vol(nres))
 allocate(res_bot_vol(nres))
 allocate(res_max_flow(nres))
 allocate(res_min_flow(nres))
 allocate(res_start_node(nres))
 allocate(res_end_node(nres))
 allocate(nodes_x(nreach,0:ns_max))
 allocate(res_num(nreach,0:ns_max))
 allocate(res_pres(nreach,0:ns_max))
 allocate(res_upstream(nreach,0:ns_max))
 allocate(rmile_node(0:ns_max))
 allocate(xres(0:ns_max))
 allocate(dx_res(0:heat_cells))


!
!--------------- read in reservoir info (if reservoirs present) --------------
!
if(reservoir) then
    read(37,*)
    do nreservoir = 1,nres
        read(37,*) dam_number(nreservoir)  &
              , res_grid_lat(nreservoir), res_grid_lon(nreservoir) &
              , res_top_vol(nreservoir) &
              , res_bot_vol(nreservoir),res_depth_feet(nreservoir) &
              , res_width_feet(nreservoir), res_length_feet(nreservoir) &
              ,res_start_node(nreservoir), res_end_node(nreservoir) &
              , res_min_flow(nreservoir)

       print *, 'nreservoir', nreservoir, 'dam_num', dam_number(nreservoir) &
    , 'end', res_end_node(nreservoir) , 'start', res_start_node(nreservoir)
    end do
end if

!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
ncell=0
!
ns_max_test=-1
!
!     Card Group IIb. Reach characteristics
!
do nr=1,nreach
!
!     Initialize NSEG, the total number of segments in this reach
!
  nseg=0
  write(*,*) ' Starting to read reach ',nr
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
  read(90,'(i5,11x,i4,10x,i5,15x,i5,15x,f10.0,i5)') no_cells(nr) &
      ,head_name,trib_cell,main_stem,rmile0
!
! Set boundary distance
!
      x_dist(nr,0)=miles_to_ft*rmile0
!
!     If this is reach that is tributary to cell TRIB_CELL, give it the
!     pointer TRIB(TRIB_CELL) the index of this reach for further use.
!     Also keep track of the total number of tributaries for this cell
!
  if (trib_cell.gt.0) then
    no_tribs(trib_cell) = no_tribs(trib_cell)+1
    trib(trib_cell,no_tribs(trib_cell)) = nr
  end if
!
!     Reading Mohseni parameters for each headwaters (UW_JRY_2011/06/18)
!
  read(90,*) alphaMu(nr),beta(nr) &
            ,gmma(nr),mu(nr),smooth_param(nr)
!
!     Reading Reach Element information
!
  first_cell=.true.
  do nc=1,no_cells(nr)
    ncell=ncell+1
!
!   Relate the cell number from the linear list to the local cell number
!
    reach_cell(nr,nc)=ncell
!
!   Read the data for point sources
!
    if (source) then
!
!  Place holder for point source input
!
    end if 
!
!     The headwaters index for each cell in this reach is given
!     in the order the cells are read
!
!     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
!
!     Variable ndelta read in here.  At present, number of elements
!     is entered manually into the network file (UW_JRY_2011/03/15)
!

        if(reservoir) then
           read(90,'(5x,i5,5x,i5,8x,i5,6x,a8,6x,a10,7x,f10.0,f5.0,i6)')  &  !nodes for each reach
              node,nrow,ncol,lat,long,rmile1,ndelta(ncell),node_res
           if(node_res .gt.0) then
             res_presx = .true.
           else
             res_presx = .false.
           end if

        !      X,Y matrix with nodes in each reach
            nodes_x(nr,ncell) = node    ! matrix of nodes for each cell
            res_num(nr,ncell) = node_res   ! matrix of reservoir index values
            res_pres(nr,ncell) = res_presx  ! matrix for if reservoir is present for each cell
            rmile_node(ncell) = rmile1

        else    ! IF no reservoir information in _Network file:
           read(90,'(5x,i5,5x,i5,8x,i5,6x,a8,6x,a10,7x,f10.0,f5.0)')  &
              node,nrow,ncol,lat,long,rmile1,ndelta(ncell)
        end if

!
!    Set the number of segments of the default, if not specified
!
    if (ndelta(ncell).lt.1) ndelta(ncell)=n_default
    if(first_cell) then
      first_cell=.false.
      head_cell(nr)=ncell
      x_dist(nr,0)=miles_to_ft*rmile0
    end if
!
! Added variable ndelta (UW_JRY_2011/03/15)
!
    dx(ncell)=miles_to_ft*(rmile0-rmile1)/ndelta(ncell)
    rmile0=rmile1
    nndlta=0
200 continue
    nndlta=nndlta+1
    nseg=nseg+1
    segment_cell(nr,nseg)=ncell
    x_dist(nr,nseg)=x_dist(nr,nseg-1)-dx(ncell)


!
!  calcualte distance to next upstream reservoir
!
    if(reservoir) then
        if(res_presx) then  !if cell in reservoir
            res_upstream(nr,ncell) = .false.

        else  if(any(res_pres(nr,:)) ) then   ! any reservoirs upstream of this cell in this reach
            xres(1:size(res_end_node)) = ncell - res_end_node
            res_upstream(nr,ncell) = .true.

        end if

    end if

!
!   Write Segment List for mapping to temperature output (UW_JRY_2008/11/19)
!
    open(22,file=TRIM(spatial_file),status='unknown') ! (changed by WUR_WF_MvV_2011/01/05)
    write(22,'(4i6,1x,a8,1x,a10,f5.0)') nr,ncell,nrow,ncol,lat,long,nndlta
!
! 
!
!  Added variable ndelta  (UW_JRY_2011/03/15)
!
    if(nndlta.lt.ndelta(ncell)) go to 200  
    no_celm(nr)=nseg
    segment_cell(nr,nseg)=ncell
    x_dist(nr,nseg)=miles_to_ft*rmile1
!
! End of cell and segment loop
!
  end do
!
! If this is a reach that is tributary to another, set the confluence cell to the previous 
! cell. This is necessary because the last cell in the reach has the same cell number
! as that of the cell it enters. This is to account for the half portion of the cell the
! the parcel traverses to the center of the grid.
!
if (trib_cell .gt. 0) then
    conflnce(trib_cell,no_tribs(trib_cell)) = ncell
end if

if(ns_max_test.lt.nseg) ns_max_test=nseg
!
! End of reach loop
!
end do
if(ns_max_test.gt.ns_max) then
  write(*,*) 'RBM is terminating because'
  write(*,*) 'NS_MAX exceeded. Change NS_MAX in Block_Network to: ',ns_max_test
  stop
end if
!
nwpd=1
xwpd=nwpd
dt_comp=86400./xwpd
!
!     ******************************************************
!                         Return to RMAIN
!     ******************************************************
!
900 continue
!
!
end subroutine BEGIN
