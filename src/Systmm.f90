Module SYSTM
!
integer:: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer:: nc,nd,ndd,nm,nr,ns, nresx, ncellx
integer:: nr_trib,ntrb,ntribs
integer:: nrec_flow,nrec_heat
integer:: n1,n2,nnd,nobs,ndays,nyear,nd_year,ntmp
integer:: npart,nseg,nwpd, nd2 
real::    dt_comp,dt_calc,dt_total,hpd,Q1,Q2,q_dot,q_surf,z
real   :: rminsmooth
real   :: T_0,T_dist
real(8):: time
real   :: x,x_bndry,xd,xdd,xd_year,x_head,xwpd,year
real,dimension(:),allocatable:: T_smth  ! ,T_trib, T_head
real,dimension(:,:,:),allocatable:: temp
integer,dimension(:,:,:),allocatable:: res_run
!
!
logical :: leap_year
logical:: DONE
!
! Indices for lagrangian interpolation
!
integer:: npndx,ntrp
! integer, dimension(3):: ndltp=(/-2,-3,-3/)
! integer, dimension(3):: nterp=(/3,4,3/)

!
real, parameter:: pi=3.14159,rfac=304.8

!       reservoir variables 
! real, dimension (:), allocatable ::  T_epil,T_hypo, volume_e_x,volume_h_x,stream_T_in
! real, dimension (:), allocatable ::  density_epil, density_hypo, density_in
!
!
contains
!
SUBROUTINE SYSTMM(temp_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_Reservoir
use Block_Flow
!
Implicit None
!
character (len=200):: temp_file
!
integer :: njb, resx2, i, j
!real,dimension(4):: ta,xa
real :: tntrp
!
! array for each reservoir for if inflow from that parcel
! to reservoir has been calculated
logical, dimension(heat_cells) :: res_inflow
res_inflow = .false.  
!FUNCTION Leap_Year(nyear)
!
!     stream reservoir 
!

!
! Allocate the arrays
!
allocate (temp(nreach,-2:ns_max,2))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_diff(heat_cells))
allocate (Q_trib(nreach))
allocate (width(heat_cells))
allocate (u(heat_cells))
allocate (dt(2*heat_cells))
allocate (dbt(heat_cells))
allocate (ea(heat_cells))
allocate (Q_ns(heat_cells))
allocate (Q_na(heat_cells))
allocate (press(heat_cells))
allocate (wind(heat_cells))
allocate (dt_res(2*heat_cells))
allocate (resx(2*heat_cells))

! 
!   reservoir allocatable
!
allocate (depth_e(nres))
allocate (depth_h(nres))
allocate (surface_area(nres))
allocate (T_epil(nres))
T_epil = 15
allocate (T_hypo(heat_cells))
T_hypo = 10
allocate (stream_T_in(nres))
stream_T_in = 15
allocate (density_epil(nres))
density_epil = 0
allocate (density_hypo(nres))
density_hypo = 0
allocate (density_in(nres))
density_in = 0
allocate (volume_e_x(nres))
allocate (volume_h_x(nres))
allocate (dV_dt_epi(nres))
allocate (dV_dt_hyp(nres))
allocate (K_z(nres))
K_z = 0.1
allocate (temp_change_ep(nres))
temp_change_ep = 0
allocate (temp_change_hyp(nres))
temp_change_hyp = 0
allocate (res_run(nres))
res_run = .false.
allocate (res_start(nres))
res_start = .false.
allocate (T_res(nres))
T_res = 15
allocate (T_res_in(nres))
T_res_in = 15
allocate (Q_trib_tot(heat_cells))
allocate (T_trib_tot(heat_cells))
allocate (Q_res_in(nres))
! allocate(temp_out(4))
temp_out = 10
allocate (trib_res(heat_cells))
![CONSIDER ADDING A SUBROUTINE THAT INTIALIZES VARIABLES LIKE:T_epil, T_hyp, K_z, etc - JRY]
!
! Initialize some arrays
!
dt_part=0.
x_part=0.
no_dt=0
nstrt_elm=0
temp=15
! Initialize headwaters temperatures
!
T_head=15
!!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
T_smth=15


!
!    initialize reservoir geometery variables
!
depth_e = res_depth_feet(:) * depth_e_frac * ft_to_m  ! ft_to_m converst from feet to m
! print *,'res_depth',res_depth_feet(:),'depth_e_frac',depth_e_frac,  'depth_e', depth_e(:)
depth_h = res_depth_feet(:) * depth_h_frac * ft_to_m  ! ft_to_m converst from feet to m
surface_area = res_width_feet(:) *  res_length_feet(:) * ft_to_m * ft_to_m  ! ft_to_m converst from feet to m
volume_e_x = surface_area(:) * depth_e(:) 
volume_h_x = surface_area(:) * depth_h(:) 

! depth_e_inital = depth_e
! volume_e_initial = volume_e_x
! depth_h_inital = depth_h
! volume_h_initial = volume_h_x
!     initial temperature

!
!     open the output file
!

open(20,file=TRIM(temp_file),status='unknown')
!
!
! Initialize dummy counters that facilitate updating simulated values
!
n1=1
n2=2
nobs=0
ndays=0
xwpd=nwpd
hpd=1./xwpd
!
!     Year loop starts
!

do nyear=start_year,end_year
  write(*,*) ' Simulation Year - ',nyear,start_year,end_year
  nd_year=365
  !if (mod(nyear,4).eq.0) nd_year=366
      FUNCTION Leap_Year(nyear)
        leap_year = .FALSE.
      IF (MOD(nyear,4) .EQ. 0)   leap_year = .TRUE.
      IF (MOD(nyear,100) .EQ. 0) leap_year = .FALSE.
      IF (MOD(nyear,400) .EQ. 0) leap_year = .TRUE.
      RETURN
      END FUNCTION Leap_Year


  if(Leap_Year) nd_year = 366
  !
  !     Day loop starts
  !
  do nd=1,nd_year
    year=nyear
    xd=nd
    xd_year=nd_year
    !
    !     Start the numbers of days-to-date counter
    !
    ndays=ndays+1
    !
    !    Daily period loop starts
    !
    do ndd=1,nwpd
      xdd = ndd
      time=year+(xd+(xdd-0.5)*hpd)/xd_year 
      res_run = .false.  ! re-initialize reservoir fun T/F array
      res_start = .false. ! re-initialize T/F for reservoir start
      Q_trib_tot = 0 ! re-set the tributary flow to 0
      T_trib_tot = 0 ! re-set the tributary flow to 0
      trib_res = .false.
      !
      ! Read the hydrologic and meteorologic forcings
      !
      call READ_FORCING
      !
      !     Begin cycling through the reaches
      !
      do nr=1,nreach
        !
        nc_head=segment_cell(nr,1)
        !
        !     Determine smoothing parameters (UW_JRY_2011/06/21)
        !
        rminsmooth=1.0-smooth_param(nr)
        print *,'nr',nr,'T_smth',T_smth(nr)

        T_smth(nr)=rminsmooth*T_smth(nr)+smooth_param(nr)*dbt(nc_head)
        !     
        !     Variable Mohseni parameters (UW_JRY_2011/06/16)
        ! 
        T_head(nr)=mu(nr)+(alphaMu(nr) &
            /(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr)))))  
        !
        temp(nr,0,n1)=T_head(nr)
        temp(nr,-1,n1)=T_head(nr)
        temp(nr,-2,n1)=T_head(nr)
        temp(nr,no_celm(nr)+1,n1)=temp(nr,no_celm(nr),n1)
        x_head=x_dist(nr,0) ! calculated distance to headw based on first x_dist value (the cell furthest upstream)
        x_bndry=x_head-50.0 

        do ns=1,no_celm(nr) ! cycle through all segments in a reach

        ! -----------------------------------------------------------------------
        !
        !                      start RIVER loop  
        !
        ! -----------------------------------------------------------------------


          !     Establish particle tracks

          call Particle_Track(nr,x_head,x_bndry)

          DONE=.FALSE.

          ! if segment in river research, or the start of reservoir
          if(res_pres(nr,segment_cell(nr,ns)) .eqv. .false. ) then !.or. any(segment_cell(nr,ns) == res_start_node(:)) 

            ncell=segment_cell(nr,ns) !cell of parcel for this time step
            nseg=nstrt_elm(ns) !segment water was at previous time step
            npndx=2
            
         !   print *, 'ncell',ncell, 'ns',ns, 'nseg',nseg

            !
            !      subroutine to establish where parcel started
            !

                   if( ns .eq. 2 )    write(92,*) time, segment_cell(nr,ns),ns,  T_0
            call upstream_subroutine(nseg,nr,ns,T_0, npndx, npart, n1, ncell, resx2)        

                   if( ns .eq. 2 )    write(93,*) time, segment_cell(nr,ns),ns,  T_0

            nncell=segment_cell(nr,nstrt_elm(ns)) ! cell of previous time step

            !    set ncell0 for purposes of tributary input
            ncell0=nncell
            dt_total=dt_calc  ! set total time to time parcel took to pass through first segment

            !    loop to set number of segments to cycle through based on if
            !   A) started in reservoir and B) finished downstream of reservoir

            if(reservoir .and. res_upstreamx .and. .not. res_pres(nr,ncell)) then
              nm_start = ns - cell_segment(nr,res_end_node(resx2))
            else
              nm_start = no_dt(ns)
            end if

            do nm=nm_start,1,-1  ! cycle through each segment parcel passed through
              z=depth(nncell)
              nd2 = nd  ! cut out later, just to print day in energy module
              call energy(T_0,q_surf,nncell, ns, nyear, nd2)
              q_dot=(q_surf/(z*rfac))
            !  print *,'nd',nd,'ncell',ncell, 'T_0', T_0
              T_0=T_0+q_dot*dt_calc !adds heat added only during time parcel passed this segment

              if(T_0.lt.0.0) T_0=0.0
              call trib_subroutine(nncell,ncell0, T_0,nr_trib, nr & 
                          ,ns, nseg, n2, DONE, dt_calc, dt_total)

            end do ! end loop cycling through all segments parcel passed through


            if (T_0.lt.0.5) T_0=0.5 ! set so lowest temp can be is 0.5
            temp(nr,ns,n2)=T_0    ! temperature will be used next simulation
            T_trib(nr)=T_0        ! temp of this reach, to calc trib inflow

            ! ---- loop to give downstream reservoir the temperature ---
            do i = 1, nres
              
              if(reservoir .and. ncell + 1 == res_start_node(i))then
                  !this loop is so T_res_in is for segment upstream of first
                  !  node in the reservoir
                  if(segment_cell(ns,nr) .eq. segment_cell(ns-1,nr) ) then 
                    T_res_in(i) = T_0
                  end if
               !    print *,'nres',nres, 'nd',nd,'i',i, 'ncell', ncell, 'T_0',T_0
              end if
            end do

          end if   ! end river if loop

        ! -----------------------------------------------------------------------
        !
        !                     Reservoir Subroutine 
        !
        ! -----------------------------------------------------------------------

          ! if start cell is in reservoir, but not first cell
          if(res_pres(nr,segment_cell(nr,ns))) then

            ncellx = segment_cell(nr,ns) ! cell (node) reservoir

            ! ------ read in tributary flow and temperature -------
              !this if loop adds previous trib flow/temp to current trib flow/temp
            if(trib_res(ncellx) .eqv. .false.)  then
              call trib_res_subroutine(ncellx, nr_trib, nr, ns) 
              trib_res(ncellx) = .true. ! set flag to true, so doesn't run twice
                                        ! for same cell
            end if

            do i = 1, nres
              
              ! -------- if start of reservoir -------------------
              if(reservoir .and. res_start_node(i) .eq. segment_cell(nr,ns) .and. .not. res_start(i) ) then
                nresx = i                
                Q_res_in(nresx) = Q_in(nncell)
               ! T_res_in(nresx) = T_res_in_x ! this will be advection from this reach to reservoir
                Q_trib_tot_x = 0   ! initialize trib flow in this reser. as 0
                T_trib_in_x = 0   ! initialize trib temp in this reser. as 0
                res_start(i) = .true.    ! logical so it won't add another T_res_in
              end if

              ! ------------ end of reservoir - calculate the reservoir temperature -----------  
              if (reservoir .and. res_end_node(i) .eq. segment_cell(nr,ns)   .and. .not. res_run(i) ) then
                nresx = i 
                
                ! ----- add tributary flow to reach inflow to reservoir ---
                do j = res_start_node(nresx), res_end_node(nresx)
                  T_trib_in_x =  (T_trib_in_x*Q_trib_tot_x + T_trib_tot(j)* Q_trib_tot(j)) &
                      /(Q_trib_tot_x + Q_trib_tot(j)) 
                  Q_trib_tot_x = Q_trib_tot_x + Q_trib_tot(j)

                end do

                ! -------- combinei all  trib flow/temp and reach flow/temp----------
                if(Q_trib_tot_x .gt. 0) then
                        T_res_in(nresx) = (T_res_in(nresx)*Q_res_in(nresx) + T_trib_in_x*Q_trib_tot_x) &
                          / (Q_res_in(nresx) +  Q_trib_tot_x)
                        Q_res_in(nresx) = Q_res_in(nresx) + Q_trib_tot_x
                end if

             !   print *, 'T_trib_in_x',  T_trib_in_x, 'Q_trib_tot_x', Q_trib_tot_x
              !  print *,'nresx',nresx,  'Q_res_in',Q_res_in(nresx),'T_res_in', T_res_in(nresx)  
                call stream_density(T_res_in(nresx), density_in(nresx))
                call stream_density(T_epil(nresx), density_epil(nresx))
                call stream_density(T_hypo(nresx), density_hypo(nresx))

                call flow_subroutine(flow_in_epi_x, flow_in_hyp_x, flow_epi_hyp_x &
                  , flow_out_epi_x, flow_out_hyp_x, ratio_sp, ratio_pen &
                  , nresx, dt_comp)

                call energy(T_epil(nresx), q_surf, res_end_node(nresx))
                call reservoir_subroutine (nresx, nd,q_surf, time)

                T_0 = T_res(nresx) !T_res is weighted average temperature

                res_run(i) = .true.  !set reservoir run to "true"


              !------ if segment on reservoir cell but reservoir already been simulated -------
              else if(reservoir .and. res_end_node(i) .eq. segment_cell(nr,ns) .and. res_run(i)) then
                nresx = i
                T_0 = T_res(nresx)

              end if ! end reservoir loop

            end do ! end individual reservoir loop-cycles thru all reservoirs

          end if ! end reservoir if statement

          !   The temperature is output at the beginning of the reach 
          call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_out(ncell))

          !
          !     End of computational element loop
          !
        end do ! end single segment loop

        !
        !     End of reach loop
        !
      end do   ! end total reach loop
      ntmp=n1
      n1=n2
      n2=ntmp
      !
      !     End of weather period loop (NDD=1,NWPD)
      !
    end do   ! end day loop

    temp_out(:) = T_res(:) !set reservoir temperature for next time step
    write(32,*),time, T_epil(1:nres), T_hypo(1:nres) ! , flow_in_epi_x, flow_out_epi_x,
    !
    !     End of main loop (ND=1,365/366)
    !
  end do ! end of year loop
  !
  !     End of year loop
  !
end do  ! end of large loop
!            dt_calc=dt(nncell) ! showed up after John's merge (5/17/16)
!            dt_total=dt_total+dt_calc ! showed up after John's merge (5/17/16)
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM
end module SYSTM
