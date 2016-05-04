Module SYSTM
!
integer:: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer:: nc,nd,ndd,nm,nr,ns, nresx
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
!
! Initialize some arrays
!
dt_part=0.
x_part=0.
no_dt=0
nstrt_elm=0
temp=0.5
! Initialize headwaters temperatures
!
T_head=4.0
!!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
T_smth=4.0


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
  if (mod(nyear,4).eq.0) nd_year=366
  !
  !     Day loop starts
  !
  do nd=1,nd_year
    year=nyear
    xd=nd
    xd_year=nd_year
 !   print *, 'xd', xd
 !       print *, 'nd', nd   
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
      !
      ! Read the hydrologic and meteorologic forcings
      !
      call READ_FORCING
      !
      !     Begin reach computations
      !
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
        !
        !     Establish particle tracks
        !
          call Particle_Track(nr,x_head,x_bndry)
        !
        !
        DONE=.FALSE.

        do ns=1,no_celm(nr) ! cycle through all segments in a reach

          ncell=segment_cell(nr,ns) !cell of parcel for this time step
          nseg=nstrt_elm(ns) !segment water was at previous time step
          npndx=2

          !
          !
          !      --------------   start RIVER loop  ----------------
          !
          !

    !    print *, 'nyear',nyear,'nd',nd,'nr',nr,'nseg',nseg  !  , 'T_0',T_0
    !     print *, 'ncell',ncell,  temp_out(:) ! 'upstream_subroutine','nr',nr,'ncell',ncell, 'nseg',nseg,'T_0',T_0
   
          call upstream_subroutine(nseg,nr,ns,T_0, npndx, npart, n1, ncell, resx2)        

          !
          !     Now do the third-order interpolation to
          !     establish the starting temperature values
          !     for each parcel
          !

   !       if(nseg.eq.1) then   ! if parcel started at headwater 
   !         T_0 = T_head(nr)
   !         res_upstreamx = .false.
   !         resx2 = 0

          ! if parcel started in reservoir
   !       else if (reservoir.and.res_pres(nr,segment_cell(nr,nseg))) then 
   !         T_0 = temp_out(res_num(nr,segment_cell(nr,nseg)))  !  
   !         res_upstreamx = .true.
   !         resx2 = res_num(nr,segment_cell(nr,nseg))

   !       else if (reservoir.and.any(res_pres(nr,segment_cell(nr,ns):segment_cell(nr,nseg)))  ) then

          ! these two lines gets reservoir number in reach
   !         resx(segment_cell(nr,nseg): segment_cell(nr,ns)) = res_num(nr,segment_cell(nr,nseg):segment_cell(nr,ns))
   !         resx2 = maxval(resx(segment_cell(nr,nseg):segment_cell(nr,ns)),1)
   !         T_0 = temp_out(resx2)  ! temperature of reservoir parcel crossed
   !         res_upstreamx = .true.
   !       else 
   !         res_upstreamx = .false.
   !         resx2 = 0

            !
            !     Interpolation at the downstream boundary
            !
   !         if(nseg.eq.no_celm(nr)) npndx=3  ! if segment at previous time step was last segment
            !
   !         do ntrp=1,nterp(npndx)
   !           npart=nseg+ntrp+ndltp(npndx)
   !           xa(ntrp)=x_dist(nr,npart) !river mile at npart segment
   !           ta(ntrp)=temp(nr,npart,n1)  !temperature at npart
   !         end do
   !         x=x_part(ns)
            !
            !     Call the interpolation function
            !
   !         T_0=tntrp(xa,ta,x,nterp(npndx))
   !       end if   ! end previous time step temperature if statements
          !
   !    300 continue
   !    350 continue

          !          dt_calc=dt_part(ns)
          nncell=segment_cell(nr,nstrt_elm(ns)) ! cell of previous time step

          !
          !    Set NCELL0 for purposes of tributary input
          !
          ncell0=nncell
          dt_total=dt_calc  ! set total time to time parcel took to pass through first segment

          !    loop to set number of segments to cycle through based on if
          !   A) started in reservoir and B) finished downstream of reservoir

          if(reservoir .and. res_upstreamx .and. .not. res_pres(nr,ncell)) then
            nm_start = ns - cell_segment(nr,res_end_node(resx2))
         !   print *, 'no_dt(ns)', no_dt(ns), 'nm_start', nm_start 
          else
            nm_start = no_dt(ns)
          end if

          do nm=nm_start,1,-1  ! cycle through each segment parcel passed through
            z=depth(nncell)
            nd2 = nd  ! cut out later, just to print day in energy module
            call energy(T_0,q_surf,nncell, ns, nyear, nd2)
            q_dot=(q_surf/(z*rfac))

 !   if(nr .lt. 14 .and. nseg==13) print *,'T_0',T_0,'q_dot',q_dot,'dt_calc',dt_calc
 !   if(nr .lt. 14 .and.  isnan(T_0)) print *,' -------------- T_0 is NaN 1  ---- ', '  nr',nr,'nseg',nseg
            T_0=T_0+q_dot*dt_calc !adds heat added only during time parcel passed this segment

 !   if(nr .lt. 14 .and.  isnan(T_0)) print *,' -------------- T_0 is NaN 2 ---- ', '  nr',nr,'nseg',nseg
  ! if(ncell .eq. 82 .and. ns .eq. 34  ) write(83,*) nyear,nd &
  !               , depth(nncell), z, nncell

  !do ncell:71,ns:12 for other

  ! if(ncell .eq. 82 .and. ns .eq. 34  )  write(71, *) nyear,nd,nseg, T_0, q_surf &
  !      , q_dot, dt_calc,z

  !      print *,'nwpd',nwpd, 'ns',ns,'nr',nr
 ! if(nncell == 0) print *, 'nr -------------------------------------', nr, 'nseg',nseg,'nncell',nncell, ncell0

            if(T_0.lt.0.0) T_0=0.0

          call trib_subroutine(nncell, T_0,nr_trib, nr & 
                ,ns, nseg, n2, DONE, ncell0)

          if(ncell0.ne.nncell) then
            ncell0=nncell
            DONE=.FALSE.
          end if
          dt_calc=dt(nncell)
          dt_total=dt_total+dt_calc


  ! if(nr .gt. 7 .and. nr .lt. 14 .and. nseg .gt. 10 .and. nseg .lt. 25) print*,'nr',nr,'ns',ns,'nseg',nseg, 'T_0', T_0 &
  !   , 'dt_calc',dt_calc,'dt_ncell',dt(nncell),'nncell',nncell,   'dt_total',dt_total

           !
            !     Look for a tributary.
            !
   !         Q1=Q_in(nncell)  !
   !         ntribs=no_tribs(nncell)
   !         if(ntribs.gt.0.and..not.DONE) then
   !           do ntrb=1,ntribs
   !             nr_trib=trib(nncell,ntrb) ! gives nr (reach) for each trib
   !             if(Q_trib(nr_trib).gt.0.0) then
   !               Q2=Q1+Q_trib(nr_trib)  !Q1 is reach inflow, Q2 is total flow
   !               T_0=(Q1*T_0+Q_trib(nr_trib)*T_trib(nr_trib))/Q2 !adjust temp based on trib temp/flow

                  ! --- calculate tributary temperature and flow (to add to reservoir) ---    
   !               T_trib_tot(nncell) = (T_trib_tot(nncell)*Q_trib_tot(nncell) + Q_trib(nr_trib)*T_trib(nr_trib)) &
   !                  /( Q_trib(nr_trib) + Q_trib_tot(nncell))
   !               Q_trib_tot(nncell) = Q_trib(nr_trib) + Q_trib_tot(nncell) 
   !             end if 
   !             !
   !             Q1=Q_out(nncell) !
   !             !
   !           end do
   !           DONE=.TRUE.
   !         end if

   !         if(ntribs.eq.0.and.Q_diff(nncell).gt.0) then !Q_diff is lateral flow in
   !           Q2=Q1+Q_diff(nncell)
   !           T_dist=T_head(nr)
   !           T_0=(Q1*T_0+Q_diff(nncell)*T_dist)/Q2
   !           Q1=Q2

              ! adding flow to tributatry flow for reservoir modeling
   !           T_trib_tot(nncell) = (T_trib_tot(nncell)*Q_trib_tot(nncell) + Q_diff(nncell)*T_dist) &
   !                 /( Q_trib(nr_trib) + Q_diff(nncell))

   !           Q_trib_tot(nncell) = Q_trib_tot(nncell) + Q_diff(nncell) 
   !         end if

   !  500 continue
   !       nseg=nseg+1  ! so nseg will be next segment down stream
   !       nncell=segment_cell(nr,nseg)  ! node for next node downstream
          !
          !     Reset tributary flag if this is a new celln_write
          !
   !       if(ncell0.ne.nncell) then
   !         ncell0=nncell
   !         DONE=.FALSE.
   !       end if
   !       dt_calc=dt(nncell)
   !       dt_total=dt_total+dt_calc

        end do ! end loop cycling through all segments parcel passed through


        if (T_0.lt.0.5) T_0=0.5
        temp(nr,ns,n2)=T_0
        T_trib(nr)=T_0
  !    if(nr .gt. 7 .and. nr .lt. 14)  print *,'nr',nr,'nseg',nseg, 'T_trib', T_0 
        !
        !             Stream Reservoir Subroutine
        !

        T_res_in_x = T_0  ! saved temperature from previous segment,when it
        ! gets to the start node of a reservoir            

        ! if(any(segment_cell(nr,ns) == res_start_node(:))) print *,'segment_cell_reservoir', segment_cell(nr,ns)
        do i = 1, nres
          if(reservoir .and. res_start_node(i) .eq. segment_cell(nr,ns) .and. .not. res_start(i) ) then
            nresx = i                
            Q_res_in(nresx) = Q_in(nncell)
            T_res_in(nresx) = T_res_in_x ! this will be advection from this reach to reservoir
            res_start(i) = .true.    ! logical so it won't add another T_res_in
            Q_trib_tot_x = 0   ! initialize trib flow in this reser. as 0
            T_trib_in_x = 0   ! initialize trib temp in this reser. as 0
                
            res_start(i) = .true.
            ! loop to calculate reservoir temperature  
          else if (reservoir .and. res_end_node(i) .eq. segment_cell(nr,ns)   .and. .not. res_run(i) ) then
            nresx = i 
                
            ! ----- add tributary flow to reach inflow to reservoir ---
            do j = res_start_node(nresx), res_end_node(nresx)
              T_trib_in_x =  (T_trib_in_x*Q_trib_tot_x + T_trib_tot(j)* Q_trib_tot(j)) &
                  /(Q_trib_tot_x + Q_trib_tot(j)) 
              Q_trib_tot_x = Q_trib_tot_x + Q_trib_tot(j)

           !  print *,'nresx',nresx,'j',j, 'T_trib_in' ,T_trib_in_x,'T_trib_tot',T_trib_tot(j)
            end do
            !    print *, 'Q-in', Q_res_in(nresx), 'Q-trib', Q_trib_tot_x
            ! --- combine trib flow/temp and reach flow/temp ---
            T_res_in(nresx) = (T_res_in(nresx)*Q_res_in(nresx) + T_trib_in_x*Q_trib_tot_x) &
                    / (Q_res_in(nresx) +  Q_trib_tot_x)
            Q_res_in(nresx) = Q_res_in(nresx) + Q_trib_tot_x
            call stream_density(nresx)

            call flow_subroutine(flow_in_epi_x, flow_in_hyp_x, flow_epi_hyp_x&
                , flow_out_epi_x, flow_out_hyp_x, ratio_sp, ratio_pen &
                , nresx, dt_comp)

            call energy(T_epil(nresx), q_surf, res_end_node(nresx))

            call reservoir_subroutine (nresx, nd,q_surf, time)

            T_0 = T_res(nresx) !T_res is weighted average temperature

            res_run(i) = .true.  !set reservoir run to "true"
   ! print *,'first segment_cell_reservoir', segment_cell(nr,ns) , 'T_0', T_0            
             ! if segment on reservoir cell but reservoir already been simulated 
          else if(reservoir .and. res_end_node(i) .eq. segment_cell(nr,ns) .and. res_run(i)) then
            nresx = i
            T_0 = T_res(nresx)
   ! print *,'later segment_cell_reservoir', segment_cell(nr,ns), 'T_0', T_0  
          end if
        end do
        !
        !   Write all temperature output UW_JRY_11/08/2013
        !   The temperature is output at the beginning of the 
        !   reach.  It is, of course, possible to get output at
        !   other points by some additional code that keys on the
        !   value of ndelta (now a vector)(UW_JRY_11/08/2013)
        !
        call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_out(ncell))
        !
        !     End of computational element loop
        !
   ! print *,'day',nnd, 'T_epil', T_epil

      end do ! end single reach loop (goes through all segs in reach)
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
      4650 format(16x,12(6x,f6.0,6x))
      4700 format(f10.4,f6.0,15(f6.1,f8.3))
      4750 format(f10.4,10(i4,f8.0))
      
      temp_out(:) = T_res(:) !set reservoir temperature for next time step

    end do

     if(nresx==4 .and. nd==180 ) print *,'-----------------' , nd,'T_epil',  T_epil(nresx), 'T_hypo',T_hypo(nresx)

    write(32,*),time, T_epil(1:nres), T_hypo(1:nres) ! , flow_in_epi_x, flow_out_epi_x,
  !     print *, 'time', time, T_epil(1:nres), T_hypo(1:nres)

    !
    !     End of main loop (ND=1,365/366)
    !
  end do ! end of year loop
  !
  !     End of year loop
  !
end do  ! end of large loop
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM
end module SYSTM
