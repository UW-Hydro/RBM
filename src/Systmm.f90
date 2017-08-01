SUBROUTINE SYSTMM(temp_file,param_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_Reservoir
!
Implicit None
! 
!
character (len=200):: temp_file
character (len=200):: param_file
! 
integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer          :: nc,nd,ndd,nm,nr,ns, i
integer          :: nr_trib,ntribs
integer          :: nrec_flow,nrec_heat
integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
integer          :: npart,nseg,nx_s,nx_part,nx_head
!
! Indices for lagrangian interpolation
!
integer              :: njb,npndx,ntrp
integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)

!
real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z
real             :: Q_dstrb,Q_inflow,Q_outflow,Q_ratio,Q_trb,Q_trb_sum
real             :: T_dstrb,T_dstrb_load,T_trb_load
real             :: rminsmooth
real             :: T_0,T_dist
real(8)          :: time
real             :: x,xd,xdd,xd_year,xwpd,year
real             :: tntrp
real             :: dt_ttotal
real,dimension(4):: ta,xa
!
real,dimension(:),allocatable     :: T_smth

logical:: DONE
!
!
! Allocate the arrays
!
allocate (temp(nreach,0:ns_max,2))
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
T_epil = 10
allocate (T_hypo(nres))
T_hypo = 10
allocate (stream_T_in(nres))
stream_T_in = 10
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
T_res = 10
allocate (T_res_in(nres))
T_res_in = 10
allocate (Q_trib_tot(heat_cells))
allocate (T_trib_tot(heat_cells))
allocate (Q_res_in(nres))
allocate(temp_out(nres))
temp_out = 10
allocate(temp_out_i(nres))
temp_out_i = 10
allocate (trib_res(heat_cells))
allocate(diffusion_tot(nres))
allocate(advec_hyp_tot(nres))
allocate(advec_epi_tot(nres))
allocate(qsurf_tot(nres))
allocate(reservoir_storage(nres))
allocate(reservoir_storage_prev(nres))


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
!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
T_smth=4.0

!
!    initialize reservoir geometery variables
!
depth_e = res_depth_feet(:) * depth_e_frac * ft_to_m  ! ft_to_m converst from feet to m
depth_h = res_depth_feet(:) * depth_h_frac * ft_to_m  ! ft_to_m converst from feet to m
surface_area = res_width_feet(:) *  res_length_feet(:) * ft_to_m * ft_to_m  !ft_to_m converst from feet to m
volume_e_x = surface_area(:) * depth_e(:)
volume_h_x = surface_area(:) * depth_h(:)


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
  DO nd=1,nd_year
    year=nyear
    xd=nd
    xd_year=nd_year
!     Start the numbers of days-to-date counter
!
    ndays=ndays+1
!
!    Daily period loop starts
!
      DO ndd=1,nwpd
      xdd = ndd
      time=year+(xd+(xdd-0.5)*hpd)/xd_year 
      !   ------------ reservoir arrays -----------
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
      temp(nr,1,n1)=T_head(nr)
!
! Begin cell computational loop
!
        do ns=1,no_celm(nr)
! 
        ! -----------------------------------------------------------------------
        !
        !                 River loop  
        !
        ! -----------------------------------------------------------------------


          ! if segment in river reach, or the first segment of reservoir
         if((res_pres(nr,segment_cell(nr,ns)) .eqv. .false.)  .or. &
            (any(segment_cell(nr,ns) == res_start_node(:)) .eqv. .false. .and. &
             res_pres(nr,segment_cell(nr,ns-1)) .eqv. .false.) ) then

            DONE=.FALSE.
  
           ! Testing new code 8/8/2016
           !
           !     Establish particle tracks
           !
           call Particle_Track(nr,ns,nx_s,nx_head)
           !
           ncell=segment_cell(nr,ns)
           !
           !     Now do the third-order interpolation to
           !     establish the starting temperature values
           !     for each parcel
           !
           nseg=nstrt_elm(ns)
           !
           !     Perform polynomial interpolation
           !
           !
           !     Interpolation inside the domain
           !
           npndx=2
           !
           !    If parcel started in or above reservoir
           !

           ! --- if parcel started upstream of reservior and finished downstream -----
           if (reservoir.and.any(res_pres(nr,segment_cell(nr,nseg):segment_cell(nr,ns)) ) &
            .and. .not. res_pres(nr,segment_cell(nr,ns)) & 
                .and. .not.res_pres(nr,segment_cell(nr,nseg))) then

              do i = nseg, ns, 1
                if(res_pres(nr,segment_cell(nr,i))) then
                 T_0 = temp_out_i(res_num(nr,segment_cell(nr,i)))  !  
                 res_upstreamx = .true.
                 resx2 = res_num(nr,segment_cell(nr,i))
                 ncell0res = res_end_node(resx2)
                end if
               end do

           ! ------- if parcel started in reservoir and finished downstream -----------
           else if (reservoir.and.res_pres(nr,segment_cell(nr,nseg))) then
             T_0 = temp_out_i(res_num(nr,segment_cell(nr,nseg)))  !  
             res_upstreamx = .true.
             resx2 = res_num(nr,segment_cell(nr,nseg))
             ncell0res = res_end_node(resx2)

           ! ---- if parcel started in river or headwater, and did not cross reservoir ---
           else
            resx2 = 0
            res_upstreamx = .false.

           end if
           !
           !     Interpolation at the upstream boundary if the
           !     parcel has reached that boundary
           !
           if(nx_head.eq.0) then
             T_0 = T_head(nr)
           else 
           !
           !
           !     Interpolation at the upstream or downstream boundary
           !
           if(nseg .eq. 1 .or. nseg .eq. no_celm(nr)) npndx=1
           !
           do ntrp=nterp(npndx),1,-1
              npart=nseg+ntrp+ndltp(npndx)
              xa(ntrp)=x_dist(nr,npart)
              ta(ntrp)=temp(nr,npart,n1)
           end do
!
! Start the cell counter for nx_s
!
            x=x_part(nx_s)
!
!     Call the interpolation function
!
            T_0=tntrp(xa,ta,x,nterp(npndx))
          end if
if(ncell .eq. 19)  write(68,*) nyear,',', ndays,',',nd,',',ncell,',',ns &
     ,',',T_0,',',nx_head,',',T_head(nr), ',',npndx, ',' ,nterp(npndx)



!
!
          nncell=segment_cell(nr,nstrt_elm(ns))
!
!    Initialize inflow
!
          Q_inflow = Q_in(nncell)
          Q_outflow = Q_out(nncell)
!
!    Set NCELL0 for purposes of tributary input
!
          ncell0=nncell
          dt_total=0.0
          do nm=no_dt(ns),1,-1
            dt_calc=dt_part(nm)
            z=depth(nncell)
            call energy(T_0,q_surf,nncell)
!
            q_dot=(q_surf/(z*rfac))
            T_0=T_0+q_dot*dt_calc
            if(T_0.lt.0.0) T_0=0.0
!
!    Add distributed flows
!    
            T_dstrb_load  = 0.0
!
            Q_dstrb = Q_diff(nncell)
!
! Temperature of distributed inflow assumed = 10.0 deg C
!
            if(Q_dstrb.gt.0.001) then
              T_dstrb  = 10.0
            else
              T_dstrb  = 10.0
            end if
              T_dstrb_load  = Q_dstrb*T_dstrb

!
!     Look for a tributary.
!
            ntribs=no_tribs(nncell)
            Q_trb_sum   = 0.0
            T_trb_load  = 0.0
            if(ntribs.gt.0.and..not.DONE) then
!
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(Q_trib(nr_trib).gt.0.0) then
                  Q_trb        = Q_trib(nr_trib)
                  Q_trb_sum    = Q_trb_sum + Q_trb
!
!  Update water temperature with tributary input
!
                  T_trb_load   = (Q_trb*T_trib(nr_trib))       &
                               +  T_trb_load
                end if
              end do
!
              DONE=.TRUE.
            end if
!
!  Update inflow and outflow
!
            Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum
            Q_ratio = Q_inflow/Q_outflow       
!
! Do the mass/energy balance
!
            T_0  = T_0*Q_ratio                              &
                 + (T_dstrb_load + T_trb_load)/Q_outflow    &
                 + q_dot*dt_calc              
!
            if (T_0.lt.0.5) T_0 =0.5
            Q_inflow = Q_outflow
!
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new cell
!
            if(ncell0.ne.nncell) then
              ncell0=nncell
              Q_inflow = Q_in(nncell)
               DONE=.FALSE.
            end if
            dt_total=dt_total+dt_calc
          end do

        end if ! end river loop

        ! -----------------------------------------------------------------------
        !
        !                     Reservoir Loop
        !
        ! -----------------------------------------------------------------------

        ! -------------  if cell is in reservoir ----------------
        if(reservoir .and. res_pres(nr,segment_cell(nr,ns))) then



        end if  ! end reservoir loop


        if (T_0.lt.0.5) T_0=0.5
            temp(nr,ns,n2)=T_0
            T_trib(nr)=T_0
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
            call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_inflow,Q_outflow)


if(ncell .eq. 82)  write(67,*) nyear,',', ndays,',',nd,',',ncell,',',ns,',',T_0

!
!     End of computational element loop
!
              end do
!     End of reach loop
!
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
!
!     End of weather period loop (NDD=1,NWPD)
!
4650 format(16x,12(6x,f6.0,6x))
4700 format(f10.4,f6.0,15(f6.1,f8.3))
4750 format(f10.4,10(i4,f8.0))
          end do
!
!     End of main loop (ND=1,365/366)
!
        end do
!
!     End of year loop
!
      end do
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM
