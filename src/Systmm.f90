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
    integer          :: nc,nd,ndd,nm,nr,ns
    integer          :: nr_trib,ntribs
    integer          :: nrec_flow,nrec_heat
    integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
    integer          :: npart,nseg,nx_s,nx_part,nx_head
    integer          :: ns_res_num, res_no, nsegment
    !
    ! Indices for lagrangian interpolation
    !
    integer              :: njb,npndx,ntrp
    integer, dimension(2):: ndltp=(/-1,-2/)
    integer, dimension(2):: nterp=(/2,3/)
    !
    real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z
    real             :: Q_dstrb,Q_inflow,Q_outflow,Q_ratio,Q_inflow_origin
    real             :: Q_inflow_out, Q_outflow_out, Q_tot
    real             :: sto_pre, sto_post
    real             :: Q_trb,Q_trb_sum,Q_sto,Q_sto_out, T_sto
    real             :: T_dstrb,T_dstrb_load,T_trb_load,T_sto_load
    real             :: Q_res_dstrb_tot, Q_res_trib_tot
    real             :: T_res_dstrb, T_res_dstrb_load_tot, T_res_trib_load_tot
    real             :: rminsmooth
    real             :: T_0,T_dist
    real(8)          :: time
    real             :: x,xd,xdd,xd_year,xwpd,year
    real             :: tntrp
    real             :: dt_ttotal
    real,dimension(4):: ta,xa
    !
    real,dimension(:),allocatable     :: T_head,T_smth,T_trib

    logical:: DONE
    logical:: ns_res_pres
    !
    !
    ! Allocate the arrays
    !
    allocate (temp(nreach,0:ns_max,2))
    allocate (sto(nreach,ns_max,2))
    allocate (temp_sto(nreach,ns_max,2))
    allocate (T_head(nreach))
    allocate (T_smth(nreach))
    allocate (T_trib(nreach))
    allocate (depth(heat_cells))
    allocate (Q_in(heat_cells))
    allocate (Q_out(heat_cells))
    allocate (Q_local(heat_cells))
    allocate (Q_diff(heat_cells))
    allocate (delta_sto_flux(heat_cells))
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
    !
    allocate (res_trib_calc(heat_cells))
    allocate (initial_storage(nres))
    initial_storage=.TRUE.
    allocate (res_start(nres))
    allocate (res_storage(nres,2))
    res_storage=0
    allocate (T_epil(nres))
    T_epil = 4
    allocate (T_hypo(nres))
    T_hypo = 4
    allocate (T_res_in(nres))
    allocate (Q_res_in(nres))
    allocate (T_res_inflow(nres))
    T_res_inflow = 4
    allocate (density_epil(nres))
    density_epil = 0
    allocate (density_hypo(nres))
    density_hypo = 0
    allocate (density_in(nres))
    density_in = 0
    allocate (water_withdrawal(nres))
    allocate (temp_change_ep(nres))
    allocate (temp_change_hyp(nres))
    allocate (temp_out(nres))
    allocate (T_res(nres))
    allocate (diffusion_tot(nres))
    allocate (advec_hyp_tot(nres))
    allocate (advec_epi_tot(nres))
    allocate (qsurf_tot(nres))
    allocate (volume_e_x(nres))
    allocate (volume_h_x(nres))
    allocate (surface_area(nres))
    allocate (K_z(nres))
    allocate (ns_res_start(nres))
    allocate (ns_res_end(nres))
    allocate (Q_res_outflow(nres))
    allocate (Q_res_inflow(nres))
    allocate (dV_dt_epi(nres))
    allocate (dV_dt_hyp(nres))
    allocate (depth_e(nres))
    allocate (depth_h(nres))
    allocate (volume_h_min(nres))
    !
    !
    ! Initialize some arrays
    !
    dt_part=0.
    x_part=0.
    no_dt=0
    nstrt_elm=0
    temp=0.5
    sto=0
    temp_sto=0.5
    ! Initialize headwaters temperatures
    !
    T_head=4.0
    !
    !
    ! Initialize smoothed air temperatures for estimating headwaters temperatures
    !
    T_smth=4.0
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
    !     Reservoir geomorphic information initialization
    !
    depth_e = res_depth_meter(:) * depth_e_frac
    depth_h = res_depth_meter(:) * depth_h_frac
    surface_area = res_area_km2(:) * km_to_m * km_to_m
    volume_e_x = surface_area(:) * depth_e(:)
    volume_h_x = surface_area(:) * depth_h(:)
    print *, 'depth_epil', depth_e
    print *, 'depth_hypo', depth_h
    !
    !     Year loop starts
    !
    DO nyear=start_year,end_year
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
                !
                ! Read the hydrologic and meteorologic forcings
                !
                call READ_FORCING
                !
                !     Begin reach computations
                !
                res_start=.FALSE.
                !
                !     Begin cycling through the reaches
                !
                DO nr=1,nreach
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
                    DO ns=1,no_celm(nr)
                        !
                        ncell=segment_cell(nr,ns)
                        !if (res_pres(ncell)) then
                            !write(*,*) ncell
                        !end if
                        !
                        !   If this segment is not located in the reservoir
                        !
                        if (.NOT. res_pres(ncell)) then
                            !
                            DONE=.FALSE.
                            !
                            ! Testing new code 8/8/2016
                            !
                            !     Establish particle tracks
                            !
                            call Particle_Track(nr,ns,nx_s,nx_head,ns_res_pres,ns_res_num)
                            !
                            !     Now do the third-order interpolation to
                            !     establish the starting temperature values
                            !     for each parcel
                            !
                            nseg=nstrt_elm(ns)
                            !
                            !     When water particle hits reservoir tracking upstream,
                            !     T_0 equals to the temperature of reservoir outlet.
                            !
                            if(ns_res_pres) then
                                T_0 = temp(nr,nseg-1,n1)
                            else
                                !
                                !     Perform polynomial interpolation
                                !
                                !
                                !     Interpolation inside the domain
                                !
                                npndx=2
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
                                        !if(nr.eq.2 .and.ns.eq.2 .and. nd.lt.30) write(*,*) &
                                        !    'nd', nd, 'ntrp', ntrp, 'npart',npart,xa(ntrp), ta(ntrp)
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
                            end if
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
                            !
                            ! Set initial river storage (Initial storage for first grid cell is 0)
                            !
                            !!if(nyear.eq.start_year .and. nd.eq.1 .and. ndd.eq.1 .and. ns.ge.3) then
                                !sto(nr,ns,n1) = width(ncell) * depth(ncell) * dx(ncell)/2
                            !!    temp_sto(nr,ns,n1) = T_head(nr)
                            !!end if
                            !
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
                                T_dstrb = smooth_param(nr)*dbt(nncell)+rminsmooth*T_smth(nr)
                                !
                                Q_dstrb = Q_local(nncell)
                                !
                                T_dstrb_load  = T_dstrb * Q_dstrb
                                !
                                !    Add flow from transient river storage
                                !
                                !T_sto_load = 0.0
                                !
                                !Q_sto = delta_sto_flux(nncell)
                                !
                                !if(nm.eq.1) Q_sto_out = delta_sto_flux(nncell)
                                !
                                !     Look for a tributary.
                                !
                                ntribs=no_tribs(nncell)
                                Q_trb_sum   = 0.0
                                T_trb_load  = 0.0
                                if(ntribs.gt.0) then
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
                                end if
                                Q_trb_sum = Q_trb_sum/2
                                T_trb_load = T_trb_load/2
                                !
                                ! Update inflow/outflow based on which segment
                                ! the water particle is located
                                ! We assume there are two segments in one grid cell
                                !
                                !if(mod(nseg,2).eq.0.and..not.DONE) then
                                !    Q_inflow=Q_outflow - Q_dstrb - Q_sto - Q_trb_sum
                                !    DONE = .TRUE.
                                !else
                                !    Q_outflow=Q_inflow + Q_dstrb + Q_sto + Q_trb_sum
                                !    DONE = .TRUE.
                                !end if
                                !
                                ! if Q_sto is smaller than 0, we need to adjust the water
                                ! balance.
                                !
                                !!T_sto = temp_sto(nr,nseg,n1)
                                !!if(Q_sto.ge.0) then
                                    !!sto_pre  = sto(nr,nseg,n1)
                                    !!sto_post = sto_pre - Q_sto*sec_day
                                    !!if(sto_post .lt. 0) then
                                    !!    sto_post = 0
                                    !!    T_sto_load=sto_pre/sec_day*T_sto
                                    !!end if
                                    !
                                    ! When Q_sto is smaller than 0, we need to adjust
                                    ! Q_dstrb, Q_trb and Q_local to make water balance.
                                    !
                                !!else
                                    !
                                    ! If local flow can fulfill the deficiency in river storage
                                    !
                                    !!sto_pre  = sto(nr,nseg,n1)
                                    !!sto_post = sto_pre - Q_sto*sec_day
                                    !!T_sto_load=0
                                    !if (Q_dstrb + Q_sto .gt. 0) then
                                    !    Q_dstrb=Q_dstrb + Q_sto
                                    !    T_dstrb_load  = Q_dstrb*T_dstrb
                                    !    if(nm.eq.1) then
                                    !        sto(nr,ns,n1)=sto(nr,ns,n2) - Q_sto*sec_day
                                    !        temp_sto(nr,ns,n1) = ( sto(nr,ns,n2)*temp_sto(nr,ns,n2) &
                                    !                             + (-Q_sto * T_dstrb))/sto(nr,ns,n1)
                                    !    end if
                                    !else
                                    !    if (Q_dstrb + Q_trb_sum + Q_sto .gt. 0) then
                                    !        Q_trb_sum_origin = Q_trb_sum
                                    !        T_trb_load_origin = T_trb_load
                                    !        Q_trb_sum = Q_dstrb + Q_trb_sum + Q_sto
                                    !        T_trb_load = (Q_trb_sum/Q_trb_sum_origin)*T_trb_load
                                    !        if(nm.eq.1) then
                                    !            sto(nr,ns,n1)=sto(nr,ns,n2) - Q_sto*sec_day
                                    !            temp_sto(nr,ns,n1) = ( sto(nr,ns,n2)*temp_sto(nr,ns,n2) &
                                    !                                 + Q_dstrb*T_dstrb                  &
                                    !                                 + T_trb_load_origin - T_trb_load)  &
                                    !                                 /sto(nr,ns,n1)
                                    !        end if
                                    !        Q_dstrb=0
                                    !        T_dstrb_load=0
                                    !    else
                                    !        Q_inflow_origin=Q_inflow
                                    !        Q_inflow=Q_inflow + Q_dstrb + Q_trb_sum + Q_sto
                                    !        if(nm.eq.1) then
                                    !            sto(nr,ns,n1)=sto(nr,ns,n2) - Q_sto*sec_day
                                    !            temp_sto(nr,ns,n1) = ( sto(nr,ns,n2)*temp_sto(nr,ns,n2) &
                                    !                                 + Q_dstrb*T_dstrb                  &
                                    !                                 + T_trb_load                       &
                                    !                                 + (Q_inflow_origin-Q_inflow)*T_0)   &
                                    !                                 /sto(nr,ns,n1)
                                    !        end if
                                    !        Q_dstrb=0
                                    !        T_dstrb_load=0
                                    !        Q_trb_sum=0
                                    !        T_trb_load=0
                                    !    end if
                                    !end if
                                    !Q_sto = 0
                                    !T_sto_load=0
                                !!end if
                                !
                                !  Update inflow and outflow
                                !
                                !!if(nr.eq.8 .and.ns.eq.3 .and. nd.lt.300 .and. nm.eq.1) write(*,*) &
                                !!    'nd', nd,  'T_0_pre', T_0, &
                                !!    Q_outflow, Q_inflow, Q_dstrb, Q_trb_sum, Q_sto, sto_post/sec_day
                                !!Q_tot = Q_inflow + Q_dstrb + Q_trb_sum + Q_sto + sto_post/sec_day
                                !!Q_ratio = Q_inflow/Q_tot
                                !
                                ! Do the mass/energy balance
                                !
                                Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum
                                if (Q_outflow .gt. 0) then
                                    T_0 = (T_0*Q_inflow + T_dstrb_load + T_trb_load)            &
                                        /(Q_inflow + Q_dstrb + Q_trb_sum)
                                else
                                    T_0 = T_0
                                end if
                                !if(nr.eq.8 .and.ns.eq.3 .and. nd.lt.300 .and. nm.eq.1) write(*,*) &
                                !    'nd', nd, 'T_0 post', T_0, 'energy', q_dot*dt_calc, &
                                !    'Inflow', Q_ratio, 'local', T_dstrb_load/Q_tot, &
                                !    'tributary', T_trb_load/Q_tot, 'storage flow', T_sto_load/Q_tot, &
                                !    'storage', sto_post/sec_day*T_sto/Q_tot
                                !
                                if (T_0.lt.0.5) T_0 =0.5
                                !
                                Q_inflow_out=Q_inflow
                                Q_outflow_out=Q_outflow
                                !
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
                            if (T_0.lt.0.5) T_0=0.5
                            temp(nr,ns,n2)=T_0
                            T_trib(nr)=T_0
                            !!sto(nr,ns,n2)=sto_post
                            !!temp_sto(nr,ns,n2)= T_0
                        !
                        !   if the segment is located in reservoir
                        !
                        call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell), &
                            Q_inflow_out, Q_outflow_out)
                        end if

                        if (res_pres(ncell)) then
                            !  res_no is the reservoir number
                            res_no = res_num(ncell)
                            !write(*,*) res_no
                            !
                            !     If it is located in the first segment in the
                            !     river reach, initialization for calculation
                            !
                            if(ncell .eq. res_start_node(res_no) .and. &
                                .NOT. res_start(res_no)) then
                                ns_res_start(res_no) = ns
                                Q_res_in(res_no) = Q_in(ncell)                  ! Inflow for reservoir
                                T_res_in(res_no) = temp(nr,ns-1,n1)               ! Inflow temperature
                                !if (res_no.eq.28 .and. nd.lt.115) write(*,*) &
                                !    nr, ns-1, n1, T_res_in(res_no)
                                Q_res_trib_tot = 0                              ! Initialize total tributary flow
                                T_res_trib_load_tot = 0                         ! Initialize tributary energy load
                                Q_res_dstrb_tot = Q_local(ncell)                ! Initialize local flow
                                T_res_dstrb = smooth_param(nr)*dbt(ncell)  &    ! Initialize local flow temperature
                                    + rminsmooth*T_smth(nr)
                                if (T_res_dstrb .lt. 0) T_res_dstrb=0
                                T_res_dstrb_load_tot = Q_res_dstrb_tot*T_res_dstrb
                                !if (res_no.eq.28 .and. nd.lt.115) write(*,*) &
                                !    'local', ncell, Q_res_dstrb_tot, T_res_dstrb
                                res_start(res_no) = .TRUE.
                                !
                                !     Calculate tributary flow for reservoir grid cells
                                !
                                ntribs=no_tribs(ncell)
                                if(ntribs.gt.0) then
                                    do ntrb=1,ntribs
                                        nr_trib=trib(ncell,ntrb)
                                        if(Q_trib(nr_trib).gt.0.0) then
                                            Q_trb             = Q_trib(nr_trib)
                                            Q_res_trib_tot    = Q_res_trib_tot + Q_trb
                                            !
                                            !  Update water temperature with tributary input
                                            !
                                            T_res_trib_load_tot   = (Q_trb*T_trib(nr_trib))       &
                                                +  T_res_trib_load_tot
                                        end if
                                    end do
                                    res_trib_calc(ncell) = .TRUE.
                                end if
                            else
                                !
                                !     Update total local flow
                                !
                                T_res_dstrb = smooth_param(nr)*dbt(ncell) &
                                    + rminsmooth*T_smth(nr)
                                if (T_res_dstrb .lt. 0) T_res_dstrb = 0
                                T_res_dstrb_load_tot = T_res_dstrb*Q_local(ncell) + &
                                    T_res_dstrb_load_tot
                                Q_res_dstrb_tot = Q_local(ncell) + Q_res_dstrb_tot
                                !
                                !     Update tributary flow
                                !
                                ntribs=no_tribs(ncell)
                                if(ntribs.gt.0 .and. .NOT. res_trib_calc(ncell)) then
                                    do ntrb=1,ntribs
                                        nr_trib=trib(ncell,ntrb)
                                        if(Q_trib(nr_trib).gt.0.0) then
                                            Q_trb             = Q_trib(nr_trib)
                                            Q_res_trib_tot    = Q_res_trib_tot + Q_trb
                                            !
                                            !  Update water temperature with tributary input
                                            !
                                            T_res_trib_load_tot   = (Q_trb*T_trib(nr_trib))       &
                                                +  T_res_trib_load_tot
                                        end if
                                    end do
                                    res_trib_calc(ncell) = .TRUE.
                                end if
                                if(ncell .eq. res_end_node(res_no) .and. &
                                    .NOT. res_pres(segment_cell(nr,ns+1))) then
                                    !write(*,*) 'test: start of loop 2'
                                    ns_res_end(res_no) = ns
                                    Q_res_outflow(res_no) = Q_out(ncell)
                                    !
                                    !      Sum up all streamflow into inflow
                                    !
                                    Q_res_inflow(res_no) = Q_res_in(res_no) + Q_res_dstrb_tot + Q_res_trib_tot
                                    if (Q_res_inflow(res_no) .gt. 0) then
                                        T_res_inflow(res_no) = (Q_res_in(res_no)*T_res_in(res_no) +  &
                                            T_res_trib_load_tot + T_res_dstrb_load_tot)/  &
                                            (Q_res_in(res_no) + Q_res_dstrb_tot + Q_res_trib_tot)
                                    else
                                        T_res_inflow(res_no) = T_res_in(res_no)
                                    end if
                                    !if (res_no.eq.28) write(*,*) &
                                    !    'inflow',Q_res_in(res_no), &
                                    !    'trb', Q_res_trib_tot, &
                                    !    'local', Q_res_dstrb_tot
                                    !
                                    !     Calculate the density of inflow, compared with epilimnion/hypolimnion
                                    !
                                    call stream_density(T_res_inflow(res_no), density_in(res_no))
                                    call stream_density(T_epil(res_no), density_epil(res_no))
                                    call stream_density(T_hypo(res_no), density_hypo(res_no))
                                    !
                                    !     Calculate flow subroutine
                                    !
                                    !call flow_subroutine(flow_in_epi_x, flow_in_hyp_x, flow_epi_hyp_x, &
                                    !    flow_out_epi_x, flow_out_hyp_x, res_no)
                                    call flow_subroutine(res_no)
                                    !
                                    !     Energy balance
                                    !
                                    call energy(T_epil(res_no), q_surf, ncell)
                                    !
                                    !     Calculate reservoir temperature
                                    !
                                    call reservoir_subroutine_implicit(res_no,q_surf)
                                    !
                                    T_0 = T_hypo(res_no) ! In reservoir, water is released from hypolimnion
                                    if (T_0.lt.0.5) T_0=0.5
                                    temp(nr,ns,n2)=T_0
                                    T_trib(nr)=T_0
                                    !
                                    !     Write output
                                    !
                                    do nsegment=ns_res_start(res_no),ns_res_end(res_no)
                                        call WRITE(time,nd,nr,ncell,nsegment,T_0, &
                                            T_head(nr),dbt(segment_cell(nr,nsegment)), &
                                            Q_res_inflow(res_no), Q_res_outflow(res_no), &
                                            res_storage_post, T_res(res_no))
                                    end do

                                else

                                end if
                            end if
                        end if
                        !
                        !   Write all temperature output UW_JRY_11/08/2013
                        !   The temperature is output at the beginning of the
                        !   reach.  It is, of course, possible to get output at
                        !   other points by some additional code that keys on the
                        !   value of ndelta (now a vector)(UW_JRY_11/08/2013)
                        !
                        !call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell), &
                        !   Q_inflow_out, Q_outflow_out,                   &
                        !   Q_sto_out, sto_post, temp_sto(nr,ns,n2))
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
4650            format(16x,12(6x,f6.0,6x))
4700            format(f10.4,f6.0,15(f6.1,f8.3))
4750            format(f10.4,10(i4,f8.0))
                write(75,*) nyear, nd, T_epil(:), T_hypo(:)
                write(76,*) nyear, nd, depth_e(:), depth_h(:)
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
