SUBROUTINE SYSTMM(temp_file,res_file,param_file)
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
    character (len=200):: res_file
    !
    integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
    integer          :: nc,nd,ndd,nm,nr,ns
    integer          :: nr_trib,ntribs
    integer          :: nrec_flow,nrec_heat
    integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
    integer          :: npart,nseg,nx_s,nx_part,nx_head
    integer          :: ns_res_num, res_no, nsegment
    integer          :: nseg_temp
    !
    ! Indices for lagrangian interpolation
    !
    integer              :: njb,npndx,ntrp
    integer, dimension(2):: ndltp=(/-1,-2/)
    integer, dimension(2):: nterp=(/2,3/)
    !
    real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z,q_dot_pre
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
    allocate (volume_e_min(nres))
    allocate (res_stratif_start(nres))
    allocate (res_turnover(nres))
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
    surface_area = res_area_km2(:) * km_to_m * km_to_m
    !
    !     Year loop starts
    !
    DO nyear=start_year,end_year
        write(*,*) ' Simulation Year - ',nyear,start_year,end_year
        nd_year=365
        if (mod(nyear,4).eq.0) nd_year=366
        !
        !     Reset stratification and turnover
        !
        res_stratif_start = .FALSE.
        res_turnover = .FALSE.
        !
        !     At the start of each year, reset the depth of epilimnion
        !     and hypolimnion.
        !
        res_depth_meter = depth_e(:) + depth_h(:)
        depth_e = res_depth_meter(:) * depth_e_frac
        depth_h = res_depth_meter(:) * depth_h_frac
        volume_e_x = surface_area(:) * depth_e(:)
        volume_h_x = surface_area(:) * depth_h(:)
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
                res_trib_calc=.FALSE.
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
                            do nm=no_dt(ns),1,-1
                                dt_calc=dt_part(nm)
                                z=depth(nncell)
                                call energy(T_0,q_surf,nncell,z,nd)
                                !
                                ! apply a different numerical method to solve
                                ! energy balance for river temperature
                                !
                                q_dot=(q_surf/(z*rfac))
                                !
                                ! Initialize dH/dt (temperature increase because
                                ! of energy balance)
                                !
                                if (nyear.eq.start_year .and. nd.eq.1) q_dot_pre = q_dot
                                !
                                ! caluclate stream temperature based on Heun
                                ! mehtod
                                !
                                T_0=T_0+q_dot*dt_calc
                                !
                                ! estimate the error by numerical method 
                                !
                                error_EE=deriv_2nd*dt_calc**2/2
                                !
                                !if (ncell.eq.4089 .and. nd.gt.220.and.nd.lt.230) write(*,*) &
                                !     nyear,nd,ns,nm,wind(nncell),q_ns(nncell)+q_na(nncell),T_0,error_EE
                                if(T_0.lt.0.0) T_0=0.0
                                q_dot_pre = q_dot
                                !
                                !    Add distributed flows
                                !
                                T_dstrb = smooth_param(nr)*dbt(nncell)+rminsmooth*T_smth(nr)
                                !
                                Q_dstrb = Q_local(nncell)
                                !
                                T_dstrb_load  = T_dstrb * Q_dstrb
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
                                Q_trb_sum = Q_trb_sum/ndelta(nncell)
                                T_trb_load = T_trb_load/ndelta(nncell)
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
                        !   
                        !   test output for a specific grid cell 
                        !
                        !
                        do nseg_temp=1,nseg_out_num
                            if (nseg_out(nr,ncell,nseg_temp).eq.ns) then
                                call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell), &
                                    Q_inflow_out, Q_outflow_out)
                            !if (ncell.eq.4089) write(*,*)nyear,nd,nr,ncell,ns,T_0
                            end if
                        end do
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
                                Q_res_trib_tot = 0                              ! Initialize total tributary flow
                                T_res_trib_load_tot = 0                         ! Initialize tributary energy load
                                Q_res_dstrb_tot = Q_local(ncell)                ! Initialize local flow
                                T_res_dstrb = smooth_param(nr)*dbt(ncell)  &    ! Initialize local flow temperature
                                    + rminsmooth*T_smth(nr)
                                if (T_res_dstrb .lt. 0) T_res_dstrb=0
                                T_res_dstrb_load_tot = Q_res_dstrb_tot*T_res_dstrb
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
                                !
                                !     If it is the last segment in the reservoir
                                !
                                if(ncell .eq. res_end_node(res_no) .and. &
                                    (.NOT. res_pres(segment_cell(nr,ns+1)) &
                                    .or. res_num(segment_cell(nr,ns+1)) .ne. res_no)) then
                                    exceed_error_bound=.FALSE.
                                    adjust_timestep=.FALSE.
                                    ns_res_end(res_no) = ns
                                    Q_res_outflow(res_no) = Q_out(ncell)     !!!!TESTINFLOW(need to uncomment)
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
                                    !
                                    !     Calculate the density of inflow, compared with epilimnion/hypolimnion
                                    !
                                    call stream_density(T_res_inflow(res_no), density_in(res_no))
                                    !call stream_density(T_epil(res_no),density_epil(res_no))
                                    !call stream_density(T_hypo(res_no),density_hypo(res_no))
                                    !
                                    !     Start the subdaily calculation
                                    !
500                                 continue
                                    numsub=1
                                    !
                                    !     calculate the number of timestep based
                                    !     on error bound
                                    !
                                    if (exceed_error_bound) then
                                        numsub=MAX(ceiling(sqrt(abs(error_e)/error_threshold)), &
                                                   ceiling(sqrt(abs(error_h)/error_threshold)))
                                        numsub=MIN(20,numsub)
                                    end if
                                    dt_res = dt_comp/numsub
                                    DO nsub=1,numsub
                                        !
                                        !     Calculate stream density 
                                        !
                                        call stream_density(T_epil(res_no), density_epil(res_no))
                                        call stream_density(T_hypo(res_no), density_hypo(res_no))
                                        !
                                        !     Calculate flow subroutine
                                        !
                                        call flow_subroutine(res_no, nyear, nd)
                                        !
                                        !     Energy balance
                                        !
                                        call energy(T_epil(res_no), q_surf, ncell)
                                        !
                                        !     Error estimation
                                        !
                                        if (.not.exceed_error_bound) then
                                            call Error_estimate(nd,res_no)
                                        end if
                                        !
                                        !     Adjust the timestep based on error
                                        !
                                        if (exceed_error_bound .and..not.adjust_timestep) then
                                             adjust_timestep=.TRUE.
                                             go to 500
                                        end if
                                        !
                                        !     Calculate reservoir temperature
                                        !
                                        !call reservoir_subroutine(res_no,q_surf)
                                        !if (res_storage(res_no,1) .gt.res_capacity_mcm(res_no)*(10**6)*0.2) then
                                            call reservoir_subroutine_implicit(res_no,q_surf,nd,dbt(ncell))
                                        !if(res_no.eq.63) write(*,*) &
                                        !        nyear,nd,nsub,dt_res,T_epil(res_no),T_hypo(res_no)
                                    end do
                                         
                                    !else
                                    !    call reservoir_single_subroutine(res_no,q_surf,nd)
                                    !end if 
                                    !
                                    T_0 = T_hypo(res_no) ! In reservoir, water is released from hypolimnion
                                    if (T_0.lt.0.5) T_0=0.5
                                    temp(nr,ns_res_start(res_no):ns_res_end(res_no),n2)=T_0
                                    T_trib(nr)=T_0
                                    !
                                    !     Write output
                                    !
                                    do nsegment=ns_res_start(res_no),ns_res_end(res_no)
                                        do nseg_temp=1,nseg_out_num
                                            if (nseg_out(nr,ncell,nseg_temp).eq.nsegment) then
                                                call WRITE(time,nd,nr,ncell,nsegment,T_0, &
                                                    T_head(nr),dbt(segment_cell(nr,nsegment)), &
                                                    Q_res_inflow(res_no), Q_res_outflow(res_no), &
                                                    res_storage_post, T_res(res_no))
                                            end if
                                        end do 
                                    end do
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
                if (reservoir) then
                    open(75,file=TRIM(res_file),status='unknown') 
                    write(75, '(i8, i6, 1680f20.10)') nyear, nd, T_epil(:), T_hypo(:), T_res_inflow(:),& 
                        depth_e(:), depth_h(:),&
                        Q_res_inflow(:), Q_res_outflow(:)
                    write(76,'(i8, i6, 240f15.10, 240f15.10)') nyear, nd, depth_e(:), depth_h(:), volume_e_x,volume_h_x
                    write(77,*) nyear, nd, diffusion_tot, advec_hyp_tot,advec_epi_tot,qsurf_tot
                end if
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
