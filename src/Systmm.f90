SUBROUTINE SYSTMM(temp_file,param_file)
    !
    use Block_Energy
    use Block_Hydro
    use Block_Network
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
                        DONE=.FALSE.
                        !
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
                        if(nyear.eq.start_year .and. nd.eq.1 .and. ndd.eq.1 .and. ns.ge.3) then
                            !sto(nr,ns,n1) = width(ncell) * depth(ncell) * dx(ncell)/2
                            temp_sto(nr,ns,n1) = T_head(nr)
                        end if
                        !
                        do nm=no_dt(ns),1,-1
                            dt_calc=dt_part(nm)
                            z=depth(nncell)
                            call energy(T_0,q_surf,nncell)
                            !
                            q_dot=(q_surf/(z*rfac))
                            !T_0=T_0+q_dot*dt_calc
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
                            T_sto_load = 0.0
                            !
                            Q_sto = delta_sto_flux(nncell)
                            !
                            if(nm.eq.1) Q_sto_out = delta_sto_flux(nncell)
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
                            if(mod(nseg,2).eq.0.and..not.DONE) then
                                Q_inflow=Q_outflow - Q_dstrb - Q_sto - Q_trb_sum
                                DONE = .TRUE.
                            else
                                Q_outflow=Q_inflow + Q_dstrb + Q_sto + Q_trb_sum
                                DONE = .TRUE.
                            end if
                            !
                            ! if Q_sto is smaller than 0, we need to adjust the water
                            ! balance.
                            !
                            T_sto = temp_sto(nr,nseg,n1)
                            if(Q_sto.ge.0) then
                                sto_pre  = sto(nr,nseg,n1)
                                sto_post = sto_pre - Q_sto*sec_day
                                if(sto_post .lt. 0) then
                                    sto_post = 0
                                    T_sto_load=sto_pre/sec_day*T_sto
                                end if
                                !
                                ! When Q_sto is smaller than 0, we need to adjust
                                ! Q_dstrb, Q_trb and Q_local to make water balance.
                                !
                            else
                                !
                                ! If local flow can fulfill the deficiency in river storage
                                !
                                sto_pre  = sto(nr,nseg,n1)
                                sto_post = sto_pre - Q_sto*sec_day
                                T_sto_load=0
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
                            end if
                            !
                            !  Update inflow and outflow
                            !
                            if(nr.eq.8 .and.ns.eq.3 .and. nd.lt.300 .and. nm.eq.1) write(*,*) &
                                'nd', nd,  'T_0_pre', T_0, &
                                Q_outflow, Q_inflow, Q_dstrb, Q_trb_sum, Q_sto, sto_post/sec_day
                            Q_tot = Q_inflow + Q_dstrb + Q_trb_sum + Q_sto + sto_post/sec_day
                            Q_ratio = Q_inflow/Q_tot
                            !
                            ! Do the mass/energy balance
                            !
                            !if(nr.eq.8 .and.ns.eq.3 .and. nd.lt.30 .and. nm.eq.1) write(*,*) &
                            !'nd', nd, 'T_0 pre', T_0
                            T_0 = T_0*Q_ratio                                           &
                                + (T_dstrb_load + T_trb_load + T_sto_load               &
                                + sto_post/sec_day*T_sto)/Q_tot                         &
                                + q_dot*dt_calc
                            if(nr.eq.8 .and.ns.eq.3 .and. nd.lt.300 .and. nm.eq.1) write(*,*) &
                                'nd', nd, 'T_0 post', T_0, 'energy', q_dot*dt_calc, &
                                'Inflow', Q_ratio, 'local', T_dstrb_load/Q_tot, &
                                'tributary', T_trb_load/Q_tot, 'storage flow', T_sto_load/Q_tot, &
                                'storage', sto_post/sec_day*T_sto/Q_tot
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
                        sto(nr,ns,n2)=sto_post
                        temp_sto(nr,ns,n2)= T_0
                        !
                        !   Write all temperature output UW_JRY_11/08/2013
                        !   The temperature is output at the beginning of the
                        !   reach.  It is, of course, possible to get output at
                        !   other points by some additional code that keys on the
                        !   value of ndelta (now a vector)(UW_JRY_11/08/2013)
                        !
                        !if(ncell.ne.nncell) write(*,*) &
                        !    ns, nd, ncell, nncell
                        call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell), &
                                   Q_inflow_out, Q_outflow_out,                   &
                                   Q_sto_out, sto_post, temp_sto(nr,ns,n2))
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
