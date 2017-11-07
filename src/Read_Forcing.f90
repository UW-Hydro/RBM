SUBROUTINE Read_Forcing
    !
    USE Block_Energy
    USE Block_Hydro
    USE Block_Network
    USE Block_Reservoir
    !
    IMPLICIT NONE
    !
    integer :: nc,ncell,nnd,no_flow,no_heat,nr,nrec_flow,nrec_heat
    integer :: nreservoir,n1,n2
    real    :: Q_avg,Q_dmmy
    real    :: z_temp,w_temp
    real    :: min_flow = 5.0
    !
    n1=1
    n2=2
    no_flow=0
    no_heat=0
    do nr=1,nreach
        do nc=1,no_cells(nr)-1
            no_flow=no_flow+1
            no_heat=no_heat+1
            !
            nrec_flow=flow_cells*(ndays-1)+no_flow
            nrec_heat=heat_cells*(ndays-1)+no_heat
            !
            read(35,*) nnd,ncell &
                ,Q_out(no_heat),Q_dmmy,Q_diff(no_heat) &
                ,depth(no_heat),width(no_heat),u(no_heat), Q_local(no_heat)
            !
            !     If the streamflow is 0, recalculate water velocity based on
            !     local streamflow
            !
            if (Q_out(no_heat) .lt. min_flow) then
                z_temp = a_z * (min_flow**b_z)
                w_temp = a_w * (min_flow**b_w)
                u(no_heat) = min_flow/(z_temp*w_temp)
                depth(no_heat) = z_temp
                width(no_heat) = w_temp
            end if
            !
            if (Q_diff(no_heat) .gt. 0) write(*,*) 'Q_diff is not equal to 0', no_heat, Q_diff(no_heat)
            !
            if(u(no_heat).lt.0.01) u(no_heat)=0.01
            if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat
            !
            read(36,*) ncell &
                ,dbt(no_heat),ea(no_heat) &
                ,Q_ns(no_heat),Q_na(no_heat),rho &
                ,press(no_heat),wind(no_heat)
            !
            if(ncell.ne.no_heat) write(*,*) 'Heat file error',ncell,no_heat
            !
            !  Added variable ndelta (UW_JRY_2011/03/15
            !
            delta_n=ndelta(ncell)
            !
            Q_avg=0.5*(Q_in(no_heat)+Q_out(no_heat))
            dt(no_heat)=dx(no_heat)/u(no_heat)
            !
            !  Added check to see if travel time of parcel exceeds the
            !  computational interval.  If so, it writes to file fort.45.
            !  One should check to see if there are NaN's at that node.
            !  (UW_JRY_2011/03/15)
            !
            if(dt(no_heat).gt.dt_comp) write(45,*) &
                'Travel time=',dt(no_heat) &
                , '> dt_comp at node -',no_heat
        end do
        !
        ! Tributary flow is Q_out from the next to the last cell
        ! However, it will be updated in Water_Balance to account
        ! for one-half the runoff from the confluence cell.
        !
        !
        !       Read the meteorology for the last cell, but not the flow
        !
        no_heat=no_heat+1
        read(36,*) ncell &
            ,dbt(no_heat),ea(no_heat) &
            ,Q_ns(no_heat),Q_na(no_heat),rho &
            ,press(no_heat),wind(no_heat)
        if (nr .eq. nreach .and. nc .eq. no_cells(nr)) then
            read(35,*) nnd,ncell &
                ,Q_out(no_heat),Q_dmmy,Q_diff(no_heat) &
                ,depth(no_heat),width(no_heat),u(no_heat), Q_local(no_heat)
            Q_trib(nr)=0
            Q_in(no_heat)=Q_out(no_heat-1)
            dt(no_heat)=dx(ncell)/u(no_heat)
        else
            Q_out(no_heat)=Q_out(no_heat-1)
            !
            ! Tributary flow from this reach equals Q_out for this cell
            !
            Q_trib(nr)=Q_out(no_heat)
            nrec_heat=heat_cells*(ndays-1)+no_heat
            !
            !  The flow and hydraulics for the last cell has to be
            !  modified so they do not
            !  take the values of the segment to which it is tributary
            !
            Q_in(no_heat)=Q_out(no_heat-1)
            u(no_heat)=u(no_heat-1)
            depth(no_heat)=depth(no_heat-1)
            width(no_heat)=width(no_heat-1)
            dt(no_heat)=dx(ncell)/u(no_heat)
        end if
    end do
    !
    ! Read in reservoir storage data
    !
    res_storage(:,n2)=res_storage(:,n1)
    read(38,*) ! Skip reading the date in storage file
    do nreservoir=1,nres
        read(38,*) res_storage(nreservoir, n1)
    end do
    !
    ! Call the water balance subroutine
    !
    call Water_Balance
!
END SUBROUTINE Read_Forcing
