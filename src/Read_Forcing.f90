SUBROUTINE Read_Forcing
!
USE Block_Energy
USE Block_Hydro
USE Block_Network
USE Block_Reservoir
! 
IMPLICIT NONE
!
integer:: nc,ncell,nnd,no_flow,no_heat,nr,nrec_flow,nrec_heat, ndx
real:: Q_avg


!       
!        upload thermal point-source input
! 
if (source) then
     read(40,*) ndx, flow_source(1:nsource)
end if


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

        read(35,*) nnd,ncell &   !  flow file
           ,Q_in(no_heat),Q_out(no_heat),Q_diff(no_heat) &  
           ,depth(no_heat),width(no_heat),u(no_heat)
        write(85,*),Q_in(no_heat),Q_out(no_heat),Q_diff(no_heat),depth(no_heat),width(no_heat),u(no_heat)
        !
        if(u(no_heat).lt.0.01) u(no_heat)=0.01
        if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat 
        !
        read(36,*) ncell &  !  heat file
            ,dbt(no_heat),ea(no_heat) &
            ,Q_ns(no_heat),Q_na(no_heat),rho &
            ,press(no_heat),wind(no_heat)
        write(86,*),dbt(no_heat),ea(no_heat),Q_ns(no_heat),Q_na(no_heat),rho,press(no_heat),wind(no_heat)
        !          
        !  if(ncell.ne.no_heat) write(*,*) 'Heat file error',ncell,no_heat 
        !
        !  Added variable ndelta (UW_JRY_2011/03/15
        !

        delta_n=ndelta(ncell)
        ! 
        Q_avg=0.5*(Q_in(no_heat)+Q_out(no_heat))
        Q_diff(no_heat)=Q_diff(no_heat)/delta_n
        dt(no_heat)=dx(no_heat)/u(no_heat)  ! time(sec) to travel between segments,  u=velocity(ft/sec)
! print *,'nd',nnd,'no_heat',no_heat, 'ncell', ncell, 'dx(ncell)', dx(ncell), 'u(no_heat)', u(no_heat), 'dt',dt(no_heat)

  !   print *, no_heat,dt(no_heat)/86400
        !
        !  Added check to see if travel time of parcel exceeds the
        !  computational interval.  If so, it writes to file fort.45.
        !  One should check to see if there are NaN's at that node.
        !  (UW_JRY_2011/03/15)
        !
        if(dt(no_heat).gt.dt_comp) write(45,*) &
             'Travel time=',dt(no_heat) &
              , '> dt_comp at node -',no_heat


        !
        !        Calculate the next upstream dam and time for water to get there 
        !
     !  if(reservoir) THEN
     !        dt_res(no_heat) =dx_res(no_heat)/u(no_heat)
     !   END IF

    end do

    !
    !      Tributary flow is Q_out from the next to the last cell
    !
    !
    !       Read the meteorology for the last cell, but not the flow
    !
    no_heat=no_heat+1 
    Q_in(no_heat)=Q_out(no_heat-1)
    Q_out(no_heat)=Q_in(no_heat)
    Q_trib(nr)=Q_out(no_heat)   ! flow out this reach is Q_trib 
    nrec_heat=heat_cells*(ndays-1)+no_heat
    read(36,*) ncell &
         ,dbt(no_heat),ea(no_heat) &
         ,Q_ns(no_heat),Q_na(no_heat),rho &
         ,press(no_heat),wind(no_heat)
    !
    !  The flow and hydraulics for the last cell has to be 
    !  modified so they do not
    !  take the values of the segment to which it is tributary
    !
    Q_in(ncell)=Q_out(ncell-1)
    !  Q_out(ncell)=Q_in(ncell-1)
    Q_diff(no_heat)=0.0
    u(no_heat)=u(no_heat-1)
    depth(no_heat)=depth(no_heat-1)
    width(no_heat)=width(no_heat-1)
    dt(no_heat)=0.5*dx(ncell)/u(no_heat)

   ! ################ This is specially for simple energy test###########!                
  ! dt(no_heat)=dx(ncell)/u(no_heat)

! print *,'nd',nnd,'no_heat',no_heat, 'ncell', ncell, 'dx(ncell)', dx(ncell), 'u(no_heat)', u(no_heat), 'dt',dt(no_heat)

!if(nnd.gt.300) stop !13505
end do

END SUBROUTINE Read_Forcing
