!
!      PROGRAM RMAIN
!
!     Dynamic river basin model for simulating water quality in
!     branching river systems with freely-flowing river segments. 
!
!     This version uses Reverse Particle Tracking in the Lagrangian
!     mode and Lagrangian interpolation in the Eulerian mode.
!
!     Topology and routing is set up to be consistent with output
!     from the Variable Infiltration Capacity (VIC) model developed by the
!     Land Surface Hydrology Group at the University of Washington.
!
!     For additional information contact:
!
!     John Yearsley
!     Land Surface Hydrology Group
!     Dept. of Civil and Environmental Engineering
!     Box 352700
!     University of Washington
!     Seattle, Washington
!     98195-2700
!
!     Modifications by Wageningen University and Research Centre (WUR)
!			Earth System Science and Climate Change Group
!			P.O. Box 47, 6700 AA  Wageningen, The Netherlands
!     for PhD research Michelle van Vliet
!			Wietse Franssen and Michelle van Vliet(WUR_WF_MvV) 
!
!
character*80 net_file,heat_file,flow_file
!
! Extra Variables (added by WUR_WF_MvV_2011/01/05)
!
character * ( 200 ) inPrefix
character * ( 200 ) outPrefix
character * ( 200 ) basin
character * ( 200 ) param_file
character * ( 200 ) temp_file
character * ( 200 ) spatial_file      
integer iargc
integer numarg

!
! Command line input (WUR_WF_MvV_2011/01/05)
!
numarg = iargc ( )
if (numarg .lt. 2) then
  write (*,*) 'Too less arguments were given'
  write (*,*) ' '
  write (*,*) 'First:  Location and prefix of input files'
  write (*,*) '        (networkfile and parameterfile)'
  write (*,*) 'Second: Location and prefix of output files'
  write (*,*) ' '
  write (*,*) 'eg: $ <program-name> ./inp/Rhine ./outp/Rhine'
  write (*,*) ' '
  stop
end if
call getarg ( 1, inPrefix )
call getarg ( 2, outPrefix )
!

! (WUR_WF_MvV_2011/01/05)
ncfile=INDEX(inPrefix,' ')-1
net_file  =inPrefix(1:ncfile)//'_Network'
param_file=inPrefix(1:ncfile)//'_Parameters'
ncfile=INDEX(net_file,' ')-1
write(*,*) 'Network file    : ',net_file(1:ncfile)
ncfile=INDEX(param_file,' ')-1
write(*,*) 'Parameter file  : ',param_file(1:ncfile)

ncfile=INDEX(outPrefix,' ')-1
temp_file  =outPrefix(1:ncfile)//'.Temp'
ncfile=INDEX(temp_file,' ')-1
write(*,*) 'Temperature file: ',temp_file(1:ncfile)

ncfile=INDEX(outPrefix,' ')-1
spatial_file  =outPrefix(1:ncfile)//'.Spat'
ncfile=INDEX(spatial_file,' ')-1
write(*,*) 'Spatial file: ',spatial_file(1:ncfile)         
!

!
!     Read name of network file from Command line (WUR_WF_MvV_2011/01/05)
!
!      write(*,*) 'Input name of network file'
!      read(*,*) net_file
ncfile=INDEX(net_file,' ')-1
OPEN(UNIT=90,FILE=net_file(1:ncfile),STATUS='OLD')
!

!
!     Read header information from control file
!
read(90,*)
read(90,'(A)') flow_file
!
!     Open file with hydrologic data
!
ncflow=INDEX(flow_file,' ')-1
write(*,*) ncflow, flow_file
open(unit=35,FILE=flow_file(1:ncflow) ,FORM='FORMATTED',ACCESS='DIRECT' ,RECL=60,STATUS='old')
!
!
read(90,'(A)') heat_file
!
!     Open file with meteorologic data
!     
ncheat=INDEX(heat_file,' ')-1
write(*,*) ncheat, heat_file
open(unit=36,FILE=heat_file(1:ncheat) ,FORM='FORMATTED',ACCESS='DIRECT' ,RECL=50,STATUS='old')
!
!     Read header on energy file to position file
!
!      read(35,*)
!
!     Call systems programs to get started
!
!     SUBROUTINE BEGIN reads control file, sets up topology and
!     important properties of reaches
!
write(*,*) 'Calling BEGIN'
CALL BEGIN(param_file, spatial_file)
!
!     SUBROUTINE SYSTMM performs the simulations
!
CALL SYSTMM(temp_file,param_file) ! (WUR_WF_MvV_2011/01/05)
!
!     Close files after simulation is complete
!
write(*,*) ' Closing files after simulation'
CLOSE(35)
CLOSE(90)
1500 FORMAT(A80)
1510 format(1x,a30)
1600 FORMAT(8F10.0)
2600 FORMAT(' NAME OF FILE CONTAINING RIVER REACH DATA')
2700 FORMAT(' NAME OF OUTPUT DATA FILE')
2800 format(' Name of file with geometric data')
2900 format(' Name of file with hydrologic data')
3000 format(' Name of file with water quality data')
STOP
END
SUBROUTINE BEGIN(param_file, spatial_file)
character*8 lat
character*10 long
integer start_year,start_month,start_day,end_year,end_month,end_day,head_name,trib_cell,first_cell
dimension ndmo(12)
character * ( 200 ) param_file ! (WUR_WF_MvV_2011/01/05)
character * ( 200 ) source_file ! (UW_JRY_2011/04/04)
character * ( 200 ) spatial_file ! (WUR_WF_MvV_2011/01/05)
logical source
INCLUDE 'rbm10_VIC.h'
data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
parameter (n_default=2)
!
!     Initialize arrays and constants
!
do nr=1,50
   do ns=1,1000
      thermal(nr,ns)=0.0
   end do
end do
!
!     Card Group I
!
!rminsmooth=1.-smooth_param
!write(*,*) ' rminsmooth  : ', rminsmooth

 read(90,1100) start_year,start_month,start_day &
 ,end_year,end_month,end_day
1100 format(2(2x,i4,2i2))
nyear1=start_year
nyear2=end_year
write(*,1100) start_year,start_month,start_day &
 ,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
!
jd=ndmo(start_month)+start_day
jul_start=julian(start_year,jd)
!
read(90,*) no_rch,flow_cells,heat_cells,source
if (source) then
!   read(90,*) source_file
   read(90,'(A)') source_file ! (WUR_WF_MvV_2011/05/23)
   print *,'source file: ', source_file ! (WUR_WF_MvV_2011/05/23)
   ncfile=INDEX(source_file,' ')-1 ! (UW_JRY_2011/04/04)
   OPEN(UNIT=40,FILE=source_file(1:ncfile),STATUS='OLD')!  (UW_JRY_2011/04/04)
end if
write(*,1200) no_rch
1200 format(10i10)
!
!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
ncell=0
nreach=0
ns_total=0
100 continue
!
!     Card Group IIb. Reach characteristics
!
nreach=nreach+1
if(nreach.gt.no_rch) go to 500
!
!     Initialize NSEG, the total number of segments in this reach
!
nseg=0
write(*,*) ' Starting to read reach ',nreach
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
read(90,1400) no_cells(nreach) &
 ,head_name,trib_cell &
 ,main_stem(nreach),rmile0
1400 format(i5,11x,i4,10x,i5,15x,i5,15x,f10.0,i5)
!
!     If this is reach that is tributary to cell TRIB_CELL, give it the
!     pointer TRIB(TRIB_CELL) the index of this reach for further use.
!     Also keep track of the total number of tributaries for this cell
!
if (trib_cell.gt.0) then
  no_tribs(trib_cell)=no_tribs(trib_cell)+1
  trib(trib_cell,no_tribs(trib_cell))=nreach
end if
!
!     Reading Mohseni parameters for each headwaters (UW_JRY_2011/06/18)
!
read(90,*) alphaMu(nreach),beta(nreach) &
  ,gmma(nreach),mu(nreach),smooth_param(nreach)

!
!     Reading Reach Element information
!
first_cell=1
do nc=1,no_cells(nreach)
  ncell=ncell+1
!
!   Read the thermal input if there is a for point sources (UW_JRY_2011/04/04)
!
  if (source) then
    read(40,*) ns,heat_inp
    if(ns.ne.ncell) then
      write(*,*) 'Input error at segment - ',ncell,ns
      stop
    end if
  end if 
!
!     The headwaters index for each cell in this reach is given
!     in the order the cells are read
!
!            head_no(nc)=nreach
!
!     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
!
!
!     Variable ndelta read in here.  At present, number of elements
!     is entered manually into the network file (UW_JRY_2011/03/15)
!
  read(90,1600) node(ncell),nrow,ncol &
   ,lat,long,rmile1,ndelta(ncell)
1600 format(5x,i5,5x,i5,8x,i5,6x,a8,6x,a10,7x,f10.0,i5)
!
!    Set the number of segments per reach to the default if 
!    ndelta(nreach) is left blank (backward compatible)
!    (UW_JRY_2011/03/15)
!
  if (ndelta(ncell).lt.1) ndelta(ncell)=n_default
  if(first_cell.eq.1) then
    first_cell=0
    head_cell(nreach)=ncell
    x_dist(nreach,0)=5280.*rmile0
  end if
!
! Added variable ndelta (UW_JRY_2011/03/15)
!
  dx(ncell)=5280.*(rmile0-rmile1)/ndelta(ncell)
  rmile0=rmile1
  nndlta=0
200 continue
  nndlta=nndlta+1
  nseg=nseg+1
!
!  Divide heat_inp by rho*Cp (rho=1000 kg/m**3, Cp=1 deg K/kg/kcal)
!
if(nndlta.eq.2) thermal(nreach,nseg)=heat_inp/1000.
  ns_total=ns_total+1
  segment_cell(nreach,nseg)=ncell
  x_dist(nreach,nseg)=x_dist(nreach,nseg-1)-dx(ncell)
!
!   Write Segment List for mapping to temperature output (UW_JRY_2008/11/19)
!
  ncfile=INDEX(spatial_file,' ')-1 ! (added by WUR_WF_MvV_2011/01/05)
  open(22,file=(spatial_file(1:ncfile)),status='unknown') ! (changed by WUR_WF_MvV_2011/01/05)
  write(22,7000) nreach,ncell,nrow,ncol,lat,long,nndlta
!
! 
7000 format(4i6,1x,a8,1x,a10,i5)
!
!  Added variable ndelta  (UW_JRY_2011/03/15)
!
  if(nndlta.lt.ndelta(ncell)) go to 200   
  no_celm(nreach)=nseg
  segment_cell(nreach,nseg)=ncell
  x_dist(nreach,nseg)=5280.*rmile1
!
end do
!
!    Define last segment for purposes of specifying tributary
!    temperatures (UW_JRY_2008/11/14)
!
go to 100
500	continue
nreach=no_rch
!
!
nwpd=1
xwpd=nwpd
dt_comp=86400./xwpd
!
!     ******************************************************
!                         Return to RMAIN
!     ******************************************************
!
RETURN
900 END

SUBROUTINE SYSTMM(temp_file,param_file) ! (modified by WUR_WF_MvV_2011/01/05)
real*4 xa(4),ta(4),T_head(5000),T_smth(5000) &
 ,dt_part(500),x_part(500)
real*8 time 
integer no_dt(1000),nstrt_elm(500) &
 ,ndltp(4),nterp(4),nptest(4)
logical DONE
INCLUDE 'rbm10_VIC.h'
data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
data pi/3.14159/,rfac/304.8/
character * ( 200 ) temp_file ! (added by WUR_WF_MvV_2011/01/05)
character * ( 200 ) param_file ! (added by WUR_WF_MvV_2011/01/05)
!
!     open the output file
!
ncfile=INDEX(temp_file,' ')-1 ! (added by WUR_WF_MvV_2011/01/05)
open(20,file=(temp_file(1:ncfile)),status='unknown') ! (changed by WUR_WF_MvV_2011/01/05)
!
do nr=1,500
  T_head(nr)=0.0
end do
cum_cuft=(3.2808*3.2808*3.2808)
cuft_cum=1./cum_cuft
n1=1
n2=2
nobs=0
!
!     Initialize the day counter used for calculating the
!     record position for direct access files
!
ndays=0
xwpd=nwpd
hpd=1./xwpd
!
!     Year loop starts
!
do nyear=nyear1,nyear2
  write(*,*) ' Simulation Year - ',nyear
  nd_year=365
  if (mod(nyear,4).eq.0) nd_year=366
!
!     Day loop starts
!
  DO ND=1,nd_year
    year=nyear
    xd=nd-1.
    xd_year=nd_year 
!
!     Start the numbers of days-to-date counter
!
    ndays=ndays+1
    nwpd=1
!
!    Daily period loop starts
!
      DO NDD=1,NWPD
      xdd = ndd
      time=year+(xd+(xdd*hpd))/xd_year 

!
!     Begin reach computations
!
      ind=nd
      ipd=ndd
      day=nd
!
!     Read advected energy and meteorology data
!
!
!     Initialize the flow and heat cell counters
!
      no_flow = 0
      no_heat  = 0
      do nr=1,nreach
!
1200 format(2i5,f10.0)
!
!
        do nc=1,no_cells(nr)-1
          no_flow=no_flow+1
          no_heat=no_heat+1
          nrec_flow=flow_cells*(ndays-1)+no_flow
          nrec_heat=heat_cells*(ndays-1)+no_heat
!
          read(35,'(2i5,2f10.1,2f6.1,f7.1,f6.2)' &
           ,rec=nrec_flow) nnd,ncell &
           ,qin(no_heat),qout(no_heat),qdiff(no_heat) &  
           ,depth(no_heat),width(no_heat),u(no_heat)
!
          if(u(no_heat).lt.0.01) u(no_heat)=0.01
          if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat 
!
          read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)' &
           ,rec=nrec_heat) ncell &
           ,dbt(no_heat),ea(no_heat) &
           ,qns(no_heat),qna(no_heat),rho &
           ,press(no_heat),wind(no_heat)
!          
          if(ncell.ne.no_heat) write(*,*) 'Heat file error',ncell,no_heat 
1400 format(i5,i6,i5,i6,4f7.0,2f6.0)
1600 format(8f10.0)
!
!  Added variable ndelta (UW_JRY_2011/03/15)
!
                     delta_n=ndelta(ncell)
! 
          qavg=0.5*(qin(no_heat)+qout(no_heat))
          qdiff(no_heat)=qdiff(no_heat)/delta_n
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
!      Tributary flow is qout from the next to the last cell
!
!
!       Read the meteorology for the last cell, but not the flow
!
        no_heat=no_heat+1 
        qin(no_heat)=qout(no_heat-1)
        qout(no_heat)=qin(no_heat)
        q_trib(nr)=qout(no_heat)       
        nrec_heat=heat_cells*(ndays-1)+no_heat
        read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)' &
         ,rec=nrec_heat) ncell &
         ,dbt(no_heat),ea(no_heat) &   
         ,qns(no_heat),qna(no_heat),rho &
         ,press(no_heat),wind(no_heat)
!
!  The flow and hydraulics for the last cell has to be 
!  modified so they do not
!  take the values of the segment to which it is tributary
!
!	            qin(ncell)=qout(ncell-1)
!	            qout(ncell)=qin(ncell-1)
        qdiff(no_heat)=0.0
        u(no_heat)=u(no_heat-1)
        depth(no_heat)=depth(no_heat-1)
        width(no_heat)=width(no_heat-1)
        dt(no_heat)=0.5*dx(ncell)/u(no_heat)
      end do
!
!     Main stem inflows and outflows for each reach first
!     Flows are cumulative and do not include tributaries if
!     tributaries are downstream of the inflow junction
!
!
!     Headwaters flow and temperature
!
!
90            continue
!
!     Begin cycling through the reaches
!
      do nr=1,nreach
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
      x_head=x_dist(nr,0)
      x_bndry=x_head-50.0
!
!     First do the reverse particle tracking
!
      ncell_write=0
      do ns=no_celm(nr),1,-1
!
!     Segment is in cell SEGMENT_CELL(NC)
!
        ncell=segment_cell(nr,ns)
        nx_s=1
        nx_part=ns
        dt_part(ns)=dt(ncell)
        dt_total=dt_part(ns)
        x_part(ns)=x_dist(nr,ns)
100 continue
!
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!
        if(dt_total.lt.dt_comp) then
          x_part(ns)=x_part(ns)+dx(segment_cell(nr,nx_part))

!
!     If the particle has started upstream from the boundary point, give it
!     the value of the boundary
!
          if(x_part(ns).ge.x_bndry) then
            x_part(ns)=x_head
            dt_part(ns)=dt(segment_cell(nr,nx_part))
            dt_total=dt_total+dt_part(ns)
!           nx_part=head_cell(nr)
            go to 200
          end if
!
!     Increment the segment counter if the total time is less than the
!     computational interval
!
          nx_s=nx_s+1
          nx_part=nx_part-1
          dt_part(ns)=dt(segment_cell(nr,nx_part))
          dt_total=dt_total+dt_part(ns)
          go to 100
        else
!
!     For the last segment of particle travel, adjust the particle location
!     such that the total particle travel time is equal to the computational
!     interval.
!
          dt_part(ns)=dt_comp-dt_total+dt_part(ns)
          x_part(ns)=x_part(ns)+u(segment_cell(nr,nx_part))*dt_part(ns)
          if(x_part(ns).ge.x_head) then
            x_part(ns)=x_head
            nx_s=nx_s-1
            dt_part(ns)=dt(head_cell(nr))
          end if
        end if
200 continue
        if(nx_part.lt.1) nx_part=1
          nstrt_elm(ns)=nx_part
          no_dt(ns)=nx_s
        end do
        n_write=0
        DONE=.FALSE.
        do ns=1,no_celm(nr)
          ncell=segment_cell(nr,ns)
!
!     Net solar radiation (kcal/meter^2/second)
!
!     qns
!
!     Net atmospheric radiation (kcal/meter^2/second)
!
!     qna
!
!     Dry bulb temperature (deg C)
!
!     dbt
!
!     Wind speed (meters/second)
!
!     wind
!
!     Factor for Bowen ratio ((deg C)^-1)
!
!     pf
!
!     Vapor pressure at given air temperature (mb)
!
!     ea
!
!     Photo period (fraction of a day.  Not used in the energy budget)
!
!     phper
!
250 continue
!
!     Now do the third-order interpolation to
!     establish the starting temperature values
!     for each parcel
!
          nseg=nstrt_elm(ns)
          npndx=1
!                       do ntrp=2,4
!                          ntest=nseg+ntrp-3
!                          if(trib(nr,ntest)) then
!                          npndx=ntrp
!                          end if
!                       end do
!
!     If starting element is the first one, then set
!     the initial temperature to the boundary value
!
          if(nseg.eq.1) then
            t0=T_head(nr)
            go to 350
          end if
!
!     Perform polynomial interpolation
!
          do ntrp=1,nterp(npndx)
            npart=nseg+ntrp+ndltp(npndx)-1
            nptest(ntrp)=npart
            xa(ntrp)=x_dist(nr,npart)
            ta(ntrp)=temp(nr,npart,n1)
          end do
          x=x_part(ns)
280 continue
!
!     Call the interpolation function
!
t0=tntrp(xa,ta,x,nterp(npndx))
ttrp=t0
300 continue
350 continue
          dt_calc=dt_part(ns)
          nncell=segment_cell(nr,nstrt_elm(ns))
!
!    Set NCELL0 for purposes of tributary input
!
          ncell0=nncell
          dt_total=dt_calc
          do nm=no_dt(ns),1,-1
            u_mps=u(nncell)/3.2808
            u_river=0.5*u_mps
            z=depth(nncell)
            call energy(t0,qsurf,A,B,nncell,u_river,qconv)
            qdot=(qsurf/(z*rfac))
            qdd=qdot*dt_calc
!
!    Added the thermal input (UW_JRY_2011/03/15)
!    Limit on delta_T is hardwired
!
            qqq=cuft_cum*qin(nncell)
            delta_T=thermal(nr,nseg)/qqq
              if(delta_T.gt.7.0) then
                delta_T=7.0
              end if
            t0=t0+qdot*dt_calc+delta_T
            if(t0.lt.0.0) t0=0.0
400 continue
!
!     Look for a tributary.
!
            q1=qin(nncell)
            ntribs=no_tribs(nncell)
            if(ntribs.gt.0.and..not.DONE) then
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(q_trib(nr_trib).lt.0.0) go to 450
                q2=q1+q_trib(nr_trib)
                t0=(q1*t0+q_trib(nr_trib)*T_trib(nr_trib))/q2
!
450 continue
                q1=qout(nncell)
              end do
              DONE=.TRUE.
            end if
            if(ntribs.eq.0.and.qdiff(nncell).gt.0) then
              q2=q1+qdiff(nncell)
              T_dist=T_head(nr)
              t0=(q1*t0+qdiff(nncell)*T_dist)/q2
              q1=q2
            end if
500 continue
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new cell
!
            if(ncell0.ne.nncell) then
              ncell0=nncell
              DONE=.FALSE.
            end if
            dt_calc=dt(nncell)
            dt_total=dt_total+dt_calc
          end do
          if (t0.lt.0.5) t0=0.5
            temp(nr,ns,n2)=t0
	    T_trib(nr)=t0
!
!   Write file 20 with all temperature output 11/19/2008
!   The temperature is now output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_2011/03/15)
!
    if(n_write.ne.ncell) then
              write (20,'(f12.6,4i6,3f8.2,f8.4,f6.1,f9.0,f8.1,f8.4)')  &                                              
                time,nd,nr,ncell,ns,t0 &
               ,T_head(nr),dbt(ncell),A,B,q1,depth(ncell),u_river
               n_write=ncell
      end if
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
    end
      SUBROUTINE ENERGY(TSURF,QSURF,A,B,ncell,u_river,qconv)
      REAL*4 LVP
      real*4 q_fit(2),T_fit(2),evrate(2)
      INCLUDE 'rbm10_VIC.h'
      data evrate/1.5e-9,0/,pf/0.640/
	  pi=3.14159
	  td=nd
	  evap_rate=evrate(1)+evrate(2)*cos(2.*pi*(td-200)/365.)
      T_fit(1)=tsurf-1.0
      T_fit(2)=tsurf+1.0
      do i=1,2
         E0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
         RB=PF*(DBT(ncell)-T_fit(i))
         LVP=597.0-0.57*T_fit(i)
!         QEVAP=1000.*LVP*EVRATE(NW,nseasn)*(WIND(ncell)+u_river)
         QEVAP=1000.*LVP*evap_rate*(WIND(ncell))
         if(qevap.lt.0.0) qevap=0.0
         QCONV=RB*QEVAP
         QEVAP=QEVAP*(E0-EA(ncell))
         QWS=6.693E-2+1.471E-3*T_fit(i)
         q_fit(i)=QNS(ncell)+QNA(ncell)-QWS-QEVAP+QCONV
         qqns=qns(ncell)
         qqna=qna(ncell)
         qqws=qws 
         qqevap=qevap
         qqconv=qconv                   
      end do
!
!     q=AT+B
!
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      qsurf=0.5*(q_fit(1)+q_fit(2))
      B=(qsurf/A)-(T_fit(1)+T_fit(2))/2.
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
      RETURN
      END

function nodays(jtime,jy0)
  dimension ndmo(12)
  data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
  jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
      return
      end
      function ndate(n,ny0)
      dimension ndmo(13)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334,365/
      ny=(n-1)/365
      njul=n-ny*365
      nm=1
 10   continue
      if(ndmo(nm+1).ge.njul) go to 50
      nm=nm+1
      go to 10
 50   continue
      nday=njul-ndmo(nm)
      nmon=100*nm
      nyear=10000*(ny0+ny)
      ndate=nyear+nmon+nday
      return
end
!
!	Third-order polynomial interpolation using Lagrange
!     polynomials.  FUNCTION is SUBROUTINE POLINT from
!     Numerial Recipes
!

FUNCTION tntrp(XA,YA,X,n)
! PARAMETER (N=4)
! DIMENSION XA(N),YA(N),C(N),D(N)
  DIMENSION XA(4),YA(4),C(4),D(4)
  NS=1
  DIF=ABS(X-XA(1))
  DO 11 I=1,N
  DIFT=ABS(X-XA(I))
  IF (DIFT.LT.DIF) THEN
    NS=I
    DIF=DIFT
  ENDIF
  C(I)=YA(I)
  D(I)=YA(I)
11 CONTINUE
  Y=YA(NS)
  NS=NS-1
  DO 13 M=1,N-1
  DO 12 I=1,N-M
  HO=XA(I)-X
  HP=XA(I+M)-X
  W=C(I+1)-D(I)
  DEN=HO-HP
  IF(DEN.EQ.0.) DEN=0.001
  DEN=W/DEN
  D(I)=HP*DEN
  C(I)=HO*DEN
12 CONTINUE
  IF (2*NS.LT.N-M)THEN
    DY=C(NS+1)
  ELSE
    DY=D(NS)
    NS=NS-1
  ENDIF
  Y=Y+DY
13 CONTINUE
  tntrp=y
  RETURN
END

function julian(iy,id)
  julian=10000*iy+id
  return
end
