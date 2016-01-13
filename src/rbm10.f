c     
c
c      PROGRAM RMAIN                                  
C
C     Dynamic river basin model for simulating water temperature in 
C     a branching river systems with freely-flowing river segments,
C     a river-run reservoirs. The numerical scheme is based on Reverse
c     Particle Tracking in the Lagrangian mode and Lagrangian
c     interpolation in the Eulerian mode.
c
c     Updated from original distribution  10/09/2002                    
C 
C     For additional information contact:
C
C     John Yearsley
C     EPA Region 10    ES-098
C     1200 Sixth Ave
C     Seattle, WA      98101
C     (206) 553-1532
C
      character*30 NAMEI
      INCLUDE 'rbm10.fi'
C
C     Open file containing river reach data
C
      WRITE(*,2600)
	read(*,1000) namei
      OPEN(UNIT=90,FILE=namei,STATUS='OLD')
c
c     Read header information from control file
c	Read advected source file name from control file
c
	read(90,*)
	read(90,*)
      read(90,1500) namei
      write(*,2900) 
      write(*,1510) namei
c
c     Open advected source file
c	This file contains hydrologic data for the 3 main river stems
c	and 12 advected point sources
c
      open(unit=35,FILE=namei,STATUS='old')
C
C     Call systems programs to get started
C
C     SUBROUTINE BEGIN reads control file, sets up topology and 
C     important properties of reaches
C
      CALL BEGIN
C
C     SUBROUTINE SYSTMM performs the simulations
C
      CALL SYSTMM
C
C     Close files after simulation is complete
C
      CLOSE(UNIT=35)
      CLOSE(UNIT=90)
C
C	FORMAT Statements for the Main Program
C
 1000 FORMAT(a30)
 1500 FORMAT(7x,A30)
 1510 FORMAT(1x,a30)
 1600 FORMAT(7x,8F10.0)
 2600 FORMAT(' NAME OF FILE CONTAINING RIVER REACH DATA')
C 2700 FORMAT(' NAME OF OUTPUT DATA FILE')
C 2800 format(' Name of file with geometric data')
 2900 format(' Name of file with hydrologic data')
C 3000 format(' Name of file with water quality data')
      STOP
      END

***********************************************************************
***********************************************************************
      SUBROUTINE BEGIN
	character*1 reach_ext
      CHARACTER*5 segtype
      character*10 end_mark
	character*13 plot_file
      character*20 dummy_file,seg_name(200)
      character*30 NAMEI
      integer trib_jnc
	integer*4 nrm_plot,nx_dist1,nx_dist2
      real*4 lat,long
      INCLUDE 'rbm10.fi'
      data nplot/200*0/,rm_trib/50*-99./,rm_plot/100*-99./
C
C     Initialize arrays and constants.  Filter parameters are set
C     for prediction mode (R_var(n)=0.0).  System variance (Q_var(n))
C     are set to values the same as estimates reported in the 
C     December 1999 Final Draft report of the Columbia River Temperature
C     Assessment
C
      Q_var(1)=0.007
      Q_var(2)=0.007
      Q_var(3)=0.002
      Q_var(4)=0.000
      Q_var(5)=0.000
      do n=1,5
         R_var(n)=0.0
      end do
      do n=1,5
         no_plots(n)=0
         do nn=-2,600
            trib(n,nn)=.FALSE.
            P_var(n,nn,1)=Q_var(n)
            P_var(n,nn,2)=Q_var(n)
         end do
         do nn=1,600
            type_res(n,nn)=.FALSE.
            type_riv(n,nn)=.FALSE.
            temp(n,nn,1)=0.0
         end do
      end do
      no_inflow=0
      no_intrib=0
      nwprov=0
C
C     Section 1 - General Model Information
C               
      read(90,*)
C
C	Reading start and end dates for the model simulation
C
      read(90,1200) start_date,end_date
	total_steps=nodays(end_date,start_date)+1
      nyear1=start_date/10000
      nyear2=end_date/10000
      nysim=nyear2-nyear1+1
      ysim=nysim
      write(*,1200) start_date,end_date
C
C	Reading the number of river reaches in the model
C
      read(90,1200) no_rch
      write(*,1200) no_rch
C                                        
C     Section 2 - Reach Information
C
	read(90,*)
      do nr=1,no_rch
         write(*,*) ' Starting to read reach ',nr
         read(90,*)
	   read(90,1480) rch_name(nr),no_elm(nr)
     .        ,no_plots(nr)
C
C     Setup locations for output.  Locations are specified in each
C     Reach by River Mile now, rather than by element number (NC)
C
         if(no_plots(nr).gt.0) then
	      write(reach_ext,'(i1)') nr
	      plot_file='Reach_#'//reach_ext
	      dummy_file=adjustl(rch_name(nr))
	      nchar=0
	      do n=1,8
	        if(dummy_file(n:n).eq.' ') go to 100
	        nchar=n
	      end do
  100       continue
  	      if(nchar.ne.0) plot_file=dummy_file(1:nchar)//'.plot'
	      nunit=40+nr
	      open(nunit,file=plot_file,status='unknown')
            nplts=no_plots(nr)
            read(90,1065) (rm_plot(nr,np),np=1,nplts)
         end if
c
c     NC is a temporary element counter
c
         nc=0
c
c     JP is a temporary plot location counter
c
         jp=1
c
c     Reading Reach Element information
c
         do ne=1,no_elm(nr)
          write(*,*) no_elm(nr),ne
C     
C     Read in the Segment description, Segment type and the beginning and 
C	ending river mile for each segment
C
		read(90,*)  
          read(90,1600) seg_name(ne),segtype,rmile1(nr,ne),rmile2(nr,ne)
          write(*,1600) seg_name(ne),segtype,rmile1(nr,ne),rmile2(nr,ne)
c     
c     Read number of computational elements, weather province, 
c     headwaters number, number of entering tributaries, Reach number
c     if the tributary is one for which temperatures are simulated
c     
            read(90,1400) no_xsctn,kwtype,khead
     .           ,ntribs,nr_trib
            write(*,1400) no_xsctn,kwtype
            x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
            nc1=nc+1
c
c     Reservoir segments
C	Records number 23 and 24
c
            if(segtype.eq.'RSRVR') then
               read(90,1700)
     .              rel_vol,rel_area,del_x
               write(*,1700) 
     .              rel_vol,rel_area,del_x
               do nx=1,no_xsctn
                  nc=nc+1
                  x_sctn=no_xsctn
                  type_res(nr,nc)=.TRUE.
c
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  if(kwtype.gt.nwprov) nwprov=kwtype
C     
C     Reservoir reaches
c     
                  dx(nr,nc)=5280.*del_x/x_sctn
c
c     Cross-sectional area and volume
c
                  res_area(nr,nc)=43560.*rel_area/x_sctn
                  res_vol(nr,nc)=43560.*rel_vol/x_sctn
c
c     Endpoints for reservoir segments
c
                  x_dist2(nr,nc)=x_dist1(nr,nc)-dx(nr,nc)
                  rmp=5280.*rm_plot(nr,jp)
	            nrm_plot=rmp
	            nx_dist1=x_dist1(nr,nc)
	            nx_dist2=x_dist2(nr,nc)
                  if(nrm_plot.ge.nx_dist2.and.   
     .                 nrm_plot.lt.nx_dist1) then
                     nc_plot(nr,jp)=nc
                     jp=jp+1
                  end if
                  x_dist1(nr,nc+1)=x_dist2(nr,nc)
               end do
            end if
c     
c     River reaches
C	Records number 21 and 22
c     
            if(segtype.eq.'RIVER') then
               if(kwtype.gt.nwprov) nwprov=kwtype
               xsctn=no_xsctn
               delta_x=5280.*(rmile1(nr,ne)-rmile2(nr,ne))
               delta_x=delta_x/xsctn
               x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
               read(90,1700) a_a,b_a,a_w,b_w
               do nx=1,no_xsctn                
                  nc=nc+1
c
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  type_riv(nr,nc)=.TRUE.
                  a_area(nr,nc)=a_a
                  b_area(nr,nc)=b_a
                  a_width(nr,nc)=a_w
                  b_width(nr,nc)=b_w
                  x_dist2(nr,nc)=x_dist1(nr,nc)-delta_x
                  rmp_ft=5280.*rm_plot(nr,jp)
                  if(rmp_ft.ge.x_dist2(nr,nc).and.   
     .                 rmp_ft.lt.x_dist1(nr,nc)) then
                     nc_plot(nr,jp)=nc
                     write(*,*) 'river',nr,nc,jp,nc_plot(nr,jp)
                     jp=jp+1
                  end if
                  x_dist1(nr,nc+1)=x_dist2(nr,nc)
                  dx(nr,nc)=delta_x
               end do
            end if
            nc2=nc
c
c     Check to see if there is a tributary in this computational element.
c     If so, define pointers for later use.
c	Records 22 and 24
c
            if(ntribs.gt.0) then
               no_inflow=no_inflow+1
               if(nr_trib.eq.0) then 
                  no_intrib=no_intrib+1
               else
                  trib_id(nr_trib)=no_inflow
               end if
               read(90,*)
               read(90,1850) trib_jnc
     .              ,rm_trib(no_inflow)
               x_trib(no_inflow)=5280.*rm_trib(no_inflow)
               x_t=5280.*rm_trib(no_inflow)+1.
               do nnc=nc1,nc2
                  xdist1=x_dist1(nr,nnc)
                  xdist2=x_dist2(nr,nnc)
                  if(x_t.le.x_dist1(nr,nnc)
     .                 .and.x_t.gt.x_dist2(nr,nnc)) then
                     trib(nr,nnc)=.TRUE.
                     trib_ndx(nr,nnc)=no_inflow
                  end if
               end do
            end if
            read(90,1851) end_mark
            write(*,2851) end_mark
            if(end_mark(1:3).ne.'End') pause 
         end do
         x_dist2(nr,0) = x_dist1(nr,1)
         x_dist2(nr,-1)= x_dist2(nr,0) +dx(nr,1)
         x_dist2(nr,-2)= x_dist2(nr,-1)+dx(nr,1)
         x_dist2(nr,nc+1)=x_dist2(nr,nc)-dx(nr,nc)
         no_celm(nr)=nc
         do nnc=1,no_celm(nr)
            rm1=x_dist1(nr,nnc)/5280.
            rm2=x_dist2(nr,nnc)/5280.
         end do
      end do    
c 
c     Section 3 - Read meteorological file names and evaporation coefficients
c
      write(*,*) ' nwprov = ',nwprov
	read(90,*)
      DO  nw=1,NWPROV
         NWTAPE=10+nw
         read(90,1030) namei,evrate(nw)
         write(*,2700)
c     
c     Open file with heat budget information
c     
         OPEN(UNIT=NWTAPE,FILE=NAMEI,STATUS='OLD')
         READ(NWTAPE,*)
	   READ(NWTAPE,*)
	   READ(NWTAPE,1025) WPNAME(nw)
         READ(NWTAPE,1027) nwpd,selev,lat,long,nwstrt,nwstop
c
c     Skip to the proper starting point, if necessary
c
         nw_skip=nodays(start_date,nwstrt)
	   do nskip=1,nw_skip
	       do nhskip=1,nwpd
	           read(nwtape,*)
	       end do
	   end do
c
c     Computational time step
c
         nw_year1=nwstrt/10000
         write(*,*) ' nwyear = ',nw_year1,nyear1,dt_comp
      end do
      xwpd=nwpd
      dt_comp=86400./xwpd
c     
c     Section 4 - Read header on advection file and advance to 
c	starting point if necessary
c     
	read(35,*)
	read(35,*)
      read(35,1042) ny_adv1,ny_adv2
      nadv_skip=nodays(start_date,ny_adv1)
      if(nadv_skip.gt.0) then
         do nskip=1,nadv_skip
               read(35,*)
               do nt_adv=1,no_intrib
                  read(35,*)
               end do
         end do
      end if
c     
c     Set DT of boundary element equal to computational DT
c
      do nr=1,no_rch
         dt(nr,0)=dt_comp
      end do
C
C	FORMAT Statements for the BEGIN Subroutine
C
 1025 FORMAT(a50)
 1027 FORMAT(5x,i5,3(5x,f5.0),2i10)
 1030 FORMAT(7X,a30,f10.0)
 1042 FORMAT(8I10)
 1065 FORMAT(7X,16F5.0)
 1200 FORMAT(/7X,8i10)
 1400 FORMAT(7X,16i5)
 1480 FORMAT(7X,a20,6i10)
 1600 FORMAT(7X,a20,5x,a5,10x,3f5.0)
 1700 FORMAT(7X,8f10.0)
 1850 FORMAT(7X,15x,i5,5x,f10.0)
 1851 FORMAT(7X,a10)
 2700 FORMAT(' energy budget file')
 2851 FORMAT(' end_mark: ',a10)
C 1020 FORMAT(A80)
C 1028 FORMAT(i5,7f10.0)
C 1035 FORMAT(1x,a30,e15.5)
C 1040 FORMAT(8F10.0)
C 1043 FORMAT(i5,2f10.0,3i5)
C 1044 FORMAT(16I5)
C 1045 FORMAT(i5,f10.0,a60)
C 1048 FORMAT(8F10.0)
C 1050 FORMAT((A30,5X,A5,10x,5F5.0))
C 1060 FORMAT(2F10.0,A20/16F5.0)
C 1063 FORMAT(F10.0,A20/16F5.0)
C 1080 FORMAT(A3)
C 1085 FORMAT(1x,a3)
C 1145 FORMAT(8F10.2)
C 1152 FORMAT(6I3)
C 1250 FORMAT(2i5,3f10.0)
C 1450 FORMAT(16(2x,a3))
C 1500 FORMAT(i9)
C 2500 FORMAT(' ENERGY BUDGET FILE FOR METEOROLOGIC PROVINCE - ',I5)
C 3000 FORMAT(1H0,'CARD SEQUENCE ERROR IN DATA FOR REACH - ',I5)
C 3500 FORMAT(' reservoir flow file')
C 
C     ******************************************************
C                         Return to RMAIN
C     ******************************************************
C
      RETURN
  900 END

***********************************************************************
***********************************************************************
      SUBROUTINE SYSTMM
      real*4 WDATA(5,7),EDATA(7),xa(4),ta(4),var(4)
     .      ,dt_part(600),x_part(600)
      integer no_dt(600),simyr,strt_elem(600)
     .     ,ndltp(4),nterp(4)
      real*4 plot_data(200,2)
      INCLUDE 'rbm10.fi'
      EQUIVALENCE (EDATA(1),QNS)
      data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
      data pi/3.14159/,rfac/304.8/
C
      time=0.0
      n1=1
      n2=2
      nobs=0
      simyr=0
 50   continue
      simyr=simyr+1
	time_step=0
      write(*,*) ' Simulation Year - ',simyr
      DO ND=1,365
         nobs=nobs+1
	   time_step=time_step+1
         DO NDD=1,NWPD
C     
C     Section 5 - Read weather data from files if time period is correct
C     
            do nw=1,nwprov
               nwr=10+nw
               READ(NWR,1028) LDUMM,(WDATA(nw,nnw),nnw=1,7)
            end do
c     
c     begin reach computations
c     
            ind=nd
            ipd=ndd    
            day=nd
c     
c     Tributary temperatures (computed)
c     
            if(no_rch.gt.1) then
               do nr=1,no_rch-1
                  nt=trib_id(nr)
                  T_trib(nt)= temp(nr,no_celm(nr),n1)
               end do
            end if
c     
            if(ndd.ne.1) go to 90
c     
c     Read flow and water quality
c     
c     Main stem inflows and outflows for each reach first
c     Flows are cumulative and do not include tributaries if
c     tributaries are downstream of the inflow junction
c     
c     
c     Headwaters flow and temperature
c   
            pi=3.14159
            ttme=nd
            read(35,1550) jy,jobs
     .           ,(qin(nr,1),T_head(nr),nr=1,no_rch)
c     
c     Check for tributaries
c     
            if(no_intrib.gt.0) then
c     
c     Tributary flow and temperature (input)
c     
               do in=1,no_intrib
                  read(35,1600)
     .                 in_trb,q_trib(in_trb),T_trib(in_trb)
c
c     Remove comment ("c" in Column 1) for Scenario 3
c
c                  if(T_trib(in_trb).gt.16.0) T_trib(in_trb)=16.0
               end do
            end if
c
c     Call the flow balance subroutine to set up the
c     system hydraulics
c
            call balance
c
 1400       format(4i5,4f10.0)
 1500       format(i5,i10,10F10.0)
 1510       format(8f10.5)
 1550       format(2i5,10f10.0)
 1600       format(i5,f10.0,f5.0)
 90         continue
            nsmpl=1
c     
c     Begin cycling through the reaches
c     
            do nr=1,no_rch
               nc=no_celm(nr)
               qsum=qin(nr,1)
               temp(nr,0,n1)=T_head(nr)
               temp(nr,-1,n1)=T_head(nr)
               temp(nr,-2,n1)=T_head(nr)
               temp(nr,nc+1,n1)=temp(nr,nc,n1)
               P_var(nr,nc+1,n1)=P_var(nr,nc,n1)
               x_head=x_dist1(nr,1)
               x_bndry=x_head-5.0
c     
c     First do the reverse particle tracking
c     
               do nc=no_celm(nr),1,-1
                  nx_s=1
                  nx_part=nc
c                  dt_part(nc,nx_s)=dt(nr,nx_part)
c                  dt_total=dt_part(nc,nx_s)
                  dt_part(nc)=dt(nr,nx_part)
                  dt_total=dt_part(nc)
                  x_part(nc)=x_dist2(nr,nx_part)
 100              continue
c
c     Determine if the total elapsed travel time is equal to the
c     computational interval
c
                  if(dt_total.lt.dt_comp) then
                     x_part(nc)=x_part(nc)+dx(nr,nx_part)
c
c     If the particle has started upstream from the boundary point, give it
c     the value of the boundary
c
                     if(x_part(nc).ge.x_bndry) then
                        x_part(nc)=x_head
                        go to 200
                     end if
c
c     Increment the segment counter if the total time is less than the
c     computational interval
c
                     nx_s=nx_s+1
                     nx_part=nx_part-1
                     dt_part(nc)=dt(nr,nx_part)
                     dt_total=dt_total+dt_part(nc)
                     go to 100
                  else
c
c     For the last segment of particle travel, adjust the particle location
c     such that the total particle travel time is equal to the computational 
c     interval.
c
                     dt_part(nc)
     .                    =dt_comp-dt_total+dt_part(nc)
                     x_part(nc)=x_part(nc)
     .                    +u(nr,nx_part)*dt_part(nc)
                     if(x_part(nc).ge.x_head) then
                        x_part(nc)=x_head
                        nx_s=nx_s-1
                        dt_part(nc)=dt(nr,1)
                     end if
 3700                format(f10.1)
 3710                format(2i5,3f10.1)
                  end if
 200              continue
                  if(nx_part.lt.1) nx_part=1
                  strt_elem(nc)=nx_part
                  no_dt(nc)=nx_s
               end do
               do nc=1,no_celm(nr)
c     
c     read meteorological data from the appropriate file
c     
                  iwr=nwtype(nr,nc)
                  if(iwr.eq.0) go to 250
c    
c     Correspondence table for energy budget terms
c     
c     Net solar radiation (kcal/meter^2/second)
c     
c     qns=wdata(iwr,1)
c     
c     Net atmospheric radiation (kcal/meter^2/second)
c     
c     qna=wdata(iwr,2)
c     
c     Dry bulb temperature (deg C)
c     
c     dbt=wdata(iwr,3)
c     
c     Wind speed (meters/second)
c     
c     wind=wdata(iwr,4)
c     
c     Factor for Bowen ratio ((deg C)^-1)
c     
c     pf=wdata(iwr,5) 
c     
c     Vapor pressure at given air temperature (mb)
c     
c     ea=wdata(iwr,6)
c     
c     Photo period (fraction of a day.  Not used in the energy budget)
c     
c     phper=wdata(iwr,7)
c     
                  do nnw=1,7
                     edata(nnw)=wdata(iwr,nnw)
                  end do
 250              continue
c     
c     Now do the third-order interpolation to
c     establish the starting temperature values
c     for each parcel
c     
                  ncell=strt_elem(nc)
                  npndx=1
                  do ntrp=2,4
                     ntest=ncell+ntrp-3
                     if(trib(nr,ntest)) then
                        npndx=ntrp
                     end if
                  end do
c                  ndltp=-2
c     
c     If starting element is the first one, then set
c     the initial temperature to the boundary value
c     
                  if(ncell.eq.1) then
                     t0=T_head(nr)
                     go to 350
                  end if
c     
c     Perform polynomial interpolation
c     
                  do ntrp=1,nterp(npndx)
                     npart=ncell+ntrp+ndltp(npndx)-1
                     xa(ntrp)=x_dist2(nr,npart)
                     ta(ntrp)=temp(nr,npart,n1)
                     var(ntrp)=P_var(nr,npart,n1)
                  end do
                  x=x_part(nc)
c     
c     Call the interpolation function
c
                  t0=tntrp(xa,ta,x,nterp(npndx))
                  var0=tntrp(xa,var,x,nterp(npndx))
 300              continue
 350              continue
                  dt_calc=dt_part(nc)
                  do nm=no_dt(nc),1,-1
                     nw=nwtype(nr,ncell)
                     z=depth(nr,ncell)
                     if(nw.gt.0) then
                        call energy(t0,qsurf,A,B,nw)
                        qdot=qsurf/(z*rfac)
                     else
c
c     Conservative properties if there are no weather stations available
c
                        qdot=0.0
                     end if
                     t0=t0+qdot*dt_calc
                     if(t0.lt.0.0) t0=0.0
                     phi=(1.+A/(z*rfac)*dt_calc)
                     var0=phi*var0*phi+Q_var(nr)
 400                 continue
c
c     Look for a tributary.  If none, skip to STATEMENT 450
c
                     if(.not.trib(nr,ncell)) go to 450
                     q1=qin(nr,ncell)
                     q2=qin(nr,ncell+1)
                     nt_trb=trib_ndx(nr,ncell)
                     t0=(q1*t0+q_trib(nt_trb)*T_trib(nt_trb))
     .                    /q2
 450                 continue
                     ncell=ncell+1
                     dt_calc=dt(nr,ncell)
                  end do
                  temp(nr,nc,n2)=t0
                  P_var(nr,nc,n2)=var0
               end do
               
 260           continue
c
c     Set up the output
c
               do np=1,no_plots(nr)
                  nc_plt=nc_plot(nr,np)
                  plot_data(np,1)=temp(nr,nc_plt,n2)
                  plot_data(np,2)=P_var(nr,nc_plt,n2)
c                  plot_data(np,2)=0.0
               end do
               time=nd
               xdd=ndd
               clock=(xdd-0.5)*dt_comp
               time=time+clock/86400.
               time=nyear1+simyr-1.0+time/365.
               if(no_plots(nr).gt.0) then
                  itplot=40+nr
c                  write(itplot,4700) time,(nc_plot(nr,np)
                  write(itplot,4700) time,(rm_plot(nr,np)
     .                 ,plot_data(np,1),plot_data(np,2)
     .                 ,np=1,no_plots(nr))
               end if
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
c     
c     End of weather period loop
c     
c 4700       format(f10.4,8(i4,f6.1,f6.3))
 4700       format(f10.4,8(f5.0,f6.1,f6.3))
 4750       format(f10.4,10(i4,f8.0))
         end do
C     
c     End of main loop
      if(time_step.ge.total_steps) go to 900
c     
      end do
c     
c     Check if there are years remaining to be simulated.  If so,
c     start at the top (STATEMENT 50)
c
c      if(simyr.lt.nysim) go to 50
      go to 50
  900 continue
c     
c     FORMAT statements
c     
 1020 FORMAT(A80)
 1025 format(a50)
 1027 format(5x,i5,3(5x,f5.0),2i10)
 1028 format(i5,7f10.0)
 1042 format(8i10)
 1045 format(i5,f10.0,a60)
c 
c     ******************************************************
c                        return to rmain
c     ******************************************************
c
  950 return
      end

***********************************************************************
***********************************************************************
      subroutine balance
      include 'rbm10.fi'
      dt_max=-10000.
      do nr=1,no_rch
         do nc=1,no_celm(nr)
            qsum=qin(nr,nc)
            if(type_res(nr,nc)) then
               s_area=res_area(nr,nc)
               volume=res_vol(nr,nc)
               depth(nr,nc)=volume/s_area
               dt(nr,nc)=volume/qin(nr,nc)
               u(nr,nc)=dx(nr,nc)/dt(nr,nc)
            end if
            if(type_riv(nr,nc)) then
               x_area=a_area(nr,nc)*(qsum**b_area(nr,nc))
               x_wide=a_width(nr,nc)*(qsum**b_width(nr,nc))
               depth(nr,nc)=x_area/x_wide
               u(nr,nc)=qsum/x_area
               dt(nr,nc)=dx(nr,nc)/u(nr,nc)
            end if
            if(trib(nr,nc)) then
               nt_trb=trib_ndx(nr,nc)
               qsum=qsum+q_trib(nt_trb)
            end if
            if(dt(nr,nc).gt.dt_comp) then
               nmax=nc
               dt_max=dt(nr,nc)
            end if
            qin(nr,nc+1)=qsum
         end do
         dt(nr,no_celm(nr)+1)=dt_comp
         nt_trb=trib_id(nr)
      if(nt_trb.gt.0) q_trib(nt_trb)=qsum
      end do
      return
      end

***********************************************************************
***********************************************************************
      SUBROUTINE ENERGY(TSURF,QSURF,A,B,NW)
      REAL*4 LVP
      real*4 q_fit(2),T_fit(2)
      INCLUDE 'rbm10.fi'
      T_fit(1)=tsurf-1.0
      T_fit(2)=tsurf+1.0
      do i=1,2
         E0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
         RB=PF*(DBT-T_fit(i))
         LVP=597.0-0.57*T_fit(i)
         QEVAP=1000.*LVP*EVRATE(NW)*WIND
         if(qevap.lt.0.0) qevap=0.0
         QCONV=RB*QEVAP
         QEVAP=QEVAP*(E0-EA)
         QWS=6.693E-2+1.471E-3*T_fit(i)
         q_fit(i)=QNS+QNA-QWS-QEVAP+QCONV
      end do
c
c     q=AT+B
c
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))
     .     /(T_fit(1)-T_fit(2))
      qsurf=0.5*(q_fit(1)+q_fit(2))
C 
C     ******************************************************
C               Return to Subroutine RIVMOD
C     ******************************************************
C
      RETURN
      END

***********************************************************************
***********************************************************************
      function nodays(jtime,jy0)
      dimension ndmo(12)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      jy1=(jtime/10000)
      jrem1=(jtime-jy1*10000)
      jm1=jrem1/100
      jd1=jrem1-jm1*100 
      nd1=365*jy1+ndmo(jm1)+jd1
      jy2=(jy0/10000)
      jrem2=(jy0-jy2*10000)
      jm2=jrem2/100
      jd2=jrem2-jm2*100 
      nd2=365*jy2+ndmo(jm2)+jd2
	nodays=nd1-nd2
      return
      end 

***********************************************************************
***********************************************************************
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

***********************************************************************
***********************************************************************
c
c	Third-order polynomial interpolation using Lagrange
c     polynomials.  FUNCTION is SUBROUTINE POLINT from
c     Numerial Recipes
c
      FUNCTION tntrp(XA,YA,X,n)
c      PARAMETER (N=4) 
c      DIMENSION XA(N),YA(N),C(N),D(N)
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
11    CONTINUE
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
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
	tntrp=y
      RETURN
      END

