c $Header: /home/CRTASS/columbia/temp_model/RCS/rbm10.f,v 2.1 2000/06/28 16:43:58 root Exp root $
c
c      PROGRAM RMAIN                                  
C
C     Dynamic river basin model for simulating water quality in 
C     branching river systems with freely-flowing river segments,
C     and river-run reservoirs.
c 
c     This version uses Reverse Particle Tracking in the Lagrangian
c     mode and Lagrangian interpolation in the Eulerian mode.
c
c     This version uses reservoir surface elevation to estimate residence times
c     dynamically rather than statically as was the case in original versions
c     of the model.  This is of interest for reservoirs such as Grand Coulee that
c     are used for flood control.
C 
C     For additional information contact:
C
C     John Yearsley
C     EPA Region 10    ES-098
C     1200 Sixth Ave
C     Seattle, WA      98101
C     (206) 553-1532
C
      character*30 advect_file,NAMEI
      INCLUDE 'rbm10.fi'
C
C     Open file containing reach data
C
      WRITE(*,2600)
      read(*,1500) namei
      OPEN(UNIT=90,FILE=namei,STATUS='OLD')
c
c     Read header information from control file
c
      do n=1,30
	   case_name(n:n)=' '
	end do
	advect_file=case_name
	read(90,1500) advect_file
      do nc=1,30
	   if(advect_file(nc:nc).ne.' ') ncadv=nc
	end do
      read(90,1500) case_name
	do nc=1,30
	   if(case_name(nc:nc).ne.' ') ncase=nc
      end do
      write(*,2900) 
      write(*,1510) case_name
c
c     Open file with hydrologic data
c
      open(unit=35,FILE=advect_file(1:ncadv)//'.advect',STATUS='old')
c
c     Open file with reservoir elevation data
c
      open(unit=36,FILE=advect_file(1:ncadv)//'.elev',STATUS='old')
c
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
      CLOSE(UNIT=30)
      CLOSE(UNIT=35)
      CLOSE(UNIT=90)
 1500 FORMAT(A30)
 1510 format(1x,a30)
 1600 FORMAT(8F10.0)
 2600 FORMAT(' NAME OF FILE CONTAINING RIVER REACH DATA')
 2700 FORMAT(' NAME OF OUTPUT DATA FILE')
 2800 format(' Name of file with geometric data')
 2900 format(' Name of file with hydrologic data')
 3000 format(' Name of file with water quality data')
      STOP
      END
      SUBROUTINE BEGIN
	character*3 plot_day
      CHARACTER*5 type
      character*10 end_mark
      character*20 seg_name(200)
      character*30 NAMEI,xplot_file
      integer start_date,end_date,trib_jnc
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
         do nn=-2,1000
            trib(n,nn)=.FALSE.
            P_var(n,nn,1)=Q_var(n)
            P_var(n,nn,2)=Q_var(n)
         end do
         do nn=1,1000
            type_res(n,nn)=0
            type_riv(n,nn)=0
            temp(n,nn,1)=0.0
         end do
      end do
	nrsrvr=0
	new_rsrvr=.TRUE.
      no_inflow=0
      no_intrib=0
      nwprov=0
C
C     Card Group I
C               
      read(90,*)
      read(90,1200) start_date,end_date
      nyear1=start_date/10000
      nyear2=end_date/10000
      nysim=nyear2-nyear1+1
      ysim=nysim
      write(*,1200) start_date,end_date
 1200 format(/8i10)
      read(90,1200) no_rch
      write(*,1200) no_rch
c
 1400    format(16i5)
 1450    format(16(2x,a3))
C                                        
C     Card Group IIb. Reach characteristics
C
      do nr=1,no_rch
         write(*,*) ' Starting to read reach ',nr
         read(90,1480) rch_name(nr),no_elm(nr)
     .        ,no_plots(nr),no_xplots(nr),nxp_year(nr)
C
C     Setup locations for compliance.  Locations are specified in each
C     Reach by River Mile now, rather than by element number (NC)
c     6/14/2001
C
         if(no_plots(nr).gt.0) then
            nplts=no_plots(nr)
            read(90,1065) (rm_plot(nr,np),np=1,nplts)
	   end if
	   do nc=1,20
	      namei(nc:nc)=' '
	   end do
	   namei=rch_name(nr)
	   do nc=1,20
            if(namei(nc:nc).eq.' ') go to 200
            nchar=nc
	   end do
  200    continue
	   write(*,9500) namei(1:nchar),nchar
 9500 format(a30,i10)
      pause
	   jchar=ncase+nchar+1
	   if(no_xplots(nr).gt.0) then
	      read(90,1400) (day_xplot(nr,nxp),nxp=1,no_xplots(nr))
	      do nxp=1,no_xplots(nr)
	         write(plot_day,'(i3)') day_xplot(nr,nxp)
	         xplot_file=namei(1:nchar)//'_'
	write(*,9500) xplot_file
	         xplot_file=xplot_file(1:nchar+1)//case_name(1:ncase)
	write(*,9500) xplot_file
			 xplot_file=xplot_file(1:jchar)//'.Day_'//plot_day
	write(*,9500) xplot_file
	pause
			 junit=100+20*(nr-1)+nxp
			 open(junit,file=xplot_file,status='unknown') 
		  end do       
         end if
 1480    format(a20,6i10)
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
C     Card Type 3. Reach description, Reach type, begin and 
C     end river mile
C  
            read(90,1600) seg_name(ne),type,rmile1(nr,ne),rmile2(nr,ne)
	.                   ,zdd
            write(*,1600) seg_name(ne),type,rmile1(nr,ne),rmile2(nr,ne)
 1600       format(a20,5x,a5,10x,3f5.0)
c     
c     Read number of computational elements, weather province, 
c     headwaters number, number of entering tributaries, Reach number
c     if the tributary is one for which temperatures are simulated
c     
            read(90,1400) no_xsctn,kwtype,kcntrl
     .           ,ntribs,nr_trib
            write(*,1400) no_xsctn,kwtype
            xsctn=no_xsctn
            delta_x=5280.*(rmile1(nr,ne)-rmile2(nr,ne))
            delta_x=delta_x/xsctn
            x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
            nc1=nc+1
c
c     Reservoir segments
c
            if(type.eq.'RSRVR') then
	         if(new_rsrvr) then
	            new_rsrvr=.FALSE.
	         end if
	         if(kcntrl.gt.0) then
	            nrsrvr=nrsrvr+1
	            no_rsrvrs=nrsrvr
	            new_rsrvr=.TRUE.

			 end if
		     read(90,1700) a_a,b_a,a_w,b_w
               do nx=1,no_xsctn
                  nc=nc+1
                  x_sctn=no_xsctn
                  type_res(nr,nc)=nrsrvr
                  a_area(nr,nc)=a_a
                  b_area(nr,nc)=b_a
                  a_width(nr,nc)=a_w
                  b_width(nr,nc)=b_w
c
c     Set the dead storage elevation for 
c
                  z_bottom(nr,nc)=zdd
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  if(kwtype.gt.nwprov) nwprov=kwtype
C     
C     Reservoir reaches
c     
                  dx(nr,nc)=delta_x
c
c     Endpoints for reservoir segments
c
                  x_dist2(nr,nc)=x_dist1(nr,nc)-dx(nr,nc)
                  rmp_ft=5280.*rm_plot(nr,jp)
                  if(rmp_ft.ge.x_dist2(nr,nc).and.   
     .                 rmp_ft.lt.x_dist1(nr,nc)) then
                     nc_plot(nr,jp)=nc
                     jp=jp+1
                  end if
                  x_dist1(nr,nc+1)=x_dist2(nr,nc)
 1700             format(8f10.0)
               end do
            end if
c     
c     River reaches
c     
            if(type.eq.'RIVER') then
               if(kwtype.gt.nwprov) nwprov=kwtype
               x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
               read(90,1700) a_a,b_a,a_w,b_w
               do nx=1,no_xsctn                
                  nc=nc+1
c
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  type_riv(nr,nc)=1
                  a_area(nr,nc)=a_a
                  b_area(nr,nc)=b_a
                  a_width(nr,nc)=a_w
                  b_width(nr,nc)=b_w
                  x_dist2(nr,nc)=x_dist1(nr,nc)-delta_x
                  rmp_ft=5280.*rm_plot(nr,jp)
                  if(rmp_ft.ge.x_dist2(nr,nc).and.   
     .                 rmp_ft.lt.x_dist1(nr,nc)) then
                     nc_plot(nr,jp)=nc
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
 1850          format(15x,i5,5x,f10.0)
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
 1851       format(a10)
            write(*,2851) end_mark
            if(end_mark(1:3).ne.'End') pause 
 2851       format(' end_mark: ',a10)
C     
         end do
         x_dist2(nr,0) = x_dist1(nr,1)
         x_dist2(nr,-1)= x_dist2(nr,0) +dx(nr,1)
         x_dist2(nr,-2)= x_dist2(nr,-1)+dx(nr,1)
         x_dist2(nr,nc+1)=x_dist2(nr,nc)-dx(nr,nc)
         no_celm(nr)=nc
	   if(nc.gt.1000) then
	      write(*,*) 'Number of computational elements exceeds 1000'
     .                 ,' in Reach - ',nr
	      pause
	   end if
         do nnc=1,no_celm(nr)
            rm1=x_dist1(nr,nnc)/5280.
            rm2=x_dist2(nr,nnc)/5280.
         end do
      end do    
c 
c     Open meteorological file
c
      write(*,*) ' nwprov = ',nwprov
      DO  nw=1,NWPROV
         NWTAPE=10+nw
         read(90,1030) namei,evrate(nw)
         write(*,*) 'evrate = ',evrate(nw),' nw = ',nw
         pause
         write(*,2700) 
c     
c     Open file with heat budget information
c     
         OPEN(UNIT=NWTAPE,FILE=NAMEI,STATUS='OLD')
         READ(NWTAPE,1025) WPNAME(nw)
         READ(NWTAPE,1027) nwpd,selev,lat,long,nwstrt,nwstop
c
c     Computational time step
c
         nw_year1=nwstrt/10000
         write(*,*) ' nwyear = ',nw_year1,nyear1,dt_comp
	   ny_skip=nyear1-nw_year1
	   do ny_w=1,ny_skip
	      do nd=1,365
	         do nwp=1,nwpd
	            read(nwtape,*) jd
	         end do
	      end do
	   end do
	   write(*,*) 'weather tape - ',nwtape,'  day - ',jd
	pause
      end do
      xwpd=nwpd
      dt_comp=86400./xwpd
c     
c     Read header on advection file and advance to starting point
c     if necessary
c     
      read(35,1042) ny_adv1,ny_adv2
      ny_adv1=ny_adv1/10000
      if(ny_adv1.lt.nyear1) then
         do ny_adv=1,nyear1-ny_adv1
            do nd_adv=1,365
               read(35,*)
               do nt_adv=1,no_intrib
                  read(35,*)
               end do
            end do
         end do
      end if
c
c     Read the reservoir elevation file header and advance to the
c     correct starting point, if necessary.
c
      read(36,1042) ny_res1,ny_res2
	write(*,*) ' Starting to read elevation file ',ny_res1,ny_res2
	pause
      ny_res1=ny_res1/10000
      if(ny_res1.lt.nyear1) then
	write(*,*) 'records to skip for elevation - ',nyear1-ny_res1
	pause
         do ny_res=1,nyear1-ny_res1
            do nd_res=1,365
               read(36,1028) jdd,(z_cntrl(nrs),nrs=1,no_rsrvrs)
            end do
         end do
      end if
      
c     
c     Set DT of boundary element equal to computational DT
c
      do nr=1,no_rch
         dt(nr,0)=dt_comp
      end do
 1020 FORMAT(A80)
 1025 format(a50)
 1027 format(5x,i5,3(5x,f5.0),2i10)
 1028 format(i5,(10f10.0))
 1030 format(a30,f10.0)
 1035 format(1x,a30,e15.5)
 1040 FORMAT(10F10.0)
 1042 FORMAT(10I10)
 1043 format(i5,2f10.0,3i5)
 1044 FORMAT(16I5)
 1045 format(i5,f10.0,a60)
 1048 FORMAT(8F10.0)
 1050 FORMAT((A30,5X,A5,10x,5F5.0))
 1060 FORMAT(2F10.0,A20/16F5.0)
 1063 FORMAT(F10.0,A20/16F5.0)
 1065 FORMAT(16F5.0)
 1080 FORMAT(A3)
 1085 format(1x,a3)
 1145 FORMAT(8F10.2)
 1152 FORMAT(6I3)
 1250 format(2i5,3f10.0)
 1500 format(i9)
 2500 FORMAT(' ENERGY BUDGET FILE FOR METEOROLOGIC PROVINCE - ',I5)
 2700 format(' energy budget file')
 3000 FORMAT(1H0,'CARD SEQUENCE ERROR IN DATA FOR REACH - ',I5)
 3500 format(' reservoir flow file')
C 
C     ******************************************************
C                         Return to RMAIN
C     ******************************************************
C
      RETURN
  900 END
      SUBROUTINE SYSTMM
      real*4 WDATA(5,7),EDATA(7),xa(4),ta(4),var(4),heat_load(5,20)
     .      ,dt_part(1000),x_part(1000),xceed(5,3,50)
      integer no_dt(1000),simyr,strt_elem(1000)
     .     ,ndltp(4),nterp(4)
      real*4 plot_data(200,2)
      INCLUDE 'rbm10.fi'
      EQUIVALENCE (EDATA(1),QNS)
      data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
      data pi/3.14159/,rfac/304.8/
	data heat_load/100*0.0/
C
      time=0.0
      n1=1
      n2=2
      nobs=0
      simyr=0
 50   continue
	nyear=nyear1+simyr
	simyr=simyr+1
      write(*,*) ' Simulation Year - ',simyr
      DO ND=1,365
         nobs=nobs+1
	   do nr=1,no_rch
	      if(no_xplots(nr).gt.0.and.nxp_year(nr).eq.nyear) then
	         do nxp=1,no_xplots(nr)
	            if(nd.eq.day_xplot(nr,nxp)) then
	                  junit=100+20*(nr-1)+nxp
	                  do nc=1,no_celm(nr)
	                     riv_mile=(x_dist1(nr,nc)+x_dist2(nr,nc))
     .					         /10560.
	                     write(junit,2500) riv_mile,temp(nr,nc,n1)
 2500 format(2f10.2)
	                  end do
	                  close(junit)
	            end if
	         end do
	      end if
	   end do	    
         DO NDD=1,NWPD
C     
C     Read weather data from files
C     
            do nw=1,nwprov
               nwr=10+nw
               READ(NWR,1028) LDUMM,(WDATA(nw,nnw),nnw=1,7)
            end do
c
     
c
c     
c     
c     begin reach computations
c     
            ind=nd
            ipd=ndd    
            day=nd
c     
c     Tributary temperatures (computed)
c     Moved here to make diurnal tributary updating possible
c     11/08/2001
c
            if(no_rch.gt.1) then
               do nr=1,no_rch-1
                  nt=trib_id(nr)
                  T_trib(nt)= temp(nr,no_celm(nr),n1)
               end do
            end if
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
            read(35,1550,end=800) jy,jobs
     .           ,(qin(nr,1),T_head(nr),nr=1,no_rch)
c     
c     Check for tributaries
c     
            if(no_intrib.gt.0) then
c     
c     Tributary flow and temperature (input)
c     
               do in=1,no_intrib
                  read(35,1600,end=800)
     .                 in_trb,q_trib(in_trb),T_trib(in_trb)
c
c     Remove comment ("c" in Column 1) for Scenario 3
c
                  if(T_trib(in_trb).gt.16.0) T_trib(in_trb)=16.0
               end do
            end if
c
c	Read the surface elevations of the reservoirs
c     
            read(36,1028,end=800) jday,(z_cntrl(nrs)
     .                           ,nrs=1,no_rsrvrs)
	        if(simyr.eq.1.and.nd.eq.1) then
	           do nr=1,no_rch
	              do nc=1,no_celm(nr)
	                 if(type_res(nr,nc).gt.0) then
	                    nrs=type_res(nr,nc)
	                    z=z_cntrl(nrs)-z_bottom(nr,nc)
	                    if(nrs.eq.7.and.z.lt.30.0) z=30.0
	                    seg_vol(nr,nc)=dx(nr,nc)
     .			    	              *a_area(nr,nc)
     .                                *exp(z*b_area(nr,nc))
	                 end if
	              end do
	           end do
	           read(36,1028) jday,(z_cntrl(nrs),nrs=1,no_rsrvrs)
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
	        nrsrvr=1
	        new_rsrvr=.TRUE.
c     
c     Begin cycling through the reaches
c     
            do nr=1,no_rch
	           ic_plot=2
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
                        dt_part(nc)=dt(nr,1)
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
     .                         +u(nr,nx_part)*dt_part(nc)
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
  280             continue
c     
c     Call the interpolation function
c
                  t0=tntrp(xa,ta,x,nterp(npndx))
	              t_prior=t0
                  var0=tntrp(xa,var,x,nterp(npndx))
 300              continue
 350              continue
                  dt_calc=dt_part(nc)
	              do nm=no_dt(nc),1,-1
                     nw=nwtype(nr,ncell)
	                 u_river=0.5*u(nr,ncell)/3.2808
                     z=depth(nr,ncell)	               
                     if(nw.gt.0) then
                        call energy(t0,qsurf,A,B,nw,u_river)
c     qdot=(a_tmp(nw)+b_tmp(nw)*t0)/(z*rfac)
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
                     q1=qin(nr,ncell)
                     q2=qin(nr,ncell+1)
                     if(.not.trib(nr,ncell)) go to 450
                     nt_trb=trib_ndx(nr,ncell)
	                 if(q_trib(nt_trb).lt.0.0) go to 450
                     t1=temp(nr,ncell-1,n1)
                     t00=t0
                     t0=(q1*t0+q_trib(nt_trb)*T_trib(nt_trb))
     .                    /q2
 450                 continue
                     ncell=ncell+1
                     dt_calc=dt(nr,ncell)
                  end do
                  temp(nr,nc,n2)=t0
                  P_var(nr,nc,n2)=var0
c
c     End of computational element loop
c
               end do             
c
c     Set up the output
c
               do np=1,no_plots(nr)
                  nc_plt=nc_plot(nr,np)
c
                  heat_load(nr,np)=heat_load(nr,np)
     .                                +qin(nr,nc_plt+1)
     .                                *temp(nr,nc_plt,n2)/1000.
                  plot_data(np,1)=temp(nr,nc_plt,n2)
                  plot_data(np,2)=P_var(nr,nc_plt,n2)
	if(plot_data(np,2).lt.0.0) plot_data(np,2)=0.0
                  t_test1=plot_data(np,1)-sqrt(plot_data(np,2))
                  t_test2=plot_data(np,1)
                  t_test3=plot_data(np,1)+sqrt(plot_data(np,2))
c
c     Test for exceedances of the benchmark
c
                  if(t_test1.gt.20.) 
     .                 xceed(nr,1,np)=xceed(nr,1,np)+1.
                  if(t_test2.gt.20.) 
     .                 xceed(nr,2,np)=xceed(nr,2,np)+1.
                  if(t_test3.gt.20.) 
     .                 xceed(nr,3,np)=xceed(nr,3,np)+1.
c
               end do
c              if(no_plots(nr).gt.0) then
c                 do np=1,no_plots(nr)	           
c                     nc_plt=nc_plot(nr,np)
c                     plot_data(np,1)=temp(nr,nc_plt,n1)
c                     plot_data(np,2)=qin(nr,nc_plt+1)/1000.
c                     plot_data(np,2)=u(nr,nc_plt)
c                  end do
c	           end if
               time=nd
	           xd=nd
               xdd=ndd
               clock=(xdd-0.5)*dt_comp
               time=time+clock/86400.
	           year1=nyear1-1900
               time=year1+simyr-1.0+time/365.
               zero=0.0
c
c     The plotting points are now also the compliance points
c
               if(no_plots(nr).gt.0) then
                  itplot=40+nr
	            write(itplot,4700) time,xd
     .            ,(rm_plot(nr,np),plot_data(np,1)
     .            ,np=1,no_plots(nr))
               end if
c
c     End of reach loop
c
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
c     
c     End of weather period loop (NDD=1,NWPD)
c     
 4650       format(16x,12(6x,f6.0,6x))
 4700       format(f10.4,f6.0,15(f6.1,f8.3))
 4750       format(f10.4,10(i4,f8.0))
         end do
C     
c     End of main loop (ND=1,365)
c     
      end do
c     
c     Check if there are years remaining to be simulated.  If so,
c     start at the top (STATEMENT 50)
c
      if(simyr.lt.nysim) go to 50
  800 continue
      ysim=nysim
      do nr=1,no_rch
         write(39,3800) nr
 3800    format(' Reach No. - ',i5)
         xobs=365*ysim
         nplots=no_plots(nr)
         do np=1,nplots
            write(39,3900) rm_plot(nr,np)
     .           ,(xceed(nr,nt,np)/xobs,nt=1,3)
	        write(39,4000) heat_load(nr,np)
 3900       format(f6.1,3f10.4)
 4000 format(6x,e15.5)      
         end do
      end do
c     
c     FORMAT statements
c     
 1020 FORMAT(A80)
 1025 format(a50)
 1027 format(5x,i5,3(5x,f5.0),2i10)
 1028 format(i5,(10f10.0))
 1042 format(8i10)
 1045 format(i5,f10.0,a60)
c 
c     ******************************************************
c                        return to rmain
c     ******************************************************
c
  950 return
      end
      subroutine balance
      include 'rbm10.fi'
      dt_max=-10000.
      do nr=1,no_rch
         do nc=1,no_celm(nr)
		    q1=qin(nr,nc)
			q2=q1
			if(trib(nr,nc)) then
               nt_trb=trib_ndx(nr,nc)
               q2=q1+q_trib(nt_trb)
            end if
c			
            if(type_res(nr,nc).gt.0) then
			   nrs=type_res(nr,nc)
			   fctr=z_cntrl(nrs)-z_bottom(nr,nc)
	           if(nrs.eq.7.and.fctr.lt.30.0) fctr=30.0
               x_area=a_area(nr,nc)*exp(fctr*b_area(nr,nc))
               x_wide=a_width(nr,nc)*exp(fctr*b_width(nr,nc))
			   depth(nr,nc)=x_area/x_wide
			   u(nr,nc)=0.5*(q1+q2)/x_area
			   dt(nr,nc)=dx(nr,nc)/u(nr,nc)
	           res_time=(dx(nr,nc)*x_area/(0.5*(q1+q2)))
		  end if
            if(type_riv(nr,nc).gt.0) then
			   fctr=q1
               x_area=a_area(nr,nc)*(fctr**b_area(nr,nc))
               x_wide=a_width(nr,nc)*(fctr**b_width(nr,nc))
               depth(nr,nc)=x_area/x_wide
               u(nr,nc)=0.5*(q1+q2)/x_area
               dt(nr,nc)=dx(nr,nc)/u(nr,nc)
	           res_time=(dx(nr,nc)*x_area/(0.5*(q1+q2)))
			end if
            if(dt(nr,nc).gt.dt_comp) then
               nmax=nc
               dt_max=dt(nr,nc)
            end if
            qin(nr,nc+1)=q2
         end do
         dt(nr,no_celm(nr)+1)=dt_comp
c
c    Assign flow to headwaters reaches that are also tributaries
c
         nt_trb=trib_id(nr)
		 q_trib(nt_trb)=q2
      end do
      return
      end
      SUBROUTINE ENERGY(TSURF,QSURF,A,B,NW,u_river)
      REAL*4 LVP
      real*4 q_fit(2),T_fit(2)
      INCLUDE 'rbm10.fi'
      T_fit(1)=tsurf-1.0
      T_fit(2)=tsurf+1.0
      do i=1,2
         E0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
         RB=PF*(DBT-T_fit(i))
         LVP=597.0-0.57*T_fit(i)
         QEVAP=1000.*LVP*EVRATE(NW)*(WIND+u_river)
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







