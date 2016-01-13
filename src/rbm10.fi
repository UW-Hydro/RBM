C
C     UNDIMENSIONED INTEGER VARIABLES
C
      COMMON/BLOCK1/
     *  NDAYS,NDPRNT,NFI,ND,NDD,NELM,NN,NNPSEG
     * ,nobs,NR,NRR,no_intrib,ntrb,nysim,nsmpl
     * ,no_inflow,NWR,NWMOD,NWPD,nyear1,nyear2
     * ,n1,n2,no_rch,no_rsrvrs
C
C     UNDIMENSIONED FLOATING POINT VARIABLES
C
      COMMON/BLOCK2/dt_comp
     * ,NCONST,PLOT
     * ,QSUM,time
     * ,XTITLE,ysim
C
C     BLOCKF.COM - CONTAINS ARRAYS FOR RIVER SYSTEM CONFIGURATION
C                  HYDROLOGY IN RNGKMOD
C                  5/14/79
C
      COMMON/BLOCK3/
     *  no_celm(5),no_elm(5),no_plots(5),nc_plot(5,20),NPLOT(200)
     * ,cntrl_pnt(50),rm_plot(5,20),rm_trib(50),trib(5,-2:1000)
     * ,trib_id(5),trib_ndx(5,1000),type_riv(5,1000),type_res(5,1000)
     * ,x_trib(50)
C
C     BLOCKG.COM - CONTAINS ARRAYS FOR HYDROLOGY IN RNGKMOD
C                  5/14/79
C
      COMMON/BLOCK4/qin(5,1000),q_trib(50)
C
C     BLOCKH.COM - CONTAINS ARRAYS FOR CONCENTRATIONS AND RATE
C                  CONSTANTS FOR RNGKMOD.
C                  5/14/79
C
      COMMON/BLOCK5/ T_1,temp(5,-2:1000,2),T_trib(50),T_head(10)
     *	,T_nnpnt(1000)

C     BLOCKI.COM - CONTAINS ARRAYS WITH DEPTH AND VELOCITY
C                  CHARACTERISTICS AND RIVER MILE INDICES
C                  FOR RNGKMOD.
C
      COMMON/BLOCK6/dx(5,1000),ELEV(5,1000),depth(5,1000)
     *	,dt(5,0:1000),seg_vol(5,1000)
	*    ,z_cntrl(100),z_dead(100),z_bottom(5,1000)
     *	,a_area(5,1000),b_area(5,1000)
     *    ,a_width(5,1000),b_width(5,1000)
     *  ,RMILE1(5,1000),RMILE2(5,1000)
     *  ,x_dist1(5,1000),x_dist2(5,-2:1000),u(5,1000)
C
C     BLOCKJ.COM - CONTAINS ARRAYS FOR REACH NAME INFORMATION
C                  FOR RNGKMOD.
C
      COMMON/BLOCK7/ HDNAME(20),PNAME(1000),RNAME(1000)
     *              ,rch_name(5),WPNAME(10),NPNAME(10)
     *              ,cp_name(50)
C
C     BLOCKK.COM - Contains arrays with meteorological data
C
      COMMON/BLOCK8/ EVRATE(10),QNS,QNA,DBT,WIND,PF,EA,PHPER
     *	,a_tmp(10),b_tmp(10),eql_tmp(10)
     *	,nwtype(5,1000),nwprov
c
c     State estimation components
c 
      COMMON/BLOCK9/H(2,5,5,365),P_var(5,-2:1000,2)
     *             ,Q_var(5),R_var(5),z_obs(2,5,5,365)
c
c     Integers characterizing the number of measurement locations
c     in each Reach and the element associated with the location
c
      COMMON/BLOCKA/meas_pnt(5,10),no_meas(5)
c
c     Declare variables
c
      COMMON/BLOCKB/ncase,case_name
     .             ,day_xplot(5,20),nxp_year(5),no_xplots(5)
c
      CHARACTER*3 cp_name
      CHARACTER*20 HDNAME,rch_name,NPNAME
      CHARACTER*30 CASE_NAME,RNAME
      CHARACTER*50 WPNAME
      CHARACTER*60 PNAME
      CHARACTER*80 XTITLE
      integer day_xplot,trib_id,trib_ndx,type_riv,type_res
      LOGICAL PLOT,trib
c
c     Version is the same as version 1.3
c































