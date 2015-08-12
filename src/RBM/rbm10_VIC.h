!
!     UNDIMENSIONED INTEGER VARIABLES
!
COMMON/BLOCK1/ &
  flow_cells,heat_cells &
 ,NDAYS,nreach,ntrb,nysim &
 ,nyear1,nyear2 &
 ,n1,n2,no_rch,nwpd,nd
!
!     UNDIMENSIONED FLOATING POINT VARIABLES
!
!
!  The following variables are now dimensioned (UW_JRY_2011/06/20)
!
!      COMMON/BLOCK2/smooth_param,rminsmooth,dt_comp &
      COMMON/BLOCK2/dt_comp &
! ,NCONST,PLOT
 ,QSUM,time,tw_init &
 ,XTITLE,ysim,delta_n 
!
!  The following variables are now dimensioned (UW_JRY_2011/06/20)
!
! ,mu,alphaMu,gamma,beta
!
!
COMMON/BLOCK3/ &
  no_celm(5000),no_cells(5000),no_tribs(5000),node(5000) &
!
!    Added a variable ndelta here (UW_JRY_2011/03/15)
!
 ,main_stem(5000),last_seg(5000),ndelta(1000)
!
!    Added last_seg(5000) to COMMON/BLOCK3/ to handle tributary
!    input temperature  11/14/2008
!
!
!    Added ndelta(1000) to accommodate variable segment lengths 
!    within a given reach.  Used primarily for reservoir reaches
!    where particles may not reach a boundary in one daily time step
!    UW_JRY_2011/03/10)
!
COMMON/BLOCK4/qin(5000),q_trib(5000),qout(5000) &
 ,qdiff(5000),depth(5000),width(5000)
!
!
COMMON/BLOCK5/ temp(500,-2:1000,2),T_trib(5000),thermal(50,1000)
!

!     BLOCK6.COM - CONTAINS ARRAYS WITH DEPTH AND VELOCITY
!                  CHARACTERISTICS AND RIVER MILE INDICES
!
COMMON/BLOCK6/dx(5000),dt(5000) &
 ,x_dist(500,0:1000),u(5000)
!
!
!     BLOCK8.COM - Contains arrays with meteorological data
!
COMMON/BLOCK8/ &
  QNS(5000),QNA(5000),DBT(5000),WIND(5000) &
 ,EA(5000),PRESS(5000),PF,PHPER,qqns,qqna,qqws,qqevap,qqconv &
!
! add arrays for Mohseni parameters.  Assigned at headwaters UW_JRY 6/18/2011
!
 ,mu(1000),alphaMu(1000),beta(1000),gmma(1000),smooth_param(1000)
!
!     Declare variables
!
COMMON/BLOCK9/segment_cell(500,1000),trib(5000,50) &
 ,head_cell(5000)
integer flow_cells,heat_cells,segment_cell,trib,head_cell
!
!  Declare mu (UW_JRY_2011/06/20)
real*4 mu
