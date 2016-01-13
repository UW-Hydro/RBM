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
COMMON/BLOCK2/dt_comp &
 ,QSUM,tw_init &
 ,XTITLE,ysim,delta_n 
!
!  The following variables are now dimensioned (UW_JRY_2011/06/20)
!
! ,mu,alphaMu,gamma,beta
!
!
COMMON/BLOCK3/ &
  no_celm(500),no_cells(500),no_tribs(500),node(500) &
!
!    Added a variable ndelta here (UW_JRY_2011/03/15)
!
 ,main_stem(500),last_seg(500),ndelta(100)
!
!    Added last_seg(5000) to COMMON/BLOCK3/ to handle tributary
!    input temperature  11/14/2008
!
!
!    Added ndelta(100) to accommodate variable segment lengths 
!    within a given reach.  Used primarily for reservoir reaches
!    where particles may not reach a boundary in one daily time step
!    UW_JRY_2011/03/10)
!
COMMON/BLOCK4/qin(500),q_trib(500),qout(500) &
 ,qdiff(500),depth(500),width(500)
!
!
COMMON/BLOCK5/ temp(50,-2:1000,2),T_trib(500),thermal(50,1000)
!

!     BLOCK6.COM - CONTAINS ARRAYS WITH DEPTH AND VELOCITY
!                  CHARACTERISTICS AND RIVER MILE INDICES
!
COMMON/BLOCK6/dx(500),dt(500) &
 ,x_dist(50,0:1000),u(500)
!
!
!     BLOCK8.COM - Contains arrays with meteorological data
!
COMMON/BLOCK8/ &
  QNS(500),QNA(500),DBT(500),WIND(500) &
 ,EA(500),PRESS(500),PF,PHPER,qqns,qqna,qqws,qqevap,qqconv &
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
