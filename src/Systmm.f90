Module SYSTM
!
integer:: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer:: nc,nd,ndd,nm,nr,ns
integer:: nr_trib,ntrb,ntribs
integer:: nrec_flow,nrec_heat
integer:: n1,n2,nnd,nobs,ndays,nyear,nd_year,ntmp
integer:: npart,nseg,nwpd 
real::    dt_comp,dt_calc,dt_total,hpd,Q1,Q2,q_dot,q_surf,z
real   :: rminsmooth
real   :: T_0,T_dist
real(8):: time
real   :: x,x_bndry,xd,xdd,xd_year,x_head,xwpd,year
real,dimension(:),allocatable:: T_head,T_smth,T_trib
real,dimension(:,:,:),allocatable:: temp
!
!
logical:: DONE
!
! Indices for lagrangian interpolation
!
integer:: npndx,ntrp
integer, dimension(3):: ndltp=(/-2,-3,-3/)
integer, dimension(3):: nterp=(/3,4,3/)

!
real, parameter:: pi=3.14159,rfac=304.8
!
!
contains
!
SUBROUTINE SYSTMM(temp_file,param_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
!
Implicit None
!
character (len=200):: temp_file
character (len=200):: param_file
!
integer::njb
!
real             :: tntrp
real,dimension(4):: ta,xa
!
! Allocate the arrays
!
allocate (temp(nreach,-2:ns_max,2))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_diff(heat_cells))
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
! Initialize headwaters temperatures
!
T_head=4.0
!!
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
do nyear=start_year,end_year
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
      do nr=1,nreach

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
      temp(nr,-1,n1)=T_head(nr)
      temp(nr,-2,n1)=T_head(nr)
      temp(nr,no_celm(nr)+1,n1)=temp(nr,no_celm(nr),n1)
      x_head=x_dist(nr,0)
      x_bndry=x_head-50.0
!
!     Establish particle tracks
!
      call Particle_Track(nr,x_head,x_bndry)
!
!
        DONE=.FALSE.
        do ns=1,no_celm(nr)
          ncell=segment_cell(nr,ns)
!
!     Now do the third-order interpolation to
!     establish the starting temperature values
!     for each parcel
!
          nseg=nstrt_elm(ns)
!
!     Perform polynomial interpolationnr_trib
!
!
!     Interpolation inside the domain
!
          npndx=2
!
!     Interpolation at the upstream boundary
!
          if(nseg.eq.1) then
            T_0 = T_head(nr)
          else 
!
!     Interpolation at the downstream boundary
!
          if(nseg.eq.no_celm(nr)) npndx=3
!
          do ntrp=1,nterp(npndx)
            npart=nseg+ntrp+ndltp(npndx)
            xa(ntrp)=x_dist(nr,npart)
            ta(ntrp)=temp(nr,npart,n1)
          end do
          x=x_part(ns)
!
!     Call the interpolation function
!
          T_0=tntrp(xa,ta,x,nterp(npndx))
          end if
!
300 continue
350 continue
!          dt_calc=dt_part(ns)
          nncell=segment_cell(nr,nstrt_elm(ns))
!
!    Set NCELL0 for purposes of tributary input
!
          ncell0=nncell
          dt_total=dt_calc
          do nm=no_dt(ns),1,-1
            z=depth(nncell)
            call energy(T_0,q_surf,nncell)
            q_dot=(q_surf/(z*rfac))
            T_0=T_0+q_dot*dt_calc
            if(T_0.lt.0.0) T_0=0.0
!
!     Look for a tributary.
!
            Q1=Q_in(nncell)
            ntribs=no_tribs(nncell)
            if(ntribs.gt.0.and..not.DONE) then
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(Q_trib(nr_trib).gt.0.0) then
                  Q2=Q1+Q_trib(nr_trib)
                  T_0=(Q1*T_0+Q_trib(nr_trib)*T_trib(nr_trib))/Q2
                end if
!
                Q1=Q_out(nncell)
!
              end do
              DONE=.TRUE.
            end if
            if(ntribs.eq.0.and.Q_diff(nncell).gt.0) then
              Q2=Q1+Q_diff(nncell)
              T_dist=T_head(nr)
              T_0=(Q1*T_0+Q_diff(nncell)*T_dist)/Q2
              Q1=Q2
            end if
500 continue
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new celln_write
!
            if(ncell0.ne.nncell) then
              ncell0=nncell
              DONE=.FALSE.
            end if
            dt_calc=dt(nncell)
            dt_total=dt_total+dt_calc
          end do
          if (T_0.lt.0.5) T_0=0.5
            temp(nr,ns,n2)=T_0
	    T_trib(nr)=T_0
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
            call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell))
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
end SUBROUTINE SYSTMM
end module SYSTM
