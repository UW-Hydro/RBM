SUBROUTINE Particle_Track(nr,ns,nx_s,nx_head)
USE Block_Hydro
USE Block_Network
IMPLICIT NONE
integer, intent(IN):: nr,ns
integer:: ncell, nx_head, nx_part, nx_s,count_steps
real:: dt_total,dt_xcess
!
!     First do the reverse particle tracking
!
!     Segment is in cell SEGMENT_CELL(NC)
   nx_s=0
   nx_part=ns
   dt_total=0.0
!
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!
    do while(dt_total.lt.dt_comp.and.nx_part.gt.0)
       nx_s=nx_s+1
       x_part(nx_s)=x_dist(nr,nx_part-1)
       ncell=segment_cell(nr,nx_part)       
       dt_part(nx_s)=dt(ncell)
       dt_total=dt_total+dt_part(nx_s)
!
!     Increment the segment counter if the total time is less than the
!     computational interval
!
        nx_part=nx_part-1
    end do
!
! If the elapsed time is greater than the computational interval and
! the parcel has not exceeded the boundary
!

    do while (dt_total.lt.dt_comp.and.nx_part.gt.0)
      nx_s=nx_s+1
      x_part(nx_s)=x_dist(nr,nx_part-1)
      ncell=segment_cell(nr,nx_part)
      dt_part(nx_s)=dt(ncell)
      nx_part=nx_part-1
      dt_total=dt_total+dt(ncell)
    end do
    if (dt_total.gt.dt_comp) then
      dt_xcess=dt_total-dt_comp
      dt_part(nx_s)=dt_part(nx_s)-dt_xcess
      x_part(nx_s)=x_dist(nr,nx_part+1)+u(ncell)*dt_part(nx_s)
    end  if
    nx_head=nx_part
    nx_part=max(1,nx_part)
    nstrt_elm(ns)=nx_part
    no_dt(ns)=nx_s 
END SUBROUTINE Particle_Track
!END MODULE P_Track
