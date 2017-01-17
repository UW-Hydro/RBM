SUBROUTINE Particle_Track(nr,x_head,x_bndry)
USE Block_Hydro
USE Block_Network
IMPLICIT NONE
integer, intent(IN):: nr
real, intent(IN):: x_head, x_bndry 
integer:: ncell, ns, nx_part, nx_s
real:: dt_total

!
!     First do the reverse particle tracking
!
DO ns=no_celm(nr),1,-1
!
!     Segment is in cell SEGMENT_CELL(NC)
!
   ncell=segment_cell(nr,ns)
   nx_s=1
   nx_part=ns
   dt_part(ns)=dt(ncell)
   dt_total=dt_part(ns)
   x_part(ns)=x_dist(nr,ns)

100 CONTINUE

!
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!
   IF(dt_total.lt.dt_comp) THEN
      x_part(ns)=x_part(ns)+dx(segment_cell(nr,nx_part))
!
!     If the particle has started upstream from the boundary point, give it
!     the value of the boundary
!
      IF(x_part(ns).GE.x_bndry) THEN
         x_part(ns)=x_head
         dt_part(ns)=dt(segment_cell(nr,nx_part))
         dt_total=dt_total+dt_part(ns)
         GO TO 200
      END IF
!
!     Increment the segment counter if the total time is less than the
!     computational interval
!

      nx_s=nx_s+1
      nx_part=nx_part-1
      dt_part(ns)=dt(segment_cell(nr,nx_part))
      dt_total=dt_total+dt_part(ns)
      GO TO 100
   ELSE
!
!     For the last segment of particle travel, adjust the particle location
!     such that the total particle travel time is equal to the computational
!     interval.
!
      dt_part(ns)=dt_comp-dt_total+dt_part(ns)
      x_part(ns)=x_part(ns)+u(segment_cell(nr,nx_part))*dt_part(ns)
      IF(x_part(ns).GE.x_head) THEN
         x_part(ns)=x_head
         nx_s=nx_s-1
         dt_part(ns)=dt(head_cell(nr))
      END IF
   END IF
200 CONTINUE
   IF(nx_part.LT.1) nx_part=1
   nstrt_elm(ns)=nx_part
   no_dt(ns)=nx_s
END DO
END SUBROUTINE Particle_Track
!END MODULE P_Track
