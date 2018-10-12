SUBROUTINE Particle_Track(nr,x_head,x_bndry)
    USE Block_Hydro
    USE Block_Network
    use Block_Reservoir
    IMPLICIT NONE
    integer, intent(IN):: nr
    real, intent(IN):: x_head, x_bndry 
    integer:: ncell, ns, nx_part, nx_s
    real:: dt_total
    !
    !    This loop goes through each reach seperately
    !
    !     First do the reverse particle tracking
    !
    DO ns=no_celm(nr),1,-1  ! loops from outlet upstream to headwater in each reach
        !
        !     Segment is in cell SEGMENT_CELL(NC)
        !
        ncell=segment_cell(nr,ns) ! cell (node #) for specific reach and segment 
        nx_s=1
        nx_part=ns                   ! index number for segment
        dt_part(ns)=dt(ncell)        ! time to travel between segments
        dt_total=dt_part(ns)         ! total elapsed travel time between segment
        x_part(ns)=x_dist(nr,ns)     ! river mile distance of that segment

100 CONTINUE
    !
    !     Determine if the total elapsed travel time is equal to the
    !     computational interval
    !
    IF(dt_total.lt.dt_comp) THEN ! if time to travel between segments is less than computation interval (usually 86400 seconds),
                                  !     i.e. water will travel between segments faster than time interval
        x_part(ns)=x_part(ns) + dx(segment_cell(nr,nx_part)) ! add distance between segments to river mile of segment 
        !
        !     If the particle has started upstream from the boundary point, give it
        !     the value of the boundary
        !
        IF(x_part(ns).GE.x_bndry) THEN ! if next upstream cell is headwater
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
        nx_part=nx_part-1 !nx_part is segment #, so now it is set one upstream
        dt_part(ns)=dt(segment_cell(nr,nx_part)) !time to travel for next upstream segment
        dt_total=dt_total+dt_part(ns)  ! total time to travel PLUS time to travel for this segment
        GO TO 100   ! loops back to top, until the if statement isn't satisfied!
    ElSE    ! if time to travel between segments is greater than computational interval
        !
        !     For the last segment of particle travel, adjust the particle location
        !     such that the total particle travel time is equal to the computational
        !     interval.
        !
        dt_part(ns)=dt_comp-dt_total+dt_part(ns) 
        x_part(ns)=x_part(ns)+u(segment_cell(nr,nx_part))*dt_part(ns) !update river mile to where water was in previous time step
                    ! nx_part is segment where water was at previous time step. 
        IF(x_part(ns).GE.x_head) THEN
            x_part(ns)=x_head
            nx_s=nx_s-1
            dt_part(ns)=dt(head_cell(nr))
        END IF
    END IF
    
200 CONTINUE

    IF(nx_part.LT.1) nx_part=1
    nstrt_elm(ns)=nx_part   !segment of water at previous time step
    no_dt(ns)=nx_s    ! number of segments water traveled through

END DO

END SUBROUTINE Particle_Track
