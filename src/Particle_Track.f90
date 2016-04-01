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
DO ns=no_celm(nr),1,-1  ! loops from no_celm(nr) to 1, backwards by 1
    !
    !     Segment is in cell SEGMENT_CELL(NC)
    !
    ncell=segment_cell(nr,ns) !cell and/or node # for specific reach and segment 
    nx_s=1
    nx_part=ns                   ! index number for segment?
    dt_part(ns)=dt(ncell)        ! time it would take to go from headw to segment
                                 ! so if less than dt_comp, parcel comes from headw
    dt_total=dt_part(ns)         ! total elapsed travel time
    x_part(ns)=x_dist(nr,ns)     ! river mile distance of that segment

    ! 
    !      Calculate time to next upstream reservoir 
    !

!    if(res_pres = 'TRUE')THEN  !'res_pres' is just TRUE/FALSE statement in 5th line of Network file
          
!      dt_part_res = dt_res(ncell)  ! time for water in cell to get to next upstream dam

 !       if('if dam upstream')THEN

  !      ELSE('if dam not upstream')

   !     END IF

  !  ELSE   ! if reservoir presence is "FALSE"

 !   END IF


    100 CONTINUE

    !
    !     Determine if the total elapsed travel time is equal to the
    !     computational interval
    !
    IF(dt_total.lt.dt_comp) THEN ! if time for water to travel from cell to headw is less than computation interval
       x_part(ns)=x_part(ns) + dx(segment_cell(nr,nx_part)) ! where particle started from: location this time step + distance water traveled  
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

  !  ELSE IF(                ) Then    ! if the travel time puts water in reservoir






   !     GO TO 100

    ElSE    ! if the total travel time is the same as the time step (computational interval)
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
