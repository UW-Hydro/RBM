SUBROUTINE Particle_Track(nr,ns,nx_s,nx_head,ns_res_pres,ns_res_num)
    USE Block_Hydro
    USE Block_Network
    USE Block_Reservoir
    IMPLICIT NONE
    integer, intent(IN):: nr,ns
    integer:: ncell, nx_head, nx_part, nx_s,count_steps, nx_part_next
    integer:: ns_res_num
    real   :: dt_total,dt_xcess
    logical:: ns_res_pres
    !
    !     First do the reverse particle tracking
    !
    !     Segment is in cell SEGMENT_CELL(NC)
    nx_s=0
    nx_part=ns
    nx_part_next=ns
    dt_total=0.0
    ns_res_pres=.FALSE.
    ns_res_num = 0
    !
    !     Determine if the total elapsed travel time is equal to the
    !     computational intervasl
    !
    do while(dt_total.lt.dt_comp.and.nx_part.gt.0.and.nx_part_next.gt.0 &
            .and. .NOT. ns_res_pres)
        nx_s=nx_s+1
        x_part(nx_s)=x_dist(nr,nx_part_next)
        ncell=segment_cell(nr,nx_part_next)
        dt_part(nx_s)=dt(ncell)
        dt_total=dt_total+dt_part(nx_s)
        !
        !     Increment the segment counter if the total time is less than the
        !     computational interval
        !
        if (nx_s.gt.1) nx_part=nx_part-1
        nx_part_next=nx_part_next-1
        ns_res_pres = res_pres(segment_cell(nr,nx_part_next))
    end do
    ! If water particle travels back to a reservoir,
    ! ns_res_num denotes the number of the grid cell.
    !
    if (ns_res_pres) then
        ns_res_num=res_num(segment_cell(nr,nx_part_next))
    end if
    !
    ! If the elapsed time is greater than the computational interval and
    ! the parcel has not exceeded the boundary
    !
    if (dt_total.gt.dt_comp) then
        dt_xcess=dt_total-dt_comp
        dt_part(nx_s)=dt_part(nx_s)-dt_xcess
        x_part(nx_s)=x_dist(nr,nx_part+1)+u(ncell)*dt_part(nx_s)
    end if
    nx_head=nx_part
    nx_part=max(1,nx_part)
    nstrt_elm(ns)=nx_part
    no_dt(ns)=nx_s
END SUBROUTINE Particle_Track
