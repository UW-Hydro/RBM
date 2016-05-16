subroutine upstream_subroutine(nseg,nr,ns,T_0, npndx, npart, n1, ncell, resx2)

use Block_Reservoir
use Block_Hydro
use Block_Network
use Block_Flow

implicit none

integer :: resx2, npndx, npart, ntrp, nseg, nr, ns, n1, ncell
real :: x, T_0, tntrp
integer, dimension(3):: ndltp=(/-2,-3,-3/)
integer, dimension(3):: nterp=(/3,4,3/)

          !
          !   Loop to establish where parcel started
          !   (e.g. in reservoir, downstream of reservoir, headwater, etc.)
          !

          ! -------------------- if parcel started at a headwater ------------
          if(nseg.eq.1) then   
            T_0 = T_head(nr)
            res_upstreamx = .false.
            resx2 = 0
            !  print *,'ncell',ncell, 'T_0', T_0
!            if(ns .eq. 3) print *,'headw',  'T_head',   T_head(nr)

          ! ------- if parcel started in reservoir but finished downstream -----------
          else if (reservoir.and.res_pres(nr,segment_cell(nr,nseg))) then
            T_0 = temp_out(res_num(nr,segment_cell(nr,nseg)))  !  
            res_upstreamx = .true.
            resx2 = res_num(nr,segment_cell(nr,nseg))

      !if(ns .eq. 3) print *,'ns',ns,'nseg',nseg, 'T_0', T_0
          ! -----------  if parcel is in reservoir (didn't finish downstream) -----------
          !              BUT not first cell in reservoir - since upstream flow is read in
          else if (reservoir.and.any(res_pres(nr,segment_cell(nr,ns):segment_cell(nr,nseg))) .and. &
               any(segment_cell(nr,ns) .ne. res_start_node(:)) ) then

            !-- these two lines gets reservoir number in reach ---
            resx(segment_cell(nr,nseg): segment_cell(nr,ns)) = res_num(nr,segment_cell(nr,nseg):segment_cell(nr,ns))
            resx2 = maxval(resx(segment_cell(nr,nseg):segment_cell(nr,ns)),1)
            T_0 = temp_out(resx2)  ! temperature of reservoir parcel crossed
            res_upstreamx = .true.

          ! ----------- if parcel started in river and ended in river  -----------
          !            (i.e. did not start in headw, did not start in reservoir)
          !            OR if start node in the reservoir
          else
            res_upstreamx = .false.
            resx2 = 0

            !
            !  third order interpolation at the downstream boundary
            !  to get starting temperature values at each parcel
            !
            if(nseg.eq.no_celm(nr)) npndx=3  ! if segment at previous time step was last segment
            !
            do ntrp=1,nterp(npndx)
              npart=nseg+ntrp+ndltp(npndx)
              xa(ntrp)=x_dist(nr,npart) !river mile at npart segment
              ta(ntrp)=temp(nr,npart,n1)  !temperature at npart
            end do
            x=x_part(ns)
            !
            !     Call the interpolation function
            !
            T_0=tntrp(xa,ta,x,nterp(npndx))

          ! ------------------------- end of the large if loop ----------------
          end if


end subroutine upstream_subroutine
