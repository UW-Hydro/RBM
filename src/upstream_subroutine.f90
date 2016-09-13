subroutine upstream_subroutine(nseg,nr,ns,T_0, npndx, npart, n1, ncell, resx2)

use Block_Reservoir
use Block_Hydro
use Block_Network
use Block_Flow

implicit none

integer :: resx2, resx3,npndx, npart, ntrp, nseg, nr, ns, n1, ncell, i
real :: x, T_0, tntrp
integer, dimension(3):: ndltp=(/-2,-3,-3/)
integer, dimension(3):: nterp=(/3,4,3/)

          !
          !   Loop to establish where parcel started
          !   (e.g. in reservoir, downstream of reservoir, headwater, etc.)
          !

          ! --- if parcel started upstream of reservior and finished downstream -----
          if (reservoir.and. any(res_pres(nr,segment_cell(nr,nseg):segment_cell(nr,ns)) ) &
            .and. .not. res_pres(nr,segment_cell(nr,ns)) .and. .not. res_pres(nr,segment_cell(nr,nseg))) then

              do i = nseg, ns, 1
                if(res_pres(nr,segment_cell(nr,i))) then
                 T_0 = temp_out_i(res_num(nr,segment_cell(nr,i)))  !  
                 res_upstreamx = .true.
                 resx2 = res_num(nr,segment_cell(nr,i))
                 ncell0res = res_end_node(resx2)
  !  if(ncell .eq. 174) write(*,*) 'start up of res, finish down of res', 'nr',nr,'nseg',nseg,  'segment_cell' &
  !      ,segment_cell(nr,nseg),'resx2',resx2

                end if
               end do
! print *, 'nseg', nseg, 'ns',ns,'T_0',T_0, 'parcel started upstream of res, finished downstream of res'

          ! -------------------- if parcel started at a headwater ------------
          else if( nseg.eq.1 .and. &
             ( (.not. res_pres(nr,segment_cell(nr,nseg))) .or. ( res_pres(nr,segment_cell(nr,ns))) ) ) then
            T_0 = T_head(nr)
            res_upstreamx = .false.
            resx2 = 0
         !   if(nr .eq. 299) write(*,*) 'parcel started at headw'
            !  print *,'ncell',ncell, 'T_0', T_0
!            if(ns .eq. 3) print *,'headw',  'T_head',   T_head(nr)
! print *, 'nseg', nseg, 'ns',ns,'T_0',T_0, 'parcel started at headw'


          ! ------- if parcel started in reservoir and finished downstream  -----------
          else if (reservoir.and.res_pres(nr,segment_cell(nr,nseg))) then
            T_0 = temp_out_i(res_num(nr,segment_cell(nr,nseg)))  !  
            res_upstreamx = .true.
            resx2 = res_num(nr,segment_cell(nr,nseg))
            ncell0res = res_end_node(resx2)
 ! print *, 'nseg', nseg, 'ns',ns,'T_0',T_0, 'parcel started in res, finished downstream'
! print *, 'res_num', res_num(nr,segment_cell(nr,nseg)),'T_0', T_0
! print *, 'temp_out', temp_out_i(:)
  ! if(ns .gt. 7) print *,'parcel started in reservoir',  '
  ! nseg',nseg,'T_0',T_0
   !if(ns .eq. 3) print *,'ns',ns,'nseg',nseg, 'T_0', T_0

          ! -----------  if parcel is in reservoir (didn't finish downstream) -----------
          !              BUT not first cell in reservoir - since upstream flow is read in
          else if (reservoir.and.any(res_pres(nr,segment_cell(nr,nseg):segment_cell(nr,ns))) &
                .and.  res_pres(nr,segment_cell(nr,ns-1)) ) then

            !-- these two lines gets reservoir number in reach ---
            resx(segment_cell(nr,nseg): segment_cell(nr,ns)) = res_num(nr,segment_cell(nr,nseg):segment_cell(nr,ns))
            resx2 = maxval(resx(segment_cell(nr,nseg):segment_cell(nr,ns)),1)
            T_0 = temp_out_i(resx2)  ! temperature of reservoir parcel crossed
            res_upstreamx = .true.

! print *, 'nseg', nseg, 'ns',ns,'T_0',T_0, 'parcel in reservoir and not first seg'
  ! if(ns .gt. 7) print *,'if parcel is in reservoir', '  nseg',nseg,'T_0',T_0

           ! ----------- if parcel started in river and ended in river  -----------
          !            (i.e. did not start in headw, did not start in reservoir)
          !            OR if start node in the reservoir
          else
            resx2 = 0
            res_upstreamx = .false.
  !          if(nr .eq. 299) write(*,*) 'parcel started and ended in river'

! print *,'res_pres', res_pres(nr,1:5),'any-res pres', any(res_pres(nr,segment_cell(nr,ns):segment_cell(nr,nseg)) )
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

            tntrp_x = tntrp(xa,ta,x,nterp(npndx))  ! non-essential - just to print output

!  print *, 'nseg', nseg, 'ns',ns,'T_0',T_0, 'parcel in river, no res, or last res cell'
   !     print *,'---------------', ' ns', ns
   !     print *, 'nterp',nterp(npndx),'ta', ta
   !     print *, 'xa',xa, 'x',x

   !    if(ns .gt. 7) print *,'if parcel started and ended in river', '  nseg',nseg,'T_0',T_0
          ! ------------------------- end of the large if loop ----------------
          end if


end subroutine upstream_subroutine
