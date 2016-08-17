subroutine trib_subroutine(nncell,ncell0, T_0,nr_trib & 
      , nr, ns, nseg, n2, DONE, dt_calc, dt_total, ncell)
use Block_Reservoir
use Block_Hydro
use Block_Network
use Block_Flow

implicit none

real ::Q1, T_0, Q2, T_dist, dt_calc,dt_total, temp_source
integer :: ntribs, nncell,nr_trib, nr, ns, nseg, n2, ncell0, ncell
logical :: DONE

            !
            !   Loop to add all tributary flow and temperature entering segment
            !
            !

  !   print *, 'trib_subroutine loop', '   no_tribs(nncell)', no_tribs(nncell),'nncell', nncell
            Q1=Q_in(nncell)  !flow entering cell from reach
            ntribs=no_tribs(nncell)

            ! ---------- loop through each tributary to add flow and temp -----------
            if(ntribs.gt.0.and..not.DONE) then

              ! -------------- cycle through each tributary ----------
!  print *, 'segment in trib subroutine', ns, 'nr_trib', nr_trib,'ntribs',ntribs
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb) ! gives nr (reach) for each trib
             !   if(nncell .eq. heat_cells) nr_trib = 0
                if(Q_trib(nr_trib).gt.0.0) then

!  print *,'cell', nncell,'ns',ns, 'stream temperature', T_0, 'trib_temperature', T_trib

                  ! --- add trib flow and temperature to total flow and temp --
                  Q2=Q1+Q_trib(nr_trib)  !Q1 is reach inflow, Q2 is total flow
                  T_0=(Q1*T_0+Q_trib(nr_trib)*T_trib(nr_trib))/Q2 !adjust temp based on trib temp/flow


      ! if(ncell == 941) print *,'T_0',T_0,'Q1',Q1,'Q_trib',Q_trib(nr_trib),'T_trib' &
      !    ,T_trib(nr_trib),'Q_2',Q2

                end if
                !
                Q1=Q2    ! Q_out(nncell) !
                !
              end do

              DONE=.TRUE. ! signify when all tribs have been cycled through
            end if

            ! ----------- loop to add in later flow (if no tribs) -------------
            if(ntribs.eq.0.and.Q_diff(nncell).gt.0) then !Q_diff is lateral flow in

              ! ----- add lateral flow to total flow ------             
              Q2=Q1+Q_diff(nncell)
              T_dist=T_head(nr)
              T_0=(Q1*T_0+Q_diff(nncell)*T_dist)/Q2
              Q1=Q2
            end if


          ! ------------- add point source thermal flow if present -----------

            if(source) then
              if(source_cell_tf(nncell)) then
                temp_source = T_0 + 5  ! just adding 5 deg C to flow
                sourcex = source_num_cell(nncell)
                T_0=(Q1*T_0 + flow_source(sourcex)  * temp_source) / (Q1 + flow_source(sourcex) )
               ! T_0=(Q1*T_0 + 1000  * temp_source) / (Q1 + 1000 ) ! John Sevier has flow of 1000 cfs
               ! print *, 'nncell', nncell, 'ncell0', ncell0
              end if
            end if

   ! write(*,*) 'trib subroutine: post thermal inflow nd ', nd
          nseg=nseg+1  ! so nseg will be next segment down stream
          nncell=segment_cell(nr,nseg)  ! node for next node downstream
! if(nr.eq.299 .and. nseg .gt. 150)  write(*,*) 'nr',nr,'nseg',nseg,'segment_cell(nr,nseg)', segment_cell(nr,nseg)
          !
          !     Reset tributary flag if this is a new celln_write
          !
          if(ncell0.ne.nncell) then
            ncell0=nncell
            DONE=.FALSE.
          end if

  !  write(*,*) 'trib subroutine: pre dt_calc '
  !  write(*,*) 'dt_calc', dt_calc,  'nncell',nncell

          ! ------ add time to pass between segments to total time ------
          dt_calc=dt(nncell)

  !  write(*,*) 'trib subroutine: post dt_calc'

          dt_total=dt_total+dt_calc
!    print *, 'segment in end of trib subroutine', ns, 'nr_trib', nr_trib

end subroutine trib_subroutine
