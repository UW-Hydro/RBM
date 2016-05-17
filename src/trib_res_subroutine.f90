subroutine trib_res_subroutine( ncellx,nr_trib, nr, ns)
use Block_Reservoir
use Block_Hydro
use Block_Network
use Block_Flow

implicit none

real :: T_dist
integer :: ntribs, nr_trib, nr, ns, ncellx

ncellx = segment_cell(nr,ns)

            !
            !   Loop to add all tributary flow and temperature entering segment
            !
            !
            ntribs=no_tribs(ncellx)

            ! ---------- loop through each tributary to add flow and temp -----------
            if(ntribs.gt.0) then

              ! -------------- cycle through each tributary ----------
              do ntrb=1,ntribs
                nr_trib=trib(ncellx,ntrb) ! gives nr (reach) for each trib
                if(Q_trib(nr_trib).gt.0.0) then
                  ! --- calculate tributary temperature and flow (to add to reservoir) ---    
                  T_trib_tot(ncellx) = (T_trib_tot(ncellx)*Q_trib_tot(ncellx) &
                     +  Q_trib(nr_trib)*T_trib(nr_trib)) &
                     /( Q_trib(nr_trib) + Q_trib_tot(ncellx))
                  Q_trib_tot(ncellx) = Q_trib(nr_trib) + Q_trib_tot(ncellx)

                end if
                !
              end do

            end if

            ! ----------- loop to add in later flow (if no tribs) -------------
            if(ntribs.eq.0.and.Q_diff(ncellx).gt.0) then !Q_diff is lateral flow in

              T_dist=T_head(nr)  ! temperature of later flow
              ! adding flow to tributatry flow for reservoir modeling
              T_trib_tot(ncellx) = (T_trib_tot(ncellx)*Q_trib_tot(ncellx) &
                + Q_diff(ncellx)*T_dist) / (Q_trib(nr_trib) + Q_diff(ncellx))
              Q_trib_tot(ncellx) = Q_trib_tot(ncellx) + Q_diff(ncellx)
            end if

            trib_res(ncellx) = .true.

end subroutine trib_res_subroutine
