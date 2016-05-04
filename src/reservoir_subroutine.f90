SUBROUTINE reservoir_subroutine(nresx, nd, q_surf,time)
! SUBROUTINE reservoir_subroutine(T_epil,T_hypo, volume_e_x,volume_h_x)
   use Block_Reservoir
   use Block_Flow
   use Block_Network

 implicit none

 real :: dayx, q_surf
 real(8):: time
 integer :: nd,  nresx, nyear
!  real, dimension (:), allocatable :: T_epil,T_hypo,volume_e_x,volume_h_x,stream_T_in


 ! ---------------- turnover loop driven only by T_epil and T_hyp ----------
        dayx = nd  ! day of year
        if ( ( T_epil(nresx) - T_hypo(nresx)) .lt. (2) .and. dayx .gt. 245) then !245 is September 1st
                if( (T_epil(nresx) - T_hypo(nresx)) .lt. (0) ) then
                         K_z(nresx) = 1 ! set high K_z when moderately unstable
                else
                         K_z(nresx) = 0.1 ! set moderate K_z when system is unstable
                end if
        else ! if T_epil greater than T_hypo
                  K_z(nresx) = 0.001  ! set the diffusion coeff. in m^2/day
                  K_z(nresx) = K_z(nresx) / (depth_e(nresx)/2) ! divide by approx thickness of thermocl.
        end if
  ! -------------------- calculate temperature terms  -------------------------
      dif_epi_x  = K_z(nresx) * surface_area(nresx) *  (T_hypo(nresx) - T_epil(nresx)) / volume_e_x(nresx)
      dif_hyp_x  = K_z(nresx) * surface_area(nresx) *  (T_epil(nresx) - T_hypo(nresx)) / volume_h_x(nresx)
!print *,'nres',nres, 'dif_epi_x', dif_epi_x,'K_z',K_z(nresx),  'surface_area(nresx)',surface_area(nresx)&
!        ,'T_hypo', T_hypo(nresx),'T_epil',T_epil(nresx), 'vol', volume_e_x(nresx)
  ! --------------------- calculate advection terms --------------------------- 
         advec_in_epix  = flow_in_epi_x * (T_res_in(nresx) - T_epil(nresx)) /volume_e_x(nresx)
         advec_epi_hyp = flow_epi_hyp_x *  (T_epil(nresx) - T_hypo(nresx)) / volume_e_x(nresx)
         advec_in_hypx = flow_in_hyp_x * (T_res_in(nresx) - T_hypo(nresx)) /volume_h_x(nresx)
! print *, 'advec_in_epix', advec_in_epix, 'flow_in_epi_x', flow_in_epi_x, 'T_res_in(nresx)', T_res_in(nresx) &
!        , 'T_hypo(nresx)',T_hypo(nresx)
  ! ------------------- calculate change in temperature  ---------------------
      ! ---------------- epilimnion -----------
         ! ------------ calculate total energy ----------
          energy_x  = (q_surf * dt_comp ) / (depth_e(nresx) * density * heat_c_kcal ) ! kcal/sec*m2 to C/day
          temp_change_ep(nresx) = advec_in_epix + dif_epi_x + energy_x

!  if(nresx==1) print *,nd,'tempchange',temp_change_ep(nresx),'advec', advec_in_epix,'dif',dif_epi_x, 'energy',energy_x ! ,q_surf,  dt_comp, depth_e(nresx), density,  heat_c_kcal

! print *,'energy',energy_x, 'q_surf',q_surf, 'dt_comp',dt_comp, 'depth_e',depth_e(nresx),'density',density, 'kcal',heat_c_kcal
!          temp_change_ep(nresx) = advec_in_epix + energy_x +  dif_epi_x !  - advec_out_epix - dV_dt_epi(nresx) ! units = C/day
! print *, 'temp_change', temp_change_ep(nresx), 'advec_in', advec_in_epix &
!        , 'energyx', energy_x, 'diffusion', dif_epi_x
 !    print *, 'T_epil', T_epil(nresx), 'T_hypo', T_hypo(nresx)

         !----- update epilimnion volume for next time step -------
          T_epil(nresx) = T_epil(nresx) + temp_change_ep(nresx)
      ! ------------------ hypolimnion ----------------
         ! ------------ calculate total energy ----------
          temp_change_hyp(nresx) = advec_in_hypx + advec_epi_hyp + dif_hyp_x !  - advec_out_hypx - dV_dt_hyp(nresx)
         !----- update epilimnion volume for next time step -------
          T_hypo(nresx) = T_hypo(nresx) +  temp_change_hyp(nresx)

!  if(nresx==1) print *, nd,energy_x ,q_surf,  dt_comp, depth_e(nresx), density,  heat_c_kcal
!  if(nresx==1) print *,temp_change_ep(nresx), T_epil(nresx)

! if(nresx==3) print *,nd, (volume_e_x(nresx) + volume_h_x(nresx))/ (flow_in_epi_x + flow_in_hyp_x)
! if(nresx==4) print *, nd,  T_epil(nresx), T_hypo(nresx)
! print *,'T_res_in',T_res_in(nresx), 'T_epil', T_epil(nresx), 'T_hypo', T_hypo(nresx)
 
  !---------- calculate combined (hypo. and epil.) temperature of outflow -----
    epix = T_epil(nresx)*(flow_out_epi_x/outflow_x)  ! portion of temperature from epilim. 
    hypox= T_hypo(nresx)*(flow_out_hyp_x/outflow_x)  ! portion of temperature from hypol.
    temp_out(nresx) = epix + hypox   ! average outflow temperature
    volume_tot = volume_e_x(nresx)  + volume_h_x(nresx)
    T_res(nresx) = (T_epil(nresx) * (volume_e_x(nresx)/volume_tot)) + &
         (T_hypo(nresx)*(volume_h_x(nresx)/volume_tot) ) ! weighted averge temp

 if(nresx.eq.4) then
  ! write(46,*),time, advec_in_epix, advec_epi_hyp, advec_in_hypx, dif_epi_x &
  !      , dif_hyp_x, energy_x, temp_change_ep(nresx), temp_change_hyp(nresx) &
  !      , T_epil(nresx), T_hypo(nresx) 

!   write(47,*),time, K_z(nresx), surface_area(nresx),  T_hypo(nresx), T_epil(nresx),  volume_e_x(nresx)
 
  write(48,*),time,advec_in_epix, flow_in_epi_x, T_res_in(nresx),  T_epil(nresx), volume_e_x(nresx) &
        , advec_in_hypx, flow_in_hyp_x, T_hypo(nresx), volume_h_x(nresx), flow_epi_hyp_x, temp_out(nresx), T_res(nresx)
 
!    write(49, *), time 
 end if
end subroutine reservoir_subroutine
