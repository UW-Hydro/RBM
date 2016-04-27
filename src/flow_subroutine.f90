SUBROUTINE flow_subroutine ( flow_in_epi_x , flow_in_hyp_x, flow_epi_hyp_x, flow_out_epi_x, flow_out_hyp_x &
                , ratio_sp, ratio_pen, nresx)! , density_epil, density_hypo, density_in)

   use Block_Reservoir
   use Block_Hydro 
   use Block_Network

implicit none

  real :: flow_in_epi_x, flow_in_hyp_x, flow_epi_hyp_x, flow_out_hyp_x,flow_out_epi_x 
  real :: ratio_sp, ratio_pen
  real :: Q1, Q2
!  real, dimension (:), allocatable :: T_epil,T_hypo, volume_e_x, volume_h_x
!  real, dimension (:), allocatable :: density_epil, density_hypo, density_in
  integer :: nresx, nnd
!  allocate (density_epil(nres))
!  allocate (density_hypo(nres))
!  allocate (density_in(nres))
 
 
      !*************************************************************************
      !      read inflow from vic/rvic simulations
      !*************************************************************************

      ! ------------------ read in in VIC flow data --------------
!print *,'nresx',nresx, 'Q_in', Q_res_in(nresx),'fttom', ftsec_to_msec, 'dtcomp',dt_comp
        Q1 = Q_res_in(nresx) * ftsec_to_msec * dt_comp
         ! converts ft^3/sec to m^3/sec, multiplys by seconds per time step

        if ( density_in(nresx) .le. density_hypo(nresx) ) then
                flow_in_hyp_x = Q1 * 0.2 
                flow_in_epi_x = Q1 * 0.8 ! majority of flow goes to epil.
        else
                flow_in_hyp_x = Q1 * 0.8
                flow_in_epi_x = Q1 * 0.2 ! majority flow goes to hypo
        end if
      !-------------- measured releases (penstock and spillway) ---------------
      !  read(56, *) datetime,Q_tot, Q_pen, Q_spill

      !  flow_out_hyp_x = Q_pen * ftsec_to_msec * dt_comp ! converts ft^3/sec to m^3/sec, multiplys by seconds per time step
      !  flow_out_epi_x = Q_spill * ftsec_to_msec * dt_comp ! converts ft^3/sec to m^3/sec, multiplys by seconds per time step
      !  flow_epi_hyp_x = flow_in_epi_x - flow_out_epi_x

      ! ---------- set outflow to inflow (same out as in)  ---
        Q2 = Q1
      !  Q_tot = 0  ! Q_tot was used to read in measured reservoir release data
        ! flow_out_hyp_x = Q_out * ftsec_to_msec * dt_comp

       ! ratio_sp = Q_spill/Q_tot
       ! ratio_pen = Q_pen/Q_tot

     !  if(Q_tot.gt.0) then
     !          flow_out_hyp_x = Q2 * ftsec_to_msec * dt_comp * (Q_pen/Q_tot)
     !          flow_out_epi_x = Q2 * ftsec_to_msec * dt_comp * (Q_spill/Q_tot)
     !          flow_epi_hyp_x = flow_in_epi_x - flow_out_epi_x
     !  else
           flow_out_hyp_x = Q2 ! * ftsec_to_msec * dt_comp
           flow_out_epi_x = 0
           flow_epi_hyp_x = flow_in_epi_x
     !  end if

 ! print *, 'flow_in_epi_x', flow_in_epi_x, 'flow_in_hyp_x', flow_in_hyp_x,'flow_out_hyp_x', flow_out_hyp_x
      ! ------------------------- calculate dV/dt terms ---------------------------
         dV_dt_epi(nresx) = (flow_in_epi_x - flow_out_epi_x - flow_epi_hyp_x) * T_epil(nresx)
         dV_dt_hyp(nresx) = (flow_in_hyp_x + flow_epi_hyp_x - flow_out_hyp_x) * T_hypo(nresx)

      !----- update epilimnion and hypolimnion volume  -------
        volume_e_x(nresx) = volume_e_x(nresx) + (flow_in_epi_x - flow_out_epi_x - flow_epi_hyp_x )
        volume_h_x(nresx) = volume_h_x(nresx) + (flow_in_hyp_x - flow_out_hyp_x + flow_epi_hyp_x)
        outflow_x = flow_out_epi_x + flow_out_hyp_x

        depth_e(nresx) = volume_e_x(nresx) / surface_area(nresx)
        depth_h(nresx) = volume_h_x(nresx) / surface_area(nresx)
  !      print *, 'depth_e', depth_e(nresx), 'depth_h', depth_h(nresx)
  !      print *,'nresx',nresx,'flow_cfs',Q1/(ftsec_to_msec * dt_comp),'vol',volume_e_x(nresx)+volume_h_x(nresx) &
  !              , 'residence time',( volume_e_x(nresx)+volume_h_x(nresx))/Q1
END SUBROUTINE flow_subroutine
