subroutine stream_density(nresx)
! subroutine stream_density(T_epil, T_hypo, stream_T_in, density_epil,density_hypo, density_in, nresx)
! subroutine stream_density(T_epil(nresx), T_hypo(nresx), stream_T_in(nresx), density_epil(nresx),density_hypo(nresx), density_in(nresx), nresx)
 use Block_Reservoir

 implicit none

 ! real, dimension (:), allocatable :: T_epil,T_hypo, volume_e_x, volume_h_x
 ! real, dimension (:), allocatable :: stream_T_in
 ! real, dimension (:), allocatable :: density_epil, density_hypo, density_in
  integer :: nresx

!  allocate (T_epil(nres))
! allocate (T_hypo(nres))
! allocate (stream_T_in(nres))
! allocate (density_epil(nres))
! allocate (density_hypo(nres))
! allocate (density_in(nres))

!-------------------calculate the density based on temperature----------------
    density_epil(nresx) = 1.000028*1e-3*((999.83952+16.945176*T_epil(nresx))-& 
        & (7.9870401e-3*(T_epil(nresx)**2)-46.170461e-6*(T_epil(nresx)**3))+&
        & (105.56302e-9*(T_epil(nresx)**4)-280.54235e-12*(T_epil(nresx)**5)))/&
        & (1+16.87985e-3*T_epil(nresx))
    density_hypo(nresx) = 1.000028*1e-3*((999.83952+16.945176*T_hypo(nresx))-&
        & (7.9870401e-3*(T_hypo(nresx)**2)-46.170461e-6*(T_hypo(nresx)**3))+ &
        & (105.56302e-9*(T_hypo(nresx)**4)-280.54235e-12*(T_hypo(nresx)**5)))/&
        & (1+16.87985e-3*T_hypo(nresx))
    density_in(nresx) = 1.000028*1e-3*((999.83952+16.945176*T_res_in(nresx))&
        & -(7.9870401e-3*(T_res_in(nresx)**2)-46.170461e-6*&
        & (T_res_in(nresx)**3))+(105.56302e-9*(T_res_in(nresx)**4)- &
        & 280.54235e-12*(T_res_in(nresx)**5)))/(1+16.87985e-3* &
        & T_res_in(nresx))

! print *, 'T_epil', T_epil(nresx), 'T_hypo', T_hypo(nresx), 'T_res_in', T_res_in(nresx)

end subroutine stream_density

