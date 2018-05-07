SUBROUTINE Energy(T_surf,q_surf,ncell,z,nd)
   use Block_Energy
   use Block_Network
   implicit none
   integer::i,ncell,nd
   integer::iter,num_timestep,tot_timestep
   real::A,B,e0,q_surf,q_conv,q_evap,q_ws,T_surf
   real::const1,const2,z
   real::deriv_1st !,deriv_conv,deriv_evap,deriv_ws
   real, dimension(2):: q_fit, T_fit
   real::q_tot,q_deriv,temp_calc
   real::rate_temp,error_estimate,timestep
!
   T_fit(1)=T_surf-1.0
   T_fit(2)=T_surf+1.0
   do i=1,2
      e0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
      rb=pf*(dbt(ncell)-T_fit(i))
      lvp=597.0-0.57*T_fit(i)
      q_evap=1000.*lvp*evap_coeff*wind(ncell)
      if(q_evap.lt.0.0) q_evap=0.0
      q_conv=rb*q_evap
      q_evap=q_evap*(e0-ea(ncell))
      q_ws=6.693E-2+1.471E-3*T_fit(i)
      q_fit(i)=q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv
!      if (ncell.eq.1827.and.nd.eq.197) write(*,*) &
!        'energy',q_ns(ncell),q_na(ncell),-q_ws,-q_evap,q_conv, &
!        'air',dbt(ncell),wind(ncell)
   end do
!
!     q=AT+B
!
!     Linear fit over the range of 2.0 deg C.
!     These results can be used to estimate the "equilibrium" 
!     temperature and linear rate constant.
!
   A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
   q_surf=0.5*(q_fit(1)+q_fit(2))
   B=(q_surf/A)-(T_fit(1)+T_fit(2))/2.
   deriv_1st=(q_surf/(z*rfac))
!   if (ncell.eq.1827) write(*,*) nd,A,B+T_surf,B,dbt(ncell)
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
!     Calculate the second derivative of Temperature
!
    const1=1000.*evap_coeff*wind(ncell)*pf
    const2=1000.*evap_coeff*wind(ncell)
    e0=2.1718E8*EXP(-4157.0/(T_surf+239.09))
    deriv_conv=const1*(-(597.0-0.57*T_surf) -0.57*(dbt(ncell)-T_surf))
    deriv_evap=const2*(-0.57*(e0-ea(ncell))+(597.0-0.57*T_surf)*e0*(4157.0/((T_surf+239.09)**2)))
    deriv_ws=1.471E-3

    deriv_2nd=deriv_1st*(deriv_conv-deriv_evap-deriv_ws)/(z*rfac) 
!
!     ******************************************************
!               Calculate equilibrium temperature
!     ******************************************************
!
if (ncell.eq.1827) then
    iter=1
    temp_equil=T_surf
    q_tot=q_surf
    q_deriv=deriv_conv-deriv_evap-deriv_ws
    temp_equil=0
    do while (iter.lt.50 .and. abs(q_tot/q_deriv).gt.0.001)
        e0=2.1718E8*EXP(-4157.0/(temp_equil+239.09))
        rb=pf*(dbt(ncell)-temp_equil)
        lvp=597.0-0.57*temp_equil
        q_evap=1000.*lvp*evap_coeff*wind(ncell)
        if(q_evap.lt.0.0) q_evap=0.0
        q_conv=rb*q_evap
        q_evap=q_evap*(e0-ea(ncell))
        q_ws=6.693E-2+1.471E-3*temp_equil
        q_tot=(q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv)
        
        rate_temp=q_tot/(z*rfac)

        const1=1000.*evap_coeff*wind(ncell)*pf
        const2=1000.*evap_coeff*wind(ncell)
        e0=2.1718E8*EXP(-4157.0/(temp_equil+239.09))
        deriv_conv=const1*(-(597.0-0.57*temp_equil) -0.57*(dbt(ncell)-temp_equil))
        deriv_evap=const2*(-0.57*(e0-ea(ncell))+(597.0-0.57*temp_equil)*e0*(4157.0/((temp_equil+239.09)**2)))
        deriv_ws=1.471E-3
        q_deriv=deriv_conv-deriv_evap-deriv_ws
        temp_equil=temp_equil-(q_tot/q_deriv)
        iter=iter+1
    enddo
    !
    !   Estimate the time to reach equilibrium tmeperature
    !
    error_estimate=deriv_2nd*dt_comp**2/2
    num_timestep=max(1,ceiling(abs(error_estimate/0.1)))
    timestep=dt_comp/num_timestep
    tot_timestep=0
    temp_calc=T_surf
    do while (tot_timestep.le.(50*num_timestep).and.abs(temp_equil-temp_calc).ge.0.01)
        e0=2.1718E8*EXP(-4157.0/(temp_calc+239.09))
        rb=pf*(dbt(ncell)-temp_calc)
        lvp=597.0-0.57*temp_calc
        q_evap=1000.*lvp*evap_coeff*wind(ncell)
        if(q_evap.lt.0.0) q_evap=0.0
        q_conv=rb*q_evap
        q_evap=q_evap*(e0-ea(ncell))
        q_ws=6.693E-2+1.471E-3*temp_calc
        q_tot=(q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv)

        rate_temp=q_tot/(z*rfac)
        temp_calc=temp_calc+rate_temp*timestep
        tot_timestep=tot_timestep+1
        if (nd.eq.197) then
            write(88,*) 'sub',tot_timestep,temp_equil,temp_calc,rate_temp,timestep
        end if
    enddo
    time_equil=tot_timestep*timestep
    if (nd.eq.197) then
        write(88,*) temp_equil,temp_calc,num_timestep,timestep,tot_timestep,error_estimate
    end if
end if
END Subroutine Energy
