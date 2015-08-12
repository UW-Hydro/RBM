SUBROUTINE Energy(T_surf,q_surf,ncell)
   use Block_Energy
   implicit none
   integer::i,ncell,nd
   real::A,B,e0,q_surf,q_conv,q_evap,q_ws,td,T_surf
   real, dimension(2):: q_fit, T_fit
!
   td=nd
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
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
END Subroutine Energy
