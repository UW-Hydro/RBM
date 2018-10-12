module Block_Energy
    !
    !   Energy budget variables
    !
    !   Incoming short wave radiation, kcal/m**2/sec
    !
    real, dimension(:), allocatable::q_ns
    !
    !   Incoming atmospheric radiation, kcal/m**2/sec
    !
    real, dimension(:), allocatable::q_na
    !
    !   Air temperature at surface, deg. C
    !
    real, dimension(:), allocatable::dbt
    !  
    !   Wind speed, m/sec
    !
    real, dimension(:), allocatable::wind
    !
    !   Vapor pressure of air at surface, mb
    !
    real, dimension(:), allocatable::ea
    !
    !   Air pressure at surface, mb
    !
    real, dimension(:), allocatable::press 
    !
    real, dimension (:), allocatable::mu,alphamu,beta,gmma,smooth_param
    !   Some important constants
    !
    real   :: lvp,rb,rho,evap_coeff=1.5e-9,pf=0.640,pi=3.14159       
    !
    real :: q_evap_tot ! to save evaporation from reservoirs
    real, dimension(:), allocatable:: q_evap_output ! to save evaporation

end module Block_Energy  
