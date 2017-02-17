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
      real             :: lvp,rb,rho
      real,parameter   :: evap_coeff=1.5e-9 !Lake Hefner coefficient, 1/meters
      real,parameter   :: pf=0.640,pi=3.14159
      real,parameter   :: rfac=304.8 !rho/Cp kg/meter**3/Kilocalories/kg/Deg K
!
end module Block_Energy
