subroutine stream_density(SOME_TEMPERATURE,THE_DENSITY)
[NO BLANK LINES - JRY]
 implicit none

 ! real :: SOME_TEMPERATURE, THE_DENSITY


!-------------------calculate the density based on temperature----------------
    THE_DENSITY = 1.000028*1e-3*((999.83952+16.945176*SOME_TEMPERATURE-& 
        & (7.9870401e-3*(SOME_TEMPERATURE**2)-46.170461e-6*(SOME_TEMPERATURE**3))+&
        & (105.56302e-9*(SOME_TEMPERATURE**4)-280.54235e-12*(SOME_TEMPERATURE**5)))/&
        & (1+16.87985e-3*SOME_TEMPERATURE)

end subroutine stream_density

