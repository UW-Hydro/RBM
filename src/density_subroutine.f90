subroutine stream_density(temp_x,density_x)
 implicit none
 real :: density_x, temp_x
!
!-------------------calculate the density based on temperature----------------
    density_x = 1.000028*1e-3*((999.83952+16.945176*temp_x-& 
        & (7.9870401e-3*(temp_x**2)-46.170461e-6*(temp_x**3))+&
        & (105.56302e-9*(temp_x**4)-280.54235e-12*(temp_x**5)))/&
        & (1+16.87985e-3*temp_x)
!
end subroutine stream_density

