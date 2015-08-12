SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_out)
!
Implicit NONE
!
integer:: nd,nr,ncell,ns 
real   :: T_0,T_dist
real(8):: time
real   :: T_head
real   :: dbt
real   :: Q_out
!
write (20,'(f12.4,4i6,3f8.2,f10.1)')           &                                              
            time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_out
end SUBROUTINE WRITE
