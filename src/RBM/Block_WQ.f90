module block_wq
!
! Dimensioned and allocated water quality variables 
!
      real, dimension(:,:), allocatable:: DO
      real, dimension(:,:), allocatable:: BOD
      real, dimension(:,:), allocatable:: PO4
      real, dimension(:,:), allocatable:: P_Org
      real, dimension(:,:), allocatable:: NO2
      real, dimension(:,:), allocatable:: NO3
      real, dimension(:,:), allocatable:: NH4
      real, dimension(:,:), allocatable:: pH
      real, dimension(:,:), allocatable:: H2CO3
      real, dimension(:,:), allocatable:: HCO3
      real, dimension(:,:), allocatable:: CO3
      real, dimension(:,:), allocatable:: ALK
      real, dimension(:,:), allocatable:: ALGAE_1
      real, dimension(:,:), allocatable:: ALGAE_2
      real, dimension(:,:), allocatable:: ZOO_1
      real, dimension(:,:), allocatable:: ZOO_2
! 
      real, dimension(:,:), allocatable:: DO_trib
      real, dimension(:,:), allocatable:: BOD_trib
      real, dimension(:,:), allocatable:: PO4_trib
      real, dimension(:,:), allocatable:: P_Org_trib
      real, dimension(:,:), allocatable:: NO2_trib
      real, dimension(:,:), allocatable:: NO3_trib
      real, dimension(:,:), allocatable:: NH4_trib
      real, dimension(:,:), allocatable:: H2CO3_trib
      real, dimension(:,:), allocatable:: HCO3_trib
      real, dimension(:,:), allocatable:: CO3_trib
      real, dimension(:,:), allocatable:: ALK_trib
      real, dimension(:,:), allocatable:: ALGAE_1_trib
      real, dimension(:,:), allocatable:: ALGAE_2_trib
      real, dimension(:,:), allocatable:: ZOO_1_trib
      real, dimension(:,:), allocatable:: ZOO_2_trib
!
      real, dimension(:,:), allocatable:: DO_nps
      real, dimension(:,:), allocatable:: BOD_nps
      real, dimension(:,:), allocatable:: PO4_nps
      real, dimension(:,:), allocatable:: P_Org_nps
      real, dimension(:,:), allocatable:: NO2_nps
      real, dimension(:,:), allocatable:: NO3_nps
      real, dimension(:,:), allocatable:: NH4_nps
      real, dimension(:,:), allocatable:: H2CO3_nps
      real, dimension(:,:), allocatable:: HCO3_nps
      real, dimension(:,:), allocatable:: CO3_nps
      real, dimension(:,:), allocatable:: ALK_nps
      real, dimension(:,:), allocatable:: ALGAE_1_nps
      real, dimension(:,:), allocatable:: ALGAE_2_nps
      real, dimension(:,:), allocatable:: ZOO_1_nps
      real, dimension(:,:), allocatable:: ZOO_2_nps
!
end module block_wq
