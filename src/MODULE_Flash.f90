!************************************************************
!  File Type:     MODULE
!  File Content:  Modules containing variables for Flash calculations
!  
!  Created:       2017.12.22 /Bkaa ScaleSim
!  MOdified:
!
!  Content:       Module Flash - parameters used in flash calculations
!************************************************************

   
   
   
   MODULE Flash
   Use Par
   Implicit None
   
   Double Precision PFL, TFL  ! Pressure (bara) and Temperature(K) used in flash
   Double Precision K(NoHC)
   Double Precision Beta
   Double Precision dFugCoef(0:2,1:NoHC,1:NoHC)
   
   Double Precision Nc(0:2,0:NoHC)
   Double Precision xy(1:2,1:NoHC)
   
   
   Integer FlashTeller
   
   
   END MODULE Flash
   
  