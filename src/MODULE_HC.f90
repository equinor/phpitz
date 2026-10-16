!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to HC-cpecies
!
!************************************************************


MODULE HC
   Use Par

! Variables

   Integer, Parameter :: PhTot   = 0
   Integer, Parameter :: PhGas   = 1
   Integer, Parameter :: PhAq    = 2
   Integer, Parameter :: PhDisoc = 3
   Integer, Parameter :: Phases  = 3


   Double Precision HCN(0:Phases,0:NoHC)       ! Total molerate of HC-components


    ! ManualUpdate update liste if new components are added
   Integer, Parameter :: iHC_CO2=1
   Integer, Parameter :: iHC_H2O=2
   Integer, Parameter :: iHC_N2=3
   Integer, Parameter :: iHC_O2=4
   Integer, Parameter :: iHC_H2S=5
   Integer, Parameter :: iHC_NO=6
   Integer, Parameter :: iHC_NO2=7
   Integer, Parameter :: iHC_SO2=8
   Integer, Parameter :: iHC_H2SO4=9
   Integer, Parameter :: iHC_HNO3=10
   Integer, Parameter :: iHC_S8=11
   Integer, Parameter :: iHC_NH3=12
   Integer, Parameter :: iHC_N2O=13
   Integer, Parameter :: iHC_N2O4=14
   Integer, Parameter :: iHC_NH4HSO4=15






   End MODULE HC


   MODULE HCprop
   Use Par

   Double Precision Zfact(0:3),FugCoef(0:3,NoHC)
   Double PRecision HCy(NoHC),Fugy(NoHC)  ! Contains y and fugcoef to be used for Henrys law calc. Gas for single gas and 2-phase. Liq when liq only.



   End MODULE HCprop
