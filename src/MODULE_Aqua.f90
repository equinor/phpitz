!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to Pitzer model
!
!************************************************************


MODULE Aqua
   Use Par


! Data for aqueous species used for model setup
   Character AqName(NoAq)*10              ! Name of aqueous specie
   Integer AqCharge(NoAq)                 ! Charge of aqueous specie
   Integer AqFeed(NoAq)                   ! Defines if specie is part of feed (given by analysis)
   Double Precision AqMw(NoAq)            ! Molecular weight of specie


! Numbering of Aq species      ManualUpdate - add new compnents here
   Integer, Parameter :: iAq_H2O=1
   Integer, Parameter :: iAq_H=2
   Integer, Parameter :: iAq_OH=3
   Integer, Parameter :: iAq_H2SO4=4
   Integer, Parameter :: iAq_HSO4=5
   Integer, Parameter :: iAq_SO4=6
   Integer, Parameter :: iAq_HNO3=7
   Integer, Parameter :: iAq_NO3=8
   Integer, Parameter :: iAq_NH3=9
   Integer, Parameter :: iAq_NH4=10
   Integer, Parameter :: iAq_CO2=11
   Integer, Parameter :: iAq_HCO3=12
   Integer, Parameter :: iAq_CO3=13



! Paremters for volumes
   Double Precision Eq_Vol_Par(NoAq,10),Vm(NoAq)


! Variables
   Double Precision AqN(NoAq)             ! Total molerate (mole/d) of aqueous species, for those were mass balance is defined
   Double Precision Aqm(NoAq)             ! Concentration of aqueous species in molality (mole/kg solvent)
   Double Precision AqP(NoAq)             ! Molerate of precipitated species (mole/d)
   Double Precision Aqg(NoAq)             ! Activity coefficients in aqueous phase. For water/MEG/EtOH etc it is activity.
   Double Precision Solvent               ! Amount of solvent in aqueous phase (kg)
   Double Precision AQN1(NoAq),AqN2(NoAq),AqNi(1:10,NoAq)
   Double Precision Aq_IS,Aq_TDS,Aq_EN,Aq_Dens
   Double Precision AqNT


!! Solvent parameters
!   Double Precision  Dens_H2O                ! Density of pure water at P and T
!   Double PRecision  Vol_H2O                 ! Molare volume of water at P and T (from Dens_H2O)
!   Double PRecision  Pvap_H2O                ! Vapour presssure of pure water at P and T
!   Double Precision  Fvap_H2O                ! Fugacity of water to be used in Raoults law
!   Double PRecision  dDdP,dlnDdP,dInvDdP,dehowdP     ! Various derivatives of dielectric property of water.



End MODULE Aqua