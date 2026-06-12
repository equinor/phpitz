
!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to aqueous equilibria
!
!************************************************************

MODULE Eq
! Parameters  ManualUpdate: Check that these constants are updated
   Integer, Parameter :: NoEq=11              ! Number of equilibrim reactions in model
   Integer, Parameter :: MaxSpecInEq=4       ! Maximum species involved in an equilibrium rx
   Integer, Parameter :: NRange=2
   Integer, Parameter :: NGasSol=4           ! Number of gas/HC components that dissolve in water. Used to calculate Henrys constant.

! Variables
   Integer Eq_Type(NoEq,0:NRange)            ! Which type of temperature equation used, and for NRange ranges
   Integer Eq_TypeRx(NoEq)                   ! What type of reaction this equilibrum constant is used for
   Integer Eq_TypeV(NoEq)                    ! Function type for volume change
   Integer Eq_TypeKappa(NoEq)                ! Function type for compressibility
   Double Precision Eq_Therm(NoEq)           ! Thermodynamic equilibrium constants, f(T,P)
   Double Precision Eq_Stoec(NoEq)           ! Stoechiometric equilibrum constants, f(T,P,m)
   Double Precision Eq_Gam(NoEq)             ! Gamma-factor, all activity coefficients for a given equilibria
   Double Precision Eq_spec(NoEq,0:MaxSpecInEq)       ! Contains number of species (index=0) and which species for each equilibria (index(1-max)
   Double Precision Eq_exp(NoEq,MaxSpecInEq) ! Contains the number of each species, i.e. the exponential therm of each species
   Character  Eq_Name(NoEq)*8
   Integer GasSolub(NGasSol,2)

! Numbering of equilibria  ManualUpdate Add new equilibrium names here
   Integer, Parameter :: Eqn_Kw=1
   Integer, Parameter :: Eqn_K1H2SO4=2
   Integer, Parameter :: Eqn_K2H2SO4=3
   Integer, Parameter :: Eqn_KHNO3=4
   Integer, Parameter :: Eqn_KH_H2SO4=5
   Integer, Parameter :: Eqn_KH_HNO3=6
   Integer, Parameter :: Eqn_KH_NH3=7
   Integer, Parameter :: Eqn_KNH4=8
   Integer, Parameter :: Eqn_KHCO2=9
   Integer, Parameter :: Eqn_K1CO2=10
   Integer, Parameter :: Eqn_K2CO2=11



   !Integer, Parameter :: Eqn_Kw           = 1
   !Integer, Parameter :: Eqn_Ksp_NaCl     = 2
   !Integer, Parameter :: Eqn_KHSO4        = 3
   !Integer, Parameter :: Eqn_KassCaSO4o   = 4
   !Integer, Parameter :: Eqn_KspCaSO4G    = 5
   !Integer, Parameter :: Eqn_KspCaSO4H    = 6
   !Integer, Parameter :: Eqn_KspCaSO4A    = 7
   !Integer, Parameter :: Eqn_Ksp_BaSO4    = 8
   !Integer, Parameter :: Eqn_Ksp_KCl      = 9
   !Integer, Parameter :: Eqn_Ksp_Na2SO4   =10
   !Integer, Parameter :: Eqn_Ksp_SrSO4    =11
   !Integer, Parameter :: Eqn_Ksp_KOH2H2O  =12
   !Integer, Parameter :: Eqn_Ksp_KOHH2O   =13
   !Integer, Parameter :: Eqn_Ksp_K2SO4    =14
   !Integer, Parameter :: Eqn_Ksp_CaOH2    =15
   !Integer, Parameter :: Eqn_K_HAc        =16
   !Integer, Parameter :: Eqn_K_HFormic    =17
   !Integer, Parameter :: Eqn_K_Hprop      =18
   !Integer, Parameter :: Eqn_K_Hbut       =19
   !Integer, Parameter :: Eqn_Ksp_NaAc_3H2O=20
   !Integer, Parameter :: Eqn_Ksp_NaAc     =21
   !Integer, Parameter :: Eqn_Ksp_NaFormic_3H2O=22
   !Integer, Parameter :: Eqn_Ksp_NaFormic_2H2O=23
   !Integer, Parameter :: Eqn_Ksp_NaFormic =24
   !Integer, Parameter :: Eqn_KH_CO2       =25
   !Integer, Parameter :: Eqn_K1_CO2       =26
   !Integer, Parameter :: Eqn_K2_CO2       =27
   !Integer, Parameter :: Eqn_KH_H2S       =28
   !Integer, Parameter :: Eqn_K1_H2S       =29
   !Integer, Parameter :: Eqn_KH_CH4       =30
   !Integer, Parameter :: Eqn_Ksp_CaCO3c   =31
   !Integer, Parameter :: Eqn_Ksp_NaHCO3   =32
   !Integer, Parameter :: Eqn_Ksp_Na2CO3A  =33
   !Integer, Parameter :: Eqn_Ksp_Na2CO3M  =34
   !Integer, Parameter :: Eqn_Ksp_Na2CO3D  =35


! Parameters - Used in the calculations of the values
   Double Precision Eq_T_Par(NoEq,8,NRange)       ! Parameters fore each T-function, for up to NRange ranges.
   Double Precision Eq_T_Range(NoEq,NRange,2)     ! Temperature range for functions
   Double Precision Eq_V_Par(NoEq,8)
   Double Precision Eq_K_Par(NoEq,8)
   Double Precision Eq_IS_Par(NoEq,2)
   Double Precision Eq_dV(NoEq),Eq_dVolm(NoEq)
   Integer          SaltINEq(NoEq)                ! Which saltNo that is in eq. zero if no salt

! Solvent parameters
   Double Precision  Dens_H2O                ! Density of pure water at P and T
   Double PRecision  Vol_H2O                 ! Molare volume of water at P and T (from Dens_H2O)
   Double PRecision  Pvap_H2O                ! Vapour presssure of pure water at P and T
   Double Precision  Fvap_H2O                ! Fugacity of water to be used in Raoults law
   Double PRecision  dDdP,dlnDdP,dInvDdP,dehowdP     ! Various derivatives of dielectric property of water.




End MODULE Eq


!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to aqueous species
!
!************************************************************


!MODULE Aqua
!   Use Par
!   SAVE
!
!! Data for aqueous species used for model setup
!   Character AqName(NoAq)*10              ! Name of aqueous specie
!   Integer AqCharge(NoAq)                 ! Charge of aqueous specie
!   Integer AqFeed(NoAq)                   ! Defines if specie is part of feed (given by analysis)
!   Double Precision AqMw(NoAq)            ! Molecular weight of specie
!
!! Numbering of Aq species
!   Integer, Parameter :: iAq_H2O=1
!   Integer, Parameter :: iAq_H=2
!   Integer, Parameter :: iAq_OH=3
!   Integer, Parameter :: iAq_Na=4
!   Integer, Parameter :: iAq_Cl=5
!   Integer, Parameter :: iAq_Alk=6
!   Integer, Parameter :: iAq_K=7
!   Integer, Parameter :: iAq_SO4=8
!   Integer, Parameter :: iAq_HSO4=9
!   Integer, Parameter :: iAq_Ca=10
!   Integer, Parameter :: iAq_CaSO4o=11
!   Integer, Parameter :: iAq_Mg=12
!   Integer, Parameter :: iAq_Ba=13
!   Integer, Parameter :: iAq_Sr=14
!   Integer, Parameter :: iAq_Fe=15
!   Integer, Parameter :: iAq_Hac=16
!   Integer, Parameter :: iAq_Ac=17
!   Integer, Parameter :: iAq_Hformic=18
!   Integer, Parameter :: iAq_Formic=19
!   Integer, Parameter :: iAq_Hprop=20
!   Integer, Parameter :: iAq_Prop=21
!   Integer, Parameter :: iAq_Hbut=22
!   Integer, Parameter :: iAq_But=23
!   Integer, Parameter :: iAq_CO2=24
!   Integer, Parameter :: iAq_HCO3=25
!   Integer, Parameter :: iAq_CO3=26
!   Integer, Parameter :: iAq_H2S=27
!   Integer, Parameter :: iAq_HS=28
!   Integer, Parameter :: iAq_CH4=29
!
!
!
!! Paremters for volumes
!   Double Precision Eq_Vol_Par(NoAq,10),Vm(NoAq)
!
!! Variables
!   Double Precision AqN(NoAq)             ! Total molerate (mole/d) of aqueous species, for those were mass balance is defined
!   Double Precision Aqm(NoAq)             ! Concentration of aqueous species in molality (mole/kg solvent)
!   Double Precision AqP(NoAq)             ! Molerate of precipitated species (mole/d)
!   Double Precision Aqg(NoAq)             ! Activity coefficients in aqueous phase. For water/MEG/EtOH etc it is activity.
!   Double Precision Solvent               ! Amount of solvent in aqueous phase (kg)
!
!   Double Precision AQN1(NoAq),AqN2(NoAq),AqNi(1:10,NoAq)
!
!   Double Precision Aq_IS,Aq_TDS,Aq_EN,Aq_Dens
!
!   Double Precision AqNT
!
!End MODULE Aqua



!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to salts
!
!************************************************************



!MODULE Salt
!   Use Par
!   SAVE
!! Parameters used to define size of variables
!   Integer, Parameter :: MaxSpecInSalt=5     ! Maximum number of species in a salt
!
!   Integer, Parameter :: Indx_NaCl     = 1
!   Integer, Parameter :: Indx_CaSO4G   = 2
!   Integer, Parameter :: Indx_CaSO4H   = 3
!   Integer, Parameter :: Indx_CaSO4A   = 4
!   Integer, Parameter :: Indx_BaSO4    = 5
!   Integer, Parameter :: Indx_KCl      = 6
!   Integer, Parameter :: Indx_Na2SO4   = 7
!   Integer, Parameter :: Indx_SrSO4    = 8
!   Integer, Parameter :: Indx_KOH2H2O  = 9
!   Integer, Parameter :: Indx_KOHH2O   =10
!   Integer, Parameter :: Indx_K2SO4    =11
!   Integer, Parameter :: Indx_CaOH2    =12
!   Integer, Parameter :: Indx_NaAc_3H2O=13
!   Integer, Parameter :: Indx_NaAc     =14
!   Integer, Parameter :: Indx_NaFormic_3H2O=15
!   Integer, Parameter :: Indx_NaFormic_2H2O=16
!   Integer, Parameter :: Indx_NaFormic =17
!   Integer, Parameter :: Indx_CaCO3c   =18
!   Integer, Parameter :: Indx_NaHCO3   =19
!   Integer, Parameter :: Indx_Na2CO3A  =20
!   Integer, Parameter :: Indx_Na2CO3M  =21
!   Integer, Parameter :: Indx_Na2CO3D  =22
!
!
!! Model parameters from file
!   Integer           Nsalt                   ! Number of salts in current model, defined by Saltdata.dat
!   Character         SaltName(MaxNsalt)*10   ! Name of salts
!   Integer           Salt_eq(MaxNsalt)       ! Contains the equation no for the salt, used to calculate SR
!   Double Precision  SaltMw(MaxNsalt)        ! Molecular weight of salts
!   Double Precision  SaltWat(MaxNsalt)       ! Number of water molecules in salt
!   Double Precision  SaltVol(MaxNsalt)        ! Partial molare volume of salt
!
!   Integer           SaltAqSpec(MaxNsalt,0:MaxSpecInSalt) ! For each salt, which aqueous species
!   Double Precision  SaltAqSpecNum(MaxNsalt,MaxSpecInSalt) ! For each salt, the number of each aquesou specie
!
!   Integer           SaltHCSpec(MaxNsalt,0:MaxSpecInSalt) ! For each salt, which aqueous species
!   Double Precision  SaltHCSpecNum(MaxNsalt,MaxSpecInSalt) ! For each salt, the number of each aquesou specie
!
!   Integer           AqInSalt(NoAq,0:MaxNsalt)
!   Double Precision  AqInSaltNum(NoAq,MaxNsalt)
!
!   Integer           HCInSalt(NoHC,0:MaxNsalt)
!   Double Precision  HCInSaltNum(NoHC,MaxNsalt)
!
!
!
!   Double Precision  SRGoal(MaxNsalt)
!
!! Variables calculated in simulation
!   Double Precision  SR(MaxNsalt)            ! SR values for the different salts
!   Integer           PrecSalt(0:MaxNsalt)    ! Which salts are included in precipitation calculation
!   Double Precision  Prec(MaxNsalt),dPrec(MaxNsalt),MaxPrec(MaxNsalt)
!
!
!
!End MODULE Salt