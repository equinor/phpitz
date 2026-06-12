!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to reactions in CO2
!
!************************************************************


MODULE GasRx

   Use Par

   Implicit None

   Double PRecision, Parameter :: RxLim=1.d-15

   Double Precision CompH(NoDefHC),CompS(NoDefHC)

! Parameters needed to define reactions
   Integer, Parameter :: MaxGasRx=20  ! ASsuming maximum 20 reactions
   Integer, Parameter :: MaxCompInRx=10 !Assuming maximum 10 components in a reaction
   Integer RxInfo(0:MaxGasRx)             ! 0=number of reactions, 1... is number of species in each reaction
   Integer RxComp(MaxGasRx, MaxCompInRx)  ! Component number in each reaction
   Integer RxCoef(MaxGasRx, MaxCompInRx)  ! Reaction coefficient
   Integer CompRx(NoDefHC,MaxGasRx)       ! Coefficient in each reaction for each specie


   Integer N_GasRxm,IncludeRxINP(MaxGasRx)

   Double Precision dH(MaxGasRx),dS(MaxGasRx),dG(MaxGasRx),lnk(MaxGasRx)


   Integer GoVal(MaxGasRx)
   Double Precision dNMax(MaxGasRx,2)

   Double Precision Q(MaxGasRx),dQQ(MaxGasRx,MaxGasRx)
!   Double Precision D(MaxGasRx)


End MODULE GasRx