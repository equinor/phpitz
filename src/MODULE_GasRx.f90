!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to reactions in CO2
!
!************************************************************




   MODULE GasRx
   Use dpDef
   Use Par  ! Number of components etc
   Implicit None


   Integer, Parameter :: ncomp  = 20
   Integer, Parameter :: nelem  = 4
   Integer, Parameter :: maxnRx = 12


! Thermo data
   Real(dp) :: CompH(NoComp),CompS(NoComp)  ! Enthalpy and entropy for components
   Real(dp) :: RxH(MaxnRx), RxS(maxnRx), RxG(maxnRx),RxlnK(maxnRx)


! Phase type constants
   Integer, parameter :: GAS   = 1
   Integer, parameter :: SOLID = 2


! Component information
   Character(len=12) :: comp_name(ncomp)
   Integer :: CompType(NoComp)      ! PHASE: GAS=1, SOLID=2
   logical  :: is_reactive(ncomp)

! Reaction information
   Integer :: nRx

   Integer :: IncludeRx(0:maxnRx)
   Real(dp) :: stoich(ncomp,maxnRx)


!  real(dp) :: Ntot



   End MODULE GasRx





Module gibbs_solver
!==============================================================================
   Use dpDef

   Integer, parameter  :: max_iter = 2000
   Real(dp), parameter :: tol_conv = 1.0e-10_dp
   Real(dp), parameter :: Nmin     = 1.0e-40_dp

  ! Her kan vi ogås legge inn en kommentar

End Module gibbs_solver





