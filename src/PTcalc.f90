!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  PCcalc.f90
!  DESCRIPTION :  Administrates calculation at given P and T
!  CREATED     :  2014.01.12  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  PTcalc      - Adm routine for calcuation of PT calculation
!                 CalcEQ      - Adm routine for calculation of equilibrium
!                 SetPTvar    - Adm routine for calcuation of all PT variables
!
!****************************************************************************



      SubRoutine PTcalc(pH)
!****************************************************************************
!  SUBROUTINE  :  PTcalc
!  DESCRIPTION :  Administrates calculation at given P and T
!  VARIABLES   :
!
!  CREATED     :  2014.01.12  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Path
   Use Aqua
   Use HC
   Use GlobVar
   Implicit None

 !  Integer PrecON
 !  Common /CPrecON/ PrecON

   Double Precision  N_aq(NoAq)
   Common /CN_Aq /N_Aq

   Double Precision pH

! Local variables
!--Start---------------------------------------------------------------------
 !  If (LOG) Write(LogFileCalc,*)'Calling SetPTVar'

   !Tc=10.d0
   !
   !Do i=1,20
   !   Tc=Tc-3.d0
   !   Temp=Tc+273.15d0
   !
   ! Call Eq_Therm_Calc
   !
   !EndDo

   Call SetPTVar           ! Calculate all variables that are function of P and T only

   Call CalcGas

   Call SolveWater(pH)


   End


   Subroutine CalcGas
!****************************************************************************
!  SUBROUTINE  :  CalcGas
!  DESCRIPTION :  TEMPORARY - Sets properties of CO2 rich phase
!  VARIABLES   :
!
!  CREATED     :  2025.09.22  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use HC
   Use HCProp
   Implicit None

   Integer Phase
   Double Precision Ni(NoHC)

! Local variables
!--Start---------------------------------------------------------------------


   Ni=HCN(PhGas,1:NoHC)
   Phase=0
   Call SinglePhase(Ni,Phase)


   HCy(1:NoHC)=HCN(PhGas,1:NoHC)/HCN(PhGas,0)



   End





   SubRoutine SetPTVar
!****************************************************************************
!  SUBROUTINE  :  SetPTVar
!  DESCRIPTION :  Calculate all variables only function of P and T
!  VARIABLES   :
!
!  CREATED     :  2015.04.15  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Implicit None

   Double Precision TempTemp

!-Start----------------------------------------------------------------------



   Call SolventProps       ! Calculate properties of solvent (density, vapour pressure, volume)

   TempTemp=Temp
   If (Temp-273.15d0.LT.-20.d0) Then
      Temp=-20.d0 + 273.15d0
   EndIf

   Call PitzerParCalc      ! Calculate Pitzer parameters as f(T,P)
   Temp=TempTemp

   Call WaterPar           ! Calculate water parameters used in Pitzer model
   Call CalcAphi           ! Calculate Debye-Hückel term

  ! TempTemp=Temp
  ! Temp=-10.d0 + 273.15d0

   Call Eq_Therm_Calc      ! Thermodynamic equilibrium constants

  ! Temp=TempTemp

   Call EOSconst           ! Sets the P-T parameters in EOS



   End



   SubRoutine SolveWater(pH)
!****************************************************************************
!  SUBROUTINE  :  SolveWater
!  DESCRIPTION :  Calculates pH and concentrations in a water for a given total composition
!  VARIABLES   :
!
!  CREATED     :  2025.08.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Use Aqua
   Use Par
   Use HC
   Use HCprop
   Use Eq
   Implicit None

   Double Precision  N_aq(NoAq)
   Common /CN_Aq /N_Aq

   External pHCalc

   Double Precision, Parameter :: pHacc=1.d-12
   Double Precision, Parameter :: pHit =1.0e-4
   Double Precision, Parameter :: pHlow=-4.d0
   Double Precision, Parameter :: pHHi=20.d0
   Double Precision, Parameter :: AqAcc=1.d-6

   Double Precision pH

! Local variables
   Integer ItNo, i, j, m, NGasSolub
   Double Precision Aqg_OLD(NoAq), f, f2, pHOld

   Double Precision MaxDiff, diff

   Double Precision ytar(10), ntar(10),nmov(10)
   Integer iHC(10)



!-Start----------------------------------------------------------------------

   F2=0.5d0
   Aqm=0.d0
   f=0.2d0
   pH=2.d0



   Do j=1,3000 ! Baard 3000
      Solvent=N_Aq(iAq_H2O)*AqMw(iAq_H2O)*1.d-3


      if (j.EQ.1) Then ! Make an ainitial guess for water phase ManualUpdate
         Aqm(iAq_H2SO4)=0.8d0*N_aq(iAq_H2SO4)/Solvent
         Aqm(iAq_HSO4) =0.2d0*N_aq(iAq_H2SO4)/Solvent
         Aqm(iAq_HNO3)=0.5d0*N_aq(iAq_HNO3)/Solvent
         Aqm(iAq_NO3) =0.5d0*N_aq(iAq_HNO3)/Solvent
         Aqm(iAq_NH3) =0.8d0*N_aq(iAq_NH3)/Solvent
      EndIf


      pHOld=7.d0

      Do Itno=1,1000

         MaxDiff=0.d0

         Call PitzerCalc

         If (ItNo.EQ.1) Then
            Aqg_OLD=Aqg
         Else
            Do i=1,NoAq
               diff=dabs(Aqg(i)-Aqg_OLD(i))
               If (diff.GT.MaxDiff) Then
                  MaxDiff=diff
               EndIf
               Aqg(i)=Aqg(i)*f + Aqg_OLD(i)*(1.d0-f)
            EndDo
            Aqg_OLD=Aqg
         EndIf

         Call Eq_Stoech


         Aqm(iAq_CO2)=Ptot*Eq_Stoec(Eqn_KHCO2)


         Call pHNewton(pHCalc,pHlow,pHHi,pHacc,pH)

         !Write(*,'(I3,F15.5,3E15.5)')Itno, pH, MaxDiff, dabs(pH-pHold), dabs(Aqg(iAq_H2O)-Aqg_OLD(iAq_H2O))

         If (dabs(pH-pHold).LT.pHit .AND. MaxDiff.LT.AqAcc .AND. Itno.GE.50) Goto 150
         pHOld=pH

      EndDo

      Write(*,*)'exceeding max 2', pH


150   Continue
      !Read(*,*)

      NGasSolub=4

   ! The remaining part will require ManualUpdate if number of components included in water-CO2 equilibrium is changed
      iHC(1)=iHC_H2O     ! NB water must be number 1
      iHC(2)=iHC_H2SO4
      iHC(3)=iHC_HNO3
      iHC(4)=iHC_NH3

   ! Caculate target mole fraction in gas phase ManualUpdate
      ytar(1)=Aqg(iAq_H2O)*Fvap_H2O/(Fugy(iHC_H2O)*Ptot)
      ytar(2)=Aqm(iAq_H2SO4)/(  Eq_Stoec(Eqn_KH_H2SO4) * Ptot )
      ytar(3)=Aqm(iAq_HNO3 )/(  Eq_Stoec(Eqn_KH_HNO3 ) * Ptot )
      ytar(4)=Aqm(iAq_NH3  )/(  Eq_Stoec(Eqn_KH_NH3  ) * Ptot )

   ! Calculate target moles in CO2 phase
      Do m=1,NGasSolub
         ntar(m)=ytar(m)*HCN(PhGas,0)
      EndDo

   ! Calculate how much to move of all components
      Do m=1,NGasSolub
         If (ntar(m).GT.HCN(PhGas,iHC(m))) Then  ! Evaporate
            If (abs(ntar(m)-HCN(PhGas,iHC(m))).GT.0.1d0*HCN(PhAq,iHC(m))) THen
               nmov(m)=0.1d0*HCN(PhAq,iHC(m))
            Else
               nmov(m)=ntar(m)-HCN(PhGas,iHC(m))
            EndIf
         Else
            If (abs(ntar(m)-HCN(PhGas,iHC(m))).GT.0.1d0*HCN(PhAq,iHC(m))) THen
               nmov(m)=-0.1d0*HCN(PhAq,iHC(m))
            Else
               nmov(m)=ntar(m)-HCN(PhGas,iHC(m))
            EndIf
         EndIf

      EndDo

      If (j.GE.50)    F2=0.25d0
      If (j.GE.500)   F2=0.1d0
      If (j.GE.1000)  F2=0.02d0
      If (j.GE.1500)  F2=0.005d0


      nmov=nmov*F2  ! Adjust to stabilize simulation and avoid oscillation

      Do m=1,NGasSolub
         HCN(PhGas,iHC(m))=HCN(PhGas,iHC(m)) + nmov(m)
         HCN(PhAq ,iHC(m))=HCN(PhAq ,iHC(m)) - nmov(m)
      EndDo
      HCN(PhAq,iHC_CO2)=(Aqm(iAq_CO2) + Aqm(iAq_HCO3) + Aqm(iAq_CO3))*Solvent
      HCN(PhGas,iHC_CO2)=HCN(PhTot,iHC_CO2) - HCN(PhAq,iHC_CO2)

      HCN(PhGas,0)=sum(HCN(Phgas,1:NoHC))
      HCN(PhAq ,0)=sum(HCN(PhAq ,1:NoHC))

   ! Update concentration in water phase ManualUpdate
      N_aq(iAq_H2O)  =HCN(PhAq,iHC_H2O)
      N_aq(iAq_H2SO4)=HCN(PhAq,iHC_H2SO4)
      N_Aq(iAq_HNO3) =HCN(PhAq,iHC_HNO3)
      N_Aq(iAq_NH3)  =HCN(PhAq,iHC_NH3)

     ! Write(*,'(I5,5F15.5)')j,nmov(1), HCN(0,iHC_H2O), HCN(PhGas,iHC_H2O), HCN(PhAq,iHC_H2O), ntar(1)


      If (j.GT.100) then   ! Minimum iteration before routine is alowed to end

         If (abs(nmov(1)) .LT. HCN(PhGas,iHC_H2O)*1.d-6  .AND. abs(nmov(1)) .LT. HCN(PhAq,iHC_H2O)*1.d-6  ) Then
            Goto 250
         EndIf

         If (ntar(1)/HCN(0,iHC_H2O).GT.1.d0) Then ! Converging to NO water phase
            HCN(PhAq,1:NoHC)=0.d0
            HCN(PhGas,1:NoHC)=HCN(0,1:NoHC)

            Do i=1,3
               HCN(i,0)=sum(HCN(i,1:NoHC))
            EndDo

            Goto 250
         EndIf
      EndIf

   EndDo
   Write(*,*)'Exceeding max 1'


250 Continue

   End


   SubRoutine pHCalc(pH,F,dF)
!****************************************************************************
!  SUBROUTINE  :  pHCalc
!  DESCRIPTION :  Used in pH iterations
!  VARIABLES   :
!
!  CREATED     :  2025.08.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use Aqua
   Use Par
   Use EQ

   Implicit None

! Call variables
   Double Precision pH,F,dF


   Double Precision  N_aq(NoAq)
   Common /CN_Aq /N_Aq

! Local variables
   Double Precision H,dAqm(NoAq)
   Double precision T1, N1, dT1, dN1 ! Help variables for H2SO4
   Double Precision T2, N2, dT2, dN2 ! Help variables for HNO3
   Double Precision T3, N3, dT3, dN3 ! Help variables for NH3
!-Start----------------------------------------------------------------------


   H=10**(-pH)
   Aqm(iAq_H)=H
   Aqm(iAq_OH)=Eq_Stoec(Eqn_Kw)/H

! H2SO4
   T1=N_Aq(iAq_H2SO4)/Solvent
   N1=H/Eq_Stoec(Eqn_K1H2SO4) + 1.d0 + Eq_Stoec(Eqn_K2H2SO4)/H

   Aqm(iAq_HSO4)=T1 /N1

   Aqm(iAq_H2SO4)=H*Aqm(iAq_HSO4)/Eq_Stoec(Eqn_K1H2SO4)
   Aqm(iAq_SO4)=Eq_Stoec(Eqn_K2H2SO4)*Aqm(iAq_HSO4)/H

! HNO3
   T2=N_Aq(iAq_HNO3)/Solvent
   N2=1.d0 + Eq_Stoec(Eqn_KHNO3)/H

   Aqm(iAq_HNO3)=T2/N2
   Aqm(iAq_NO3)=Eq_Stoec(Eqn_KHNO3)*Aqm(iAq_HNO3)/H


! NH3
   T3=N_Aq(iAq_NH3)/Solvent
   N3=1.d0 + Eq_Stoec(Eqn_KNH4)/H

   Aqm(iAq_NH4)=T3/N3
   Aqm(iAq_NH3)=Eq_Stoec(Eqn_KNH4)*Aqm(iAq_NH4)/H

! CO2
   Aqm(iAq_HCO3)=Eq_Stoec(Eqn_K1CO2)*Aqm(iAq_CO2)/H
   Aqm(iAq_CO3)=Eq_Stoec(Eqn_K2CO2)*Aqm(iAq_HCO3)/H


! ManualUpdate is needed if more acid/bases are added

   F=Aqm(iAq_OH) + Aqm(iAq_HSO4) + 2.d0*Aqm(iAq_SO4) + Aqm(iAq_NO3) - Aqm(iAq_H) - Aqm(iAq_NH4) + Aqm(iAq_HCO3) + 2.d0*Aqm(iAq_CO3)

! And the derivatives
   dAqm=0.d0

   dAqm(iAq_H)=1.d0
   dAqm(iAq_OH)=-Aqm(iAq_OH)/H

! H2SO4
   dT1=0.d0
   dN1=1.d0/Eq_Stoec(Eqn_K1H2SO4) - Eq_Stoec(Eqn_K2H2SO4)/H**2

   dAqm(iAq_HSO4)=(dT1*N1 - T1*dN1)/N1**2
   dAqm(iAq_SO4)= Eq_Stoec(Eqn_K2H2SO4)*(dAqm(iAq_HSO4)*H - Aqm(iAq_HSO4))/H**2

! HNO3
   dT2=0.d0
   dN2=-Eq_Stoec(Eqn_KHNO3)/H**2

   dAqm(iAq_HNO3)=(dT2*N2 - T2*dN2)/N2**2
   dAqm(iAq_NO3) = Eq_Stoec(Eqn_KHNO3)*(dAqm(iAq_HNO3)*H - Aqm(iAq_HNO3))/H**2

! NH3
   dT3=0.d0
   dN3=-Eq_Stoec(Eqn_KNH4)/H**2

   dAqm(iAq_NH4)=(dT3*N3 - T3*dN3)/N3**2
   dAqm(iAq_NH3) = Eq_Stoec(Eqn_KNH4)*(dAqm(iAq_NH4)*H - Aqm(iAq_NH4))/H**2

   dAqm(iAq_HCO3)=-Aqm(iAq_HCO3)/H
   dAqm(iAq_CO3) =-2.d0*Aqm(iAq_CO3)/H


! ManualUpdate is needed if more acid/bases are added
   dF=dAqm(iAq_OH)  + dAqm(iAq_HSO4) + 2.d0*dAqm(iAq_SO4) + dAqm(iAq_NO3) - dAqm(iAq_H) - dAqm(iAq_NH4) + dAqm(iAq_HCO3) + 2.d0*dAqm(iAq_CO3)


   dF =-dlog(10.D0)*H*dF

  ! Write(*,*)pH,F,dF

   End

