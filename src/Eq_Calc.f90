
!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Eq_Calc.f90
!  DESCRIPTION :  Calculate equilibrium constants and related variables
!  CREATED     :  2013.01.29  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  Eq_Therm_Calc   -  Calculates the thermodynamic equilibrium constants
!                 CalcVoles       -  Calculates specie volumes for further calculations of dV (pressure dependence)
!                 CalcdV          -  Calculate volume changes using volumes calculated in CalcVoles
!                 SolventProps    -  Calculates vapour pressures of solvent(s) and density
!                 Eq_Therm_Func   -  Temperature function routines for equilibrium constants
!
!****************************************************************************


   SubRoutine Eq_Therm_Calc
!****************************************************************************
!  SUBROUTINE  :  Eq_Therm_Calc
!  DESCRIPTION :  Calculates thermodynamic equilibrium constants
!  VARIABLES   :
!
!  CREATED     :  2013.01.29  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
      Use Eq
      Use GlobVar
      Use Pitzer
      Implicit None

! Functions
      Double Precision Eq_Therm_Func

      Double Precision Poynt(NoEq)
      COmmon /CPoynt/Poynt

! Local variables
      Integer i,EqN,Range
      Double Precision T,dV,dK,Rcm,RT,Pcorr,Pg,Pgsat,T_Low,T_Hi,tc,j

!----Start-------------------------------------------------------------------
      T=Temp
      tc=Temp-273.15d0

      Rcm=R*10.d0
      RT=Rcm*Temp

      Poynt=0.d0

      If (T.LT.373.15d0) Then                               ! Set pgauge - excess pressure. Pref=1 at t<100C, Pref = Pvap at t>100C
         Pg=Ptot-1.d0
      Else
         Pg=Ptot-Pvap_H2O
      EndIf
      Pgsat=Ptot-Pvap_H2O

      Call CalcVoles                                        ! Calculate partial molar volumes as function of P, T and I

      Do EqN=1,NoEq
         Range=1
         If (Eq_Type(EqN,0).GT.1) Then
            Do i=1,Eq_Type(EqN,0)
               T_Low=Eq_T_Range(EqN,i,1)+273.15d0
               T_Hi=Eq_T_Range(EqN,i,2) +273.15d0
               IF (T.GT.T_low .AND. T.LE.T_Hi) Then
                  Range=i
                  Goto 100
               EndIf
            EndDo
         EndIf
100      Continue

         Eq_Therm(EqN)=Eq_Therm_Func(Eq_Type(EqN,Range),Eq_T_Par(EqN,1:8,Range),T)        ! Therm eq const as f(T)


         dV=0.d0 ; dK=0.d0

         If(Eq_TypeV(EqN).EQ.3) Then                                            ! ScaleSim dV temperature function
            dV=Eq_Therm_Func(Eq_TypeV(EqN),Eq_V_Par(EqN,1:8),T)                  ! Pressure ariables volume and compressibility
            dK=Eq_Therm_Func(Eq_TypeKappa(EqN),Eq_K_Par(EqN,1:8),T)*1.d-3

            dV=dV + (Eq_IS_Par(EqN,1) + Eq_IS_Par(EqN,2)*tc)*IS                  ! Correction of volume for ionic strength

         ElseIf (Eq_TypeV(EqN).EQ.4) Then                                        ! Appelo dV functio
            Call CalcdV(EqN)
            dV=Eq_dvolm(EqN)

         ElseIf (Eq_TypeV(EqN).EQ.9) Then  ! Duan eq for calculation of dV and dK
            dV=Eq_Therm_Func(Eq_TypeV(EqN),Eq_V_Par(EqN,1:8),T)                  ! Pressure ariables volume and compressibility
            dK=Eq_Therm_Func(Eq_TypeKappa(EqN),Eq_K_Par(EqN,1:8),T)

            dV=dV + (Eq_IS_Par(EqN,1) + Eq_IS_Par(EqN,2)*tc)*IS                  ! Correction of volume for ionic strength
         EndIf

         If (Eq_TypeRx(EqN).EQ.1 .OR. Eq_TypeRx(EqN).EQ.2) Then                  ! Pressure correction aqueous reaction or salt solubility Pref=1 bar or Psat
            Pcorr=dexp(Pg/RT*(-dV + dK*Pg/2.d0))

         ElseIf (Eq_TypeRx(EqN).EQ.3) Then                                       ! Pressure correction Henrys law. Pref=Psat
            Pcorr=dexp(Pgsat/RT*(-dV + dK*Pgsat/2.d0))
         Else
            Pcorr=1.d0
         EndIf
         Poynt(EqN)=Pcorr

         Eq_Therm(EqN)=Eq_Therm(EqN)*Pcorr                                    ! Final thermodynamic K at P and T


      EndDo


   End

   SubRoutine CalcVoles
!****************************************************************************
!  SUBROUTINE  :  CalcVoles
!  DESCRIPTION :  Calculate partial molare volumes - Function from Appolo
!  VARIABLES   :
!
!  CREATED     :  2016.04.11  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Use Eq
   Use Pitzer
   Use Aqua
   Use SolventPar
   Implicit None

!-Local variables
   Integer i
   Double Precision T,P,p1,p2,RIS

!-Start----------------------------------------------------------------------

   T=Temp
   P=Ptot
   RIS=dsqrt(IS)

   Do i=1,NoAq
      Vm(i)=41.84d0*(0.1d0*Eq_Vol_Par(i,1)  + 1.d2*Eq_Vol_Par(i,2)/(2600.d0+P)  + Eq_Vol_Par(i,3)/(T-228.d0) + 1.d4*Eq_Vol_Par(i,4)/((2600.d0+P)*(T-228.d0)) + Eq_Vol_Par(i,5)*1.d5*dInvDdP)

      p1=Av*0.5d0*AqCharge(i)**2*RIS/(1.d0 + Eq_Vol_Par(i,6)*By*RIS)
      P2=(Eq_Vol_Par(i,7) + Eq_Vol_Par(i,8)/(T-228.d0) + Eq_Vol_Par(i,9)*(T-228.d0))*IS**Eq_Vol_Par(i,10)

      Vm(i)=Vm(i) + p1 + p2
   EndDo
   Vm(iAq_H2O)=Vol_H2O

   End

   SubRoutine CalcdV(EqN)
!****************************************************************************
!  SUBROUTINE  :  CalcdV
!  DESCRIPTION :  Calculate dV for reactions from molar volumes
!  VARIABLES   :
!
!  CREATED     :  2016.04.11  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use Eq
   Use Aqua
!   Use Salt
   Implicit None


! Local variables
   Integer EqN,j,ion,nyx,SaltNo

!-Start----------------------------------------------------------------------

   Eq_dVolm(EqN)=0.d0                              ! Start by setting it to zero

   !Do j=1,Eq_spec(EqN,0)                           ! Loop over all species in eq
   !   Ion=Eq_spec(EqN,j)                           ! Set ion no
   !   nyx=Eq_exp(EqN,j)                            ! Set exponent for that ion
   !
   !   Eq_dVolm(EqN)=Eq_dVolm(EqN) + Vm(Ion)*(nyx)   ! Add partial molar volumes
   !EndDo
   !
   !If (SaltInEq(EqN).GT.0) Then                       ! If equilibira contain salt
   !   SaltNo=SaltInEq(EqN)                            ! Find which salt it contains
   !   Eq_Dvolm(EqN)=Eq_Dvolm(EqN) - SaltVol(SaltNo)   ! Correct for salt volume
   !EndIf

   End

   SubRoutine SolventProps
!****************************************************************************
!  SUBROUTINE  :  SolventProps
!  DESCRIPTION :  Calculate properties of solvent. Water + others
!  VARIABLES   :
!
!  CREATED     :  2015.04.14  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Use SolventPar
   Use Eq
   Implicit None
! Functions
      Double Precision Eq_Therm_Func

! Local variables
   Double Precision TK,tc,FiRefH2O,PyH2O,Rcm
   Double Precision, Parameter :: MwH2O=18.016d0

!-Start----------------------------------------------------------------------
   TK=Temp
   tC=Temp-273.15d0
   Rcm=R*10.d0

!Calculate density
   Dens_H2O=(DensPar(1) + DensPar(2)*tc + DensPar(3)*tc**2 + DensPar(4)*tc**3 + DensPar(5)*tc**4 + DensPar(6)*tc**5) / (1.d0 + DensPar(7)*tc)
   Dens_H2O=Dens_H2O + (DensPar(8) + DensPar(9)*tc + DensPar(10)*tc**2 + DensPar(11)*tc**3)*(Ptot-1.d0)

   Vol_H2O=MwH2O/Dens_H2O*1.d3

! Calculate vapour pressure
   Pvap_H2O=Eq_Therm_Func(1,VapPar(1:8),TK)                                               ! Vapour pressure of water (bar)

! Reference fugacity of water vapour at Pvap
   FiRefH2O=FiPar(1) + FiPar(2)*TK + FiPar(3)*TK**2 + FiPar(4)*TK**3 + FiPar(5)*TK**4     ! Reference fucacity

! Poynting correction
   PyH2O=DEXP(Vol_H2O*(Ptot-Pvap_H2O)/(Rcm*TK))

! Calculate fugacity
   Fvap_H2O=Pvap_H2O*FiRefH2O*PyH2O                                                       ! Water fugacity

   End

   Double Precision Function Eq_Therm_Func(EqType,Par,T)
!****************************************************************************
!  FUNCTION    :  Eq_Therm_Func
!  DESCRIPTION :  Calculates the thermodynamic equilibrium constant
!  VARIABLES   :  EqType - Type of equation
!                 Par   - Parameters in the eqation
!                 T   - Temperature in Kelvin
!
!  CREATED     :  2013.01.29  Baard Kaasa www.scalesim.no
!  CHANGED     :  2016.01.11  BKaa: EqType 4 added
!
!****************************************************************************

! Call variables
      Integer EqType
      Double Precision Par(8),T,tc,Tref

!----Start-------------------------------------------------------------------
      Tref=25.d0 + 273.15d0

      If (abs(EqType).EQ.1) Then
         Eq_Therm_Func=dexp(Par(1)/T + Par(2) + Par(3)*dlog(T) + Par(4)*T + Par(5)/T**2 + Par(6)*T**2)

      ElseIf (abs(EqType).EQ.2) Then
         Eq_Therm_Func=Par(1) + T*(Par(2) + T*(Par(3) + T*(Par(4) + T*(Par(5) + T*Par(6)))))

      ElseIf (abs(EqType).EQ.3) Then
         tc=T-273.15d0
         Eq_Therm_Func=Par(1) + Par(2)*tc + Par(3)*tc**2 + Par(4)*tc**3 + Par(5)*tc**4 + Par(6)*tc**5 &
                     + Par(7)*tc**6 + par(8)*tc**7

      ElseIf (abs(EqType).EQ.4) Then
         Eq_Therm_Func=dexp(Par(1)/T + Par(2) + Par(3)*dlog(T) + Par(4)*T + Par(5)/(680.d0-T) + Par(6)*T**2)

      ElseIf (abs(EqType).Eq.5) Then  ! Møller 1988
         Eq_Therm_Func = Par(1)                         &
                        + Par(2)*T                      &
                        + Par(3)/T                      &
                        + Par(4)*dlog(T)                &
                        + Par(5)/(T-263.d0)             &
                        + Par(6)*T**2                   &
                        + Par(7)/(680.d0-T)             &
                        + Par(8)/(T-227.d0)

         Eq_Therm_Func=dexp(Eq_Therm_Func)

      ElseIf (abs(EqType).EQ.6) Then
         tc=T-273.15d0
         Eq_Therm_Func=exp(Par(1) + tc*(Par(2) + tc*(Par(3) + tc*(Par(4) + tc*(Par(5) + tc*Par(6))))))

      ElseIf (abs(EqType).EQ.7) Then     ! Christov 2004
         Eq_Therm_Func = Par(1)                      &
                     + Par(2)*T                      &
                     + Par(3)*T**2                   &
                     + Par(4)*T**3                   &
                     + Par(5)/T                      &
                     + Par(6)*dlog(T)                &
                     + Par(7)/(T-263.d0)             &
                     + Par(8)/(680.d0-T)
         Eq_Therm_Func=dexp(Eq_Therm_Func)

      ElseIf(Abs(EqType).EQ.8) Then  ! Dorn
         Eq_Therm_Func = Par(1)                           &
                   + Par(2)*(1.d0/T - 1.d0/Tref)   &
                   + Par(3)*dlog(T/Tref)           &
                   + Par(4)*(T-Tref)               &
                   + Par(5)*(T**2 - Tref**2)       &
                   + Par(6)*Dlog(T-225.d0)

         Eq_Therm_Func=dexp(Eq_Therm_Func)

      ElseIf (Abs(EqType).EQ.9) Then ! Duan dV and dK function
         Eq_Therm_Func=Par(1)/T + Par(2) + Par(3)*dlog(T) + Par(4)*T + Par(5)/T**2 + Par(6)*T**2

      ElseIf (Abs(EqTYpe).EQ.10) Then  ! plynomal of Tc then exp
                  tc=T-273.15d0
         Eq_Therm_Func=Par(1) + Par(2)*tc + Par(3)*tc**2 + Par(4)*tc**3 + Par(5)*tc**4 + Par(6)*tc**5 &
                     + Par(7)*tc**6 + par(8)*tc**7

         Eq_Therm_Func=dexp(Eq_Therm_Func)

      EndIf




      If (EqType.LT.0) Then               ! Need the inverse
         Eq_Therm_Func=1.d0/Eq_Therm_Func
      EndIf

   End



