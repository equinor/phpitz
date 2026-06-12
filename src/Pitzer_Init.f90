!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Pitzer_Init.f90
!  DESCRIPTION :  Initialisation of Pitzer model and some general calculations
!  CREATED     :  2013.05.05  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  PitzerSetup - Sets up Pitzer model
!                 PitsParCalc - Calculates the Pitzer parameters from temperature function
!                 WaterPar    - Calculates parameters for water used in Pitzer model, density, dielectricy
!                 CalcAphi    - Calculates the Debue-Hückel parameter Aphi
!                 PitzTfunc   - General temperature function for Pitzer parameters.
!
!****************************************************************************


      SubRoutine PitzerSetup
!****************************************************************************
!  SUBROUTINE  :  PitzerSetup
!  DESCRIPTION :  Sets up Pitzer model
!  VARIABLES   :
!
!  CREATED     :  2013.05.04  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None
! Call variables


! Local variables
      Integer i
!-----Start-------------------------------------------------------------------

! Sort cations, anions and neutrals
      CatM=0
      AnX=0
      NeuN=0


      Do i=1,NoAq


         if (AqInPitz(i).EQ.1) Then
            If (AqCharge(i).GT.0) Then
               CatM(0)=CatM(0)+1
               CatM(CatM(0))=i
            ElseIf (AqCharge(i).LT.0) Then
               AnX(0)=AnX(0)+1
               AnX(AnX(0))=i
            Else
               NeuN(0)=NeuN(0)+1
               NeuN(NeuN(0))=i
            EndIf
         EndIf
      EndDo

   End


      SubRoutine PitzerParCalc
!****************************************************************************
!  SUBROUTINE  :  PitzerParCalc
!  DESCRIPTION :  Calculate the Pitzer parameters as temperature functions
!  VARIABLES   :
!
!  CREATED     :  2013.05.05  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Use GlobVar
      Implicit None

! Functions
      Double Precision PitzTfunc
! Local variables
      Integer M,M2,X,X2,N,Indx1,Indx2,Indx3
      Integer c,cp,a,ap
      Double Precision T,P
!-----Start-------------------------------------------------------------------

      Beta0  = 0.d0 ; Beta1  = 0.d0 ; Beta2  = 0.d0 ; Cphi   = 0.d0
      Theta  = 0.d0 ; Psi    = 0.d0 ; Lambda = 0.d0 ; Ksi    = 0.d0

      T=Temp
      P=Ptot


      Do M=1,CatM(0) ! Calculate Beta0, Beta1, Beta2 and C
         Do X=1,AnX(0)
            Indx1=CatM(M)
            Indx2=AnX(X)
            Beta0(Indx1,Indx2)=PitzTfunc(BetaCC(0,Indx1,Indx2,1:Npar),BetaCCEq(0,Indx1,Indx2))
            Beta1(Indx1,Indx2)=PitzTfunc(BetaCC(1,Indx1,Indx2,1:Npar),BetaCCEq(1,Indx1,Indx2))
            Beta2(Indx1,Indx2)=PitzTfunc(BetaCC(2,Indx1,Indx2,1:Npar),BetaCCEq(2,Indx1,Indx2))
            Cphi (Indx1,Indx2)=PitzTfunc(CphiCC(  Indx1,Indx2,1:Npar),CphiCCEq( Indx1,Indx2))
         EndDo
      EndDo

      Do M=1,CatM(0)-1        ! Psi for cc'a
         Do M2=M+1,CatM(0)
            Do X=1,AnX(0)
               c=CatM(M)
               cp=CatM(M2)
               a=AnX(X)

               Psi(c,cp,a)=PitzTfunc(PsiCC(c,cp,a,1:Npar),PsiCCEq(c,cp,a))
               Psi(cp,c,a)=Psi(c,cp,a)
            EndDo
         EndDo
      EndDo

      Do M=1,CatM(0)          ! psi for caa'
         Do x=1,AnX(0)-1
            Do X2=X+1,AnX(0)
               c=CatM(M)
               a=AnX(X)
               ap=AnX(X2)

               Psi(c,a,ap)=PitzTfunc(PsiCC(c,a,ap,1:Npar),PsiCCEq(c,a,ap))
               Psi(c,ap,a)=Psi(c,a,ap)
            EndDo
         EndDo
      EndDo


      Do M=1,CatM(0)-1        ! Theta for cc'
         Do M2=M+1,CatM(0)
            Indx1=CatM(M)
            Indx2=CatM(M2)

            Theta(Indx1,Indx2)=PitzTfunc(ThetCC(Indx1,Indx2,1:Npar),ThetCCEQ(Indx1,Indx2))
            Theta(Indx2,Indx1)=Theta(Indx1,Indx2)
         EndDo
      EndDo
      Do X=1,AnX(0)-1         ! Theta for aa'
         Do X2=X+1,AnX(0)
            Indx1=AnX(X)
            Indx2=AnX(X2)

            Theta(Indx1,Indx2)=PitzTfunc(ThetCC(Indx1,Indx2,1:Npar),ThetCCEq(Indx1,Indx2))
            Theta(Indx2,Indx1)=Theta(Indx1,Indx2)
         EndDo
      EndDo

      Do N=1,NeuN(0)          ! Lambda for nc
         Do M=1,CatM(0)
            Indx1=NeuN(N)
            Indx2=CatM(M)

            Lambda(Indx1,Indx2)=PitzTfunc(LambCC(Indx1,Indx2,1:Npar),LambCCEq(Indx1,Indx2))
         EndDo
      EndDo

      Do N=1,NeuN(0)          ! Lambda for na
         Do X=1,AnX(0)
            Indx1=NeuN(N)
            Indx2=AnX(X)

            Lambda(Indx1,Indx2)=PitzTfunc(LambCC(Indx1,Indx2,1:Npar),LambCCEq(Indx1,Indx2))
         EndDo
      EndDo

      Do N=1,NeuN(0)          ! Ksi for n-c-a
         Do M=1,CatM(0)
            Do X=1,AnX(0)
               Indx1=NeuN(N)
               Indx2=CatM(M)
               Indx3=AnX(X)

               Ksi(Indx1,Indx2,Indx3)=PitzTfunc(KsiCC(Indx1,Indx2,Indx3,1:Npar),KsiCCEq(Indx1,Indx2,Indx3))
            EndDo
         EndDo
      EndDo


   End


      SubRoutine WaterPar
!****************************************************************************
!  SUBROUTINE  :  WaterPar
!  DESCRIPTION :  Calculates water properties (density and dielectic constant)
!  VARIABLES   :
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use GlobVar
      Use Pitzer
      Use Aqua
      Use Eq

      Implicit None

! Local variables
      Double Precision D1000, C, B, tc

! Data
      Double Precision U(9)  ! Dielectric constant of water
      Data U / 3.4279d2, -5.0866d-3 , 9.469d-7 , -2.0525d0 , 3.1159d3 , -1.8289d2 , -8.0325d3 , 4.2142d6 , 2.1417d0 /

      Double Precision K(7)
      Data K/50.884917d0, 0.62590623d0, 0.0013848668d0, 0.000021603427d0, 0.000000072087667d0, 4.6545054d-10, 0.019859983d0/
!-----Start------------------------------------------------------------

! Dielectric constant of water calculated from:
!     Bradley and Pitzer, J. Phys. Chem., 83, 1979, p. 1599


      D1000=U(1)*dexp(U(2)*Temp + U(3)*Temp**2)
      C = U(4) + U(5)/(U(6) + Temp)
      B = U(7) + U(8)/Temp + U(9)*Temp
      Dielec = D1000 + C*dlog((B + Ptot)/(B + 1.d3))

      dDdP=C/(B+Ptot)
      dlnDdP=dDdP/Dielec
      dInvDdP=-dDdP/Dielec**2

      rhow=Dens_H2O

      tc=Temp-273.15d0
      dehowdP=(k(1) + k(2)*tc + k(3)*tc**2 + k(4)*tc**3 + k(5)*tc**4 + k(6)*tc**5) / (1.d0 + k(7)*tc)*1.d-6

   End


      SubRoutine CalcAphi
!****************************************************************************
!  SUBROUTINE  :  CalcAphi
!  DESCRIPTION :  Calculates the Debue Hückel Aphi parameter
!  VARIABLES   :
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :  2015.12.16  BKaa corrected parameters in accordanc with SI chem data
!****************************************************************************
      Use GlobVar
      Use Pitzer
      Use Eq
      Implicit None

      Double Precision, Parameter :: Pi=3.141592654
      Double Precision, Parameter :: NAv=6.02214179d23
      Double Precision, Parameter :: el=1.602176487d-19
      Double Precision, Parameter :: eps0=8.854187817d-12
      Double Precision, Parameter :: k=1.3806504d-23

! Local variabels
      Double Precision P1, P2
!-----Start------------------------------------------------------------
      P1=2.d0*Pi*Nav*rhow
      P2=el**2/(4.d0*Pi*Dielec*eps0*k*Temp)

      Aphi=1.d0/3.d0*dsqrt(P1)*P2**(1.5d0)

      By= Nav*el/1.d10*dsqrt(2.d0*rhow/(eps0*Dielec*R*Temp))

      Av=2.d0*Aphi*R*10.d0*Temp*(3.d0*dlnDdP - dehowdP )*1.013

   End


      Double Precision Function PitzTfunc(Par,EqN)
!****************************************************************************
!  SUBROUTINE  :  PitzTfunc
!  DESCRIPTION :  Temperature function for Pitzer parameters
!  VARIABLES   :
!
!  CREATED     :  2013.05.05  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use GlobVar
      Use Pitzer
      Implicit None
! Call variables
      Double Precision Par(Npar)
      Integer EqN

! Local variables
      Double Precision Tref
!-----Start------------------------------------------------------------
      Tref = 25.d0 + 273.15D0

      If (EqN.EQ.1) Then      ! Standard T-function
         PitzTfunc = Par(1)                           &
                   + Par(2)*(Temp - Tref)					&
                   + Par(3)*(Temp**2 - Tref**2)			&
                   + Par(4)*(1.d0/Temp - 1.d0/Tref)         &
                   + Par(5)*dlog(Temp/Tref)           &
                   + Par(6)*(1.d0/Temp**2 - 1.d0/Tref**2)

      ElseIf (EqN.EQ.2) Then  ! Møller 1988 T-function
         PitzTfunc = Par(1)                           &
                   + Par(2)*Temp                      &
                   + Par(3)/Temp                      &
                   + Par(4)*dlog(Temp)                &
                   + Par(5)/(Temp-263.d0)             &
                   + Par(6)*Temp**2                   &
                   + Par(7)/(680.d0-Temp)             &
                   + Par(8)/(Temp-227.d0)
      ElseIf (EqN.EQ.3) Then  ! Haarberg T-function
         PitzTfunc = Par(1)                           &
                  + Par(2)*(Temp-Tref)                &
                  + 0.5d0*Par(3)*(Temp-Tref)**2

      ElseIf (EqN.EQ.4) Then  ! Christov
         PitzTfunc = Par(1)                           &
                   + Par(2)*Temp                      &
                   + Par(3)*Temp**2                   &
                   + Par(4)*Temp**3                   &
                   + Par(5)/Temp                      &
                   + Par(6)*dlog(Temp)                &
                   + Par(7)/(Temp-263.d0)             &
                   + Par(8)/(680.d0-Temp)
      ElseIf (EqN.EQ.5) Then   ! Dorn
         PitzTfunc = Par(1)                           &
                   + Par(2)*(1.d0/Temp - 1.d0/Tref)   &
                   + Par(3)*dlog(Temp/Tref)           &
                   + Par(4)*(Temp-Tref)               &
                   + Par(5)*(Temp**2 - Tref**2)       &
                   + Par(6)*Dlog(Temp-225.d0)
      Else
         PitzTfunc=0.d0
      EndIf
   End


