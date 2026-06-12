!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Pitzer_Calc.f90
!  DESCRIPTION :  Calculates activity coefficients with Pitzer model
!  CREATED     :  2013.05.05  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  PitzerCalc - Administrative routine for calculation of activity coefficients
!                 DH         - Calculates the Debye-Hückel term fy
!                 CalcF      - Calculates the F, includes the Debye-Hückel
!                 CalcISansZ - Calculates ionic strengt and the Z parameter in the Pitzer model
!
!****************************************************************************


      SubRoutine PitzerCalc
!****************************************************************************
!  SUBROUTINE  :  Pitzer
!  DESCRIPTION :  Calculates activity coefficients with Pitzer model
!  VARIABLES   :
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

!-----Start-------------------------------------------------------------------


      Aqg=1.d0

      Call CalcISandZ

      If (IS.GT.0.d0) Then
         Call DH
         Call CalcFandSumC

         Call CatIons
         Call AnIons
         Call Netural
         Call Osmo
      EndIf

   End

      SubRoutine DH
!****************************************************************************
!  SUBROUTINE  :  DH
!  DESCRIPTION :  Calculates the Debue-Hückel term in Pitzer model
!  VARIABLES   :  IS - Ionic strength
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Local variables
      Double Precision sqrI,N
!-----Start-------------------------------------------------------------------

      N=1.d0 + Debye_b*sqrIS   ! Precalc of expression usedseveral times

      fy=-Aphi*(sqrIS/N + 2.d0/Debye_b*dlog(N))


      End

      SubRoutine CalcFandSumC
!****************************************************************************
!  SUBROUTINE  :  CalcFandSumC
!  DESCRIPTION :  Calculates F and Sum(mc*ma*Cca) in Pitzer model
!  VARIABLES   :
!
!  CREATED     :  2013.05.13  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Functions
      Double PRecision BMX,Phi,CMX,E_Theta,E_ThetaP,J

   Logical OnlyDebye
   Common /COnlyDebye / OnlyDebye

! Local variables
      Double Precision mmB,mmPhic,mmPhia
      Integer c,cp,a,ap,Indxc,Indxa,Indxcp,Indxap
!-----Start-------------------------------------------------------------------
      F=0.d0
      mmB=0.d0
      mmPhic=0.d0
      mmPhia=0.d0
      SumCca=0.d0

      Do c=1,CatM(0)
         Do a=1,AnX(0)
            Indxc=CatM(c)
            Indxa=AnX(a)
            mmB   =mmB    + Aqm(Indxc)*Aqm(Indxa)*BMX(Indxc,Indxa,1) ! Derivative of BMX
            SumCca=SumCca + Aqm(Indxc)*Aqm(Indxa)*CMX(Indxc,Indxa)
         EndDo
      EndDo

      Do c=1,CatM(0)-1
         Do cp=c+1,CatM(0)
            Indxc =CatM(c)
            Indxcp=CatM(cp)
            mmPhic=mmPhic + Aqm(Indxc)*Aqm(Indxcp)*Phi(Indxc,Indxcp,1) ! Derivative of Phi
         EndDo
      EndDo

      Do a=1,AnX(0)-1
         Do ap=a+1,AnX(0)
            Indxa =AnX(a)
            Indxap=AnX(ap)
            mmPhia=mmPhia + Aqm(Indxa)*Aqm(Indxap)*Phi(Indxa,Indxap,1) ! Derivative of Phi
         EndDo
      EndDo


      F=fy + mmB + mmPhic + mmPhia

      If (OnlyDebye)  F=fy



     ! Write(*,'(A,5F15.8)')'F',fy,mmB,mmPhic,mmPhia,F
      End

      SubRoutine CatIons
!****************************************************************************
!  SUBROUTINE  :  CatIons
!  DESCRIPTION :  Calculates activity coefficienst for cations
!  VARIABLES   :
!
!  CREATED     :  2013.05.13  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Functions
      Double Precision BMX, CMX,Phi,E_ThetaP

   Logical OnlyDebye
   Common /COnlyDebye / OnlyDebye

! Local variables
      Integer M,IndxM,a,Indxa,c,Indxc,ap,Indxap,n,Indxn
      Double PRecision SumBZC,SumPhi,SumPsi,SumL,SumLL
!-----Start-------------------------------------------------------------------


      Do M=1,CatM(0)
         IndxM=CatM(M)

         SumBZC=0.d0
         Do a=1,AnX(0)
            Indxa=AnX(a)
            SumBZC=SumBZC + Aqm(Indxa)*(2.d0*BMX(IndxM,Indxa,0) + Z*CMX(IndxM,Indxa))
         EndDo


         SumPhi=0.d0
         Do c=1,CatM(0)
            Indxc=CatM(c)

            If (Indxc.NE.IndxM) Then
               SumPsi=0.d0
               Do a=1,AnX(0)
                  Indxa=AnX(a)
                  SumPsi=SumPsi + Aqm(Indxa)*Psi(IndxM,Indxc,Indxa)
               EndDo

               SumPhi=SumPhi + Aqm(Indxc)*(2.d0*Phi(IndxM,Indxc,0) + SumPsi)

            EndIf
         EndDo

         SumPsi=0.d0
         Do a=1,AnX(0)-1
            Indxa =AnX(a)
            Do ap=a+1,AnX(0)
               Indxap=AnX(ap)
               SumPsi=SumPsi + Aqm(Indxa)*Aqm(Indxap)*Psi(IndxM,Indxa,Indxap)
            EndDo
         EndDo

         SumL=0.d0
         Do n=1,NeuN(0)
            Indxn=NeuN(n)
            SumL=SumL + Aqm(Indxn)*Lambda(Indxn,IndxM)
         EndDo
         SumL=2.d0*SumL

         SumLL=0.d0
         Do n=1,NeuN(0)
            Indxn=NeuN(n)
            Do a=1,AnX(0)
               Indxa =AnX(a)
               SumLL=SumLL + Aqm(Indxn)*Aqm(Indxa)*Ksi(Indxn,IndxM,Indxa)
            EndDo
         EndDo
         SumLL=SumLL*6.d0


!         Write(*,'(I3,6E12.4)')IndxM,  AqCharge(IndxM)**2*F , SumBZC , SumPhi , SumPsi , AqCharge(IndxM)*SumCca , SumL

         Aqg(IndxM)=AqCharge(IndxM)**2*F + SumBZC + SumPhi + SumPsi + AqCharge(IndxM)*SumCca + SumL+SumLL

         If (OnlyDebye) Aqg(IndxM)=AqCharge(IndxM)**2*F


         Aqg(IndxM)=dexp(Aqg(IndxM))
 !        write(*,*)IndxM,Aqm(IndxM),Aqg(IndxM)
      EndDo

   End

      SubRoutine AnIons
!****************************************************************************
!  SUBROUTINE  :  AnIons
!  DESCRIPTION :  Calculates activity coefficienst for Anions
!  VARIABLES   :
!
!  CREATED     :  2013.05.13  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Use Aqua
      Implicit None

! Functions
      Double Precision BMX, CMX, Phi

   Logical OnlyDebye
   Common /COnlyDebye / OnlyDebye

! Local variables
      Integer X,IndxX,c,Indxc,a,Indxa,cp,Indxcp,n,Indxn
      Double PRecision SumBZC,SumPhi,SumPsi,SumL,SumLL
!-----Start-------------------------------------------------------------------


      Do X=1,AnX(0)
         IndxX=AnX(X)

         SumBZC=0.d0
         Do c=1,CatM(0)
            Indxc=CatM(c)
            SumBZC=SumBZC +       Aqm(Indxc)*(2.d0*BMX(Indxc,IndxX,0) + Z*CMX(Indxc,IndxX))
         EndDo

         SumPhi=0.d0
         Do a=1,AnX(0)
            Indxa=AnX(a)

            If (Indxa.NE.IndxX) Then
               SumPsi=0.d0
               Do c=1,CatM(0)
                  Indxc=CatM(c)
                  SumPsi=SumPsi + Aqm(Indxc)*Psi(Indxc,IndxX,Indxa)
               EndDo
               SumPhi=SumPhi + Aqm(Indxa)*(2.d0*Phi(IndxX,Indxa,0) + SumPsi)
            EndIf
         EndDo


         SumPsi=0.d0
         Do c=1,CatM(0)-1
            Indxc =CatM(c)
            Do cp=c+1,CatM(0)
               Indxcp=CatM(cp)
               SumPsi=SumPsi + Aqm(Indxc)*Aqm(Indxcp)*Psi(Indxc,Indxcp,IndxX)
            EndDo
         EndDo

         SumL=0.d0
         Do n=1,NeuN(0)
            Indxn=NeuN(n)
            SumL=SumL + Aqm(Indxn)*Lambda(Indxn,IndxX)
         EndDo
         SumL=2.d0*SumL

         SumLL=0.d0
         Do n=1,NeuN(0)
            Indxn=NeuN(n)
            Do c=1,CatM(0)
               Indxc =CatM(c)
               SumLL=SumLL + Aqm(Indxn)*Aqm(Indxc)*Ksi(Indxn,Indxc,IndxX)
            EndDo
         EndDo
         SumLL=SumLL*6.d0


!         If (IndxX.EQ.8) Then
!            write(*,*)'SumBZC:',SumBZC
!            write(*,*)'SumPhi:',SumPhi
!            write(*,*)'SumPsi:',SumPsi
!            write(*,*)'SumCca:',AqCharge(IndxX)*SumCca
!         EndIf

         Aqg(IndxX)=AqCharge(IndxX)**2*F + SumBZC + SumPhi + SumPsi - AqCharge(IndxX)*SumCca + SumL + SumLL

         If (OnlyDebye) Aqg(IndxX)=AqCharge(IndxX)**2*F

         Aqg(IndxX)=dexp(Aqg(IndxX))

!         write(*,'(I3,4F15.5)')IndxX,AqCharge(IndxX)**2*F , SumBZC , SumPhi , AqCharge(IndxX)*SumCca

      EndDo
   End



         SubRoutine Netural
!****************************************************************************
!  SUBROUTINE  :  Neutral
!  DESCRIPTION :  Calculates activity coefficienst for neutral components
!  VARIABLES   :
!
!  CREATED     :  2013.05.13  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Use Aqua
      Implicit None

! Local variables
      Integer N,IndxN,c,Indxc,a,Indxa
      Double PRecision SumNc,SumNa,SumNca
!-----Start-------------------------------------------------------------------

      Do N=1,NeuN(0)
         IndxN=NeuN(N)                 ! Find the index for the neutral component

         SumNc=0.d0 ; SUmNa=0.d0; SumNca=0.d0       ! Set sums to zero before start


         Do c=1,CatM(0)
            Indxc=CatM(c)
            SumNc=SumNc + Aqm(Indxc)*Lambda(IndxN,Indxc)
         EndDo

         Do a=1,AnX(0)
            Indxa=AnX(a)
            SumNa=SumNa + Aqm(Indxa)*Lambda(IndxN,Indxa)
         EndDo

         Do c=1,CatM(0)
            Indxc=CatM(c)
            Do a=1,AnX(0)
               Indxa=AnX(a)
               SumNca=SumNca + Aqm(Indxc)*Aqm(Indxa)*Ksi(IndxN,Indxc,Indxa)
            EndDo
         EndDo

         Aqg(IndxN)=2.d0*(SumNc + SumNa) + SumNca
         Aqg(IndxN)=dexp(Aqg(IndxN))

      EndDo
      End



      SubRoutine Osmo
!****************************************************************************
!  SUBROUTINE  :  Osmo
!  DESCRIPTION :  Calculates osmotic coefficient
!  VARIABLES   :
!
!  CREATED     :  2013.05.13  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Functions
      Double Precision BMX, CMX, Phi

      Logical OnlyDebye
      Common /COnlyDebye / OnlyDebye

! Local variables
      Integer c,Indxc,cp,Indxcp,a,Indxa,ap,Indxap,n,Indxn
      Double PRecision SumBZC,SumPhicc,SumPhiaa,SumPsi,SumL,SumM
!-----Start-------------------------------------------------------------------

      SumM=0.d0
      Do c=1,CatM(0)
         SumM=SumM + Aqm(CatM(c))
      EndDo
      Do a=1,AnX(0)
         SumM=SumM + Aqm(Anx(a))
      EndDo
      Do n=1,NeuN(0)
         SumM=SumM + Aqm(NeuN(n))
      EndDo

      SumBZC=0.d0
      Do c=1,CatM(0)
         Indxc=CatM(c)
         Do a=1,AnX(0)
            Indxa=AnX(a)

            SumBZC=SumBZC + Aqm(Indxc)*Aqm(Indxa)*(BMX(Indxc,Indxa,2) + Z*CMX(Indxc,Indxa))
         EndDo
      EndDo

      SumPhicc=0.d0
      Do c=1,CatM(0)-1
         Indxc=CatM(c)
         Do cp=c+1,CatM(0)
            Indxcp=CatM(cp)
            SumPsi=0.d0
            Do a=1,AnX(0)
               Indxa=AnX(a)
               SumPsi=SumPsi + Aqm(Indxa)*Psi(Indxc,Indxcp,Indxa)
            EndDo

            SumPhicc=SumPhicc + Aqm(Indxc)*Aqm(Indxcp)*(Phi(Indxc,Indxcp,2) + SumPsi)

         EndDo
      EndDo

      SumPhiaa=0.d0
      Do a=1,AnX(0)-1
         Indxa=AnX(a)
         Do ap=a+1,AnX(0)
            Indxap=Anx(ap)
            SumPsi=0.d0
            Do c=1,CatM(0)
               Indxc=CatM(c)
               SumPsi=SumPsi + Aqm(Indxc)*Psi(Indxc,Indxa,Indxap)
            EndDo

            SumPhiaa=SumPhiaa + Aqm(Indxa)*Aqm(Indxap)*(Phi(Indxa,Indxap,2) + SumPsi)

         EndDo
      EndDo

      SumL=0.d0
      Do n=1,NeuN(0)
         Indxn=NeuN(n)
         Do c=1,CatM(0)
            Indxc=CatM(c)
            SumL=SumL + Aqm(Indxn)*Aqm(Indxc)*Lambda(Indxn,Indxc)
         EndDo

         Do a=1,AnX(0)
            Indxa=AnX(a)
            SumL=SumL + Aqm(Indxn)*Aqm(Indxa)*Lambda(Indxn,Indxa)
         EndDo
      EndDo

      OsmoC=(-Aphi*IS**(1.5d0)/(1.d0 + Debye_b*sqrIS) + SumBZC + SumPhicc + SumPhiaa + SumL)*2.d0/SumM + 1.d0

      If (OnlyDebye) OsmoC=(-Aphi*IS**(1.5d0)/(1.d0 + Debye_b*sqrIS) )*2.d0/SumM + 1.d0


      Aqg(iAq_H2O)=dexp(-OsmoC*AqMw(iAq_H2O)*1.d-3*SumM)

    !  Aqg(iAq_H2O)=OsmoC

      End

      SubRoutine CalcISandZ
!****************************************************************************
!  SUBROUTINE  :  CalcIS
!  DESCRIPTION :  Calculates Ionic strength and Z
!  VARIABLES   :
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Local variables
      Integer i
!-----Start-------------------------------------------------------------------
      Z=0.d0
      IS=0.d0


      Do i=1,NoAq
         if (AqInPitz(i).EQ.1) Then
            Z = Z + Aqm(i)*abs(AqCharge(i))
            IS=IS + Aqm(i)*AqCharge(i)**2
         EndIF
      EndDo
      IS=IS*0.5d0
      sqrIS=dsqrt(IS)


      End