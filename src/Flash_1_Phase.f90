

      SubRoutine SinglePhase(Ni,Phase)
!****************************************************************************
!  SUBROUTINE  :  SinglePhase
!  DESCRIPTION :  Load parameters for equilibrium constants
!  VARIABLES   :
!
!  CREATED     :  2014.01.12  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use EOS
      Use HC
      Use FLash
      Use HCprop
      Use GlobVar
      Implicit None

! Call variables
      Integer Lev
      Double Precision Ni(NoHC),ZF

      Double Precision GasDens
      Common /CGasDens / GasDens

! Local variables
      Integer i,Phase,ResPhase
      Double Precision SumAq,FugSingle(NoHC),Mwgas

!-----Start--------------------------------------------------
      Lev=2

      If (Sum(Ni).LE.0.d0) Then
         Zfact=1.d0
         FugCoef=1.d0
         HCy=0.d0
         FugY=1.d0
         Return
      EndIf

      Call SolveEOS(Ni,ZF,Phase)
      Call CalcFug(Lev,Ni,ZF,FugSingle)
!      Call CalcFug_Deriv(nc(0,1:Nmax),ZF,dFug)
!      Call CPA_Temp_deriv(nc(0,1:Nmax),ZF)

      Phase=1

! Set the properties from the calculation
      Zfact=0.d0
      FugCOef=0.d0

      ZF=ZF - EOS_C/Sum(Ni)*Ptot/(Reos*Temp)
      Zfact(Phase)=ZF

      FugCoef(Phase,1:NoHC)=FugSingle(1:NoHC)


 !     HCN(1:2,1:NoHC)=0.d0
 !     HCN(phase,1:NoHC)=Ni

! Set HCy and Fugy to be used for calc of Henrys law
      HCy=Ni/Sum(Ni)

      Fugy=FugCoef(Phase,1:NoHC)

      Mwgas=sum(Ni(1:NoHC)*Mw(1:NoHC))/sum(Ni(1:NoHC))

      GasDens=Ptot*Mwgas/(ZF*Reos*Temp)*1.d3

   End


!         SubRoutine CalcPhase(Ni,Phase,Deriv)
!!****************************************************************************
!!  SUBROUTINE  :  CalcPhase
!!  DESCRIPTION :  Calculate phase properties from EOS for a given composition and phase
!!  VARIABLES   :
!!
!!  CREATED     :  2017.12.22  Baard Kaasa www.scalesim.no
!!  CHANGED     :
!!****************************************************************************
!      Use EOS
!      Use HC
!      Use FLash
!      Use HCprop
!      Use Path
!      Implicit None
!
!! Call variables
!      Integer Phase,Deriv,Lev
!      Double Precision Ni(NoHC),ZF
!
!
!! Local variables
!      Integer i
!      Double Precision SumAq,FugSingle(NoHC)
!
!!-----Start--------------------------------------------------
!      Lev=2 ; If (Deriv.EQ.1) Lev=3
!
!      Call SolveEOS(Ni,ZF,Phase)
!
!      If (Lev.EQ.2) Then
!         Call CalcFug(Lev,Ni,ZF,FugSingle)
!      Else
!         Call CalcFugDeriv(Ni,ZF,FugSingle,dFugCoef(Phase,1:NoHC,1:NoHC))
!      EndIf
!
!!      If (Deriv.EQ.1) Call CalcFug_Deriv(Nc,ZF,dFug)
!
!!      Call CPA_Temp_deriv(nc(0,1:Nmax),ZF)
!
!! Set properties for correct phase
!      Zfact(Phase)=ZF
!      FugCoef(Phase,1:NoHC)=FugSingle(1:NoHC)
!
!      End