
   Program CO2Transport_exe
!****************************************************************************
!  PROGRM      :  CO2 Transport_exe
!  DESCRIPTION :  Test program for calculation of water phase properties
!  VARIABLES   :
!
!  CREATED     :  2025.08.12  Baard Kaasa www.scaleconsult.com
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Use Path
   Use Aqua
   Use HC
   Use HCprop
   Use GasRx
   Implicit None


   Double Precision  N_aq(NoAq)
   Common /CN_Aq /N_Aq

   Integer CalcOptions(0:10)
   Common /CCalcOptions/ CalcOptions

   Double Precision GasDens
   Common /CGasDens / GasDens

   Double Precision pH
! Local variables
   Integer i, InputFileNum,j,NRxInput,N,IncludeRx(0:MaxGasRx),iHC,k,NCases

   Character InputFileName*256,MainPath*256
   Double Precision InpN(11),InpFug(11),InpLiqN(11),Paq, Tcalc(10),NGas(11),NTOT,GasbeforepH(NoHC),Tcinp

   Double Precision Results(5,20),PT(2,20)
!-Start----------------------------------------------------------------------

   ! This is the path to your working folder
   MainPath="../datafiles/"

   InputFileName=trim(MainPath) // "CalcInput.dat"
   InputFileNum=37

   Open (UNIT=InputFileNum, FILE=InputFileName)
      Read(InputFileNum,*)Ptot,Tcinp
      Do i=1,11
         Read(InputFileNum,*)HCnames(i)
         Read(InputFileNum,*)InpN(i),InpLiqN(i)
      EndDO

      IncludeRxINP=0
      Read(InputFileNum,*)NRxInput
      Do i=1,NRxInput
         Read(InputFileNum,*)IncludeRxINP(i)
      EndDo

      CalcOptions=0
      Read(InputFileNum,*)CalcOptions(0)
      Do i=1,CalcOptions(0)
         Read(InputFileNum,*)CalcOptions(i)
      EndDo

      PT=0.d0
      Read(InputFileNum,*)Ncases
      Do i=1,Ncases
         Read(InputFileNum,*)PT(1,i),PT(2,i)
      EndDo



   Close(InputFileNum)

 !  Temp=Tcinp + 273.15d0


   Call InitVar


   FileName=trim(MainPath) // "PathInfo.inp"  ! This will read paths to datafolder and other folders
   Call LoadPaths       ! Sets the corect paths

   Call PureCompLoad
   Call BinLoad

   Call AqSetup         ! Set up Aq-model from setupfile
   Call PitzerSetup     ! Defines the species in the Pitzer model

   Call EqSetup         ! Set up equilibrium constants

! Load parameters
   Call PitzerLoad      ! Read all Pitzer paramters (T and P parameters)
   Call Eq_Load         ! Read all equilibrium constants (T and P effects)

   Call GasThermoLoad
   Call GasRxLoad


! You can store various results in variable Results(N,k) where N is type of result and k is case number.
   Results=0.d0

   Do k=1,NCases

      Ptot=PT(1,k)
      Temp=PT(2,k)+273.15d0


      Call SetPTVar

      NGas=InpN

      ! ManualUpdate This file (below here) will require some manual updae if new components are added.


      If (CalcOptions(1).EQ.1) Then
         Call SinglePhase(NGas,1)
         NTOT=sum(NGas)

         Results(1,k)=GasDens  ! This is gas density

         Write(*,*)'Gas phase properties'
         Write(*,'(A14,F10.2,A)')' Pressure   : ',Ptot,' bar'
         Write(*,'(A14,F10.2,A)')' Temperature: ',Temp-273.15d0,' C'
         Write(*,'(A14,F10.2,A)')' Density    : ',GasDens,' kg/m3'
         Write(*,*)"-------------------------------------------------------------"
         Write(*,*)'Component    Moles         ppm          FugCoef'
         Write(*,*)"-------------------------------------------------------------"
         Do N=1,NoHC
            Write(*,'(A1,A10,2ES12.3,F12.4)')" ",HCNames(N),NGas(N), NGas(N)*1.d6/NTOT, Fugy(N)
         EndDo
         Write(*,*)"-------------------------------------------------------------"
      EndIf

      If (CalcOptions(2).EQ.1) Then
         Call CalcGasThermo
         Call SinglePhase(NGas,1)
         IncludeRx=0
         IncludeRx(1:NRxInput)=IncludeRxINP(1:NRxInput)

         Call CalcQ(NGas,IncludeRx)

         Write(*,*)'Gas thermodynamics'
         Write(*,*)"-------------------------------------------------------------"
         Write(*,*)'Reaction    lnQ         dH          dS          lnK '
         Write(*,*)"-------------------------------------------------------------"
         Do i=1,NRxInput
            Write(*,'(I3,F15.4, 3F13.4)')i,Q(i),dH(i),dS(i),lnK(i)
         EndDo
         Write(*,*)"-------------------------------------------------------------"
      EndIf

      If (CalcOptions(3).EQ.1) Then
         Call CO2PhaseReactions(NGas)
      EndIf

      GasbeforepH=NGas

      If (CalcOptions(4).EQ.0) Then
         Goto 500
      EndIf


      HCN=0.d0
      Do i=1,11
         HCN(0,i)=NGas(i) + InpLiqN(i)
      EndDo
      HCN(0,0)=sum(HCN(0,1:NoHC))


   ! Put water in water phase
      HCN(PhAq,iHC_H2O)  = HCN(PhTot,iHC_H2O)

   ! Add H2SO4 to water, but not more moles than water
      iHC=iHC_H2SO4
      If (HCN(PhTot,iHC).GT.HCN(PhTot,iHC_H2O)) Then
         HCN(PhAq,iHC)=HCN(PhAq,iHC_H2O)
      Else
         HCN(PhAq,iHC)= HCN(PhTot,iHC)
      EndIf

   ! Add HNO3 to water, but not more moles than water
      iHC=iHC_HNO3
      If (HCN(PhTot,iHC).GT.HCN(PhTot,iHC_H2O)) Then
         HCN(PhAq,iHC)=HCN(PhAq,iHC_H2O)
      Else
         HCN(PhAq,iHC)= HCN(PhTot,iHC)
      EndIf

   ! The rest is kept in the gas
      Do i=1,NoHC
         HCN(PhGas,i)=HCN(PhTot,i) - HCN(PhAq,i)
      EndDo
      HCN(PhGas,0)=sum(HCN(PhGas,1:NoHC)) ! Summarise totals
      HCN(PhAq,0) =sum(HCN(PhAq,1:NoHC))


      NGas(1:NoHC)=HCN(PhGas,1:NoHC)
      Call SinglePhase(NGas,1)


      N_aq=0.d0
      N_aq(iAq_H2O)=HCN(PhAq,iHC_H2O)
      N_aq(iAq_H2SO4)=HCN(PhAq,iHC_H2SO4)
      N_Aq(iAq_HNO3)=HCN(PhAq,iHC_HNO3)


      Call PTcalc(pH)

      Write(*,*)
      Write(*,*)'******************* RESULTS ********************************************'
      Write(*,*)'(ppm)     Gas before     Gas after      FugCoef'
      Write(*,*)'------------------------------------------------------------------------'
      Do i=1, NoHC
         Write(*,'(A,3F13.4)')HCnames(i),GasbeforepH(i)/sum(GasbeforepH)*1.d6, HCN(PhGas,i)/HCN(PhGas,0)*1.d6,FugY(i)
      EndDO
     Write(*,*)'------------------------------------------------------------------------'

     Results(2,k)=HCN(PhGas,iHC_H2O)/HCN(PhGas,0)*1.d6
     Results(3,k)=HCN(PhGas,iHC_H2SO4)/HCN(PhGas,0)*1.d9
     Results(4,k)=HCN(PhGas,iHC_HNO3)/HCN(PhGas,0)*1.d6

      Write(*,*)'                   Mole rates '
      Write(*,*)'------------------------------------------------------------------------'
      Write(*,*)'               Total            Gas             Water             %aq '
      Do i=1, NoHC
         IF (HCN(PhTot,i).GT.0.d0) Then
            Paq=HCN(PhAq,i)/HCN(PhTot,i)*1.d2
         Else
            Paq=0.d0
         EndIf

         Write(*,'(A,3e17.5,F12.4)')HCnames(i),HCN(PhTot,i),HCN(PhGas,i),HCN(PhAq,i),pAq
      EndDO
      Write(*,*)'------------------------------------------------------------------------'
      Write(*,*)
      Write(*,*)'            Water phase results'
      Write(*,*)'------------------------------------------------------------------------'
      IF (HCN(PhAq,0).GT.0.d0) Then
         Write(*,'(A,F8.3,A)')'Water pH  : ', pH
         Write(*,'(A,F8.4,A)')'Water mass: ',HCN(PhAq,iHC_H2O)*18.016d0/1.d3, ' kg water'

        Write(*,*)
        Write(*,*)'Ion           Conc              Gamma'
        Write(*,*)'------------------------------------------------------------------------'
        Do i=1,NoAq
            Write(*,'(A,F17.5,F17.5)')AqName(i),Aqm(i),Aqg(i)
        EndDo
        Write(*,*)'------------------------------------------------------------------------'

      Else
         Write(*,*)'            No water phase present'
         Write(*,*)'------------------------------------------------------------------------'
      EndIf

   EndDo


500 Continue

   Write(*,'(A)')'        P             T             Density        H2O(ppm)      H2SO4(ppb)        HNO3(ppm)'
   Do k=1,NCases
      Write(*,'(6F15.5)')PT(1,k),PT(2,k),Results(1,k),Results(2,k),Results(3,k),Results(4,k)
   EndDo


   End program



   SubRoutine InitVar
!****************************************************************************
!  SUBROUTINE  :  InitVar
!  DESCRIPTION :  Initializes some variables and settings
!  VARIABLES   :
!
!  CREATED     :  2025.08.12  Baard Kaasa www.scaleconsult.com
!  CHANGED     :
!
!****************************************************************************
   Use Files
   Use EOS
   Implicit None

! Local variables
   Logical OnlyDebye
   Common /COnlyDebye / OnlyDebye

!-Start----------------------------------------------------------------------

   ReadFromFile=.True.
   OnlyDebye=.False.

   EOS_Type="SRK"

   End






!   SubRoutine TempMal
!!****************************************************************************
!!  SUBROUTINE  :  Name of routne
!!  DESCRIPTION :  Short description
!!  VARIABLES   :
!!
!!  CREATED     :  2025.xx.x  Baard Kaasa www.scaleconsult.com
!!  CHANGED     :
!!
!!****************************************************************************
!
!   Implicit None
!
!! Local variables
!
!
!!-Start----------------------------------------------------------------------
!
!
!   End