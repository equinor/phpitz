
   Program CO2Transport_exe
!****************************************************************************
!  PROGRM      :  CO2 Transport_exe
!  DESCRIPTION :  Test program for calculation of CO2 transport
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
   Use EOS
   Use HCprop
   Use GasRx
   Use Gibbs_data
   Use ResHandling
   Use ReactEquil
   Implicit None


   Double Precision :: N_aq(NoAq)
   Common /CN_Aq /N_Aq

   Integer :: CalcOptions(0:10)
   Common /CCalcOptions/ CalcOptions

   Double Precision :: GasDens
   Common /CGasDens / GasDens

   Integer :: ModelOptions(0:10)
   Common /CModelOptions / ModelOptions


   Double Precision pH, Gref(NoComp)
! Local variables
   Integer :: i,j, k, N, ProfileType , NCases, iHC, NRxInput, ResFileNum, Phase

   Character :: InputFileName*256,MainPath*256
   Double Precision :: InpN(NoComp,500),InpFug(NoComp),InpLiqN(NoComp),Paq, NGas(NoCOmp),NTOT, Tcalc(10),  GasbeforepH(NoHC), Tcinp, AllComp(NoComp), ProfileVar(20)

   Double Precision :: MultiCaseRes(NoComp,500),MultiCaseWat(10,500)

   Double Precision :: Results(50,500),PT(2,500), Inb(10), Aqmass , GasMass, mH2SO4, mHNO3, mNH3, mCO2
!-Start----------------------------------------------------------------------
   AllComp=0.d0 ; MultiCaseRes=0.d0 ; MultiCaseWat=0.d0

   ! Some general model setup
   Call InitVar


! The file PathInfo.inp contains paths to where various folders/files. You can set this in two ways.

   ! 1. If the Pathinfo.inp file is in the same folder as the executable file, you dont need to specify the path
   ! it is sufficient to use: FileName= "PathInfo.inp"

   ! 2. If the Pathinfo.inp file is located elsewhere, you can specify the path to the folder and set it like this:
   ! First set the path:   MainPath="C:\<path to folder>\"
   ! Then se ful filename: FileName=trim(MainPath) // "PathInfo.inp"


   ! Option 1: Set the filename only without the path
     ! FileName= "PathInfo.inp"  ! This will read paths to datafolder and other folders

   ! Option 2: Include full paht
   ! This is the path to your working folder - this method is used by Baard locally
   MainPath="../datafiles/"
   FileName=trim(MainPath) // "PathInfo.inp"  ! This will read paths to datafolder and other folders


   Call LoadPaths       ! Sets the corect paths


! Set up components and component types in model
   CAll CompThermoLoad
   Call ReactionsLoad

! Load other files
   Call BinLoad
   Call AqSetup         ! Set up Aq-model from setupfile
   Call PitzerSetup     ! Defines the species in the Pitzer model
   Call EqSetup         ! Set up equilibrium constants
   Call PitzerLoad      ! Read all Pitzer paramters (T and P parameters)
   Call Eq_Load         ! Read all equilibrium constants (T and P effects)
   Call EOSParLoad
  ! Call ReactionsLoad




   ! This is path to simulation input file. Either specify full path, or just filename if the file is in the same folder as the executable file
   InputFileName=trim(MainPath) // "CalcInput.dat"
!   InputFileName= "CalcInput.dat"


! Read the input fiel that defines the calcualtion
   !Call ReadInput(InputFileName,InpN,InpLiqN,Ncases,PT,ProfileVar)

   Call ReadInput(InputFileName,Ncases,PT,InpN)

! You can store various results in variable Results(N,k) where N is type of result and k is case number.
   Results=0.d0

   Do k=1,NCases

      If (NCases.GT.1) Write(*,*)'Case: ',k


      Ptot=PT(1,k)
      Temp=PT(2,k)+273.15d0

      Call SetPTVar

      AllComp =InpN(1:NoComp,k)

      ! Set moles in gas phase to calcualte fugacities. Solids are set to zero.
      Do i=1, ncomp
         If (CompType(i).EQ.1) Then
            NGas(i)    =InpN(i,k)
         EndIf
      EndDo


      ! ManualUpdate This file (below here) will require some manual updae if new components are added.


      If (CalcOptions(1).EQ.1) Then
         Phase=0
         Call SinglePhase(NGas,Phase)
         NTOT=sum(NGas)

         Results(1,k)=GasDens  ! This is gas density

         GasMass=Dot_product(NGas, Mw)*1.d-3

         Write(*,*)'CO2 phase properties'
         Write(*,*)"-------------------------------------------------------------"
         Write(*,'(A14,F10.2,A)')' Pressure   : ',Ptot,' bar'
         Write(*,'(A14,F10.2,A)')' Temperature: ',Temp-273.15d0,' C'
         Write(*,'(A14,F10.2,A)')' Density    : ',GasDens,' kg/m3'
         Write(*,'(A14,F10.2,A)')' Tot mass   : ',GasMass,' kg'
         Write(*,'(A14,F10.2,A)')' Volume     : ',GasMass/GasDens,' m3'
         Write(*,*)"-------------------------------------------------------------"
         Write(*,*)'Component       Moles         ppm          FugCoef'
         Write(*,*)"-------------------------------------------------------------"
         Do N=1,NoHC
            Write(*,'(A1,A10,ES15.3,F12.2,F15.6)')" ",comp_name(N),NGas(N), NGas(N)*1.d6/NTOT, Fugy(N)
         EndDo
         Write(*,*)"-------------------------------------------------------------"
      EndIf



      If (CalcOptions(2).EQ.1) Then

         If (CalcOptions(1).EQ.0) Call SinglePhase(NGas,1)


         If (ModelOptions(1).EQ.1) Then
            Call CalcGandK

            Write(*,*)
            Write(*,*)'Reactions and Thermodynamic Data'
            Write(*,*)"-------------------------------------------------------------"
            Write(*,*)'RxNo   dH (kJ/mole)   dS(J/moleK)   Gibbs(kJ/mole)   lnK'
            Write(*,*)"-------------------------------------------------------------"
            Do i=1,nRx
               Write(*,'(I3,4F14.2)')i,RxH(i),RxS(i),RxG(i)*1.0E-3_dp, RxlnK(i)
            EndDo
            Write(*,*)"-------------------------------------------------------------"

         ElseIf (ModelOptions(1).EQ.2) Then
            Call CalcThermo(Gref)
            Call CalcCompG(AllComp,Gref,Fugy,Temp)

            write(*,*)
            Write(*,*)'Gas thermodynamics'
            Write(*,*)"-------------------------------------------------------------"
            Write(*,*)'Component      Gref            Gibbs            Ni x Gibbs'
            Write(*,*)"-------------------------------------------------------------"
            Do i=1,ncomp
               Write(*,'(A1, A10,2F15.2,F20.0)')' ',Comp_name(i), Gref(i),CompGibbs(i), AllComp(i)*CompGibbs(i)
            EndDo
            Write(*,*)"-------------------------------------------------------------"
         EndIf

      EndIf




      If (CalcOptions(3).EQ.1) Then
         If (CalcOptions(1).EQ.0) Then
            Phase=0
            Call SinglePhase(NGas,Phase)  ! Calculate gas fugacities if it hasnt been calculated before
         EndIf

         If (ModelOptions(1).EQ.1) Then
            IF (CalcOptions(2).EQ.0) Call CalcGandK         ! Calculate gas phase reaction thermo if not calcualted before
            Call GasRxSolver(AllComp)

            MultiCaseRes(1:NoComp,k)=NGasRx(1:NoComp)

         ElseIf (ModelOptions(1).EQ.2) Then
            IF (CalcOptions(2).EQ.0) Call CalcThermo(Gref)  ! Calculate gibbs for components if not calculated before
            Call CO2PhaseReactionsNew(AllComp)

            MultiCaseRes(1:NoComp,k)=AllComp(1:NoComp)

         EndIf

      EndIf

      NGas=AllComp(1:NoHC)


      If (ProfileTYpe.EQ.2) Then
         Do i=1,13
            Results(19+i,k)=AllComp(i)/Sum(AllComp)*1.d6
         EndDo
      EndIf

      GasbeforepH=NGas

      If (CalcOptions(4).EQ.1) Then

         HCN=0.d0
         Do i=1,NoHC
            HCN(0,i)=NGas(i) + InpLiqN(i)
         EndDo
         HCN(0,0)=sum(HCN(0,1:NoHC))


         IF (HCN(PhTot,iHC_H2O).GT.0.d0 .AND. (HCN(PhTot,iHC_H2SO4) + HCN(PhTot,iHC_HNO3) + HCN(PhTot,iHC_NH3)).GT.0.d0) Then


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

            ! Add NH3 to water, but not more moles than water
               iHC=iHC_NH3
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
               Phase=1
               Call SinglePhase(NGas,Phase)


               N_aq=0.d0
               N_aq(iAq_H2O)=HCN(PhAq,iHC_H2O)
               N_aq(iAq_H2SO4)=HCN(PhAq,iHC_H2SO4)
               N_Aq(iAq_HNO3)=HCN(PhAq,iHC_HNO3)
               N_Aq(iAq_NH3)=HCN(PhAq,iHC_NH3)


               Call PTcalc(pH)

         EndIf

         IF (HCN(PhAq,0).GT.0.d0) Then


            Write(*,*)
            Write(*,*)'AQUEOUS PHASE FORMATION'
            Write(*,*)'------------------------------------------------------------------------'
            Write(*,*)'Component     Gas before(ppm)     Gas after (ppm)      FugCoef'
            Write(*,*)'------------------------------------------------------------------------'
            Do i=1, NoHC
               Write(*,'(2A,3F17.4)')' ',comp_name(i),GasbeforepH(i)/sum(GasbeforepH)*1.d6, HCN(PhGas,i)/HCN(PhGas,0)*1.d6,FugY(i)
            EndDO
           Write(*,*)'------------------------------------------------------------------------'

           Results(2,k)=HCN(PhGas,iHC_H2O)/HCN(PhGas,0)*1.d6
           Results(3,k)=HCN(PhGas,iHC_H2SO4)/HCN(PhGas,0)*1.d9
           Results(4,k)=HCN(PhGas,iHC_HNO3)/HCN(PhGas,0)*1.d6
           Results(5,k)=HCN(PhGas,iHC_NH3)/HCN(PhGas,0)*1.d6

            MultiCaseRes(1:NoComp,k)=HCN(PhGas,1:NoComp)

            If (pH.LT.5.0_dp) Then
               Write(*,*)' MOLE RATES AFTER ACID PHASE FORMATION '
            Else
               Write(*,*)' MOLE RATES AFTER AQUEOUS PHASE FORMATION '
            EndIf


            Write(*,*)'------------------------------------------------------------------------'
            Write(*,*)' Component       Total          Gas            Water         % in aq '
            Write(*,*)'------------------------------------------------------------------------'
            Do i=1, NoHC
               IF (HCN(PhTot,i).GT.0.d0) Then
                  Paq=HCN(PhAq,i)/HCN(PhTot,i)*1.d2
               Else
                  Paq=0.d0
               EndIf

               Write(*,'(2A,3ES15.3,F12.4)')' ',comp_name(i),HCN(PhTot,i),HCN(PhGas,i),HCN(PhAq,i),pAq
            EndDO
            Write(*,*)'------------------------------------------------------------------------'
            Write(*,*)

            If (pH.LT.5.0_dp) Then
               Write(*,*)' ACID PHASE RESULTS '
            Else
               Write(*,*)' AQUEOUS PHASE RESULTS '
            EndIf

            Write(*,*)'------------------------------------------------------------------------'



               Aqmass=0.d0 ; mH2SO4=0.d0 ; mHNO3=0.d0 ; mNH3=0.d0
               Do i=1,ncomp
                  Aqmass=Aqmass + HCN(PhAq,i)*Mw(i)/1.d3
               EndDo

               mH2SO4=HCN(PhAq,iHC_H2SO4)*Mw(iHC_H2SO4)/1.d3 / Aqmass * 1.d2
               mHNO3 =HCN(PhAq,iHC_HNO3)*Mw(iHC_HNO3)/1.d3 / Aqmass * 1.d2
               mNH3  =HCN(PhAq,iHC_NH3)*Mw(iHC_NH3)/1.d3 / Aqmass * 1.d2
               mCO2  =HCN(PhAq,iHC_CO2)*Mw(iHC_CO2)/1.d3 / Aqmass * 1.d2

               Write(*,'(A,F8.3,A)')' Water pH      : ', pH
               Write(*,'(A,F8.4,A)')' Water phase   : ',Aqmass, ' kg'
               Write(*,'(A,ES9.3)')  ' Water activity: ',Aqg(1)
               Write(*,'(A,F8.4,A)')' H2SO4 in water:',mH2SO4, ' wt%'
               Write(*,'(A,F8.4,A)')' HNO3  in water:',mHNO3, ' wt%'
               Write(*,'(A,F8.4,A)')' NH3   in water:',mNH3, ' wt%'
               Write(*,'(A,F8.4,A)')' CO2   in water:',mCO2, ' wt%'

               MultiCaseWat(1,k)=1.d0
               MultiCaseWat(2,k)=pH
               MultiCaseWat(3,k)=Aqmass
               MultiCaseWat(4,k)=mH2SO4
               MultiCaseWat(5,k)=mHNO3
               MultiCaseWat(6,k)=mNH3
               MultiCaseWat(7,k)=mCO2


              Write(*,*)'------------------------------------------------------------------------'
              Write(*,*)'Ion            Conc(mole/kg)          Gamma'
              Write(*,*)'------------------------------------------------------------------------'
              Do i=2,NoAq
                  Write(*,'(2A,F17.7,F17.5)')' ',AqName(i),Aqm(i),Aqg(i)
              EndDo
              Write(*,*)'------------------------------------------------------------------------'

         Else
            Write(*,*)



            Write(*,*)'AQUEOUS PHASE FORMATION'


            Write(*,*)'------------------------------------------------------------------------'
            Write(*,*)'            No aquous phase present'
            Write(*,*)'------------------------------------------------------------------------'

         EndIf
      EndIf

   EndDo


500 Continue


    If (CalcOptions(5)==1) Then

      FileName=trim(MainPath) // "MultiCase.res"

      ResFileNum=62

      Open (UNIT=ResFileNum, FILE=FileName)

      Write(ResFileNum,*)NCases

      Do k=1,NCases
         Write(ResFileNum,*)PT(1,k)
         Write(ResFileNum,*)PT(2,k)
         Do i=1,NoComp
            Write(ResFileNum,*)MultiCaseRes(i,k)
         EndDO

         Do i=1,7
            Write(ResFileNum,*) MultiCaseWat(i,k)
         EndDo


      EndDo

      CLose(ResFileNum)
    Endif

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
   Logical :: OnlyDebye
   Common /COnlyDebye / OnlyDebye

!-Start----------------------------------------------------------------------

   ! Some parameters used at different places - set to default values
   ReadFromFile=.True.
   OnlyDebye=.False.
   EOS_Type="SRK"

   End


   SubRoutine ReadInput(InputFileName,Ncases,PT,InpN)
!****************************************************************************
!  SUBROUTINE  :  InitVar
!  DESCRIPTION :  Initializes some variables and settings
!  VARIABLES   :
!
!  CREATED     :  2025.08.12  Baard Kaasa www.scaleconsult.com
!  CHANGED     :
!
!****************************************************************************
   Use GlobVar
   Use ResHandling
   Use dPDef
   Use Path
   Use GasRx

   Implicit None

! Call variables
   Character :: InputFileName*256

! common variables
   Integer :: ModelOptions(0:10)
   Common /CModelOptions / ModelOptions

   Integer :: CalcOptions(0:10)
   Common /CCalcOptions/ CalcOptions

!Local variables
   Integer :: i,j, InputFileNum, NCases
   Real(dp) :: InpN(NoComp,500), PT(2,500)
!-Start----------------------------------------------------------------------

   InputFileNum=37

   PT=0.d0 ; InpN=0.d0

   Open (UNIT=InputFileNum, FILE=InputFileName)

      CalcOptions=0
      Read(InputFileNum,*)CalcOptions(0)
      Do i=1,CalcOptions(0)
         Read(InputFileNum,*)CalcOptions(i)
      EndDo

! Reac model selections
      ModelOptions=0
      Read(InputFileNum,*)ModelOptions(0)
      Do i=1,ModelOptions(0)
         Read(InputFileNum,*)ModelOptions(i)
      EndDo


! Read in which reactions that should be included or not
      Read(InputFileNum,*)IncludeRx(0)
      Do i=1,IncludeRx(0)
         Read(InputFileNum,*)IncludeRx(i)
      EndDo

! Read various user choices
      UserChoices=0
      Read(InputFileNum,*)UserChoices(0)
      Do i=1,UserChoices(0)
         Read(InputFileNum,*)UserChoices(i)
      EndDo

      Read(InputFileNum,*)Ncases

      Do j=1,Ncases
         Read(InputFileNum,*)PT(1,j),PT(2,j)

         Do i=1, NoComp
            Read(InputFileNum,*)InpN(i,j)
         EndDo
      EndDo


   Close(InputFileNum)




   End


      SubRoutine CalcGandK
      Use GlobVar
      Use GasRx
      Implicit None

      Integer :: i,j

      RxH=0.0_dp ; RxS=0.0_dp
      Do j=1,nRx
         Do i=1,ncomp
            RxH(j)=RxH(j) + stoich(i,j)*CompH(i)
            RxS(j)=RxS(j) + stoich(i,j)*CompS(i)
         EndDo

         RxG(j)=RxH(j)*1.0E3_dp - Temp*RxS(j)
         RxlnK(j)=-RxG(j)/(Rgas*Temp)
      EndDo

   End
