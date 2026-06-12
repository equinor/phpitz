!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  GasLoad.f90
!  DESCRIPTION :  Loads parameters for gas reactions
!  CREATED     :  2025.09.30  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  GasLoad
!
!****************************************************************************


   SubRoutine GasThermoLoad
!****************************************************************************
!  SUBROUTINE  :  GasThermoLoad
!  DESCRIPTION :  Load gas thermo parameters
!  VARIABLES   :
!
!  CREATED     :  2016.06.26  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Path
   Use Files
   Use Par
   Use GasRx
   Implicit None

! Local variables
   Integer i, Info

!-----Start-------------------------------------------------------------------

   FileName=Trim(Paths(DataPath)) // GasThermoFile
   Call CheckFile
   If (Info.GT.0) Return

   Open (UNIT=InFile, FILE=FileName)
   Call FindStart(InFile)

   Read(Infile,*)

   Do i=1,NoDefHC
      Read(Infile,*)
      Read(Infile,*)CompH(i),CompS(i)
   EndDo

   Close (InFile)

   End

      SubRoutine GasRxLoad
!****************************************************************************
!  SUBROUTINE  :  GasRxLoad
!  DESCRIPTION :  Load gas gas reactions
!  VARIABLES   :
!
!  CREATED     :  2016.06.26  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Path
   Use Files
   Use Par
   Use GasRx
   Implicit None

! Local variables
   Integer RxNo, Info, Dummy, CompNo,Dummy2

!-----Start-------------------------------------------------------------------

   FileName=Trim(Paths(DataPath)) // GasRxFile
   Call CheckFile
   If (Info.GT.0) Return

   Open (UNIT=InFile, FILE=FileName)
   Call FindStart(InFile)

   Read(Infile,*)RxInfo(0)

   Do RxNo=1,RxInfo(0)
      Read(InFile,*)Dummy,RxInfo(RxNo) ! Read reactionnumber (not used therefore dummy) and number of components
      Do CompNo=1,RxInfo(RxNo)
         Read(InFile,*)Dummy, RxComp(RxNo,CompNo),RxCoef(RxNo,CompNo)
      EndDo

   EndDo

   Read(Infile,*)

   Do CompNo=1,NoHC
      Do RxNo=1,RxInfo(0)
         Read(InFile,*)Dummy,Dummy2,CompRx(CompNo,RxNo)
      EndDo
   EndDo

   Close (InFile)

   End
