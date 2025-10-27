!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  ModelSetup.f90
!  DESCRIPTION :  Set up different tables and parameters for Scalesim model
!  CREATED     :  2013.05.04  Baard Kaasa  www.scalesim.no
!                 2015.04.09  BKaa  InitVar added
!
!  ROUTINES    :  LoadPaths   -  Reads paths from ScaleSimPath.inp
!                 AqSetup     -  Read and set up Aqueous model parameters
!                 EqSetup     -  Read equilibrium relations from file
!                 SaltSetup   -  Read saltdata from file and sets up salt model
!                 InitVar     -  Initializes variables
!                 LogFileOpen -  Creates a logfile
!
!****************************************************************************


   SubRoutine LoadPaths
!****************************************************************************
!  SUBROUTINE  :  LoadPaths
!  DESCRIPTION :  Set necessary paths
!  VARIABLES   :  
!
!  CREATED     :  2013.05.04  Baard Kaasa www.scalesim.no
!  CHANGED     :  2013.08.30  BKaa: Temppath added
!                 2015.04.05  BKaa: Read Paths from file
!****************************************************************************    
   Use Path
   Use ErrMSGS
   Implicit None
      
! Local variable
   Integer i,NoInputPaths
!--Start---------------------------------------------------------------------
     
   Call CheckFile
   If (Info.GT.0) Return
      
   Open (UNIT=InFile, FILE=FileName)
   Call FindStart(InFile)
      
   Read (InFile,*)NoInputPaths
      
   Do i=1,NoInputPaths
      Read(Infile,'(A)')Paths(i)
   EndDo
   
   Close (InFile)
      
   End

   
   
   SubRoutine AqSetup
!****************************************************************************
!  SUBROUTINE  :  AqSetup
!  DESCRIPTION :  Set up aqueous model
!  VARIABLES   :  
!
!  CREATED     :  2013.05.04  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************      
      Use Path
      Use Files
      Use Aqua
      Use ErrMSGS
      Use Pitzer
      Implicit None
     
! Local variables

      Character Dummy*25
      Integer i,NAqSpecies,SpecNo
!-----Start-------------------------------------------------------------------   

      If (ReadFromFile) Then
         FileName=Trim(Paths(DataPath)) // AqSetupFile
         If (LOG) write(LogFilePath,*)"Rading file: ", AqSetupFile
      
         Call CheckFile
         
         If (Info.GT.0) Return
      
         Open (UNIT=InFile, FILE=FileName)
         Call FindStart(InFile)

         Read(Infile,*)NAqSpecies
         Read(Infile,*)                   ! Skip heading
         Do i=1,NAqSpecies
            Read(InFile,*) SpecNo, AqName(SpecNo),AqCharge(SpecNo),AqMw(SpecNo),AqFeed(SpecNo),AqInPitz(SpecNo)
         EndDo
      
         Close(InFile)
      Else
!         Call SET_AqSetUp
      EndIf
   End
   

   SubRoutine EqSetup
!****************************************************************************
!  SUBROUTINE  :  EqSetup
!  DESCRIPTION :  Set up equilibria
!  VARIABLES   :  
!
!  CREATED     :  2014.01.12  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************     
   Use Eq   
   Use Files
   Use Path
   Use ErrMSGS
   Use SolventPar
   Implicit None
      
! Local variables
   Integer EqNo,j,i,Npar,HCindex
!-----Start------------------------------------------------------------------- 
   !If (ReadFromFile) Then
      FileName=Trim(Paths(DataPath)) // EqSetupFile
      If (LOG) write(LogFilePath,*)"Rading file: ", EqSetupFile
      Eq_spec=0
      Eq_exp=0.d0
      

      Call CheckFile
      If (Info.GT.0) Return
      
      Open (UNIT=InFile, FILE=FileName)
      Call FindStart(InFile)
      
      Read(InFile,*) 
      
      Do EqNo=1,NoEq
         Read(InFile,*)Eq_Name(EqNo)
         Read(InFIle,*)Eq_TypeRx(EqNo)
         Read(InFile,*)Eq_spec(EqNo,0)
         REad(InFile,*)SaltINEq(EqNo)
         Do j=1,Eq_spec(EqNo,0)
            Read (InFile,*)Eq_spec(EqNo,j),Eq_exp(EqNo,j)
         EndDo
      EndDo
      
      Close(InFile)
   !Else
   !   Call SET_EqSetUp
   !EndIf
   
! Read Water data
   If (ReadFromFile) Then
      FileName=Trim(Paths(DataPath)) // WaterDataFile
      Call CheckFile
      If (Info.GT.0) Return
      
      Open (UNIT=InFile, FILE=FileName)
      Call FindStart(InFile)
   
      DensPar=0.d0
      Read(Infile,*)
      Do i=1,11
         Read(Infile,*)DensPar(i)
      EndDo
   
      VapPar=0.d0
      Read(Infile,*)
      Do i=1,6
         Read(Infile,*)VapPar(i)
      EndDo
   
      FiPar=0.d0
      Read(InFile,*)
      Do i=1,5
         Read(InFile,*)FiPar(i)
      EndDo
      Close (Infile)
   Else
  !    Call SET_WaterDATA
   EndIf
   
   If (ReadFromFile) Then
      FileName=Trim(Paths(DataPath)) // GasSolubFile
      If (LOG) write(LogFilePath,*)"Rading file: ", GasSolubFile
      
      Call CheckFile
      If (Info.GT.0) Return
      
      Open (UNIT=InFile, FILE=FileName)
      Call FindStart(InFile)
      Read(Infile,*)Npar
      Read(Infile,*)
   
      GasSolub=0
      Do i=1,Npar
         Read(Infile,*)j,EqNo,HCindex
         GasSolub(i,1)=EqNo
         GasSolub(i,2)=HCindex
      EndDo
      Close (Infile)
   !Else
   !   Call SET_GasSolub
   EndIf
   
   End
