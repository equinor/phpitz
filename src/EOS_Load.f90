!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  EOS_Load.f90
!  DESCRIPTION :  Loads parameters for EOS
!  CREATED     :  2016.06.27  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  PureCompLoad - Loads pure component parameters
!
!****************************************************************************


   SubRoutine PureCompLoad
!****************************************************************************
!  SUBROUTINE  :  PureCompLoad
!  DESCRIPTION :  Load pure component paramters for EOS
!  VARIABLES   :  
!
!  CREATED     :  2016.06.26  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************   
      Use Par
      Use EOS
      Use Path
      Use Files
      Use HC
      Use ErrMSGS
      Implicit None
      
! Local variables
      Integer i
!-----Start------------------------------------------------------------------- 

      CompName="" ; MW=0.d0 ; Dens=0.d0 ; Tc=0.d0 ; Pc=0.d0; Acent=0.d0 ; 

      If (ReadFromFile) Then
         FileName=Trim(Paths(DataPath)) // EOS_PureComp
      
         Call CheckFile
         If (Info.GT.0) Return
      
         If (LOG) write(LogFilePath,*)"Rading file: ", EOS_PureComp
      
         Open (UNIT=InFile, FILE=FileName)
         Call FindStart(InFile)
      
         Read(Infile,*)
         Do i=1,NoDefHC
            Read(Infile,*) CompName(i),Mw(i),Tc(i),Pc(i),Acent(i),CPen(i)
            CompName(i)=Trim(CompName(i))
         EndDo
   
         Close (InFile)
      !Else
      !   Call SET_PureComp
      EndIf
      
   End
   
   
   SubRoutine BinLoad
!****************************************************************************
!  SUBROUTINE  :  BinLoad
!  DESCRIPTION :  Reads binary interaction parameters
!  VARIABLES   :  
!
!  CREATED     :  2016.07.05  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
      Use Par
      Use EOS
      Use Path
      Use Files
      Use HC
      Use ErrMSGS
   Implicit None

! Local variables
      Integer i,j
      Double Precision binij
!-Start----------------------------------------------------------------------
      kij=0.d0
      
      If (ReadFromFile) Then
         FileName=Trim(Paths(DataPath)) // EOS_BIN_File
      
         Call CheckFile
         If (Info.GT.0) Return
      
         Open (UNIT=InFile, FILE=FileName)
         Call FindStart(InFile)
      
         Read(Infile,*)
10       Read(Infile,*)i,j,binij
      
         If (i.GT.0 .AND. j.GT.0) Then
            Kij(i,j)=binij  ; kij(j,i)=binij
            Goto 10
         EndIF
         
         Close (InFile)
   !Else
   !   Call SET_PVTbin 
   EndIf
   End