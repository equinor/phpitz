!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  DIV.f90
!  DESCRIPTION :  Contains various help subroutines used different places in program
!  CREATED     :  2013.05.05  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  CheckFile - Checks if a file is present
!                 FindStart - Reads through info part in beginning of files to fint "start"
!                 UnitConvert - Convert units
!
!****************************************************************************


   SubRoutine CheckFile
!****************************************************************************
!  SUBROUTINE  :  CheckFile
!  DESCRIPTION :  Checks if a file is present
!  VARIABLES   :  FileName - File to be checked
!
!  CREATED     :  2013.05.04  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Path
   Use ErrMSGS
   Implicit None

! Local variables
   Logical Exists
!--Start----------------------------------------------------''---------------
   Info=0
   Inquire (FILE=FileName,Exist=Exists)   ! Checks if File exists. Variables Exists set to true or false

   If (.not.Exists) Then
      Info=1
      ErrMSG='File does not exist: '
      ErrMSG=ErrMSG(1:Len_Trim(ErrMSG)) // FileName
   EndIf

   End




   SubRoutine FindStart(InFile)
!****************************************************************************
!  SUBROUTINE  :  FindStart
!  DESCRIPTION :  Reads thtrough info-part in beginning of files
!  VARIABLES   :  FileName - File to be read
!
!  CREATED     :  2013.05.04  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************

! Call variables
   Integer InFile

! Local variables
   Character Dummy*25
   Logical FoundStart
!--Start-------------------------------------------------------------------
   Rewind(InFile)

   FoundStart=.False.

   Do while (.NOT.FoundStart)
      Read(InFile,*)Dummy
      Dummy=trim(Dummy)
      If (Dummy.EQ.'START' .OR. Dummy.EQ.'Start') Then
         FoundStart=.True.
      EndIf
   EndDo


   End



