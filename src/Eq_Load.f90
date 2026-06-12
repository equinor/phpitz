!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Therm_Load.f90
!  DESCRIPTION :  Loads thermodynamic data
!  CREATED     :  2014.01.12  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  Eq_Load   Load parameters for equilibrium constants
!
!****************************************************************************


      SubRoutine Eq_Load
!****************************************************************************
!  SUBROUTINE  :  Eq_Load
!  DESCRIPTION :  Load parameters for equilibrium constants
!  VARIABLES   :
!
!  CREATED     :  2014.01.12  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Path
      Use Files
      Use Eq
      Use ErrMSGS
      Use Aqua

! Local variables
      Integer N,k,j,z
!-----Start-------------------------------------------------------------------

      Eq_Type   =0
      Eq_TypeV  =0
      Eq_T_Par  =0.d0
      Eq_T_Range=0.d0
      Eq_V_Par  =0.d0
      Eq_K_Par  =0.d0
      Eq_IS_Par =0.d0

      !If (ReadFromFile) Then
         FileName=Trim(Paths(DataPath)) // EqThermFile

         Call CheckFile
         If (Info.GT.0) Return

         Open (UNIT=InFile, FILE=FileName)
         Call FindStart(InFile)

         Read(InFile,*)

         Do N=1,NoEq
            Read(InFile,*) Eq_Name(N)
            REad(InFile,*) Eq_Type(N,0)
            Do j=1,Eq_Type(N,0)
               Read(InFile,*) Eq_Type(N,j)
               Read(InFile,*) Eq_T_Range(N,j,1),Eq_T_Range(N,j,2)
               Read(InFile,*)(Eq_T_Par(N,k,j),k=1,8)
            EndDo
         EndDo

         Do N=1,NoEq
            Read(InFile,*)
            Read(InFile,*) Eq_TypeV(N),Eq_TypeKappa(N)
            Read(InFile,*)(Eq_V_Par(N,k),k=1,6)
            Read(InFile,*)(Eq_K_Par(N,k),k=1,6)
            Read(Infile,*)(Eq_IS_Par(N,k),k=1,2)
         EndDo

         Close (InFile)
      !Else
      !   Call SET_Eq_Load
      !EndIf


!''Read parameters to calculate volumes

   If (ReadFromFile) Then
      FileName=Trim(Paths(DataPath)) // VolFileName

      Call CheckFile
      If (Info.GT.0) Return

      Open (UNIT=InFile, FILE=FileName)
      Call FindStart(InFile)

      Read(Infile,*)N

      Do i=1,N
         Read (InFile,*)(Eq_Vol_Par(i,k),k=1,10),z
      EndDo

      Close (Infile)
   Else
 !     Call SET_AqSpecVol
   EndIf

      End