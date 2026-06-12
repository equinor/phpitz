!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Pitzer_Load.f90
!  DESCRIPTION :  Loads all datafiles
!  CREATED     :  2013.05.05  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  PitzerSetup - Sets up Pitzer model
!
!****************************************************************************


      SubRoutine PitzerLoad
!****************************************************************************
!  SUBROUTINE  :  PitzerLoad
!  DESCRIPTION :  Loads all Pitzer parameters
!  VARIABLES   :
!
!  CREATED     :  2013.05.05  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Use Path
      Use ErrMSGS
      Use Files
      Implicit None

! Local variabelse
      Integer i,k,Indx1,Indx2,Indx3,FileNo,Nspec,Eqn
      Double Precision Tpar(Npar)

      Integer, Parameter :: NDatafiles=9
      Character DataFiles(NDatafiles)*20
      Data DataFiles/'PitzBeta0.dat','PitzBeta1.dat','PitzBeta2.dat','PitzC.dat', &
      'PitzPsi.dat','PitzTeta.dat','PitzLambda.dat','PitzKsi.dat','PitzAlpha.dat'/

!-----Start-------------------------------------------------------------------

! Start by setting all parameters to zero. Only read parameters will have value.
      PitzAlpha=0.d0  ! Double Precision parameters
      BetaCC   =0.d0
      CphiCC   =0.d0
      PsiCC    =0.d0
      ThetCC   =0.d0
      LambCC   =0.d0
      KsiCC    =0.d0

      BetaCCEq   =0     ! Integer parameters
      CphiCCEq   =0
      PsiCCEq    =0
      ThetCCEq   =0
      LambCCEq   =0
      KsiCCEq    =0



      PitzAlpha(1,1:NoAq,1:NoAq)=2.0d0 ! Default values - only those who are different are read in
      PitzAlpha(2,1:NoAq,1:NoAq)=0.0d0

      Do i=1,NoAq
         Do k=1,NoAq
            If ((AqCharge(i).EQ.2 .AND. AqCharge(k).EQ.-2) .OR. (AqCharge(i).EQ.-2 .AND. AqCharge(k).EQ.2)) Then
               PitzAlpha(1,i,k) =  1.4D0
               PitzAlpha(2,i,k) = 12.0D0
            EndIf
         EndDo
      EndDo


      If (ReadFromFile) Then
         Do FileNo=1,NDatafiles
            FileName=Paths(DataPath)
            FileName=FileName(1:Len_Trim(FileName)) // DataFiles(FileNo)
            If (LOG) write(LogFilePath,*)"Rading file: ", DataFiles(FileNo)

            Call CheckFile          ! FileName is automatically transferred as it is stored in MODULE Path
            If (Info.GT.0) Return

            Open (UNIT=InFile,File=FileName) ! Open the file
            Call FindStart(InFile)           ! Find the starting point of the data

            Read(InFile,*)Nspec

            If (FileNo.LE.3) Then   ! Read the Beta0, 1 and 2 files
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,(Tpar(k),k=1,Npar),Eqn
                  BetaCC(FileNo-1,Indx1,Indx2,1:Npar)=Tpar(1:Npar)
                  BetaCCEq(FileNo-1,Indx1,Indx2)=Eqn
               EndDo


            ElseIf (FileNo.EQ.4) Then     ! Read C parameters
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,(Tpar(k),k=1,Npar),Eqn
                  CphiCC(Indx1,Indx2,1:Npar)=Tpar(1:Npar)
                  CphiCCEq(Indx1,Indx2)=Eqn
               EndDo



            ElseIf (FileNo.EQ.5) Then     ! Read 3-ion interactions cc'a and caa'
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,Indx3,(Tpar(k),k=1,Npar),Eqn
                  PsiCC(Indx1,Indx2,Indx3,1:Npar)=Tpar(1:Npar)
                  PsiCCEq(Indx1,Indx2,Indx3)=Eqn
               EndDo

            ElseIf (FileNo.EQ.6) Then     ! Read 2-ion interaction of same sign cc' and aa'
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,(Tpar(k),k=1,Npar),Eqn
                  ThetCC(Indx1,Indx2,1:Npar)=Tpar(1:Npar)
                  ThetCCEq(Indx1,Indx2)=Eqn
               EndDo

            ElseIf (FileNo.EQ.7) Then     ! Read 2-specie interaction for neutral species n-c and n-a
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,(Tpar(k),k=1,Npar),Eqn
                  LambCC(Indx1,Indx2,1:Npar)=Tpar(1:Npar)
                  LambCCEq(Indx1,Indx2)=Eqn
               EndDo

            ElseIf (FileNo.EQ.8) Then     ! Read 3-specie interaction for netural species n-c-a
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,Indx3,(Tpar(k),k=1,Npar),Eqn
                  KsiCC(Indx1,Indx2,Indx3,1:Npar)=Tpar(1:Npar)
                  KsiCCEq(Indx1,Indx2,Indx3)=Eqn
               EndDo

            ElseIf (FileNo.EQ.9) Then     ! REad non-standard alpha parameters
               Do i=1,Nspec
                  Read(InFile,*)Indx1,Indx2,(Tpar(k),k=1,2) ! Its not actually temperature parameters but use Tpar instead of
                  PitzAlpha(1:2,Indx1,Indx2)=Tpar(1:2)      ! - declaring a new variables
               EndDo
            Endif
            Close(InFile)
         EndDo
      Else
 !        Call SET_PitzPar
      EndIf
      End