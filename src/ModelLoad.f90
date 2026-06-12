!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  ModelLoad.f90
!  DESCRIPTION :  Loads parameters for gas phase calculations
!  CREATED     :  2026.03.19  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  CompThermoLoad
!
!****************************************************************************


   SubRoutine CompThermoLoad
!****************************************************************************
!  SUBROUTINE  :  CompThermoLoad
!  DESCRIPTION :  Load thermo parameters for gases and solids
!  VARIABLES   :
!
!  CREATED     :  2016.06.26  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Path
   Use Files
   Use Par
   Use GasRx
   Use HC
   Use gibbs_data
   Use ErrMSGS
   Implicit None

! Local variables
   Integer :: i, j, Dummy, cnt

!-----Start-------------------------------------------------------------------

   NGasComp=0 ; NSolComp=0

   FileName=Trim(Paths(DataPath)) // CompThermoFile

   Call CheckFile
   If (Info.GT.0) Return


   Open (UNIT=InFile, FILE=FileName)
   Call FindStart(InFile)

   Read(Infile,*) Dummy,Dummy

   Do i=1,NoComp
      Read(Infile,*) Dummy, comp_name(i)
      Read(Infile,*)CompH(i), CompS(i), CompType(i), Dummy

      If (Dummy.EQ.1) Then
         is_reactive(i)=.True.
      Else
         is_reactive(i)=.False.
      EndIf
   EndDo

   Do i=1,NoComp
      Do j=1, Nelem
         Read(Infile,*)atom_matrix(j,i)
      EndDo
   EndDo

   Close (InFile)


! Calculate reaction index - which components that will be considered for Gibb minimum calcualtion
    cnt = 0
    do i = 1, ncomp
      if (is_reactive(i)) then
        cnt = cnt + 1
        react_idx(cnt) = i
      end if
    end do
    nreact = cnt

   End



   SubRoutine ReactionsLoad
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
   Use ErrMSGS
   Implicit None


! Local variables
   Integer RxNo, i, Dummy, nspec, specie, coeff

!-----Start-------------------------------------------------------------------

   FileName=Trim(Paths(DataPath)) // ReactionsFile
   Call CheckFile
   If (Info.GT.0) Return

   Open (UNIT=InFile, FILE=FileName)
   Call FindStart(InFile)

   Read(Infile,*)nRx

   stoich=0.0_dp

   Do RxNo=1,nRx
      Read(InFile,*)Dummy, nspec
      Do i=1,nspec
         Read(InFile,*)Dummy, specie, coeff
         stoich(specie,RxNo)=coeff
      EndDo
   EndDo



   Close (InFile)

   End