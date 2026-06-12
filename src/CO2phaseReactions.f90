

   Subroutine CO2PhaseReactions(NGas)
!****************************************************************************
!  SUBROUTINE  :  CO2PhaseReactions
!  DESCRIPTION :  Administration of CO2 phase reactons
!  VARIABLES   :
!
!  CREATED     :  2025.09.30  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use HCprop
   Use Par
   Use GasRx

   Implicit None

! Call variables
    Double Precision NGas(NoHC)

! Common variables
   Integer CalcOptions(0:10)
   Common /CCalcOptions/ CalcOptions


! Matrix solver and some parameters
   External DGESV
   Integer nrhs,LDA,LDB,info
   Integer  ipiv(MaxGasRx)


   Integer Phase
   Double Precision a(MaxGasRx,MaxGasRx),b(MaxGasRx),bb(MaxGasRx), Sol(MaxGasRx),D(MaxGasRx),NGasStart(NoHC)
   Integer i,ii,k,RxNo,j,jj,Itno,Teller,N,IncludeRx(0:MaxGasRx)
   Integer CompN,Coef,np
   Logical FoundIt


!-----Start-------------------------------------------------------------------
   NGasStart=NGas
   IncludeRx=0
   IncludeRx(1:RxInfo(0))=IncludeRxINP(1:RxInfo(0))  ! These are the reaction selected by the user as input
   IncludeRx(0)=Sum(IncludeRx(1:MaxGasRx))

   Phase=1

! Settings for matrix solver
   nrhs=1
   LDA=MaxGasRx
   LDB=MaxGasRx


   Call CalcGasThermo

! Check rx direction and perform tiny steps to have both reactants and products for all reactions
   Do k=1,RxInfo(0)
      Call RxCanGo(NGas)

      Do RxNo=1,RxInfo(0)
         If (IncludeRx(RxNo)) Then
            Call  FindNmax(RxNo,NGas)
            If (GoVal(RxNo).EQ.2 .OR. GOVal(RxNo).EQ.3) Then
               Call MakeTinyStep(RxNo,NGas)
            EndIf
         EndIf
      EndDo
   EndDo

   Do RxNo=1,RxInfo(0)
      If (GoVal(RxNo).EQ.4) IncludeRx(RxNo)=0
   EndDo
   IncludeRx(0)=Sum(IncludeRx(1:MaxGasRx))

   If (IncludeRx(0).EQ.0) Goto 250

! This is where the iteration starts
   Do ItNo=1,1000
      Call SinglePhase(NGas,Phase)
      Call CalcQ(NGas,IncludeRx)

      D=0.d0 ; Sol=-8.d0
      Do RxNo=1,RxInfo(0)
         If (IncludeRx(RxNo).EQ.1) then
            D(RxNo)=1.d0 - Q(RxNo)/lnK(RxNo)
            Sol(RxNo)=Q(RxNo)/lnK(RxNo)
         EndIf
      EndDo

      Write(*,'(I4,4F15.5)')ItNo,(Sol(j) ,j=1,4)

      Call CheckSolution(IncludeRx,Sol,FoundIt)
      If (FoundIt) GOto 250

      Call CheckInclude(NGas,IncludeRx)
      Call ReOrgMatrix(IncludeRx,a,b,D)
      Call dgesv(IncludeRx(0),nrhs,a,LDA,ipiv,b,LDB,info)

      If (Info.NE.0) Then
         Write(*,*)'Info=',Info,' Solving matrix in dgesv failed'
         Read(*,*)
      EndIf

      bb=0.d0 ; ii=0
      Do RxNo=1,RxInfo(0)
         If (IncludeRx(RxNo).EQ.1) Then
            ii=ii+1
            bb(RxNo)=b(ii)
         EndIf
      EndDo
      Call TakeStep(Ngas,bb)

   EndDo

250 continue

    Call PrintResults(IncludeRx,bb,NGasStart,NGas)

   End



   SubRoutine PrintResults(IncludeRx,bb,NGasStart,NGas)
!****************************************************************************
!  SUBROUTINE  :  PrintResults
!  DESCRIPTION :  Prints results from calculation
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!*******************************************************************************
      Use GasRx
      Use HC
      Use HCprop
      Implicit None
! Call variables
      Integer IncludeRx(0:MaxGasRx)
      Double Precision NGasSTart(NoHC),NGas(NoHC),bb(MaxGasRx),NTOT

! Local variables
      Integer RxNo,N, j,i , CompN, Coef
      Integer RxComps(2),RxCompN(2,NoHC)
      Character RxNames(2,NoHC)*10, InfoNames(2)*10
      Character RxDsc*100

      Data InfoNames /"Reactants", "Products"/
!-----Start-------------------------------------------------------------------
   NTOT=sum(NGas)

   Write(*,*)
   Write(*,*)'FINISH'
   Write(*,*)"-------------------------------------------------------------"
   Write(*,*)'Component   Start Moles End Moles   End ppm     FugCoef'
   Write(*,*)"-------------------------------------------------------------"
   Do N=1,NoHC
      Write(*,'(A10,3ES12.3,F12.4)')HCNames(N),NGasStart(N),NGas(N), NGas(N)*1.d6/NTOT, Fugy(N)
   EndDo
   Write(*,*)"-------------------------------------------------------------"

   !Write(*,*)'delta'
   !Write(*,'(4E15.5)')(bb(j),j=1,RxInfo(0))
   !
   !Write(*,*)'This is the matrix:'
   !Do RxNo=1,RxInfo(0)
   !   WRite(*,'(4E15.7)')(dQQ(RxNo,j),j=1,4)
   !EndDo

   WRite(*,*)
   WRite(*,*)


   Do RxNo=1,RxInfo(0)
      IF (IncludeRxINP(RxNo).EQ.1) Then
         RxComps=0
         RxNames=""
         Do N=1,RxInfo(RxNo)
            CompN = RxComp(RxNo,N)
            Coef  = RxCoef(RxNo,N)

            If (Coef.LT.0) Then
               RxComps(1)=RxComps(1)+1
               RxNames(1,RxComps(1))=HCNames(CompN)
               RxCompN(1,RxComps(1))=CompN
            Else
               RxComps(2)=RxComps(2)+1
               RxNames(2,RxComps(2))=HCNames(CompN)
               RxCompN(2,RxComps(2))=CompN
            EndIf
         EndDO

         RxDsc=""
         Do N=1,RxComps(1)
            CompN=RxCompN(1,N)
            RxDsc=Trim(RxDsc) // HCNames(CompN)
            If (N.LT.RxComps(1))Then
               RxDsc=Trim(RxDsc) //"+"
            Else
               RxDsc=Trim(RxDsc) //" = "
            EndIF
         EndDo
         Do N=1,RxComps(2)
            CompN=RxCompN(2,N)
            IF (N.EQ.1) Then
               RxDsc=Trim(RxDsc) // " " //HCNames(CompN)
            Else
               RxDsc=Trim(RxDsc) // HCNames(CompN)
            EndIf
            If (N.LT.RxComps(2))Then
               RxDsc=Trim(RxDsc) //"+"
            EndIF
         EndDo


         Write(*,*)"-------------------------------------------------------------"
         Write(*,'(A12,i3,A3,A)')' Reaction no',RxNo,"   ", Trim(RxDsc)
         Write(*,'(A,2F15.5)')' lnQ and lnK:',Q(RxNo),lnK(RxNo)
         Write(*,*)"-------------------------------------------------------------"

         Do j=1,2
            Write(*,*)InfoNames(j) !,"    Molerates           ppm"
            Do N=1,RxComps(j)
               CompN=RxCompN(j,N)
               Write(*,'(2A,ES15.5,F15.5)')" ",HCNames(CompN),NGas(CompN),Ngas(CompN)*1.d6/NTOT
            EndDo
         EndDo
      EndIf
   EndDo
   Write(*,*)"-------------------------------------------------------------"

   End



   SubRoutine CheckSolution(IncludeRx,Sol,FoundIt)
!****************************************************************************
!  SUBROUTINE  :  CheckSolution
!  DESCRIPTION :  Check if solution is found
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!*******************************************************************************
   Use GasRx
   Implicit None

! Call variables
   Integer IncludeRx(0:MaxGasRx)
   Double Precision Sol(MaxGasRx)
   Logical FoundIt
! Local variables
   Integer RxNo
!-----Start-------------------------------------------------------------------

   FoundIt=.True.
   Do RxNo=1,RxInfo(0)
      If (IncludeRx(RxNo).EQ.1) Then
         If (dabs(1.d0-Sol(RxNo)).GT.1.d-2) Then
            FoundIt=.False.
         EndIf
      EndIf
   EndDo


   End


   SubRoutine ReOrgMatrix(IncludeRx,a,b,D)
!****************************************************************************
!  SUBROUTINE  :  ReOrgMatrix
!  DESCRIPTION :  Reorganiaces matrix when some reactions are not included
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!*******************************************************************************
   Use GasRx
   Implicit None

! Call variables
   Integer IncludeRx(0:MaxGasRx)
   Double Precision a(MaxGasRx,MaxGasRx),b(MaxGasRx),D(MaxGasRx)

! Local variables
   Integer i, ii, j, jj
!-----Start-------------------------------------------------------------------

      a=0.d0 ; b=0.d0
      ii=0
      Do i=1,RxInfo(0)
         If (IncludeRx(i).EQ.1) Then
            ii=ii+1
            b(ii)=D(i)
            jj=0
            Do j=1,RxInfo(0)
               If (IncludeRx(j).EQ.1) Then
                  jj=jj+1
                  a(ii,jj)=dQQ(i,j)
               EndIf
            EndDo
         EndIf
      EndDo

   End









   SubRoutine CheckInclude(NGas,IncludeRx)
!****************************************************************************
!  SUBROUTINE  :  CheckInclude
!  DESCRIPTION :  Checks which reactions that should be included
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!*******************************************************************************
   Use GasRx
   Implicit None



! Call variables
   Integer IncludeRx(0:MaxGasRx)
   Double Precision NGas(NoHC)

!Local variables
   Integer RxNo, N,  CompN, Coef
   Double Precision Y(NoHC)
!-----Start-------------------------------------------------------------------

      Y=NGas/Sum(NGas)

      IncludeRx(1:RxInfo(0))=IncludeRxINP(1:RxInfo(0))
      Do RxNo=1,RxInfo(0)
         If (IncludeRxINP(RxNo).EQ.1) Then   ! Dont check reactions that user have said not to include
            If (Q(RxNo)/lnK(RxNo).LT.1.d0) Then  ! Rx should go forward
               Do N=1,RxInfo(RxNo)
                  CompN = RxComp(RxNo,N)
                  Coef  = RxCoef(RxNo,N)

                  If (Coef.LT.0 .AND. Y(CompN).LT.RxLim) Then
                     IncludeRx(RxNo)=0
                  EndIF
               EndDo
            Else ! Rx should go backwards
               Do N=1,RxInfo(RxNo)
                  CompN = RxComp(RxNo,N)
                  Coef  = RxCoef(RxNo,N)

                  If (Coef.GT.0 .AND. Y(CompN).LT.RxLim) Then
                     IncludeRx(RxNo)=0
                  EndIF
               EndDo
            EndIf
         EndIf
      EndDo
      IncludeRx(0)=sum(IncludeRx(1:MaxGasRx))
   End


   SubRoutine TakeStep(Ngas,b)
!****************************************************************************
!  SUBROUTINE  :  CheckAdjustSTep
!  DESCRIPTION :  Checks calculated Newton step and calculate necesary damping
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GasRx
   Implicit None
! Call variables
   Double Precision NGas(NoHC),b(MaxGasRx),MinDamp, dN(NoHC)

! Local variables
   Integer RxNo,CompN,Coef,N, InRx(NoHC)
   Double PRecision f,Fact
!-----Start-------------------------------------------------------------------
   MinDamp=1.d0

   f=1.d0 ; Fact=1.d0

! Map components included in a reaction
   InRx=0
   Do RxNo=1,Rxinfo(0)
      Do N=1,RxInfo(RxNo)
         CompN = RxComp(RxNo,N)
         InRx(CompN)=1
      EndDo
   EndDo


   dN=0.d0
   Do RxNo=1,Rxinfo(0)
      Do N=1,RxInfo(RxNo)
         CompN = RxComp(RxNo,N)
         Coef  = RxCoef(RxNo,N)

         dN(CompN)=dN(CompN) + b(RxNo)*Coef
      EndDo
   EndDo

   Do N=1,NoHC
      If (dN(N).LT.0.d0 .AND. InRx(N).EQ.1) Then
         If (NGas(N)+dN(N).LE.0.d0) Then
            Fact=-NGas(N)/dN(N)
         EndIf
      EndIf
      If (Fact.LT.f) Then
         f=Fact*0.9d0
      EndIf

   EndDo

   dN=dN*f

   Do N=1,NoHC
      NGas(N)=NGas(N)+dN(N)
   EndDo


   End



   SubRoutine CalcQ(NGas,IncludeRx)
!****************************************************************************
!  SUBROUTINE  :  CalcQ
!  DESCRIPTION :  Calculates Q value for current composition
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GasRx
   Use HCProp
   Use Globvar
   Implicit None
! Call variables
   Integer IncludeRx(0:MaxGasRx)
   Double Precision NGas(NoHC)

! Local variables
   Logical SetToZero
   Integer i,RxNo, N, CompN, Coef, L, CoefQ, CoefL
   Double Precision NT, f(NoHC)
!-----Start-------------------------------------------------------------------
   NT=sum(NGas)
   Do i=1,NoHC
      f(i)=NGas(i)/NT*Fugy(i)*Ptot
   EndDo

   Q=0.d0

   Do RxNo=1,RxInfo(0)
      SetToZero=.False.
   !   If (IncludeRx(RxNo).EQ.1) Then
         Do N=1,RxInfo(RxNo)
            CompN = RxComp(RxNo,N)
            Coef  = RxCoef(RxNo,N)
            If (f(CompN).GT.0.d0) Then
               Q(RxNo)=Q(RxNo) + Coef*log(f(CompN))
            Else
               SetToZero=.True.
            EndIf
         EndDo
         If (SetToZero) Q(RxNo)=0.d0
   !   EndIf
   EndDO


   dQQ=0.d0               ! Set derivative for RxNo with respect to change L to 0
   Do RxNo = 1,RxInfo(0)  ! Outer loop over reactions  Q1, Q2 etc
      Do L = 1,RxInfo(0)  ! Inner loop for cross derivatives
         Do N = 1,NoHC    ! Loop over all compnents
            CoefQ = CompRx(N, RxNo) ! Get compnent number
            CoefL = CompRx(N, L)    ! Get reaction coefficeint for this component in this reaction

            If (CoefQ.NE.0) Then    ! Add contribution to derivative of coeff <> 0
               dQQ(RxNo, L) = dQQ(RxNo, L) + CoefQ * CoefL / NGas(N)
            End If
         EndDo
      EndDo
   EndDo

   Do RxNo=1,RxInfo(0)
      Do L=1,RxInfo(0)
         dQQ(RxNo,L)=dQQ(RxNo,L)/lnK(RxNo)
      EndDo
   EndDo

   End


   SubRoutine MakeTinyStep(RxNo,NGas)
!****************************************************************************
!  SUBROUTINE  :  MakeTinyStep
!  DESCRIPTION :  Determines if a reaction can go 1) both ways, 2) backwards 3) forward 4) no at all
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GasRx
   Use Par
   Implicit None

   Integer RxNO
   Double Precision NGas(NoHC)

! Local variables
   Logical ReactPos, ProdPos
   Integer N, CompN, Coef

   Double precision StepFraction, dN
!-----Start-------------------------------------------------------------------
   StepFraction=20.d0

   IF (GOVal(RxNo).EQ.3) Then
      dN=dNmax(RxNo,1)/StepFraction
   ElseIf (GoVal(RxNo).EQ.2) Then
      dN=-dNmax(RxNo,2)/StepFraction
   EndIf

   Do N = 1,RxInfo(RxNo)
      CompN = RxComp(RxNo,N)
      Coef  = RxCoef(RxNo,N)

      NGas(CompN) = NGas(CompN) + dN * Coef
   EndDo


   End



   SubRoutine RxCanGo(NGas)
!****************************************************************************
!  SUBROUTINE  :  RxCanGo
!  DESCRIPTION :  Determines if a reaction can go 1) both ways, 2) backwards 3) forward 4) no at all
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GasRx
   Use Par
   Implicit None

   Double Precision NGas(NoHC)

   Logical ReactPos, ProdPos
   Integer RxNo, N, CompN, Coef
!-----Start-------------------------------------------------------------------
   GoVal=0
   Do RxNo=1,RxInfo(0)
      ReactPos=.True. ; ProdPos=.True.

      Do N=1,RxInfo(RxNo)
         CompN = RxCOmp(RxNo,N)
         Coef  = RxCoef(RxNo,N)

         If (Coef<0) Then
            If (NGas(CompN).EQ.0.d0) Then
               ReactPos=.False.
            EndIf

         Else
            If (NGas(CompN).EQ.0.d0) Then
               ProdPos=.False.
            EndIf
         EndIf
      EndDo

      If (ReactPos .AND. ProdPos) Then !Have all, rx can go both ways
         GoVal(RxNo)=1
      ElseIf (.NOT. ReactPos .AND. .NOT. ProdPos) Then !Missing bot reactant and production - cannot go
         GoVal(RxNo)=4
      ElseIf (.Not.ReactPos) Then !Missing a reactant, can go backwards
         GoVal(RxNo) = 2
      ElseIf (.Not.ProdPos) Then ! Missing a product, can go forward
         GoVal(RxNo) = 3
      Else
         GoVal(rxNo)=5
      EndIf
   EndDo

   End


   SubRoutine FindNmax(RxNo,NGas)
!****************************************************************************
!  SUBROUTINE  :  FindNmax
!  DESCRIPTION :  Determines maximum reaction moles for a reaction
!  VARIABLES   :
!
!  CREATED     :  2025.10.01  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GasRx
   Use Par
   Implicit None

! Call variables
   Integer RxNo
   Double Precision NGas(NoHC)

! Local variables
   Integer N, CompN, Coef
!-----Start-------------------------------------------------------------------

   dNmax(RxNo,1) = Sum(NGas(1:NoHC))
   dNmax(RxNo,2) = Sum(NGas(1:NoHC))

   Do N = 1,RxInfo(RxNo)
      CompN = RxComp(RxNo,N)
      Coef = RxCoef(RxNo,N)

      If (Coef.LT.0) Then ! Dette er et product som skal brukes
         If (-dNmax(RxNo,1)*Coef.GT.NGas(CompN)) Then
            dNmax(RxNo,1) = -NGas(CompN)/Coef
         End If
      Else
         If (dNmax(RxNo,2)*Coef.GT.NGas(CompN)) Then
             dNmax(RxNo,2) = NGas(CompN)/Coef
         End If
      End If

   EndDo


   End



   SubRoutine CalcGasThermo
!****************************************************************************
!  SUBROUTINE  :  CalcGasThermo
!  DESCRIPTION :  Calculates gas reaction thermo
!  VARIABLES   :
!
!  CREATED     :  2025.09.30  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use GlobVar
   Use GasRx
   Implicit None


! Local variables
   Integer RxNo, NC, N, CompN, Coef

!-----Start-------------------------------------------------------------------

   Do RxNo=1,RxInfo(0)
      NC=RxInfo(RxNo)

      dH(RxNo)=0.d0
      dS(RxNo)=0.d0
      dG(RxNo)=0.d0
      Do N=1,NC
         CompN=RxComp(RxNo,N)
         Coef=RxCoef(RxNo,N)
         dH(RxNo)=dH(RxNo) + Coef*CompH(CompN)
         dS(RxNo)=dS(RxNo) + Coef*CompS(CompN)
      EndDo

      dG(RxNo)=dH(RxNo)*1.d3 - Temp*dS(RxNo)
      lnK(RxNo)=(-dG(RxNo)/(Temp*R))

  !    Write(*,*)dH(RxNo), dS(RxNo), lnK(RxNo)

   EndDo

   End