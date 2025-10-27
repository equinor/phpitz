
!************************************************************
!  File Type:     F90 Subroutines
!  File Content:  SolveEOS    Solves the volume
!                 EOSabc      Calculates parameters a,b and c in the Equiation of State
!                 FindVol     Administrative routine to calculate volume (compressibility)
!                 VolSolve    Solves volume of phase or possible phases and finds correct solution
!                 Volit       Volume iteration to find volume
!                 SolveQubic  Solves cubic equaation when volume can be found from quibic EOS
!                 Gibbs       Calculates Gibbs energy to determine which volume is stabel
!************************************************************


      SubRoutine SolveEOS(Ni,ZF,Phase)
!**********************************************************************
! DESC:        Solves EOS
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  
!
! Created:     20.09.2010 BKaa
! Update log:  
!**********************************************************************
      Use EOS
      Implicit None

! Call variables
      Integer Phase
      Double Precision Ni(NoHC),ZF

      Integer i
!-----Start--------------------------------------------------
      
!      Call CPA_setup(nz)
      Call EOSabc(Ni,phase)
      Call FindVol(Ni,ZF,phase)

   End

   
      SubRoutine EOSabc(Ni,phase)
!**********************************************************************
! DESC:        Calculates the a,b and c constants used in the EOS
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  EOS_a, EOS_b, EOS_c
!
! Created:     20.09.2010 BKaa
! Update log:  
!**********************************************************************
      Use EOS
      Implicit None

! Call variables
      Integer phase
      Double Precision Ni(NoHC)


! Common variables
      Double Precision dBfunc,GE,tel(NoHC),Nev(NoHC)
      Common /CGE         / GE,tel,Nev

!--Local variables
      Integer i,j
      Double Precision ledd

!--Start-------------------------------------------------------
      EOS_A=0.d0 ; EOS_B=0.d0 ; EOS_C=0.d0   ; EOS_Cvpt=0.d0

      EOS_B=Dot_Product(Ni(1:NoHC),bi(1:NoHC))
      EOS_C=Dot_Product(Ni(1:NoHC),cPen(1:NoHC))
      
      
      
      Do i=1,NoHC
         Do j=1,NoHC
            EOS_A=EOS_A + dsqrt(ai(i)*ai(j))*Ni(i)*Ni(j)*(1.d0-kij(i,j))
         EndDo
      EndDo
   
      
      
      If (EOS_Type.EQ."VPT") EOS_Cvpt=Dot_Product(Ni(1:NoHC),Cvpti(1:NoHC))
      
      Phase_A(phase)=EOS_A
      Phase_B(phase)=EOS_B
      Phase_C(phase)=EOS_C
      Phase_Cvpt(phase)=EOS_Cvpt
      
      End SubRoutine

   
      SubRoutine FindVol (Ni,ZF,Phase)
!**********************************************************************
! DESC:        Solves the volume of the EOS 
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  
!
! Created:     10.03.2012 BKaa
! Update log:  
!**********************************************************************
   Use EOS
   Use GlobVar
   Implicit None

! Call variables
   Integer Phase,Lev,J
   Double Precision Ni(NoHC),X(3),ZF

! Functions
   Double Precision GibbsE

! Local variables
   Double Precision n,RT,V
   Double Precision Fug_liq(NoHC),Fug_gas(NoHC)
!--Start--------------------------------------------------
   
   n=Sum(Ni) ; RT=Reos*Temp
  
   If (.True.) Then
      Call SolveQubic(Ni,X,J)
   Else
      Call VolSolve(Ni,Phase,ZF)
      Return
   EndIf

   
   If (J.EQ.1) Then  ! Only one solution
      ZF=X(1)
      
      V=ZF*n*Reos*Temp/Ptot
      
      If (phase.eq.-1) Then
         If (V-1.75*EOS_b.LT.0.d0) Then   ! Determine if it is a "liquid" or "gas" type phase - this is just an estimate
            Phase=1
         Else
            Phase=2
         EndIf
      EndIf
      
   Else              ! Two solutions
      If (Phase.EQ.1) Then ! Want the liquid root
         ZF=x(1)
      ElseIf (Phase.EQ.2) Then ! Want the gas root
         ZF=X(3)
      Else  ! Determine most stable root
         Lev=2   ! Calculation of fugacity, but not derivatives
         Call CalcFug(Lev,Ni,X(1),Fug_liq) 
         Call CalcFug(Lev,Ni,X(3),Fug_gas)
      
         If (GibbsE(Ni,Fug_liq).LE.GibbsE(Ni,Fug_gas)) Then
            ZF=X(1)
!            Phase=1
         Else
            ZF=X(3)
!            Phase=2
         EndIF     
      EndIf
   EndIf
   
   ZF=ZF - EOS_C/n*Ptot/(Reos*Temp)
   
   
   End
   
   
   SubRoutine VolSolve(Ni,Phase,ZF)
!**********************************************************************
! DESC:        Solves the volume of an  EOS
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  
!
! Created:     2016.06.27 Baard Kaasa www.scalesim.no
! Update log:  
!**********************************************************************
   Use EOS
   Use GlobVar
   Implicit None
   
! Call variables   
   Integer Phase,Lev
   Double Precision Ni(NoHC),ZF
! Functions
   Double Precision GibbsE
   
! Local variables
   Double Precision n,nRT,Tred,Vgas,Vliq
   Double Precision Gibbs_Liq,Gibbs_Gas,ZF_liq,ZF_gas,Fug_liq(NoHC),Fug_gas(NoHC)
   integer i

!--Start------------------------------------------------------
   n=sum(Ni)
   nRT=n*Reos*Temp
   
   Tred=n*Temp/Dot_Product(Ni(1:NoHC),Tc(1:NoHC))

   Vgas=nRT/Ptot                   ! Gas estimate - ideal solution
   Vliq=EOS_B*(2.d0 + Tred)/2.d0    ! Liquid estimate - a little larger than EOS_B
  
   
! If phase is known - calculte the solution
   If (Phase.EQ.1) Then             ! Know that it is liquid and want liquid root
      Call VOLit(Ni,Vliq,1)
      ZF=Ptot*Vliq/(nRT)
      Return
      
   ElseIf (Phase.EQ.2) Then         ! Know that it is gas and want gas root
      Call VOLit(Ni,Vgas,2)

      If (Vgas.EQ.-1.d0) Then
         Call VOLit(Ni,Vliq,1)
         Vgas=Vliq
      EndIf
      ZF=Ptot*Vgas/(nRT)
      Return
   EndIf
   

! Phase is not known - need to check both
   Call VOLit(Ni,Vliq,1)         ! Calculate liquid solution
   Call VOLit(Ni,Vgas,2)         ! Calculate gas solution   
    

      
   If (dabs((Vgas-Vliq)/Vgas).LT.1.d-4) Then        ! Two identical solutions found - i.e. there was only one solution
      Vgas=(Vgas+Vliq)/2.d0
      ZF=Ptot*Vgas/(nRT)
      
      If (Vgas-1.75*EOS_b.LT.0.d0) Then   ! Determine if it is a "liquid" or "gas" type phase - this is just an estimate
         Phase=1
      Else
         Phase=2
      EndIf
      Return
   EndIf   

! Have found two solutions - check if both are valid


   If (Vgas*Vliq.LE.0.d0) Then  ! One of the solutions are not valid
      If (Vgas.LT.0.d0) Then    ! Did not find a gas solution - liquid solution is correct
         ZF=Ptot*Vliq/(nRT)
         Phase=1
      ElseIf (Vliq.LT.0.d0) Then    ! Did not find a liq solution - gas solution is correct
         ZF=Ptot*Vgas/(nRT)
         Phase=2
      EndIf       
      Return
   EndIf
            
! Have two valid solutions, must determine the correct (Lowst Gibbs energy
   ZF_liq=Ptot*Vliq/(nRT)
   ZF_gas=Ptot*Vgas/(nRT)
   Lev=2 ! Calculate fugacity, but not derivative
   Call CalcFug(Lev,Ni,ZF_liq,Fug_liq) 
   Call CalcFug(Lev,Ni,ZF_gas,Fug_gas) 
   Gibbs_Liq=GibbsE(Ni,Fug_liq)
   Gibbs_Gas=GibbsE(Ni,Fug_gas)   

   If (Gibbs_Liq.LE.Gibbs_Gas) Then    ! The liquid solution was correct
      ZF=ZF_liq
      Phase=1
   Else                                ! The gas solution was correct
      ZF=ZF_gas
      Phase=2
   EndIf   
   
   
   End
   
   
   SubRoutine Volit(Ni,V,Phase)
   Use EOS
   Use GlobVar
   Use F_funcs
   Implicit None
   
! Call variables
   Double Precision Ni(NoHC),V
   Integer Phase
! Local variables
   Integer, Parameter::Maxit=50
   Double Precision, Parameter::dVkrit=1.d-4,Hkrit=1.d-4
   Double Precision Vmin,Vmax,n,RT,nRT
   Integer i,Lev
   
   Double Precision H,dHdV,dV,dFdV,d2FDVV
!--Start------------------------------------------------------
   Vmin=EOS_B
   Vmax=-1.d0
   Lev=1       ! Only volume derivatives
   
   n=Sum(Ni(1:NoHC)) ;   RT=Reos*Temp ;     nRT=n*RT
   
   Do i=1,Maxit
      Call EOS_F_func(Lev,Ni,V)
      
      dFdV=dF_dV                   ! Total pressure contribution
      H=nRT/V - RT*dFdV - Ptot
      
      d2FdVV=d2F_dVV
      dHdV=-nRT/V**2 - RT*d2FdVV
         
      dV=-H/dHdV                          ! Calculate Newton step
      
      If (Vmin*Vmax.GT.0.d0) Then         ! This can only occure after both upper and lower limit is set
         IF (dV.GT.2.d0*(Vmax-Vmin)) Then ! Have found a top point
            V=-1.d0
            Return
         EndiF
      EndIf
         
      If (dHdV.GT.0.d0) Then              ! In 2 phase reqion corresponding to intermediate solution
         If (Phase.EQ.1 .OR. PHase.EQ.3) then
            Vmax=V
            V=(V+Vmin)/2.d0
         ElseIf (Phase.EQ.2) Then
            Vmin=V
            If (Vmax.GT.0.d0) Then
               V=(V+Vmax)/2.d0
            Else
               V=1.5d0*V
            EndIf
         EndIf
         dV=0.d0
      Else                                ! dHdV has correct slope
         If (H.GT.0.d0) Then              ! Volume is too small
            Vmin=V
         ElseIf (H.LT.0.d0) Then          ! Volume is too high
            Vmax=V
         Else
            Return                        ! Found the solution by coinsidense - will probably never happen
         EndIf

         If (V+dV.LE.Vmin) Then           ! Newton step goes below lower limit
            dV=-0.5d0*(V-Vmin)
         ElseIf (Vmax.GT.0.d0 .AND. V+dV.GE.Vmax) Then       ! Newton step above upper limit - if upper limit is set
            dV=0.5d0*(Vmax-V)
         EndIf

         V=V+dV                  ! Calculate new Volume
!         write(*,'(I3,5F15.5)')i,H,V,dV
      EndIf
      
      
      If (dabs(dV/V).LT.dVkrit .AND. dabs(H).LT. Hkrit) Then
         Return
      EndIf
   EndDo
   V=-2.d0
   End
   
   
   SubRoutine SolveQubic(Ni,X,J)
!**********************************************************************
! DESC:        Solves the volume of a cubic EOS
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  
!
! Created:     2016.06.27 Baard Kaasa www.scalesim.no
! Update log:  
!**********************************************************************
   Use EOS
   Use GlobVar
   Implicit None

! Call variables
   Integer Phase,Lev
   Double Precision Ni(NoHC),n,ZF,nRT,BB,AA,RR,SS,TT,X(3)
   Integer J
!--Start------------------------------------------------------
   n=sum(Ni)
   nRT=n*Reos*Temp
   
   BB=EOS_B*Ptot/nRT
   AA=EOS_A*ptot/nRT**2
      
   RR=BB*(Delta1+Delta2) - BB - 1.d0
   SS=Delta1*Delta2*BB**2 - (Delta1+Delta2)*(BB+BB**2) + AA
   TT=-Delta1*Delta2*(BB**2 + BB**3) - AA*BB
   
   Call CubEqSolve(RR,SS,TT,X,J)
   
   End
   
   
   Double Precision Function GibbsE(Ni,Fug)
!****************************************************************************
!  SUBROUTINE  :  GibbsE
!  DESCRIPTION :  Calculate Gibbs excess energy G=Sum[ xi*ln(xi*fi) ]
!  VARIABLES   :  Nc  - Mole vector
!                 Fug - Fucacity coefficients
!
!  CREATED     :  2016.06.30  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Par
   Implicit None
! Call variables
   Double Precision Ni(NoHC),Fug(NoHC)
   
! Local variables
   Double Precision Sum
   Integer i

!-Start----------------------------------------------------------------------
   Sum=0.d0
   DO i=1,NoHC
      If (Ni(i).GT.0.d0) Then
         Sum=Sum + Ni(i)*dlog(Ni(i)*Fug(i))
      EndIf
   EndDo
   GibbsE=Sum
   End
   
   