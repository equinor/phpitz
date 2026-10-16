
!**********************************************************************
!  File Type:     F90 Subroutines
!  File Content:  EOSconst       Subroutine for calculation of EOS constants (non-compositional)
!                 EOSTypeConst   Finds the correct type of EOS and calculates constants omega_A and Omega_B
!                 EOSmfunc       Calculates the m in the EOS
!
!**********************************************************************


      SubRoutine EOSconst
!**********************************************************************
! DESC:        Calculates constants for the EOS only dependent on P and T, not composition
!
! CALLED BY:   -
!
! CALLS:       None
!
! PARAMETERS:  ai,bi    critical parameters for the EOS
!
! Created:     2016.06.26  Baard Kaasa www.scalesim.no
! Updated:
!**********************************************************************
      Use GlobVar
      Use EOS
 !     Use LOGRes
      Implicit None

! Functions
      Double Precision EOSmfunc

!--Local variables
      Integer i
      Double Precision RTc
      Double Precision aci(NoHC)

!--Start---------------------------------------------------------------------

      Call EOSTypeConst

      aci=0.d0 ; bi=0.d0 ; ai=0.d0 ; Cvpti=0.d0


      If (EOS_Type.EQ."SRK" .OR. EOS_Type.EQ."PR") Then     ! Calculate from acentric factor and critical properties
         Do i=1,NoHC
            RTc=Reos*Tc(i)
            EOSmi(i)  = EOSmfunc(Acent(i))

            IF (Pc(i).GT.0.d0) Then
               aci(i) = Omega_A*(RTc)**2/Pc(i)
               bi(i)  = Omega_B*(RTc)/Pc(i)
            EndIf
         EndDo


      EndIf

      Do i=1,NoHC
         IF (Tc(i).GT.0.d0) Then
		      ai(i) = aci(i)*(1.d0 + EOSmi(i)*(1.d0-dsqrt(Temp/Tc(i))))**2
         EndIf
      EndDo

	   End SubRoutine EOSconst


      SubRoutine EOSTypeConst
!**********************************************************************
! DESC:        Sets the parameters specific for the selected EOS
!
! CALLS:       None
!
! PARAMETERS:  Omega_A, Omega_B  Constants used in EOS
!
!**********************************************************************
      Use EOS
      Implicit None

! Local variables
      Integer cpar

      Double Precision D1,D2,F1,F2,Apluss,Y,Omega_Rat,Zcrt
!--Start---------------------------------------------------------------------

      If (EOS_Type.EQ."SRK") Then
         cpar=0
!         CPen=CpenSRK
!         HV_Lambda=dlog(2.d0)

      ElseIf (EOS_Type.EQ."PR") Then
         cpar=1
!         CPen=CpenPR
!         HV_Lambda=1.d0/(2.d0*dsqrt(2.d0))*dlog((dsqrt(2.d0)+1.d0)/(dsqrt(2.d0)-1.d0))


      Else
         cpar=0      ! Sets SRK as a default
!         CPen=CpenSRK
!         HV_Lambda=dlog(2.d0)

      EndIf

! Calculate delta 1 and 2. From Michelsen/Mollerup, Thermodynamic models: Fundamentals and computational aspects, 2nd ed, 2007, page75
      Delta1=0.5d0*(cpar + 1.d0 + dsqrt((cpar+1.d0)**2 + 4.d0*cpar))
      Delta2=0.5d0*(cpar + 1.d0 - dsqrt((cpar+1.d0)**2 + 4.d0*cpar))

      D1=1.d0 + cpar
      D2=-cpar

      F1=(1.d0 + delta1)**(1.d0/3.d0)
      F2=(1.d0 + delta2)**(1.d0/3.d0)

      Apluss=F1*F2*(F1+F2)
      Y=1.d0 + Apluss
      Omega_Rat=(3.d0*Y*(Y+D1) + D1*D1 - D2)/(3.d0*Y + cpar)

      Zcrt=Y/(3.d0*Y + cpar)

      Omega_B=Zcrt/Y
      Omega_A=Omega_Rat*Omega_B

      End SubRoutine EOSTypeConst


      Double Precision Function EOSmfunc(w)
!**********************************************************************
! DESC:        Calculates m as function of Acentric factor
!
! CALLED BY:   EOSconst
!
! CALLS:       None
!
! PARAMETERS:  w = acentric factor
!
! Created:     2016.06.26  Baard Kaasa www.scalesim.no
! Update log:
!**********************************************************************
      Use EOS

      Implicit None

! Call variables
      Double precision w

!--Local variables
      Double Precision mSRK(3),mPR(3)

	   Data mSRK /0.48d0   ,1.574d0  ,-0.176d0  /
	   Data mPR  /0.37464d0,1.54226d0,-0.26992d0/

!-----Start------------------------------------------------------------

	   if (EOS_Type.EQ."SRK" .OR. EOS_Type.EQ."CPA") Then
         EOSmfunc=mSRK(1) + mSRK(2)*w + mSRK(3)*w**2
      ElseIF (EOS_Type.EQ."PR") Then
         EOSmfunc=mPR(1)  + mPR(2)*w  + mPR(3)*w**2
	   EndIf

	   End Function EOSmfunc
