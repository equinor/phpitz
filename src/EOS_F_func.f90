      SubRoutine EOS_F_func(Lev,Ni,V)
!**********************************************************************
! DESC:        Calculates all the F-functions from the EOS
!
! PARAMETERS:
!
!
! Created:     20.03.2012 BKaa
! Update log:
!**********************************************************************
      Use GlobVar
      Use EOS
      Use F_funcs
      Implicit None

! Call variables
      Integer Lev
      Double Precision Ni(NoHC),V

! Function variables
      Double Precision hfunc1,hVfunc1,hBfunc1,hVVfunc1
      Double Precision gfunc1,gVfunc1,gBfunc1,gVVfunc1
      Double Precision dBfunc,dDfunc


! Local variables
      Integer i
      Double Precision n,D,B

      Double Precision h,hV,hB,hVV
      Double Precision g,gV,gB,gVV

      Double Precision Fn,FV,FB,FD
      Double Precision FVV


      Double Precision dBdi(NoHC),dDdi(NoHC)
!-----Start--------------------------------------------------
      n=sum(Ni)

      If (EOS_Type.EQ."SRK" .OR. EOS_Type.EQ."PR") Then
         D=EOS_a  ; B=EOS_b

   ! Calculate F only
	      h   = hfunc1  (V,B)
         g   = gfunc1  (V,B)
         F  = -n*g  - D/Temp*h

         IF (Lev.EQ.0) Return

   ! Calculate volume derivatives
         hV  = hVfunc1 (V,B)
         gV  = gVfunc1 (V,B)

         hVV = hVVfunc1(V,B)
         gVV = gVVfunc1(V,B)

         FV = -n*gV - D/Temp*hV
         FVV = -n*gVV - D/Temp*hVV

         dF_dV = FV
         d2F_dVV = FVV

         IF (Lev.EQ.1) Return

   ! Calculate 1. order compositional derivatives.

         Do i=1,NoHC
            dBdi(i)=dBfunc(i)
            dDdi(i)=dDfunc(i,Ni)
         EndDo

         gB  = gBfunc1 (V,B)
         hB  = hBfunc1 (V,B)
         Fn = -g
         FB = -n*gB - D/Temp*hB
         FD = -h/Temp


         Do i=1,NoHC
            dF_dni(i)=Fn + FB*dBdi(i) + FD*dDdi(i)
         EndDo

      EndIf

      End