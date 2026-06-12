   SubRoutine CalcFug(Lev,Ni,ZF,Fug)
!****************************************************************************
!  SUBROUTINE  :  CalcFug
!  DESCRIPTION :  Calculates fugacity coefficients from EOS
!  VARIABLES   :
!
!  CREATED     :  2016.06.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use Par
   Use GlobVar
   Use F_funcs
   Use EOS
   Implicit None

! Call in variables
   Integer Lev
   Double Precision Ni(NoHC),ZF,Fug(NoHC)

! Call out variale
   Double Precision V

! Local variables
   Integer i
!-Start----------------------------------------------------------------------
   Fug=0.d0
   V  = Sum(Ni)*ZF*Reos*Temp/Ptot

   Call EOS_F_func(Lev,Ni,V)

   If (Lev.EQ.0) Return ! Calucation of F only
   If (Lev.EQ.1) Return ! Calculation of volume derivative only

 ! Calcualte fugacity coefficients
   Do i=1,NoHC
!      Fug(i)=dexp(dF_dni(i) + dQdni(i) - dlog(ZF)) ! When CPA is included
      Fug(i)=dexp(dF_dni(i) - dlog(ZF))
   EndDo

   If (Lev.EQ.2) Return ! Calculation of fugacity (1. derivative)

   If (Lev.EQ.3) Return ! Calculated fugacity derivatives


! Calculate 1. order comp derivatives of F (fugacity from EOS)
! Calculate 1. order comp derivatives of Q (fugacity from CPA)
!   If (nassoc(0).GT.0) Call CPA_comp_deriv(nc,V)


   End

     SubRoutine CalcFugDeriv(Ni,ZF,Fug,dFug)
!****************************************************************************
!  SUBROUTINE  :  CalcFug
!  DESCRIPTION :  Calculates fugacity coefficients from EOS
!  VARIABLES   :
!
!  CREATED     :  2016.06.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Use Par
   Use GlobVar
   Use F_funcs
   Use EOS
   Implicit None

! Call in variables
   Integer Lev
   Double Precision Ni(NoHC),ZF,Fug(NoHC),dFug(NoHC,NoHC)

! Call out variale
   Double Precision V

! Local variables
   Integer i,j
   Double Precision SumN,dPdV,dPdni(NoHC),d2F_dVdn

!-Start----------------------------------------------------------------------
   Lev=4
   dFug=0.d0
   SumN=Sum(Ni)

   V  = SumN*ZF*Reos*Temp/Ptot

   Call EOS_F_func(Lev,ni,V)

   If (Lev.EQ.0) Return ! Calucation of F only
   If (Lev.EQ.1) Return ! Calculation of volume derivative only


 ! Calcualte fugacity coefficients
   Do i=1,NoHC
!      Fug(i)=dexp(dF_dni(i) + dQdni(i) - dlog(ZF)) ! When CPA is included
      Fug(i)=dexp(dF_dni(i) - dlog(ZF))
   EndDo

   If (Lev.EQ.2) Return ! Calculation of fugacity (1. derivative)

!     dPdV=-R*Temp*(d2F_dVV + d2Q_dV2) - Sumn*R*Temp/V**2      ! Includes EOS and CPA
      dPdV=-Reos*Temp*(d2F_dVV ) - SumN*Reos*Temp/V**2      ! Includes EOS

      Do i=1,NoHC
         d2F_dVdn=d2F_dnidV(i) !+ D2Q_DVdni(i)                  ! Sum contribution from EOS and CPA
         dPdni(i)=Reos*Temp*(-d2F_dVdn + 1.d0/V)                  ! Includes EOS and CPA
      EndDO

      Do i=1,NoHC
         Do j=1,NoHC
           dFug(i,j)=d2F_dninj(i,j)  + 1.d0/Sumn + 1.d0/(Reos*Temp)*dPdni(i)*dPdni(j)/dPdV !+ D2Q_Dninj(i,j)
!           write(*,*)i,j,d2F_dninj(i,j)
         EndDo
      EndDo



   If (Lev.EQ.3) Return ! Calculated fugacity derivatives


! Calculate 1. order comp derivatives of F (fugacity from EOS)
! Calculate 1. order comp derivatives of Q (fugacity from CPA)
!   If (nassoc(0).GT.0) Call CPA_comp_deriv(nc,V)


   End