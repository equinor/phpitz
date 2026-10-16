!************************************************************
!  File Type:     F90 Subroutines
!  File Content:  Derivatives of mixing rules
!
!
!  Created:       2016.06.28 Baard Kaasa www.scalesim.no
!  Modified:
!************************************************************


   Double Precision Function dBfunc(i)
!************************************************************
!  Descr: Compositional derivative of B
!  Input: i - component number i
!
!************************************************************
   Use EOS
   Implicit None
   Integer i
!--Start----------------------------------------------------
   dBfunc=bi(i)
   End

   Double Precision Function d2Bfunc(i,j)
!************************************************************
!  Descr: Second order compositional derivative of B
!  Input: i - component number i
!
!************************************************************
   Use EOS
   Implicit None
   Integer i,j
!--Start----------------------------------------------------
   d2Bfunc=0.d0      ! If B=Sum[b(i)*x(i)] the second derivative is zero. Included for consistency
   End


      Double Precision Function dCvptfunc(i)
!************************************************************
!  Descr: Compositional derivative of C in VPT EOS
!  Input: i - component number i
!
!************************************************************
   Use EOS
   Implicit None
   Integer i
!--Start----------------------------------------------------
   dCvptfunc=Cvpti(i)
   End


   Double Precision Function d2Cvptfunc(i,j)
!************************************************************
!  Descr: Second order compositional derivative of C in VPT EOS
!  Input: i - component number i
!
!************************************************************
   Use EOS
   Implicit None
   Integer i,j
!--Start----------------------------------------------------
   d2Cvptfunc=0.d0      ! If BC=Sum[C(i)*x(i)] the second derivative is zero. Included for consistency
   End



   Double Precision Function dDfunc(i,Ni)
!************************************************************
!  Descr: Compositional derivative of D with respect to i
!  Input: i  - Component number
!         Nc - Molevector
!
!************************************************************
   Use EOS
   Implicit None
   Double Precision Ni(NoHC)
   Integer i

! Local variables
   Double Precision Sum
   Integer k
!--Start----------------------------------------------------

   Sum=0.d0
   Do k=1,NoHC
      Sum=Sum + Ni(k)*(1.d0 - Kij(i,k))*dsqrt(ai(i)*ai(k))
   EndDo
   dDfunc=2.d0*Sum
   End


   Double Precision Function d2Dfunc(i,j)
!************************************************************
!  Descr: Second compositional derivative of D with respect to i and j
!  Input: i,j - component numbers
!
!************************************************************
   Use EOS
   Implicit None

   Integer i,j
! Local variables
!--Start----------------------------------------------------
   d2Dfunc=2.d0*(1.d0 - Kij(i,j))*dsqrt(ai(i)*ai(j))
   End


   SubRoutine dDdTfunc(Ni,DT,DTT,DiT)
!************************************************************
!  Descr: Derivatives of D with T, second order T and comp+T
!  Input: Nc - component vector
!
!************************************************************
   Use EOS
   Use GlobVar
   Implicit None

   Double Precision Ni(NoHC),DT,DTT,DiT(NoHC)

! Local variables
   Integer i,j
   Double Precision Cip,qi,dqidT,d2qidTT,dlnaidT(NoHC),d2lnaidTT(NoHC)
   Double Precision aij,daij_dT(NoHC,NoHC),d2aij_dT(NoHC,NoHC)

!--Start----------------------------------------------------

      Do i=1,NoHC                                 ! Calculate temperature derivatives of ai
         Cip    = -EOSmi(i)/dsqrt(Tc(i))
         qi     = 1.d0 + EOSmi(i) + Cip*dsqrt(Temp)
         dqidT  = 0.5d0*Cip/dsqrt(Temp)
         d2qidTT= -0.25d0*Cip*Temp**(-1.5d0)

         dlnaidT(i)=2.d0/qi*dqidT
         d2lnaidTT(i)=-2.d0/qi**2*dqidT**2 + 2.d0/qi*d2qidTT
      EndDo


      Do i=1,NoHC                                 ! Calculate temperature derivatives of aij
         Do j=1,NoHC
            aij=dsqrt(ai(i)*ai(j))*(1.d0-Kij(i,j))
            daij_dT(i,j)=0.5d0*aij*(dlnaidT(i) + dlnaidT(j))
            d2aij_dT(i,j)=0.5d0*daij_dT(i,j)*(dlnaidT(i) + dlnaidT(j)) + 0.5d0*aij*(d2lnaidTT(i) + d2lnaidTT(j))
         EndDo
      EndDo


      DT=0.d0                                      ! Calculate temperature derivatives of D
      DTT=0.d0
      Do i=1,NoHC
         Do j=1,NoHC
            DT =DT  + Ni(i)*Ni(j)*daij_dT(i,j)
            DTT=DTT + Ni(i)*Ni(j)*d2aij_dT(i,j)
         EndDo
      EndDo


      DiT=0.d0                                     ! Calculate T and comp derivative of D
      Do i=1,NoHC
         Do j=1,NoHC
            DiT(i)=DiT(i) + Ni(j)*daij_dT(i,j)
         EndDo
         DiT(i)=2.d0*DiT(i)
      EndDo

   End