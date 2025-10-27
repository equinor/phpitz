!************************************************************
!  File Type:     MODULE
!  File Content:  Modules containing variables for EOS calculations
!  
!  Created:       2016.06.27 /Bkaa ScaleSim
!  MOdified:
!
!  Content:       Module EOS  General EOS parameters
!************************************************************

MODULE EOS
   Use Par
!   Implicit None
   Save

   Double Precision, Parameter :: Reos=83.14473d0    ! Gas constant for EOS. Units = [cm3 bar/K mole]. Gives volume in cm3.
!  Double Precision, Parameter :: Reos=83.147295d0    ! Gas constant for EOS. Units = [cm3 bar/K mole]. Gives volume in cm3.
   
! EOS parameters
      Character EOS_Type*5       !"SRK" or "PR" or "VPT" or "CPA"
      Double Precision Delta1, Delta2  ! Variables dependent on type of EOS

! Pure component parameters
      Character CompName(NoHC)*10
      Double Precision Pc(NoHC),Tc(NoHC),Acent(NoHC),Mw(NoHC),Dens(NoHC),CPen(NoHC),Zcrit(NoHC),Vcrit(NoHC)

! Binary interaction parameters
      Double Precision Kij(NoHC,NoHC)

! Calculated EOS variables
      Double Precision Omega_A,Omega_B,Omega_Cvpt
      Double Precision ai(NoHC),bi(NoHC),EOSmi(NoHC),Cvpti(NoHC)
      
! EOS variables
      Double Precision EOS_a,EOS_b,EOS_c,EOS_Cvpt
      Double Precision Phase_A(0:2),Phase_B(0:2),Phase_C(0:2),Phase_Cvpt(0:2)
      
   END MODULE EOS
   
   
   

   
   
   MODULE F_funcs
      Use Par
      Implicit None

      Double Precision F
      Double Precision dF_dni(NoHC),d2F_dninj(NoHC,NoHC)
      
      Double Precision dF_dT, d2F_dTT
      Double Precision dF_dV, d2F_dVV
      
      Double Precision d2F_dVdT
      
      Double Precision d2F_dnidT(NoHC), d2F_dnidV(NoHC)
      
End MODULE F_funcs