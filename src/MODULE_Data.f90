
!************************************************************
!  File Type:     MODULE
!  File Content:  All parameters and data used in models
!
!************************************************************


MODULE Vap_P      ! Parameters used for calculation of vapour pressures

      SAVE

      ! Water
      Double Precision K_H2O(4),V_H2O(3),Fi_Ref_H2O(5)
      DATA K_H2O/ -7862.5D0, 86.707D0, -11.754D0, 0.010672D0 /                ! Temperature function
      DATA V_H2O/ 22.898D0, -3.6442D-2, 6.8107D-5/                            ! Volume change for Poyting correction
      Data Fi_Ref_H2O /2.24d-11, -4.23d-08, 2.59d-05, -6.59d-03, 1.60d+00 /   ! Reference fugacity of pure water

   End Module Vap_P



   MODULE SolventPar

      SAVE
      Double Precision DensPar(11),VapPar(8),FiPar(5)

   End Module SolventPar