!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Eq_ActCoeff.f90
!  DESCRIPTION :  Calculate activity coefficients
!  CREATED     :  2013.01.29  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  Eq_Stoech   - Calculates stochiometric equilibrium constants from thermodynamic and activity coefficients
!
!****************************************************************************

      SubRoutine Eq_Stoech
!****************************************************************************
!  SUBROUTINE  :  Eq_Stoech
!  DESCRIPTION :  Calculates stoichiometric equilibrium constants
!  VARIABLES   :  Eq_gam - Activity coefficients fraction for a equilibria
!                 Eq_spec - Contains the species that are involved in each equilibria
!                 Eq_exp  - Contains the exponent (number of each specie) involved in the equilibra
!                 Eq_stoec - Stoichiometric equilibrium constants
!
!  CREATED     :  2013.01.29  Baard Kaasa www.scalesim.no
!  CHANGED     :  2013.05.03  Baard Kaasa, Updated name
!
!****************************************************************************
      Use Eq
      Use Aqua
      Use HCprop
      Implicit None


! Local variables
      Integer EqN, i, Ion,HCindex
      Double Precision nyx

!----Start-------------------------------------------------------------------

      Eq_Gam=1.d0                                    ! Start by setting all to 1

  !    Call Eq_Therm_Calc

    !  write(*,*)'stoec'
      Do EqN=1,NoEq                                  ! Run through all equilibria
         Do i=1,Eq_spec(EqN,0)                       ! Run through all species in each equilibra
            Ion=Eq_spec(EqN,i)                       ! Specie no
            nyx=Eq_exp(EqN,i)                        ! No of species in rx, i.e. exponent in expression  (negative means reactant)


            Eq_Gam(EqN)=Eq_Gam(EqN)*Aqg(Ion)**(nyx)  ! Activity coefficient of ion times exponent

         EndDo

         Eq_Stoec(EqN)=Eq_Therm(EqN)/Eq_Gam(EqN)     ! Calculate stoichiometric equilibrium constants from thermodynamic and gammas

      EndDo


      Do i=1,NGasSol                                  ! Correct Henrys law equilibria with gas-fugacity
         EqN=GasSolub(i,1)
         HCindex=GasSolub(i,2)
         Eq_Stoec(EqN)=Eq_Stoec(EqN)*Fugy(HCindex)
      EndDo

      End