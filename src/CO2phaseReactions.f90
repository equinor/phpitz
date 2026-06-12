

   Subroutine CO2PhaseReactionsNew(AllComp)
!****************************************************************************
!  SUBROUTINE  :  CO2PhaseReactions
!  DESCRIPTION :  Administration of CO2 phase reactons
!  VARIABLES   :
!
!  CREATED     :  2026.03.17  Baard Kaasa www.scaleconsult.com
!  CHANGED     :
!****************************************************************************

   Use Par
   use gibbs_data
   use gibbs_solver
   Use HCprop
   Use Globvar
   Use ResHandling
   Use GasRx
   Implicit None

   Real(dP) :: AllComp(NoComp)

   Real(dp) :: N(ncomp), N_feed(ncomp)
   Real(dp) :: Gref(ncomp), fug_coeff(ncomp)
   Real(dp) :: P, T
   Real(dp) :: Gtot_init, Gtot_final
   Real(dp) :: b_feed(nelem), b_final(nelem)
   Integer  :: ierr, i,k
   character(len=1) :: elem_label(nelem)

   Real(dP) :: Gas_phase_ntot
   Real(dP) :: total_gibbs


  Real(dp) :: lam(nelem)            ! Lagrange multipliers from solver
  logical  :: sp_active(ncomp)      ! species active flags from solver
  Real(dp) :: delta_solid, Gref_crit

!-----Start-------------------------------------------------------------------

   Call CalcThermo(Gref)
   P=Ptot
   T=Temp

   elem_label(1)='H'
   elem_label(2)='O'
   elem_label(3)='N'
   elem_label(4)='S'

   N_feed=AllComp

   fug_coeff=1.0_dp
   fug_coeff(1:NoHC)=Fugy(1:NoHC)

   N = N_feed

  write(*,'(A)') ''
  write(*,'(A)') '================================================================='
  write(*,'(A)') '       FEED COMPOSITION'
  write(*,'(A)') '================================================================='
  write(*,'(A20, A18, A14, A8)') 'Species', 'Moles', 'Reactive', 'Phase'
  write(*,'(A)') '----------------------------------------------------------------------'
  do i = 1, ncomp
    if (.not. is_reactive(i)) then
      write(*,'(A20, ES18.6, A14, A8)') comp_name(i), N_feed(i), '  No (inert)', '  -'
    else if (CompType(i) == GAS) then
      write(*,'(A20, ES18.6, A14, A8)') comp_name(i), N_feed(i), '  Yes', '  Gas'
    else
      write(*,'(A20, ES18.6, A14, A8)') comp_name(i), N_feed(i), '  Yes', '  Solid'
    end if
  end do

  call element_totals(N_feed, b_feed)
  write(*,'(A)') ''
  write(*,'(A)') ' Element totals (reactive species):'
  do i = 1, nelem
    write(*,'(A4,A1,A2,ES14.6)') '   ', elem_label(i), ' =', b_feed(i)
  end do

  Gtot_init = total_gibbs(N, Gref, fug_coeff, T)
  write(*,'(A,ES18.8,A)') ' Initial G_total = ', Gtot_init, ' J'
  write(*,'(A)') ''

  call solve_gibbs_min(N, Gref, fug_coeff, T, ierr, lam, sp_active)

  Gtot_final = total_gibbs(N, Gref, fug_coeff, T)

  write(*,'(A)') ''
  write(*,'(A)') '================================================================='
  write(*,'(A)') '       EQUILIBRIUM COMPOSITION'
  write(*,'(A)') '================================================================='
  write(*,'(A20, A18, A18, A18, A8)') 'Species', 'Feed (mol)', 'Equil (mol)', 'Mole frac', 'Phase'
  write(*,'(A)') '------------------------------------------------------------------------------------'

  do i = 1, ncomp
    if (CompType(i) == SOLID) then
      ! Compute Delta for this solid
      delta_solid = Gref(i) / (Rgas * T)
      do k = 1, nelem
        delta_solid = delta_solid + lam(k) * atom_matrix(k,i) / (Rgas * T)
      end do
      if (sp_active(i)) then
        write(*,'(A20, ES18.6, ES18.6, A18, A8)') &
          comp_name(i), N_feed(i), N(i), '     (solid)', '  Solid'
      else
        write(*,'(A20, ES18.6, ES18.6, A18, A8)') &
          comp_name(i), N_feed(i), N(i), '  (not stable)', '  Solid'
      end if
      ! Print Delta and critical Gref on next line
      Gref_crit = 0.0_dp
      do k = 1, nelem
        Gref_crit = Gref_crit - lam(k) * atom_matrix(k,i)
      end do
      write(*,'(A20, A, F10.4, A, F14.2, A, F14.2, A)') &
        '', '  Delta=', delta_solid, '  Gref_crit=', Gref_crit, &
        '  Gap=', Gref(i) - Gref_crit, ' J/mol'
    else
      write(*,'(A20, ES18.6, ES18.6, ES18.6, A8)') &
        comp_name(i), N_feed(i), N(i), N(i)/gas_phase_ntot(N), '  Gas'
    end if
  end do

! Write results back to AllComp
  AllComp=N

  call element_totals(N, b_final)
  write(*,'(A)') ''
  write(*,'(A)') ' Element balance verification:'
  write(*,'(A20, A18, A18, A18)') 'Element', 'Feed', 'Equilibrium', 'Rel. Error'
  write(*,'(A)') '---------------------------------------------------------------------------'
  do i = 1, nelem
    if (abs(b_feed(i)) > 0.0_dp) then
      write(*,'(A20, ES18.6, ES18.6, ES18.6)') &
        elem_label(i), b_feed(i), b_final(i), &
        (b_final(i)-b_feed(i))/b_feed(i)
    end if
  end do

  write(*,'(A)') ''
  write(*,'(A,ES18.8,A)') ' Initial  G_total = ', Gtot_init,  ' J'
  write(*,'(A,ES18.8,A)') ' Final    G_total = ', Gtot_final, ' J'
  write(*,'(A,ES18.8,A)') ' Delta G          = ', Gtot_final - Gtot_init, ' J'
  write(*,'(A)') ''

  if (ierr == 0) then
    write(*,'(A)') ' Status: CONVERGED SUCCESSFULLY'
  else
    write(*,'(A)') ' Status: WARNING - DID NOT CONVERGE'
  end if
  write(*,'(A)') '================================================================='

   End




   SubRoutine CalcThermo(Gref)
!****************************************************************************
!  SUBROUTINE  :  CalcThermo
!  DESCRIPTION :  Calculates Gibbs energry for each component G=H-TS
!  VARIABLES   :
!
!  CREATED     :  2026.03.19  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Use Par
   Use Globvar
   Use GasRx
   Implicit None

! Call variables
   Double Precision :: Gref(NoCOmp)

! Local variables
   Integer :: i
!-----Start-------------------------------------------------------------------

   ! Calculate reference Gibbs energy Gref=Hf - TxSf  Hf, Sf = Enthalpy/Entropy of formation
   ! Calculation can be improved by including temperature functions for H and S using the heat capacity

   Do i=1,NoComp
      Gref(i)=CompH(i)*1.d3 - Temp*CompS(i) ! Units are: H is kJ/mole, S is J/mole K, Gref= J/mole
   EndDO

   End


