


  !----------------------------------------------------------------------------
  ! Compute Ntot for gas phase only (solids don't contribute)
  !----------------------------------------------------------------------------
  Function gas_phase_ntot(N) result(Ntot)
   Use gibbs_data
   Use gibbs_solver
   Implicit None

    Real(dp), intent(in) :: N(ncomp)
    Real(dp) :: Ntot
    Integer  :: i

    Ntot = 0.0_dp
    do i = 1, ncomp
      if (CompType(i) == GAS) Ntot = Ntot + N(i)
    end do
   end function






  !----------------------------------------------------------------------------
  function total_gibbs(N, Gref, fug_coeff,  T) result(Gtot)
   Use GlobVar
   Use gibbs_data
   Use gibbs_solver

   Implicit None

   Real(dP) :: Gas_phase_ntot

    Real(dp), intent(in) :: N(ncomp), Gref(ncomp), fug_coeff(ncomp),  T
    Real(dp) :: Gtot, Ntot_gas, yi, mu_i
    Integer  :: i

    Ntot_gas = gas_phase_ntot(N)

    Gtot = 0.0_dp
    do i = 1, ncomp
      if (N(i) > Nmin * 10.0_dp) then
        if (CompType(i) == SOLID) then
          ! Solid: activity = 1, so mu = Gref
          mu_i = Gref(i)
        else
          ! Gas: mu = Gref + RT * ln(yi * fi * P)
          yi = N(i) / Ntot_gas
          mu_i = Gref(i) + Rgas * T * log(yi * fug_coeff(i) * Ptot)
        end if
        Gtot = Gtot + N(i) * mu_i
      end if
    end do
   end function



  !----------------------------------------------------------------------------
   SubRoutine CalcCompG(N,Gref,Fug_coeff,T)
   Use GlobVar
   Use gibbs_data
   Use gibbs_solver

   Implicit None

   Real(dP) :: Gas_phase_ntot

   Real(dp), intent(in) :: N(ncomp), Gref(ncomp), fug_coeff(ncomp),  T

   Real(dp) :: Ntot_gas, yi
   Integer :: i

   Ntot_gas = gas_phase_ntot(N)
   Do i=1, ncomp
      if (N(i) > Nmin * 10.0_dp) then
        if (CompType(i) == SOLID) then
          ! Solid: activity = 1, so mu = Gref
          CompGibbs(i) = Gref(i)
        else
          ! Gas: mu = Gref + RT * ln(yi * fi * P)
          yi = N(i) / Ntot_gas
          CompGibbs(i) = Gref(i) + Rgas * T * log(yi * fug_coeff(i) * Ptot)
        end if
      end if

   EndDo


   End




  !----------------------------------------------------------------------------
   subroutine element_totals(N, b)
      Use gibbs_data
      Use gibbs_solver
      Real(dp), intent(in)  :: N(ncomp)
      Real(dp), intent(out) :: b(nelem)
      Integer :: k, i


      b = 0.0_dp
      do k = 1, nelem
         do i = 1, ncomp
            if (is_reactive(i)) b(k) = b(k) + atom_matrix(k,i) * N(i)
         end do
      end do
   end subroutine



  !----------------------------------------------------------------------------
  ! Solve Gibbs minimization using Newton-Raphson on the full KKT system
  ! with phase stability checking for solids AND trace gas removal.
  !
  ! PRIMARY VARIABLES: w_i = ln(N_i) for reactive species
  !                    lam_k = Lagrange multipliers
  !
  ! Species can be ACTIVE or INACTIVE:
  !
  ! ACTIVE GAS species i:
  !   f_i = Gref/(RT) + w_i - ln(Ntot_gas) + ln(fi*P) + sum_k(lam_k*a_ki/RT) = 0
  !   Contributes to Ntot_gas and element balances.
  !
  ! ACTIVE SOLID species i (phase present):
  !   f_i = Gref/(RT) + sum_k(lam_k*a_ki/RT) = 0
  !   Contributes to element balances (not Ntot_gas).
  !
  ! INACTIVE species i (gas trace or absent solid):
  !   f_i = w_i - w_min = 0   (force N ~ 0)
  !   Does NOT contribute to element balances or Ntot_gas.
  !
  ! A gas species is deactivated when it hits the w_min floor and its
  ! residual indicates it wants to go even lower (trace gas removal).
  ! After convergence, inactive species are checked for reactivation.
  !
  ! A solid is deactivated when N -> 0 and reactivated when its
  ! stability criterion Delta = Gref/(RT) + sum(lam*a/RT) < 0.
  !----------------------------------------------------------------------------
  subroutine solve_gibbs_min(N, Gref, fug_coeff, T, ierr, lam_out, active_out)
      Use gibbs_data
      Use gibbs_solver

    Real(dp), intent(inout) :: N(ncomp)
    Real(dp), intent(in)    :: Gref(ncomp), fug_coeff(ncomp),  T
    Integer, intent(out)    :: ierr
    Real(dp), intent(out)   :: lam_out(nelem)       ! Lagrange multipliers
    logical, intent(out)    :: active_out(ncomp)    ! species active flags

    Integer :: ndim
    Real(dp), allocatable :: Jac(:,:), F(:), dx(:), Jac_work(:,:), rhs(:)

    Real(dp) :: b_target(nelem), b_curr(nelem)
    Real(dp) :: w(ncomp)             ! w_i = ln(N_i) for reactive species
    Real(dp) :: lam(nelem)           ! Lagrange multipliers (J/mol)
    Real(dp) :: Ntot_gas, RT
    Real(dp) :: Gtot_curr
    Real(dp) :: max_F, max_F_old, damping
    Real(dp) :: max_dw
    Integer  :: iter, i, j, k, gi, gj, info
    Real(dp) :: ln_Ntot_gas
    Real(dp) :: delta_i              ! stability test value
    Real(dp) :: fi_check             ! residual check for reactivation

    ! Activity tracking (unified for gas and solid)
    logical  :: species_active(ncomp)  ! is this species active in the solve?
    logical  :: phase_changed          ! did we activate/deactivate anything?
    Integer  :: outer_iter
    Integer, parameter :: max_outer = 30

    ! Lower bound for ln(N): corresponds to ~1e-200
    Real(dp), parameter :: w_min = -460.0_dp
    ! Threshold for deactivating (in ln-space)
    Real(dp), parameter :: w_deactivate = -200.0_dp
    ! Residual threshold: if species is at w_min and residual is this positive,
    ! it means the species wants to go even lower -> deactivate
    Real(dp), parameter :: resid_deactivate = 10.0_dp

    ierr = 0
    RT = Rgas * T
    ndim = nreact + nelem

    !Legger inn en kommentar her også

    allocate(Jac(ndim, ndim), F(ndim), dx(ndim), Jac_work(ndim,ndim), rhs(ndim))

    ! Element totals from feed
    call element_totals(N, b_target)

    ! Initially all species are active
    species_active = .true.

    ! Deactivate species that contain elements with zero total atoms in feed.
    ! These species can never form and would only slow convergence.
    do i = 1, ncomp
      if (is_reactive(i)) then
        do k = 1, nelem
          if (atom_matrix(k,i) > 0.0_dp .and. abs(b_target(k)) < 1.0e-30_dp) then
            species_active(i) = .false.
            exit
          end if
        end do
      end if
    end do

    ! Smart initialization:
    ! Species with zero feed but whose elements exist should start inactive.
    ! The reactivation mechanism will bring them in if needed.
    ! Species with nonzero feed start active at their feed values.
    do i = 1, ncomp
      if (is_reactive(i)) then
        if (.not. species_active(i)) then
          ! Already deactivated (missing element)
          N(i) = 0.0_dp
          w(i) = w_min
        else if (N(i) < 1.0e-20_dp) then
          ! Zero or near-zero feed: start inactive, let reactivation handle it
          species_active(i) = .false.
          N(i) = 0.0_dp
          w(i) = w_min
        else
          w(i) = log(N(i))
        end if
      end if
    end do

    ! Initialize Lagrange multipliers to zero
    lam = 0.0_dp
    max_F_old = 1.0e30_dp
    damping = 0.0_dp

    write(*,'(A)') '================================================================='
    write(*,'(A)') ' Gibbs Energy Minimization - Full KKT Newton (log-space)'
    write(*,'(A)') ' with phase stability and trace species handling'
    write(*,'(A)') '================================================================='
    write(*,'(A,ES14.6,A)') ' Temperature  = ', T, ' K'
    write(*,'(A,ES14.6,A)') ' Pressure     = ', Ptot, ' bar'
    write(*,'(A,I4,A,I4)')  ' System size  = ', nreact, ' reactive + ', nelem, ' elements'
    write(*,'(A)') ''
    write(*,'(A)') ' Phase assignments:'
    do i = 1, ncomp
      if (is_reactive(i)) then
        if (CompType(i) == GAS) then
          write(*,'(A4,A12,A)') '   ', comp_name(i), '  Gas'
        else
          write(*,'(A4,A12,A)') '   ', comp_name(i), '  Solid'
        end if
      end if
    end do
    write(*,'(A)') '-----------------------------------------------------------------'
    write(*,'(A)') ' Iter    Max|F|          G_total (J)       Damping'
    write(*,'(A)') '-----------------------------------------------------------------'

    !=========================================================================
    ! Outer loop: Newton until convergence, then check species stability
    !=========================================================================
    do outer_iter = 1, max_outer

    do iter = 1, max_iter

      ! Reconstruct N from w
      do i = 1, ncomp
        if (is_reactive(i)) then
          N(i) = exp(max(w(i), w_min))
        end if
      end do

      ! Gas-phase total moles (only active gas species)
      Ntot_gas = 0.0_dp
      do i = 1, ncomp
        if (CompType(i) == GAS .and. species_active(i)) then
          Ntot_gas = Ntot_gas + N(i)
        end if
      end do
      ! Inert gases always contribute
      do i = 1, ncomp
        if (.not. is_reactive(i) .and. CompType(i) == GAS) then
          Ntot_gas = Ntot_gas + N(i)
        end if
      end do
      ln_Ntot_gas = log(Ntot_gas)

      !----------------------------------------------------------------------
      ! Build residual vector F
      !----------------------------------------------------------------------
      do i = 1, nreact
        gi = react_idx(i)

        if (.not. species_active(gi)) then
          ! INACTIVE: force w_i = w_min
          F(i) = w(gi) - w_min

        else if (CompType(gi) == SOLID) then
          ! ACTIVE SOLID
          F(i) = Gref(gi) / RT
          do k = 1, nelem
            F(i) = F(i) + lam(k) * atom_matrix(k, gi) / RT
          end do

        else
          ! ACTIVE GAS
          F(i) = Gref(gi)/RT + w(gi) - ln_Ntot_gas + log(fug_coeff(gi) * Ptot)
          do k = 1, nelem
            F(i) = F(i) + lam(k) * atom_matrix(k, gi) / RT
          end do
        end if
      end do

      ! Element balance constraints (only active species)
      b_curr = 0.0_dp
      do k = 1, nelem
        do i = 1, ncomp
          if (is_reactive(i) .and. species_active(i)) then
            b_curr(k) = b_curr(k) + atom_matrix(k,i) * N(i)
          end if
        end do
      end do

      do k = 1, nelem
        if (abs(b_target(k)) > 0.0_dp) then
          F(nreact + k) = (b_curr(k) - b_target(k)) / b_target(k)
        else
          ! No atoms of this element: replace with lam_k = 0
          F(nreact + k) = lam(k) / RT
        end if
      end do

      max_F = maxval(abs(F))
      Gtot_curr = total_gibbs(N, Gref, fug_coeff, T)

      if (mod(iter,50) == 1 .or. iter <= 20 .or. max_F < 1.0e-3_dp) then
        write(*,'(I5, ES16.6, ES20.8, ES12.3)') iter, max_F, Gtot_curr, damping
      end if

      ! Check convergence
      if (max_F < tol_conv) then
        write(*,'(I5, ES16.6, ES20.8)') iter, max_F, Gtot_curr
        exit  ! -> phase stability check
      end if

      !----------------------------------------------------------------------
      ! Build Jacobian
      !----------------------------------------------------------------------
      Jac = 0.0_dp

      do i = 1, nreact
        gi = react_idx(i)

        if (.not. species_active(gi)) then
          ! INACTIVE: equation is w_i - w_min = 0
          Jac(i, i) = 1.0_dp

        else if (CompType(gi) == SOLID) then
          ! ACTIVE SOLID: no w dependence (already zero)
          continue

        else
          ! ACTIVE GAS
          do j = 1, nreact
            gj = react_idx(j)
            if (.not. species_active(gj)) cycle
            if (CompType(gj) == GAS) then
              if (i == j) then
                Jac(i, j) = 1.0_dp - N(gj) / Ntot_gas
              else
                Jac(i, j) = -N(gj) / Ntot_gas
              end if
            end if
          end do
        end if
      end do

      ! df_i/dlam_k (skip inactive)
      do i = 1, nreact
        gi = react_idx(i)
        if (.not. species_active(gi)) cycle
        do k = 1, nelem
          Jac(i, nreact + k) = atom_matrix(k, gi) / RT
        end do
      end do

      ! dg_k/dw_j (skip inactive species and empty element rows)
      do k = 1, nelem
        if (abs(b_target(k)) < 1.0e-30_dp) then
          ! No atoms of this element: row is lam_k/RT = 0, d/dlam_k = 1/RT
          Jac(nreact + k, nreact + k) = 1.0_dp / RT
        else
          do j = 1, nreact
            gj = react_idx(j)
            if (.not. species_active(gj)) cycle
            Jac(nreact + k, j) = atom_matrix(k, gj) * N(gj) / b_target(k)
          end do
        end if
      end do

      !----------------------------------------------------------------------
      ! Solve: Jac * dx = -F
      !----------------------------------------------------------------------
      Jac_work = Jac
      rhs = -F
      call solve_linear(ndim, Jac_work, rhs, dx, info)

      if (info /= 0) then
        write(*,'(A,I5,A,I5)') ' WARNING: solve_linear singular at iter ', iter, ' info=', info
        do i = 1, nreact
          dx(i) = -F(i) * 0.1_dp
        end do
        do k = 1, nelem
          dx(nreact + k) = 0.0_dp
        end do
      end if

      !----------------------------------------------------------------------
      ! Damping
      !----------------------------------------------------------------------
      damping = 1.0_dp

      max_dw = 0.0_dp
      do i = 1, nreact
        if (abs(dx(i)) > max_dw) max_dw = abs(dx(i))
      end do
      if (max_dw > 5.0_dp) then
        damping = 5.0_dp / max_dw
      end if

      do k = 1, nelem
        if (abs(dx(nreact+k) * damping) > 0.5_dp * RT) then
          damping = min(damping, 0.5_dp * RT / abs(dx(nreact+k)))
        end if
      end do

      damping = max(damping, 1.0e-6_dp)

      ! Apply update
      do i = 1, nreact
        gi = react_idx(i)
        w(gi) = w(gi) + damping * dx(i)
        if (w(gi) < w_min) w(gi) = w_min
      end do
      do k = 1, nelem
        lam(k) = lam(k) + damping * dx(nreact + k)
      end do

      !----------------------------------------------------------------------
      ! Mid-iteration deactivation checks
      !----------------------------------------------------------------------
      do i = 1, nreact
        gi = react_idx(i)
        if (.not. species_active(gi)) cycle

        if (CompType(gi) == SOLID) then
          ! Solid: deactivate if N -> 0
          if (w(gi) < w_deactivate) then
            write(*,'(A,A12,A)') ' >> Deactivating solid: ', comp_name(gi), &
              ' (N -> 0 during iteration)'
            species_active(gi) = .false.
            w(gi) = w_min
          end if
        else
          ! Gas: deactivate if stuck at w_min with large positive residual
          ! (meaning the species wants to go even lower)
          if (w(gi) <= w_min + 1.0_dp) then
            ! Recompute this species' residual
            fi_check = Gref(gi)/RT + w(gi) - ln_Ntot_gas + log(fug_coeff(gi)*Ptot)
            do k = 1, nelem
              fi_check = fi_check + lam(k) * atom_matrix(k, gi) / RT
            end do
            if (fi_check > resid_deactivate) then
              write(*,'(A,A12,A,ES10.2,A)') ' >> Deactivating trace gas: ', &
                comp_name(gi), ' (residual=', fi_check, ', at w_min floor)'
              species_active(gi) = .false.
              w(gi) = w_min
            end if
          end if
        end if
      end do

      max_F_old = max_F

    end do  ! end inner Newton loop

    !========================================================================
    ! Phase stability / reactivation check
    !========================================================================
    if (iter > max_iter) then
      write(*,'(A)') ' WARNING: Newton did not converge within max iterations.'
      ierr = 1
      exit
    end if

    phase_changed = .false.

    ! Check all inactive species for possible reactivation.
    ! IMPORTANT: only reactivate ONE species per outer iteration to maintain
    ! solver stability. Pick the most favourable candidate.
    ! Variables to track best candidate:
    ! best_react_idx = reactive index of best candidate (-1 = none)
    ! best_react_val = stability/residual value (most negative wins)
    block
      Integer  :: best_ri
      Real(dp) :: best_val, check_val
      logical  :: best_is_solid

      best_ri  = -1
      best_val = 0.0_dp
      best_is_solid = .false.

      do i = 1, nreact
        gi = react_idx(i)
        if (species_active(gi)) cycle

        ! Skip species whose elements are absent from feed
        delta_i = 0.0_dp
        do k = 1, nelem
          if (atom_matrix(k, gi) > 0.0_dp .and. abs(b_target(k)) < 1.0e-30_dp) then
            delta_i = 1.0_dp
            exit
          end if
        end do
        if (delta_i > 0.5_dp) cycle

        if (CompType(gi) == SOLID) then
          check_val = Gref(gi) / RT
          do k = 1, nelem
            check_val = check_val + lam(k) * atom_matrix(k, gi) / RT
          end do
          write(*,'(A,A12,A,ES12.4,A)') '    Check solid: ', comp_name(gi), &
            '  (stability = ', check_val, ')'
        else
          check_val = Gref(gi)/RT + w_min - ln_Ntot_gas + log(fug_coeff(gi)*Ptot)
          do k = 1, nelem
            check_val = check_val + lam(k) * atom_matrix(k, gi) / RT
          end do
          write(*,'(A,A12,A,ES12.4,A)') '    Check gas:   ', comp_name(gi), &
            '  (residual = ', check_val, ')'
        end if

        ! Track most negative (most favourable for reactivation)
        if (check_val < best_val) then
          best_val = check_val
          best_ri  = i
          best_is_solid = (CompType(gi) == SOLID)
        end if
      end do

      ! Reactivate the single best candidate (if any)
      if (best_ri > 0 .and. best_val < -tol_conv) then
        gi = react_idx(best_ri)
        if (best_is_solid) then
          write(*,'(A,A12,A,ES12.4,A)') ' >> Activating solid: ', comp_name(gi), &
            '  (stability = ', best_val, ')'
          species_active(gi) = .true.
          w(gi) = log(1.0e-5_dp)
        else
          write(*,'(A,A12,A,ES12.4,A)') ' >> Reactivating gas: ', comp_name(gi), &
            '  (residual = ', best_val, ')'
          species_active(gi) = .true.
          w(gi) = w_min + 10.0_dp
        end if
        phase_changed = .true.
      end if
    end block

    if (.not. phase_changed) then
      write(*,'(A)') '-----------------------------------------------------------------'
      write(*,'(A)') ' Converged (all stability checks passed).'
      write(*,'(A)') ''

      ! Detailed stability analysis for inactive solid phases
      write(*,'(A)') ' Solid phase stability analysis:'
      write(*,'(A)') ' (Delta < 0 means solid should form; Delta > 0 means not stable)'
      write(*,'(A)') ''


      do i = 1, nreact
        gi = react_idx(i)
        if (CompType(gi) /= SOLID) cycle

        ! Compute stability
        delta_i = Gref(gi) / RT
        fi_check = 0.0_dp  ! sum of lambda*a terms
        do k = 1, nelem
          delta_i  = delta_i + lam(k) * atom_matrix(k, gi) / RT
          fi_check = fi_check + lam(k) * atom_matrix(k, gi)
        end do

        write(*,'(A,A12)') '   Species: ', comp_name(gi)
        if (species_active(gi)) then
          write(*,'(A,ES14.6)') '     Status:         ACTIVE (present), moles = ', N(gi)
        else
          write(*,'(A)')        '     Status:         INACTIVE (not present)'
        end if
        write(*,'(A,F12.4)')    '     Delta (stab.):  ', delta_i
        write(*,'(A,F14.2,A)')  '     Current Gref:   ', Gref(gi), ' J/mol'
        ! Critical Gref = -sum(lam_k * a_ki) = -fi_check
        write(*,'(A,F14.2,A)')  '     Critical Gref:  ', -fi_check, ' J/mol  (Gref needed for Delta=0)'
        write(*,'(A,F14.2,A)')  '     Gap:            ', Gref(gi) - (-fi_check), ' J/mol  (current - critical)'
        write(*,'(A)')          ''
      end do

      write(*,'(A)') '================================================================='

      ! Reconstruct final N
      do i = 1, ncomp
        if (is_reactive(i)) then
          if (species_active(i)) then
            N(i) = exp(max(w(i), w_min))
          else
            N(i) = 0.0_dp
          end if
        end if
      end do

      lam_out = lam
      active_out = species_active
      deallocate(Jac, F, dx, Jac_work, rhs)
      return
    end if

    write(*,'(A)') ' >> Species set changed, re-solving...'
    write(*,'(A)') '-----------------------------------------------------------------'

    end do  ! end outer loop

    ! If outer loop exhausted
    do i = 1, ncomp
      if (is_reactive(i)) then
        if (species_active(i)) then
          N(i) = exp(max(w(i), w_min))
        else
          N(i) = 0.0_dp
        end if
      end if
    end do

    if (ierr == 0) then
      write(*,'(A)') ' WARNING: Species stability loop did not settle.'
      ierr = 2
    end if
    lam_out = lam
    active_out = species_active
    deallocate(Jac, F, dx, Jac_work, rhs)

  end subroutine solve_gibbs_min



