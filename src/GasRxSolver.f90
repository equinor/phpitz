





   !======================================================================
! Module: ReactEquil
! Solve gas-phase reaction equilibrium using reaction extents.
!
! With extreme K values (lnK=40..660), many reactions go to completion.
! The solver:
!   1. Drives completion reactions first
!   2. Uses coupled Newton-Raphson for remaining active reactions
!      with proper line search on sum-of-squared-residuals
!======================================================================
module ReactEquil
   Use Globvar
   Use par
   Use GasRx
   Use HCprop
   Use ResHandling   ! UserChoices(0:10)
  implicit none

  integer, parameter :: MAXCOMP = 20, mxnRx = 20

!  integer :: nComp
  character(len=20) :: compNameRx(MAXCOMP)
  real(dp) :: NGasRx(MAXCOMP), NGas0(MAXCOMP), fugCoefRx(MAXCOMP)

 ! integer  :: nRx
 ! real(dp) :: stoich(MAXCOMP, mxnRx)
  real(dp) :: lnK(mxnRx)
  logical  :: rxUserDisabled(mxnRx)
  ! Staged-inclusion flag: when staging reactions in ascending |lnK| order
  ! during Phase 1, reactions not yet introduced have rxActive=.false. and
  ! are held out of the active set. Phase 2 sets all enabled reactions
  ! rxActive=.true. so the full system converges together.
  logical  :: rxActive(mxnRx)
  ! Persistent flag tracking whether a reaction was at-completion the
  ! last time buildActiveList examined it. Used for hysteresis: once
  ! a reaction goes at-completion, it stays at-completion until its
  ! limiting species climbs to yResumeFactor * yMinPract * Ntot.
  logical  :: rxAtCompletion(mxnRx)
  integer  :: rxStatus(mxnRx)  ! 0=active,1=converged,2=completion,3=disabled,4=no Q

  real(dp) ::  NtotRx

  ! Strict convergence tolerance (log-space residual). 1e-6 in lnQ-lnK
  ! means Q/K matches within about 1 ppm - plenty tight for practical work.
  real(dp) :: tolConv     = 1.0e-6_dp

  ! Absolute-moles floor, used by initReactions / takeStep to avoid
  ! genuinely negative populations. Kept very small - this is a numerical
  ! floor, not a physical one.
  real(dp) :: minConc     = 1.0e-30_dp

  ! Practical-zero threshold on MOLE FRACTION. A reaction that still
  ! wants to deplete a species with y_i < yMinPract is declared "at
  ! completion" - going further has no physical meaning because the
  ! species is already effectively absent. With a CO2 diluent of ~1e6
  ! moles, yMinPract=1e-10 corresponds to ~1e-4 absolute moles.
  real(dp) :: yMinPract   = 1.0e-10_dp

  ! Hysteresis factor for at-completion / re-activation. A reaction
  ! that has been marked at-completion only re-activates when its
  ! limiting species climbs to yResumeFactor * yMinPract. With a
  ! wide band (e.g. 100) the active set is much less prone to
  ! toggling when GS sweeps or LM-damped Newton transiently push a
  ! borderline species across the deactivation threshold.
  real(dp) :: yResumeFactor = 100.0_dp

  ! Reporting threshold (mole fraction) for the final status table.
  ! Reactions whose limiting species sits below yReportZero are
  ! flagged AT COMPLETION even if they didn't formally hit
  ! rxAtCompletion during iteration. With yMinPract=1e-10 and
  ! a 10x factor, any species below ~1 ppb is "practically gone"
  ! and the reaction is at completion for engineering purposes.
  real(dp) :: yReportZero = 1.0e-9_dp

  ! Practical convergence tolerance. If Newton stalls but all active
  ! residuals are below this, accept the solution rather than loop
  ! forever (log-space: exp(1e-2) ~ 1.01, i.e. 1% match in Q/K).
  real(dp) :: pracTol     = 1.0e-2_dp

  real(dp) :: initFrac    = 1.0e-4_dp
  integer  :: maxIter     = 10000 !5000
  real(dp) :: maxStepFrac = 0.3_dp  ! 0.3dp
  ! Controls all diagnostic output from the solver. Set from
  ! UserChoices(1) in setupProblem (0 = silent, 1 = print everything).
  logical  :: PrintDetails     = .false.

   contains



!======================================================================
! Get all parameters and define reactions
!======================================================================
   subroutine setupProblem(AllComp)
      Real(dP) :: AllComp(ncomp)
      Integer i

      ! Diagnostic-output flag from input file: UserChoices(1)=1 -> print
      ! everything (default), UserChoices(1)=0 -> stay silent.
      PrintDetails = (UserChoices(1) == 1)

      compNameRx(1:ncomp)=comp_name

      NGas0 = 0.0_dp
      NGas0(1:ncomp)=AllComp
      NGasRx = NGas0
      NtotRx = sum(NGas0(1:nComp))


      fugCoefRx = 1.0_dp
      fugCoefRx(1:ncomp)=Fugy

      lnK(1:maxnRx)=RxlnK

     rxUserDisabled = .false.
    Do i=1,nRx
       If (IncludeRx(i).EQ.0) Then
          rxUserDisabled(i)=.True.
      EndIf
    EndDo

    ! All non-user-disabled reactions start active. Phase 1 will gate
    ! some off temporarily for staged inclusion, then turn them all back on.
    Do i=1,nRx
       rxActive(i) = .not. rxUserDisabled(i)
    EndDo

    ! Hysteresis state: nothing at completion at the start.
    rxAtCompletion = .false.

    rxStatus = 0
  end subroutine setupProblem





!======================================================================
! Calculate Q values
!======================================================================
  function calcLnQ(j, lnQval) result(ok)
    integer, intent(in) :: j
    real(dp), intent(out) :: lnQval
    logical :: ok
    integer :: i
    real(dp) :: fi
!--Start---------------------------------------------------------------
    ok = .true.; lnQval = 0.0_dp
    do i = 1, nComp
      if (stoich(i,j) /= 0.0_dp) then
        if (NGasRx(i) <= 0.0_dp) then; ok=.false.; return; end if
           If (CompType(i).EQ.1) Then
               fi = (NGasRx(i)/NtotRx)*fugCoefRx(i)*Ptot
           Else
               fi=1.0_dp
            EndIf
        lnQval = lnQval + stoich(i,j)*log(fi)
      end if
    end do
    end function calcLnQ


!======================================================================
! Take a step in L direction
!======================================================================
  subroutine takeStep(j, L)
    integer, intent(in) :: j
    real(dp), intent(in) :: L
    integer :: i
!--Start---------------------------------------------------------------
    do i = 1, nComp
      if (stoich(i,j)/=0.0_dp) then
        NGasRx(i) = NGasRx(i) + stoich(i,j)*L
        if (NGasRx(i) < 0.0_dp) NGasRx(i) = 0.0_dp
      end if
    end do
  end subroutine takeStep


!======================================================================
! Calculate maximum step leght either forward or backward for each reaction
!======================================================================
  function maxAllowedStep(j, dir) result(mxL)
    integer, intent(in) :: j
    real(dp), intent(in) :: dir
    real(dp) :: mxL, c
    integer :: i
!--Start---------------------------------------------------------------
    mxL = huge(1.0_dp)
    do i = 1, nComp
      c = stoich(i,j)*dir
      if (c < 0.0_dp) then
        if (NGasRx(i)>0.0_dp) then
          mxL = min(mxL, NGasRx(i)/abs(c))
        else
          mxL = 0.0_dp; return
        end if
      end if
    end do
  end function maxAllowedStep

!======================================================================
! Take a step in L direction
!======================================================================
  subroutine initReactions()
    integer :: pass, j, i
    real(dp) :: minR, minP, step
    logical :: needF, needB, canF, canB
!--Start---------------------------------------------------------------
    do pass = 1, nRx*3
      do j = 1, nRx
        if (rxUserDisabled(j)) cycle
        minR=huge(1.0_dp); minP=huge(1.0_dp)
        needF=.false.; needB=.false.; canF=.true.; canB=.true.
        do i=1,nComp
          if (stoich(i,j)<0.0_dp) then
            if (NGasRx(i)<=minConc) then; needB=.true.; canF=.false.
            else; minR=min(minR,NGasRx(i)/abs(stoich(i,j))); end if
          else if (stoich(i,j)>0.0_dp) then
            if (NGasRx(i)<=minConc) then; needF=.true.; canB=.false.
            else; minP=min(minP,NGasRx(i)/abs(stoich(i,j))); end if
          end if
        end do
        if (needF.and.needB) cycle
        if (needF.and.canF) then
          step=minR*initFrac; call takeStep(j,step)
        else if (needB.and.canB) then
          step=minP*initFrac; call takeStep(j,-step)
        end if
      end do
    end do
    NtotRx = sum(NGasRx(1:nComp))
      end subroutine initReactions

!======================================================================
! Move reactions with extreem K values to one end
!======================================================================
  subroutine driveToCompletion(j, frac)
    integer, intent(in) :: j
    real(dp), intent(in) :: frac
    integer :: it, stall
    real(dp) :: lnQv, res, dir, mxL, resOld
    logical :: ok
!--Start---------------------------------------------------------------
    ! Geometric decay: after N steps the limiting species is at (1-frac)^N
    ! of its start value. With frac=0.9 we hit 1e-10 in ~10 steps.
    resOld = huge(1.0_dp)
    stall  = 0
    do it = 1, 2000
      ok = calcLnQ(j, lnQv); if (.not.ok) return
      res = lnQv - lnK(j)
      if (abs(res) < tolConv) return
      if (res<0.0_dp) then; dir=1.0_dp; else; dir=-1.0_dp; end if
      mxL = maxAllowedStep(j, dir)
      if (mxL <= 0.0_dp) return
      call takeStep(j, dir*frac*mxL)
      NtotRx = sum(NGasRx(1:nComp))

      ! Stall detection: if residual has stopped moving meaningfully,
      ! we've hit the practical limit (reactant already at ~minConc).
      if (abs(res - resOld) < 1.0e-4_dp * max(abs(res), 1.0_dp)) then
        stall = stall + 1
        if (stall > 3) return
      else
        stall = 0
      end if
      resOld = res
    end do
  end subroutine driveToCompletion

  subroutine solveLinear(A, b, x, n, singular)
!======================================================================
! Solve set of linear equations
!======================================================================
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(mxnRx,mxnRx), b(mxnRx)
    real(dp), intent(out) :: x(mxnRx)
    logical, intent(out) :: singular
    integer :: i, j, k, pr
    real(dp) :: pv, fac, tmp, AA(mxnRx,mxnRx), bb(mxnRx)
!--Start---------------------------------------------------------------
    singular=.false.; x=0.0_dp
    AA(1:n,1:n)=A(1:n,1:n); bb(1:n)=b(1:n)
    do k=1,n
      pr=k; pv=abs(AA(k,k))
      do i=k+1,n
        if(abs(AA(i,k))>pv) then; pv=abs(AA(i,k)); pr=i; end if
      end do
      if(pv<1e-30_dp) then; singular=.true.; return; end if
      if(pr/=k) then
        do j=1,n; tmp=AA(k,j); AA(k,j)=AA(pr,j); AA(pr,j)=tmp; end do
        tmp=bb(k); bb(k)=bb(pr); bb(pr)=tmp
      end if
      do i=k+1,n
        fac=AA(i,k)/AA(k,k)
        do j=k+1,n; AA(i,j)=AA(i,j)-fac*AA(k,j); end do
        bb(i)=bb(i)-fac*bb(k)
      end do
    end do
    do i=n,1,-1
      x(i)=bb(i)
      do j=i+1,n; x(i)=x(i)-AA(i,j)*x(j); end do
      x(i)=x(i)/AA(i,i)
    end do
    end subroutine solveLinear

!======================================================================
! Check which reactions that are still active or if rx/prod are below minimum
!======================================================================
  subroutine buildActiveList(nAct, actList, lnQval, res)
    integer, intent(out) :: nAct, actList(mxnRx)
    real(dp), intent(out) :: lnQval(mxnRx), res(mxnRx)
    integer :: j, i
    logical :: ok, canGo
    real(dp) :: nThreshDeact, nThreshReact, nThreshUse
!--Start---------------------------------------------------------------
    ! Two-level (hysteresis) threshold. A reaction newly going
    ! at-completion uses the strict threshold (nThreshDeact); a reaction
    ! that was already at completion needs its limiting species above
    ! the looser threshold (nThreshReact = yResumeFactor * nThreshDeact)
    ! before it re-activates. This stops the active set from toggling
    ! when species sit right at the practical-zero boundary.
    nThreshDeact = yMinPract * max(NtotRx, 1.0_dp)
    nThreshReact = yResumeFactor * nThreshDeact

    nAct=0; actList=0; lnQval=0.0_dp; res=0.0_dp
    do j=1,nRx
      if(rxUserDisabled(j)) then; rxStatus(j)=3; cycle; end if
      ! Reactions not yet staged in (Phase 1 homotopy) are held out
      if(.not. rxActive(j)) then; rxStatus(j)=3; cycle; end if
      ok=calcLnQ(j,lnQval(j))
      if(.not.ok) then; rxStatus(j)=4; rxAtCompletion(j)=.true.; cycle; end if
      res(j) = lnQval(j) - lnK(j)

      ! Pick threshold based on previous state (hysteresis)
      if (rxAtCompletion(j)) then
        nThreshUse = nThreshReact   ! must climb above resume threshold to re-activate
      else
        nThreshUse = nThreshDeact   ! drop below strict threshold to mark at completion
      end if

      if (res(j) < -tolConv) then
        ! Wants forward: blocked if any reactant (stoich<0) is depleted
        canGo=.true.
        do i=1,nComp
          if(stoich(i,j)<0.0_dp .and. NGasRx(i)<nThreshUse) then; canGo=.false.; exit; end if
        end do
        if(.not.canGo) then; rxStatus(j)=2; rxAtCompletion(j)=.true.; cycle; end if
      end if
      if (res(j) > tolConv) then
        ! Wants backward: blocked if any product (stoich>0) is depleted
        canGo=.true.
        do i=1,nComp
          if(stoich(i,j)>0.0_dp .and. NGasRx(i)<nThreshUse) then; canGo=.false.; exit; end if
        end do
        if(.not.canGo) then; rxStatus(j)=2; rxAtCompletion(j)=.true.; cycle; end if
      end if
      rxStatus(j) = 0
      rxAtCompletion(j) = .false.
      nAct=nAct+1; actList(nAct)=j
    end do
      end subroutine buildActiveList


!======================================================================
! Evaluate SSR for a trial set of steps
!======================================================================
  function trialSSR(NGasSave, nAct, actList, deltaL, alpha) result(ssr)
    real(dp), intent(in) :: NGasSave(MAXCOMP), deltaL(mxnRx), alpha
    integer, intent(in) :: nAct, actList(mxnRx)
    real(dp) :: ssr, lnQv
    integer :: ja, j, i
    logical :: ok
!--Start---------------------------------------------------------------
    ! Set NGasRx to saved + alpha*step
    NGasRx = NGasSave
    do ja=1,nAct
      j=actList(ja)
      do i=1,nComp
        if(stoich(i,j)/=0.0_dp) then
          NGasRx(i) = NGasRx(i) + stoich(i,j)*alpha*deltaL(ja)
        end if
      end do
    end do
    do i=1,nComp; if(NGasRx(i)<0.0_dp) NGasRx(i)=0.0_dp; end do
    NtotRx=sum(NGasRx(1:nComp))
    ssr=0.0_dp
    do ja=1,nAct
      j=actList(ja); ok=calcLnQ(j,lnQv)
      if(ok) then; ssr=ssr+(lnQv-lnK(j))**2
      else; ssr=ssr+1e30_dp; end if
    end do
    end function trialSSR



  !--------------------------------------------------------------------
  ! MAIN SOLVER
  !--------------------------------------------------------------------
  subroutine solveEquilibrium(ierr)
    integer, intent(out) :: ierr

    integer :: iter, j, i, ja, ka, sortIdx
    integer :: nAct, actList(mxnRx), rxOrder(mxnRx)
    real(dp) :: lnQval(mxnRx), res(mxnRx)
    real(dp) :: dQQ(mxnRx,mxnRx), deltaL(mxnRx), b_rhs(mxnRx)
    real(dp) :: maxRes, dir, mxL, sLim
    real(dp) :: crossDamp, NGasSave(MAXCOMP), ssr0, ssr1, alpha
    logical :: singular
    integer :: nFail
    real(dp) :: diagBoost
    ! Best-so-far tracking for robust recovery
    real(dp) :: bestRes, NGasBest(MAXCOMP), maxResPre
    integer  :: bestIter, stallIter, gsIt, gsBlock

    ierr = 1

    if (PrintDetails) then
      write(*,'(A)') '================================================'
      write(*,'(A)') ' Reaction Equilibrium Solver'
      write(*,'(A,F8.2,A)') ' Pressure = ', Ptot, ' bar'
      write(*,'(A)') '================================================'
    end if

    call initReactions()

    ! ============================================================
    ! Phase 1: Homotopy / staged reaction inclusion
    !
    ! Reactions are introduced one at a time in ascending |lnK| order.
    ! Reactions with lnK near zero equilibrate at "normal" mole fractions
    ! and behave well numerically; reactions with extreme |lnK| effectively
    ! act as constraints that drive one species toward zero.
    !
    ! By introducing the well-conditioned reactions first, the system
    ! finds a balanced state where moderate equilibria are satisfied.
    ! Adding extreme-K reactions LAST means that by the time they enter,
    ! their preferred reactants may already be partially consumed by
    ! more thermodynamically favorable pathways involving multiple
    ! moderate reactions in series. This avoids the local-minimum trap
    ! where a single high-|lnK| reaction grabs a shared reactant before
    ! competing pathways can use it.
    !
    ! Concrete example: 8 H2S + 4 O2 -> 8 H2O + S8 (lnK=748) versus
    ! the H2SO4 pathway H2S + 3 NO2 -> ... -> H2SO4 (combined lnK=234
    ! with TWO reactions, but per atom of H2S: 234 vs 748/8 = 93.5).
    ! The H2SO4 pathway wins thermodynamically per mole of S, but only
    ! if it's allowed to operate before Rx 4 consumes the H2S.
    ! ============================================================
    if (PrintDetails) write(*,'(A)') &
      ' Phase 1: Homotopy (stage in by ascending |lnK|)...'

    ! Sort enabled reactions by |lnK| ascending.
    do j=1,nRx; rxOrder(j)=j; end do
    do i=1,nRx-1
      do j=i+1,nRx
        if (abs(lnK(rxOrder(j))) < abs(lnK(rxOrder(i)))) then
          sortIdx=rxOrder(i); rxOrder(i)=rxOrder(j); rxOrder(j)=sortIdx
        end if
      end do
    end do

    ! Hold all reactions out, then stage them in one at a time.
    rxActive = .false.

    do i = 1, nRx
      j = rxOrder(i)
      if (rxUserDisabled(j)) cycle

      rxActive(j) = .true.

      if (PrintDetails) write(*,'(A,I2,A,I2,A,F10.3)') &
        '  Stage ', i, ': activate Rx ', j, '   lnK = ', lnK(j)

      ! Inner loop: Gauss-Seidel-Newton on the currently active subset.
      ! Each call to singleRxStep takes a bounded Newton step on one
      ! reaction; cycling over the whole active list converges the
      ! coupled system. We use a tight tolerance per stage so each
      ! incremental problem is solved well before the next reaction
      ! is added.
      do iter = 1, 500
        call buildActiveList(nAct, actList, lnQval, res)
        if (nAct == 0) exit
        maxRes = 0.0_dp
        do ja=1,nAct; maxRes=max(maxRes,abs(res(actList(ja)))); end do
        if (maxRes < 1.0e-3_dp) exit

        do ja = 1, nAct
          call singleRxStep(actList(ja))
        end do
      end do

      if (PrintDetails) &
        write(*,'(A,ES11.4,A,I4,A,I2)') &
          '           maxRes=', maxRes, '  iters=', iter, '  nAct=', nAct
    end do

    ! All enabled reactions back on for Phase 2 final convergence.
    do i=1,nRx
      rxActive(i) = .not. rxUserDisabled(i)
    end do

    ! ============================================================
    ! Phase 2: Coupled Newton-Raphson
    ! Key improvement: use Levenberg-Marquardt style diagonal boost
    ! when line search fails, and golden section for step size
    ! ============================================================
    if (PrintDetails) write(*,'(A)') ' Phase 2: Newton-Raphson with LM damping...'

    nFail = 0
    crossDamp = 1.0_dp
    diagBoost = 0.0_dp   ! Levenberg-Marquardt parameter
    bestRes   = huge(1.0_dp)
    NGasBest  = NGasRx
    bestIter  = 0
    stallIter = 0

    do iter = 1, maxIter
      call buildActiveList(nAct, actList, lnQval, res)
      if (nAct==0) then; ierr=0; exit; end if

      maxRes = 0.0_dp
      do ja=1,nAct; maxRes=max(maxRes,abs(res(actList(ja)))); end do

      ! Best-so-far tracking: if we're making things worse, we can revert.
      if (maxRes < bestRes) then
        bestRes   = maxRes
        NGasBest  = NGasRx
        bestIter  = iter
        stallIter = 0
      else
        stallIter = stallIter + 1
      end if

      if (PrintDetails .and. (mod(iter,500)==1 .or. iter<=10 .or. maxRes<tolConv*10)) &
        write(*,'(A,I6,A,ES11.4,A,I2,A,F5.2,A,ES8.1)') &
          '  NR',iter,'  maxRes=',maxRes,'  nAct=',nAct, &
          '  xD=',crossDamp,'  lm=',diagBoost

      if (maxRes < tolConv) then
        if (PrintDetails) write(*,'(A,I6)') ' *** Converged at iter ',iter
        ierr=0; exit
      end if

      ! Practical-convergence exit: if we're within pracTol and Newton
      ! has stopped making progress for a while, accept. The user's
      ! physical tolerance is much looser than tolConv anyway.
      if (bestRes < pracTol .and. stallIter > 300) then
        if (PrintDetails) &
          write(*,'(A,ES11.4,A,I6)') &
            ' *** Practical convergence (bestRes=',bestRes,') at iter ',iter
        NGasRx = NGasBest
        NtotRx = sum(NGasRx(1:nComp))
        ierr   = 0
        exit
      end if

      ! Hard-stall guard: no progress at all for a long stretch, no matter
      ! the residual level. Accept best-so-far so we don't spin forever.
      if (stallIter > 1500) then
        if (PrintDetails) &
          write(*,'(A,ES11.4,A,I6)') &
            ' *** Hard stall, accepting bestRes=',bestRes,' at iter ',iter
        NGasRx = NGasBest
        NtotRx = sum(NGasRx(1:nComp))
        if (bestRes < 1.0_dp) then
          ierr = 0
        else
          ierr = 1
        end if
        exit
      end if

      NGasSave=NGasRx
      ssr0 = 0.0_dp
      do ja=1,nAct; ssr0=ssr0+res(actList(ja))**2; end do

      ! Build Jacobian: dQQ(ja,ka) = d(lnQ_j)/d(L_k)
      dQQ = 0.0_dp
      do ja = 1, nAct
        do ka = 1, nAct
          do i = 1, nComp
            if (stoich(i,actList(ja))/=0.0_dp .and. NGasRx(i)>0.0_dp) &
              dQQ(ja,ka) = dQQ(ja,ka) + &
                stoich(i,actList(ja))*stoich(i,actList(ka))/NGasRx(i)
          end do
          if (ja/=ka) dQQ(ja,ka) = dQQ(ja,ka)*crossDamp
        end do
        ! Levenberg-Marquardt diagonal boost
        dQQ(ja,ja) = dQQ(ja,ja) + diagBoost
      end do

      ! RHS = -(lnQ - lnK)
      do ja=1,nAct; b_rhs(ja)=-res(actList(ja)); end do
      call solveLinear(dQQ, b_rhs, deltaL, nAct, singular)

      if (singular) then
        ! Pure diagonal fallback
        do ja=1,nAct
          j=actList(ja); dQQ(ja,ja)=0.0_dp
          do i=1,nComp
            if(stoich(i,j)/=0.0_dp .and. NGasRx(i)>0.0_dp) &
              dQQ(ja,ja)=dQQ(ja,ja)+stoich(i,j)**2/NGasRx(i)
          end do
          dQQ(ja,ja) = dQQ(ja,ja) + diagBoost
          if(abs(dQQ(ja,ja))>1e-30_dp) then
            deltaL(ja)=-res(j)/dQQ(ja,ja)
          else; deltaL(ja)=0.0_dp; end if
        end do
      end if

      ! Per-reaction step limiting (prevent going negative)
      do ja=1,nAct
        j=actList(ja)
        if(abs(deltaL(ja))<1e-30_dp) cycle
        dir=sign(1.0_dp,deltaL(ja))
        mxL=maxAllowedStep(j,dir)
        sLim=maxStepFrac*mxL
        if(abs(deltaL(ja))>sLim .and. sLim>0.0_dp) deltaL(ja)=sign(sLim,deltaL(ja))
      end do

      ! Backtracking line search
      alpha = 1.0_dp
      do i=1,40
        ssr1 = trialSSR(NGasSave, nAct, actList, deltaL, alpha)
        if (ssr1 < ssr0*(1.0_dp - 1.0e-4_dp*alpha)) exit  ! Armijo condition
        alpha = alpha * 0.5_dp
      end do

      ! Accept or reject
      if (ssr1 < ssr0) then
        ! Accept the step (NGasRx already set by trialSSR)
        nFail = 0
        ! Decrease LM parameter (more Newton-like)
        diagBoost = diagBoost * 0.5_dp
        if (diagBoost < 1.0e-12_dp) diagBoost = 0.0_dp
        ! Restore cross-derivatives if they were reduced
        if (crossDamp < 1.0_dp) crossDamp = min(1.0_dp, crossDamp + 0.1_dp)
      else
        ! Reject - restore and increase LM parameter
        NGasRx = NGasSave
        NtotRx = sum(NGasRx(1:nComp))
        nFail = nFail + 1

        ! Increase LM damping (more steepest-descent-like)
        if (diagBoost < 1.0e-10_dp) then
          diagBoost = 1.0_dp
        else
          diagBoost = diagBoost * 10.0_dp
        end if

        ! If many consecutive failures, reduce cross-derivatives
        if (nFail > 5) then
          crossDamp = crossDamp * 0.5_dp
          if (crossDamp < 0.01_dp) crossDamp = 0.0_dp
          nFail = 0
        end if

        ! If really stuck with huge LM, do bounded GS sweeps with
        ! progress-based termination. We limit total sweeps and require
        ! measurable residual reduction per block, otherwise we give up
        ! and let the outer loop hand back bestRes.
        if (diagBoost > 1.0e6_dp) then
          if (PrintDetails) write(*,'(A)',advance='no') '  GS sweeps... '
          diagBoost = 0.0_dp
          crossDamp = 1.0_dp

          call buildActiveList(nAct, actList, lnQval, res)
          maxResPre = 0.0_dp
          do ja=1,nAct; maxResPre=max(maxResPre,abs(res(actList(ja)))); end do

          gsBlock = 0
          gsIt    = 0
          do
            ! one block of 500 sweeps
            do gsIt = 1, 500
              do ja=1,nAct
                call singleRxStep(actList(ja))
              end do
            end do
            gsBlock = gsBlock + 1

            call buildActiveList(nAct, actList, lnQval, res)
            if (nAct==0) exit
            maxRes = 0.0_dp
            do ja=1,nAct; maxRes=max(maxRes,abs(res(actList(ja)))); end do
            if (maxRes < tolConv) exit

            ! Progress check: need to drop residual by at least 10% per block
            if (maxRes > 0.9_dp * maxResPre) exit
            maxResPre = maxRes

            ! Hard cap: at most 8 blocks (4000 sweeps)
            if (gsBlock >= 8) exit
          end do

          if (PrintDetails) write(*,'(A,ES11.4)') 'maxRes=',maxRes

          ! Update best-so-far if GS helped
          if (maxRes < bestRes) then
            bestRes  = maxRes
            NGasBest = NGasRx
            bestIter = iter
            stallIter = 0
          end if

          if (maxRes < tolConv) then
            ierr=0; exit
          end if

          ! If GS made no progress, step out and let stall-guard finish us
          if (maxRes > 0.99_dp * maxResPre) then
            if (PrintDetails) write(*,'(A)') '   GS stalled, reverting to best-so-far'
            NGasRx = NGasBest
            NtotRx = sum(NGasRx(1:nComp))
          end if
        end if
      end if
    end do

    ! If we fell out of the loop without an explicit exit, restore best.
    if (iter > maxIter .and. bestRes < huge(1.0_dp)) then
      NGasRx = NGasBest
      NtotRx = sum(NGasRx(1:nComp))
      if (bestRes < pracTol) then
        if (PrintDetails) &
          write(*,'(A,ES11.4)') ' *** Max iterations, accepting bestRes=', bestRes
        ierr = 0
      else if (bestRes < 1.0_dp) then
        if (PrintDetails) &
          write(*,'(A,ES11.4)') ' *** Max iterations, weak match bestRes=', bestRes
        ierr = 0
      else
        if (PrintDetails) &
          write(*,'(A,ES11.4)') ' *** Max iterations, no convergence. bestRes=', bestRes
        ierr = 1
      end if
    end if

    call printResults()
      end subroutine solveEquilibrium


!======================================================================
! Single-reaction Newton step (for GS fallback)
!======================================================================
  subroutine singleRxStep(j)
    integer, intent(in) :: j
    integer :: i
    real(dp) :: lnQv, r, diag, dL, dir, mxL, sLim
    logical :: ok
!--Start---------------------------------------------------------------
    ok=calcLnQ(j,lnQv); if(.not.ok) return
    r=lnQv-lnK(j)
    if(abs(r)<tolConv) return
    diag=0.0_dp
    do i=1,nComp
      if(stoich(i,j)/=0.0_dp .and. NGasRx(i)>0.0_dp) &
        diag=diag+stoich(i,j)**2/NGasRx(i)
    end do
    if(abs(diag)<1e-30_dp) return
    dL=-r/diag
    dir=sign(1.0_dp,dL)
    mxL=maxAllowedStep(j,dir)
    sLim=maxStepFrac*mxL
    if(abs(dL)>sLim .and. sLim>0.0_dp) dL=sign(sLim,dL)
    call takeStep(j,dL)
    NtotRx=sum(NGasRx(1:nComp))
  end subroutine singleRxStep


!======================================================================
! Print out results
!======================================================================
  subroutine printResults()
    integer :: i, j
    real(dp) :: lnQv, SumN, nThreshReport
    logical :: ok
!--Start---------------------------------------------------------------
    ! Reporting threshold (looser than the iteration's strict
    ! practical-zero threshold). A reaction whose limiting species
    ! is below this is flagged at completion even if it ended the
    ! iteration without rxAtCompletion=.true. This catches borderline
    ! species sitting just above the strict iteration threshold but
    ! physically essentially gone.
    nThreshReport = yReportZero * max(NtotRx, 1.0_dp)

    do j = 1, nRx
      if (rxUserDisabled(j)) then
        rxStatus(j) = 3
        cycle
      end if

      ok = calcLnQ(j, lnQv)
      if (.not.ok) then
        rxStatus(j) = 4
        cycle
      end if

      ! Iteration set rxAtCompletion: report as AT COMPLETION.
      if (rxAtCompletion(j)) then
        rxStatus(j) = 2
        cycle
      end if

      ! Strict residual converged: CONVERGED.
      if (abs(lnQv - lnK(j)) < tolConv) then
        rxStatus(j) = 1
        cycle
      end if

      ! Loose reporting: check whether the residual sign points to a
      ! species that is below the report threshold. If yes, the
      ! reaction is practically at completion regardless of the
      ! iteration's strict view.
      rxStatus(j) = 0
      if (lnQv < lnK(j)) then
        ! Wants forward: check reactants
        do i=1,nComp
          if (stoich(i,j)<0.0_dp .and. NGasRx(i)<nThreshReport) then
            rxStatus(j) = 2; exit
          end if
        end do
      else
        ! Wants backward: check products
        do i=1,nComp
          if (stoich(i,j)>0.0_dp .and. NGasRx(i)<nThreshReport) then
            rxStatus(j) = 2; exit
          end if
        end do
      end if
    end do

   ! if (PrintDetails) then
      write(*,'(A)') ''
      write(*,'(A)') '================================================================================='
      write(*,'(A)') ' CO2 PHASE REACTIONS '
      write(*,'(A)') '================================================================================='
      write(*,*) ' Component    Feed (mol)     Equil (mol)       Equil (ppm)'
      write(*,'(A)') '---------------------------------------------------------------------------------'
   ! end if

    SumN=0.0_dp
    Do i=1,ncomp
       If (CompType(i).EQ.1) Then
          SumN=SumN+NgasRx(i)
       EndIf
    EndDo

   ! if (PrintDetails) then
      do i=1,nComp
        write(*,'(A,A10,ES15.6,ES15.6,F15.4)')'  ',compNameRx(i), NGas0(i), NGasRx(i), NgasRx(i)*1.0E6_dp/sumN
      end do
      write(*,'(A)') '================================================================================='



      write(*,'(/,A)') ' Reaction status:'
      write(*,'(A)') '---------------------------------------------------------------------------------'
      do j=1,nRx
        ok = calcLnQ(j, lnQv)
        select case(rxStatus(j))
        case(1)
          write(*,'(A,I2,A,F12.3,A,F12.3,A)') &
            '  Rx',j,':  lnQ=',lnQv,'  lnK=',lnK(j),'  CONVERGED'
        case(2)
          write(*,'(A,I2,A,F12.3,A,F12.3,A)') &
            '  Rx',j,':  lnQ=',lnQv,'  lnK=',lnK(j), &
            '  AT COMPLETION (reactant depleted)'
        case(3)
          write(*,'(A,I2,A)') '  Rx',j,':  USER DISABLED'
        case(4)
          write(*,'(A,I2,A)') '  Rx',j,':  Cannot compute Q'
        case default
          if (ok) then
            write(*,'(A,I2,A,F12.3,A,F12.3,A,ES10.3)') &
              '  Rx',j,':  lnQ=',lnQv,'  lnK=',lnK(j), &
              '  NOT CONVERGED |res|=',abs(lnQv-lnK(j))
          else
            write(*,'(A,I2,A)') '  Rx',j,':  Cannot compute Q'
          end if
        end select
      end do
      write(*,'(A)') '================================================================================='
    ! end if
  end subroutine printResults

   end module ReactEquil


   SubRoutine GasRxSolver(AllComp)
   use ReactEquil
   Use GasRx
   implicit none
   integer :: ierr

   Real(dp) :: AllComp(ncomp)

   call setupProblem(AllComp)
  ! rxUserDisabled(4) = .true.   ! Example: disable a reaction

   call solveEquilibrium(ierr)

   AllComp(1:Ncomp)=NGasRx(1:Ncomp)

   if (PrintDetails) then
      if(ierr==0) then
         write(*,'(/,A)') 'All active reactions converged or at completion.'
      else
         write(*,'(/,A)') 'WARNING: Some reactions did not converge.'
      end if
   end if

end SubRoutine GasRxSolver


