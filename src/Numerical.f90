!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Numerical.f90
!  DESCRIPTION :  Contains various general numerical routines
!  CREATED     :  2015.04.10  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  Newton   - General Newton Raphson iteration
!                 pHNewton - Special version of Newton routine for pH iteration using alkalinity
!
!****************************************************************************






   SubRoutine pHNewton(funcd,x1,x2,xacc,root)
!****************************************************************************
!  SUBROUTINE  :  pHNewton
!  DESCRIPTION :  Special version of Newton routine for pH iteration using alkalinity
!  VARIABLES   :  funcd - External function supplying function value and first derivative
!                 xL,xH - Lower and upper limit (braket) for root
!                 FL,FH - Value of F corresponding to xL and xH. FL will be negative, FH will be positive
!                 root  - The solution!
!
!  CREATED     :  2015.04.10  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Implicit None

   External funcd

! Call variables
   Double Precision x1,x2,xacc,root

! Parameters
   Integer, Parameter :: Maxit=100
   Double Precision, Parameter :: MaxStep=2.d0

! Local variables
   Integer j
   Double Precision F,dF,dx,dxold,xL,xH
!-Start----------------------------------------------------------------------
   xL=x1
   xH=x2

   If (root.LE.xL .OR. root.GE.xH) root=(xL+xH)/2.d0  ! Unless root is within bracketed limits - give it a start estimate

   dxold=xH-xL
   dx=dxold

   Call funcd(root,F,dF)

   Do j=1,Maxit
      If (((root-XH)*dF-F)*((root-XL)*dF-F).GE.0.d0 .OR. dAbs(2.d0*F).GT.dAbs(dxold*dF)) Then   ! Check if step is outside limits or progress is to slow

         dxold=dx

         If (F.GT.0.d0) Then              ! Make sure that bisection step goes in correct direction
            dx=-0.5d0*(root-XL)
         Else
            dx=0.5d0*(XH-root)
         EndIf
      Else
         dxold=dx
         dx=-F/dF
      EndIf

      If (dabs(dx).GT.MaxStep) Then       ! Reduce step if it is too large
         dx=dx*MaxStep/dabs(dx)
      EndIf

      root=root+dx                        ! Calculate new solution

      If (dAbs(dx).LT.xacc) Return     ! Found a solution that is good enough

      CALL Funcd(root,F,dF)            ! Calculate new solution


      If(F.LT.0.d0) Then               ! Set corrected upper/lower limit
         XL=root
      ElseIf (F.GT.0.d0) Then
         XH=root
      Else
         Return                        ! BKaa: Found the solution - very seldom
      EndIf

   EndDo

   End


      SubRoutine SolveMat(a,b,n,np)
!**********************************************************************
!  Name:       SolveMat
!----------------------------------------------------------------------
!  Descr:      Solves the matrix of SR-derivatives
!  Called by:
!  Calls:
!----------------------------------------------------------------------
!  Variables:
!
!**********************************************************************
      Use GasRx
      Implicit None

! Call variables
      Integer np
      Double Precision a(np,np),b(np),d
      Integer n,indx(np),Code

 !-----Start------------------------------------------------------------

      CALL LUdecomp(a,n,np,indx,d,Code)
      CALL LUbacksub(a,n,np,indx,b)

   End SubRoutine SolveMat





   SubRoutine LUdecomp(a,n,np,indx,d,Code)
!****************************************************************************
!  SUBROUTINE  :  LUdecomp
!  DESCRIPTION :  General routine for LU decomposition of a matrix
!  VARIABLES   :  a  - matrix - will be changed when returned
!                 n  - Actual size of matrix a
!
!  CREATED     :  2015.04.16  Baard Kaasa www.scalesim.no From Numerical recipes
!  CHANGED     :
!
!****************************************************************************
   Implicit None

! Call variables
   Integer n,np,indx(n),Code
   Double Precision d,a(np,np)

! Local parameters
   Integer, Parameter :: nMax=50
   Double Precision, Parameter :: Tiny=1.0d-20

! Local variables
   Integer i,imax,j,k
   Double Precision aamax,dum,sum,vv(nMAX)

!--Start---------------------------------------------------------------
   d=1.d0 ; Code=0

      Do i=1,n
         aamax=0.d0
         Do j=1,n
            If (DABS(a(i,j)).GT.aamax) aamax=DABS(a(i,j))
         EndDo

         If (aamax.EQ.0.d0) Then
            Code=1
    !        Return
             Write(*,*) 'singular matrix in ludcmp'
         EndIf
         vv(i)=1.d0/aamax
      EndDo

      Do j=1,n
         Do i=1,j-1
            sum=a(i,j)
            Do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            EndDo
            a(i,j)=sum
         EndDo
         aamax=0.d0
         Do i=j,n
            sum=a(i,j)
            Do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            EndDo
            a(i,j)=sum
            dum=vv(i)*DABS(sum)
            If (dum.GE.aamax) THEN
               imax=i
               aamax=dum
            EndIf
         EndDo
         If (j.NE.imax) THEN
            Do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            EndDo
            d=-d
            vv(imax)=vv(j)
         EndIf
         indx(j)=imax
         If (a(j,j).EQ.0.d0) a(j,j)=TINY
         If (j.NE.n) THEN
            dum=1.d0/a(j,j)
            Do i=j+1,n
               a(i,j)=a(i,j)*dum
            EndDo
         EndIf
      EndDo


   End

   SubRoutine LUbacksub(a,n,np,indx,b)
!****************************************************************************
!  SUBROUTINE  :  LUbacksub
!  DESCRIPTION :  General routine for LU backsubsitution to solve matric
!  VARIABLES   :  a  - Modified matrix from LUDecomp
!                 b  - Input vector - will contain solution after LUbacksub
!
!  CREATED     :  2015.04.16  Baard Kaasa www.scalesim.no From Numerical recipes
!  CHANGED     :
!
!****************************************************************************
   Implicit None

! Call variables
   Integer n,np,indx(n)
   Double Precision a(np,np),b(np)

! Local variables
   Integer i,ii,j,ll
   Double Precision sum
!--Start---------------------------------------------------------------
      ii=0
      Do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         If (ii.NE.0) THEN
            Do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            EndDo
         ElseIf (sum.NE.0.d0) THEN
            ii=i
         EndIf
         b(i)=sum
      EndDo

      Do i=n,1,-1
         sum=b(i)
         Do j=i+1,n
            sum=sum-a(i,j)*b(j)
         EndDo
         b(i)=sum/a(i,i)
      EndDo

   End

   SubRoutine CubEqSolve(B,C,D,X,j)
!****************************************************************************
!  SUBROUTINE  :  CubRootSolve
!  DESCRIPTION :  Solves cubic equation of type X^3 + BX^2 + CX + D = 0
!                 Uses Cardano's Goniometric Finite method
!  VARIABLES   :  B,C,D - constants in cubic equation
!                 x1,x2,x3 - roods
!                 j=1 Only one Real root X1 (X2=X3=0)
!                 J=3 Three Real roots X1 <= X2 <= X3
!
!  CREATED     :  2016.06.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!
!****************************************************************************
   Implicit None

! Call variables
   Double Precision B,C,D
   Double Precision X(3)
   Integer j

   Double Precision CubRt

! Local variables
   Double Precision P,Q,R,Y1,Y2,Y3,R1
!--Start------------------------------
   X=0.d0 ; J=3

   P=(3.D0*C-B**2)/9.D0
   Q=B**3/27.D0-B*C/6.D0+D/2.D0
   R=Q**2 + P**3

   If (R.GT.0.d0) Then                 ! Only one Real solutions J=1
      J=1
      R1=dsqrt(R)
      Y1=CubRt(-Q+R1) + CubRt(-q-R1)
      Goto 20
   EndIf

   If (R.EQ.0.d0) Then                 ! Three roots - two of them equal - very seldom
      Y2=CubRt(Q)
      Y1=-2.d0*Y2
      Y3=Y2
      If (Q.LT.0.d0) Call SwapNum(Y1,Y3)
      Goto 20
   EndIf

   Y3=2.d0*dsqrt(-P)*dcos(dacos(-Q/dsqrt(-P**3))/3.d0)   ! Three different Real roots
   R1=dsqrt(Y3**2 - 12.d0*P - 4.d0*Y3**2)
   Y2=(-Y3 + R1)/2.d0
   Y1=(-Y3 - R1)/2.d0

20 X(1)=Y1 - B/3.d0                      ! First root is calculated
   If (J.EQ.3) Then                    ! Other roots calculated if Real
      X(2)=Y2 - B/3.d0
      X(3)=Y3 - B/3.d0
   EndIf

   End


   Double Precision Function CubRt(X)
!****************************************************************************
!  SUBROUTINE  :  3rdRoot
!  DESCRIPTION :  Calculates 3rd root of a Real number
!
!  VARIABLES   :  X Input value - Real number
!
!  CREATED     :  2016.06.28  Baard Kaasa www.scalesim.no
!  CHANGED     :  2017.09.29  BKaa - added check for zero
!
!****************************************************************************
   Implicit None

! Call variables
   Double Precision X
! Local variables
   Double Precision abX,Sgn
!--Sart----------------------------------------------------------------------
   If (X.NE.0.d0) Then
      abX=dabs(X)
      Sgn=X/abX

      CubRT=Sgn*abX**(1.d0/3.d0)
   Else
      CubRT=0.d0
   EndIf
   End

   Subroutine SwapNum(X,Y)
!****************************************************************************
!  SUBROUTINE  :  SwapNum
!  DESCRIPTION :  Swaps value of X and Y
!
!  VARIABLES   :  X, Y - Input values to be swapped
!
!  CREATED     :  2016.06.28  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
   Implicit None
   Double Precision X,Y,Help
!--Start--------------------------------------------------------------------
   Help=X
   X=Y
   Y=Help

   End


