!****************************************************************************
!                 FILE INFO
!
!  FILENAME    :  Pitzer_Funcs.f90
!  DESCRIPTION :  The various functions used by Pitzer model
!  CREATED     :  2013.05.11  Baard Kaasa  www.scalesim.no
!
!  ROUTINES    :  BMX      - Calculates B from Beta0, Beta1 and Beta2
!                 g, gP    - g function and its derivative
!                 E_Theta  - Calculates E_Theta for unsymetrical mixing
!                 E_ThetaP - Derivative of E_Theta
!                 J        - Calculates J for unsymetrical mixing
!
!
!
!****************************************************************************

      Double Precision Function BMX(M,X,p)
!****************************************************************************
!  SUBROUTINE  :  BMX
!  DESCRIPTION :  Calculates BMX parameter
!  VARIABLES   :  M and X are cation and anion indexes in MultiScale ref
!                 p indicates if it is BMX(p=0), derivative (p=1) or for osmotic (p=2)
!
!  CREATED     :  2013.05.11  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Call variables
      Integer M,X,p

! Functions
      Double Precision g,gp

! Local variables
      Double Precision x1,x2
!-----Start-------------------------------------------------------------------

      x1=PitzAlpha(1,M,X)*sqrIS
      x2=PitzAlpha(2,M,X)*sqrIS

      If (p.EQ.0) Then
         BMX=Beta0(M,X) + Beta1(M,X)*g(x1) + Beta2(M,X)*g(x2)
      ElseIf (p.EQ.1) Then
         BMX=(Beta1(M,X)*gp(x1) + Beta2(M,X)*gp(x2))/IS
      ElseIf (p.EQ.2) Then
         BMX=Beta0(M,X) + Beta1(M,X)*dexp(-x1) + Beta2(M,X)*dexp(-x2)
      EndIf

   End

      Double Precision Function g(x)
!****************************************************************************
!  SUBROUTINE  :  g
!  DESCRIPTION :  Calculates g
!  VARIABLES   :
!
!  CREATED     :  2013.05.011  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Double Precision x
!-----Start-------------------------------------------------------------------
      If (x.EQ.0.d0) Then
         g=0.d0
      Else
         g=2.d0*(1.d0 - (1.d0 + x)*dexp(-x))/x**2
      EndIf

      End

      Double Precision Function gP(x)
!****************************************************************************
!  SUBROUTINE  :  gP
!  DESCRIPTION :  Calculates derivative if g
!  VARIABLES   :
!
!  CREATED     :  2013.05.011  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Double Precision x
!-----Start-------------------------------------------------------------------
      If (x.EQ.0.d0) Then
         gP=0.d0
      Else
         gP=-2.d0*(1.d0 - (1.d0 + x + 0.5d0*x**2)*dexp(-x))/x**2
      EndIf
   End


      Double Precision Function CMX(M,X)
!****************************************************************************
!  SUBROUTINE  :  CMX
!  DESCRIPTION :  Calculates CMX
!  VARIABLES   :  M and X are cation and anion indexes in MultiScale ref
!
!  CREATED     :  2013.05.12  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Use Aqua
      Implicit None

! Call variables
      Integer M,X
! Local variables
      Double Precision ChargeProd
!-----Start-------------------------------------------------------------------

      ChargeProd=dble(AqCharge(M)*AqCharge(X))
      ChargeProd=dabs(ChargeProd)

      CMX=Cphi(M,X) / (  2.d0*dsqrt( ChargeProd ) )

   End


      Recursive Double Precision Function Phi(i,j,p)
!****************************************************************************
!  SUBROUTINE  :  Phi
!  DESCRIPTION :  Calculates Phi
!  VARIABLES   :  i and j are two ions of same sign in MultiScale numbering
!                 p indicates if it is Phi(p=0), derivative (p=1) or for osmotic (p=2)
!
!  CREATED     :  2013.05.11  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Call variables
      Integer i,j,p

! Functions
      Double Precision E_Theta,E_ThetaP

! Local variables
      Double Precision alpha1,alpha2,x1,x2
      Integer zi,zj
!-----Start-------------------------------------------------------------------

      zi=AqCharge(i)
      zj=AqCharge(j)

      If (p.EQ.0) Then
         Phi=Theta(i,j) + E_Theta(zi,zj)
      ElseIf (p.EQ.1) Then
         Phi=E_ThetaP(zi,zj)
      ElseIf (p.EQ.2) Then
         Phi=Theta(i,j) + E_Theta(zi,zj) + IS*E_ThetaP(zi,zj)
     EndIf

   End

      Double Precision Function E_Theta(zi,zj)
!****************************************************************************
!  SUBROUTINE  :  Etheta
!  DESCRIPTION :  Calculates E_Theta function in Pitzer model
!  VARIABLES   :  zM,zN - charge of ions (have same sign, but only different value)
!
!  CREATED     :  2013.05.011  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Function
      Double Precision J

! Call variables
      Integer zi,zj

! Local variables
      Double Precision help,Xij,Xii,Xjj
!-----Start-------------------------------------------------------------------

      If (zi.NE.zj) Then
         help=6.d0*Aphi*sqrIS
         Xij=zi*zj*help
         Xii=zi*zi*help
         Xjj=zj*zj*help

         E_Theta=zi*zj/(4.d0*IS)*(J(Xij,0) - 0.5d0*J(Xii,0) - 0.5d0*J(Xjj,0))
      Else
         E_theta=0.d0
      EndIf

   End

      Double Precision Function E_ThetaP(zi,zj)
!****************************************************************************
!  SUBROUTINE  :  EthetaP
!  DESCRIPTION :  Calculates derivative of E_Theta function in Pitzer model
!  VARIABLES   :  zM,zN - charge of ions (have same sign, but different value)
!
!  CREATED     :  2013.05.011  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Use Pitzer
      Implicit None

! Function
      Double Precision J,E_Theta

! Call variables
      Integer zi,zj

! Local variables
      Double Precision help,Xij,Xii,Xjj,Pij,Pii,Pjj
!-----Start-------------------------------------------------------------------

      If (zi.NE.zj) Then
         help=6.d0*Aphi*sqrIS
         Xij=zi*zj*help
         Xii=zi*zi*help
         Xjj=zj*zj*help

         Pij=Xij*J(Xij,1)
         Pii=Xii*J(Xii,1)
         Pjj=Xjj*J(Xjj,1)

         E_ThetaP= -E_Theta(zi,zj)/IS + zi*zj/(8.d0*IS**2)*(Pij - 0.5d0*Pii - 0.5d0*Pjj)
      Else
         E_thetaP=0.d0
      EndIf

   End

      Double Precision Function J(x,p)
!****************************************************************************
!  SUBROUTINE  :  J
!  DESCRIPTION :  Calculates J function for unsymetrical mixing
!  VARIABLES   :  x - 6zizjAphiI^1/2
!                 p - indicates derivative or not  0=function 1=derivaitive
!
!  CREATED     :  2013.05.09  Baard Kaasa www.scalesim.no
!  CHANGED     :
!****************************************************************************
      Implicit None

! Call variables
      Double Precision x
      Integer p

! Local variables
      Double Precision c(4),T,N,dT,dN
      Data c/4.581d0, -0.7237d0, -0.012d0, 0.528d0 /

!-----Start-------------------------------------------------------------------

! J is calculated from empirical function given by Pitzer, J. Solution Chem, 4, 249, (1975)

      T = x
      N = 4.d0 + c(1)*(X**c(2))*dexp(c(3)*x**c(4))

      If (p.EQ.0) Then
         J= T/N
      ElseIf (p.EQ.1) Then
         dT = 1

         dN = c(2)*x**(c(2)-1.d0)*dexp(c(3)*x**c(4))
         dN = dN + x**c(2)*dexp(c(3)*x**c(4))*(c(3)*c(4)*x**(c(4) - 1.d0))
         dN = dN*c(1)

         J = (dT*N - T*dN) / (N**2)
      EndIf

   End








