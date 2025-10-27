!************************************************************
!  File Type:     MODULE
!  File Content:  Parameters related to Pitzer model
!
!************************************************************


MODULE Pitzer
   
   Use Aqua
   
   SAVE
   Integer CatM(0:NoAq),AnX(0:NoAq),NeuN(0:NoAq)
   
   Integer AqInPitz(NoAq)                                                  ! Defines if specie is included in Pitzer calculation

! Pitzer temperature parameters - these are the temperature coefficients   
      Integer, Parameter :: Npar=8
      Double Precision BetaCC(0:2,NoAq,NoAq,Npar)
      Double Precision CphiCC(    NoAq,NoAq,Npar)
      Double Precision PsiCC (    NoAq,NoAq,NoAq,Npar)
      Double Precision ThetCC(    NoAq,NoAq,Npar)
      Double Precision LambCC(    NoAq,NoAq,Npar)
      Double PRecision KsiCC (    NoAq,NoAq,NoAq,Npar)

      Integer BetaCCEq(0:2,NoAq,NoAq)
      Integer CphiCCEq(    NoAq,NoAq)
      Integer PsiCCEq (    NoAq,NoAq,NoAq)
      Integer ThetCCEq(    NoAq,NoAq)
      Integer LambCCEq(    NoAq,NoAq)
      Integer KsiCCEq(     NoAq,NoAq,NoAq)
      
   
! Pitzer parameters
      Double Precision Beta0(NoAq,NoAq),Beta1(NoAq,NoAq),Beta2(NoAq,NoAq)  ! Main interactions Beta ca
      Double Precision Cphi(NoAq,NoAq)                                     ! Main interaction C ca
      Double Precision Theta(NoAq,NoAq)                                    ! Pair of same sign cc and aa
      Double Precision Psi(NoAq,NoAq,NoAq)                                 ! Tripple interactions cc'a and caa'
      Double Precision Lambda(NoAq,NoAq)                                   ! Neutral-ion pairs nc and na
      Double Precision Ksi(NoAq,NoAq,NoAq),Tau(NoAq,NoAq,NoAq)             ! Neutral-ion tripplets nca
      Double Precision PitzAlpha(2,NoAq,NoAq)      
      
! Pitzer constants
      Double Precision, Parameter :: Debye_b=1.2d0
      
! Water properties
      Double Precision Dielec,rhow,Aphi,By,Av,Ay
      Double Precision IS,Z,sqrIS

! Debue-Hückel, F and SumCca
      Double PRecision F,fy,SumCca,OsmoC

End MODULE Pitzer