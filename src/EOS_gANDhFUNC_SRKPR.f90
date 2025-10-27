
!************************************************************
!  File Type:     F90 Subroutines
!  File Content:  g and f functions given by Michelsen
!
!                 NB-NB-NB   VALID FOR SKR AND PR  NB-NB-NB
!
!                 NB f-functions are nambed h-functions to avoid confusion wiht captial F functions
!
!  Created:       2016.06.28 Baard Kaasa www.scalesim.no
!  Modified:
!************************************************************


      Double Precision Function gfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start----------------------------------
      gfunc1=dlog((V-B)/V)
	   End Function gfunc1


      Double Precision Function hfunc1(V,B)
      Use EOS
      Use GlobVar
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      hfunc1=1.d0/(Reos*B*(delta1-delta2))*dlog((V+delta1*B)/(V+delta2*B))
	   End Function hfunc1


      Double Precision Function gVfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      gvfunc1=B/(V*(V-B))
	   End Function gvfunc1


      Double Precision Function gBfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      gBfunc1=-1.d0/(V-B)
	   End Function gBfunc1


      Double Precision Function hVfunc1(V,B)
      Use EOS
      Use GlobVar
      Implicit None
      Double Precision V,B
!--Start------------------------------------
	   hvfunc1=-1.d0/(Reos*(V+delta1*B)*(V+delta2*B))
	   End Function hvfunc1


      Double Precision Function hBfunc1(V,B)
      Implicit None
      Double Precision V,B
! Function variables
      Double Precision hfunc1, hvfunc1
!--Start------------------------------------
      hbfunc1=-(hfunc1(V,B) + V*hvfunc1(V,b))/B
	   End Function hbfunc1


      Double Precision Function gVVfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      gvvfunc1 = - 1.d0/(V - B)**2 + 1.d0/V**2
      End Function gvvfunc1


      Double Precision Function hVVfunc1(V,B)
      Use EOS
      Use GlobVar
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      hvvfunc1 = 1.d0/(Reos*B*(delta1-delta2))*(-1.d0/(V+delta1*B)**2 + 1.d0/(V+delta2*B)**2)
      End Function hvvfunc1


      Double Precision Function gBVfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      gbvfunc1 = 1.d0/(V-B)**2
      End Function gbvfunc1


      Double Precision Function hBVfunc1(V,B)
      Implicit None
      Double Precision V,B,hvfunc1,hvvfunc1
!--Start------------------------------------
      hbvfunc1 = -(2.d0*hvfunc1(V,B) + V*hvvfunc1(V,B))/B
      End Function hbvfunc1


      Double Precision Function hbbfunc1(V,B)
      Implicit None
      Double Precision V,B,hbfunc1,hbvfunc1
!--Start------------------------------------
      hbbfunc1 = -(2.d0*hbfunc1(V,B) + V*hbvfunc1(V,B))/B
      End Function hbbfunc1


      Double Precision Function gbbfunc1(V,B)
      Implicit None
      Double Precision V,B
!--Start------------------------------------
      gbbfunc1 = - 1.d0/(V-B)**2
      End Function gbbfunc1













