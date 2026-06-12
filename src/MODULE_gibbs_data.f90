




   Module gibbs_data
   Use GlobVar
   Use GasRx


   ! Total number of reactive species
   Integer :: nreact
   Integer :: react_idx(ncomp)   ! maps reactive index -> global index

   Real(dP) :: CompGibbs(ncomp)
   Real(dp) :: atom_matrix(nelem, ncomp)


End Module gibbs_data