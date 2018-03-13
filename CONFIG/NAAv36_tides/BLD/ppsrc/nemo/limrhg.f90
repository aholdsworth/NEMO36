MODULE limrhg
   !!======================================================================
   !!                     ***  MODULE  limrhg  ***
   !!   Ice rheology : sea ice rheology
   !!======================================================================
   !! History :   -   !  2007-03  (M.A. Morales Maqueda, S. Bouillon) Original code
   !!            3.0  !  2008-03  (M. Vancoppenolle) LIM3
   !!             -   !  2008-11  (M. Vancoppenolle, S. Bouillon, Y. Aksenov) add surface tilt in ice rheolohy 
   !!            3.3  !  2009-05  (G.Garric) addition of the lim2_evp cas
   !!            3.4  !  2011-01  (A. Porter)  dynamical allocation 
   !!            3.5  !  2012-08  (R. Benshila)  AGRIF 
   !!            3.6  !  2016-06  (C. Rousset) Rewriting (conserves energy)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_rhg( k1 , k2 )         ! Dummy routine
      WRITE(*,*) 'lim_rhg: You should not have seen this print! error?', k1, k2
   END SUBROUTINE lim_rhg

   !!==============================================================================
END MODULE limrhg
