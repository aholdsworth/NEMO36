MODULE trcdta
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  reads passive tracer data 
   !!=====================================================================
   !! History :   1.0  !  2002-04  (O. Aumont)  original code
   !!              -   !  2004-03  (C. Ethe)  module
   !!              -   !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!            3.4   !  2010-11  (C. Ethe, G. Madec)  use of fldread + dynamical allocation 
   !!            3.5   !  2013-08  (M. Vichi)  generalization for other BGC models
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dta( kt, sf_dta, ptrfac, ptrc)        ! Empty routine
      WRITE(*,*) 'trc_dta: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dta
   !!======================================================================
END MODULE trcdta
