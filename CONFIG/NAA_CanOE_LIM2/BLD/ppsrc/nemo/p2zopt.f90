MODULE p2zopt
   !!======================================================================
   !!                         ***  MODULE p2zopt  ***
   !! TOP :   LOBSTER Compute the light availability in the water column
   !!======================================================================
   !! History :    -   !  1995-05  (M. Levy) Original code
   !!              -   !  1999-09  (J.-M. Andre, M. Levy) 
   !!              -   !  1999-11  (C. Menkes, M.-A. Foujols) itabe initial
   !!              -   !  2000-02  (M.A. Foujols) change x**y par exp(y*log(x))
   !!   NEMO      2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!             3.2  !  2009-04  (C. Ethe, G. Madec)  minor optimisation + style
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_opt( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_opt: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_opt

   !!======================================================================
END MODULE  p2zopt
