MODULE p2zexp
   !!======================================================================
   !!                         ***  MODULE p2zsed  ***
   !! TOP :   LOBSTER Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :    -   !  1999    (O. Aumont, C. Le Quere)  original code
   !!              -   !  2001-05 (O. Aumont, E. Kestenare) add sediment computations
   !!             1.0  !  2005-06 (A.-S. Kremeur) new temporal integration for sedpoc
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!             3.5  !  2012-03  (C. Ethe)  Merge PISCES-LOBSTER
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_exp( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_exp: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_exp

   !!======================================================================
END MODULE  p2zexp
