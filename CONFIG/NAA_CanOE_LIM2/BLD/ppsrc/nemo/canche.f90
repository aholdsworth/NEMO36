MODULE canche
   !!======================================================================
   !!                         ***  MODULE canche  ***
   !! TOP :   CanOE Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_che( kt )                   ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'can_che: You should not have seen this print! error?', kt
   END SUBROUTINE can_che

   !!======================================================================
END MODULE  canche
