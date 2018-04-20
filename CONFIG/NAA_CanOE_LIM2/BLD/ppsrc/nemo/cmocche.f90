MODULE cmocche
   !!======================================================================
   !!                         ***  MODULE cmocche  ***
   !! TOP :   CMOC Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_che( kt )                   ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'cmoc_che: You should not have seen this print! error?', kt
   END SUBROUTINE cmoc_che

   !!======================================================================
END MODULE  cmocche
