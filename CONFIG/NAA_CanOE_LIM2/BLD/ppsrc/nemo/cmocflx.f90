MODULE cmocflx
   !!======================================================================
   !!                         ***  MODULE cmocflx  ***
   !! TOP :   CMOC CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_flx( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'cmoc_flx: You should not have seen this print! error?', kt
   END SUBROUTINE cmoc_flx

   !!======================================================================
END MODULE  cmocflx
