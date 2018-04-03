MODULE canflx
   !!======================================================================
   !!                         ***  MODULE canflx  ***
   !! TOP :   CanOE CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_flx( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'can_flx: You should not have seen this print! error?', kt
   END SUBROUTINE can_flx

   !!======================================================================
END MODULE  canflx
