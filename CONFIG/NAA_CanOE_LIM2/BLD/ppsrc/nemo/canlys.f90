MODULE canlys
   !!======================================================================
   !!                         ***  MODULE canlys  ***
   !! TOP :   CanOE 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'can_lys: You should not have seen this print! error?', kt
   END SUBROUTINE can_lys
   !!======================================================================
END MODULE  canlys
