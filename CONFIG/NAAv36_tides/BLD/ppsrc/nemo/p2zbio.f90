MODULE p2zbio
   !!======================================================================
   !!                         ***  MODULE p2zbio  ***
   !! TOP :   LOBSTER
   !!======================================================================
   !! History :    -   !  1999-07  (M. Levy) Original code
   !!              -   !  2000-12  (E. Kestenare) assign a parameter to name individual tracers
   !!              -   !  2001-03  (M. Levy)  LNO3 + dia2d 
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_bio( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_bio: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_bio

   !!======================================================================
END MODULE p2zbio
