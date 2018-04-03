MODULE trcsms_canoe
   !!======================================================================
   !!                         ***  MODULE trcsms_canoe  ***
   !! TOP :   CanOE Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_canoe( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_canoe: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_canoe

   !!======================================================================
END MODULE trcsms_canoe 
