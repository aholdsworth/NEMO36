MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_my_trc( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_my_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_my_trc

   !!======================================================================
END MODULE trcsms_my_trc
