MODULE trcsms_age
   !!======================================================================
   !!                         ***  MODULE trcsms_age  ***
   !! TOP :   Main module of the AGE tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No AGE model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_age( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_age: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_age

   !!======================================================================
END MODULE trcsms_age
