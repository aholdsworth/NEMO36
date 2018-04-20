MODULE trcsms_cmoc
   !!======================================================================
   !!                         ***  MODULE trcsms_cmoc  ***
   !! TOP :   CMOC Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_cmoc( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_cmoc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_cmoc

   !!======================================================================
END MODULE trcsms_cmoc 
