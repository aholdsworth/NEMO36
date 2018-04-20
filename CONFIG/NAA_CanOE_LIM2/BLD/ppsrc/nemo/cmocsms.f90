MODULE cmocsms
   !!======================================================================
   !!                         ***  MODULE cmocsms  ***
   !! TOP :   CMOC Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!           CMOC1  !  2017-08 (A. Holdsworth) adapted CMOC for NEMO3.6
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'cmoc_sms: You should not have seen this print! error?', kt
   END SUBROUTINE cmoc_sms

   !!======================================================================
END MODULE cmocsms 
