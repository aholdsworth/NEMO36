MODULE p2zsms
   !!======================================================================
   !!                         ***  MODULE p2zsms  ***
   !! TOP :   Time loop of LOBSTER model
   !!======================================================================
   !! History :   1.0  !            M. Levy
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                     No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_sms: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_sms

   !!======================================================================
END MODULE p2zsms
