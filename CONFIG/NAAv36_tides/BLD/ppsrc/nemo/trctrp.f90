MODULE trctrp
   !!======================================================================
   !!                       ***  MODULE trctrp  ***
   !! Ocean Physics    : manage the passive tracer transport
   !!======================================================================
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             3.3  !  2010-07 (C. Ethe) Merge TRA-TRC
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        No TOP models
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp( kstp )              ! Empty routine
      INTEGER, INTENT(in) ::   kstp
      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kstp
   END SUBROUTINE trc_trp
   
   !!======================================================================
END MODULE trctrp
