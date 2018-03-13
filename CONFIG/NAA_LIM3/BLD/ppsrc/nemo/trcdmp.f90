MODULE trcdmp
   !!======================================================================
   !!                       ***  MODULE  trcdmp  ***
   !! Ocean physics: internal restoring trend on passive tracers
   !!======================================================================
   !! History :  OPA  !  1991-03  (O. Marti, G. Madec)  Original code
   !!                 !  1996-01  (G. Madec) statement function for e3
   !!                 !  1997-05  (H. Loukos)  adapted for passive tracers
   !!    NEMO    9.0  !  2004-03  (C. Ethe)    free form + modules
   !!            3.2  !  2007-02  (C. Deltel)  Diagnose ML trends for passive tracers
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dmp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dmp


   !!======================================================================
END MODULE trcdmp
