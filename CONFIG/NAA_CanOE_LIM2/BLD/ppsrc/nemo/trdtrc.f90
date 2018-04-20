MODULE trdtrc
   !!======================================================================
   !!                       ***  MODULE  trdtrc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  3.0  !  2010-07  (C. Ethe)  Original code (from trdtrc.F90)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------

   INTERFACE trd_trc
      MODULE PROCEDURE trd_trc_trp, trd_trc_bio
   END INTERFACE

CONTAINS

   SUBROUTINE trd_trc_trp( ptrtrd, kjn, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   kjn     ! tracer index
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_trp

   SUBROUTINE trd_trc_bio( ptrbio, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrbio  ! Temperature or U trend
      WRITE(*,*) 'trd_trc_trp : You should not have seen this print! error?', ptrbio(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_trc_bio

   !!======================================================================
END MODULE trdtrc
