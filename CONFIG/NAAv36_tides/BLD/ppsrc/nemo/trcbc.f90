MODULE trcbc
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  module for passive tracer boundary conditions
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_bc_read( kt )        ! Empty routine
      WRITE(*,*) 'trc_bc_read: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bc_read

   !!======================================================================
END MODULE trcbc
