MODULE trcdia
   !!======================================================================
   !!                       *** MODULE trcdia ***
   !! TOP :   Output of passive tracers
   !!======================================================================
   !! History :   OPA  !  1995-01 (M. Levy)  Original code
   !!              -   !  1998-01 (C. Levy) NETCDF format using ioipsl interface
   !!              -   !  1999-01 (M.A. Foujols) adapted for passive tracer
   !!              -   !  1999-09 (M.A. Foujols) split into three parts
   !!   NEMO      1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!                  !  2008-05 (C. Ethe re-organization)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dia( kt )                      ! Empty routine   
      INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_dia   

   !!======================================================================
END MODULE trcdia
