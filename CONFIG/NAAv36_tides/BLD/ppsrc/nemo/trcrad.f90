MODULE trcrad
   !!======================================================================
   !!                       ***  MODULE  trcrad  ***
   !! Ocean passive tracers:  correction of negative concentrations
   !!======================================================================
   !! History :   -   !  01-01  (O. Aumont & E. Kestenare)  Original code
   !!            1.0  !  04-03  (C. Ethe)  free form F90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                         NO TOP model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rad( kt )              ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rad
   
   !!======================================================================
END MODULE trcrad
