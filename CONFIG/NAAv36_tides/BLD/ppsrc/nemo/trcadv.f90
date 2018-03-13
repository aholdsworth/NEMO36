MODULE trcadv
   !!==============================================================================
   !!                       ***  MODULE  trcadv  ***
   !! Ocean passive tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  05-11  (G. Madec)  Original code
   !!            3.0  !  10-06  (C. Ethe)   Adapted to passive tracers
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_adv: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv

  !!======================================================================
END MODULE trcadv
