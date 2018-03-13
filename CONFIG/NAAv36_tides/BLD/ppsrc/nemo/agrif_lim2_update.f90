

MODULE agrif_lim2_update
   !!======================================================================
   !!                       ***  MODULE agrif_lim2_update ***
   !! Nesting module :  update surface ocean boundary condition over ice
   !!                   from a child grif
   !! Sea-Ice model  :  LIM 2.0 Sea ice model time-stepping
   !!======================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)  initial version
   !!            3.4   !  08-2012  (R. Benshila, C. Herbaut) update and EVP
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE agrif_lim2_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_lim2_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_lim2_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_lim2_update_empty
END MODULE agrif_lim2_update
