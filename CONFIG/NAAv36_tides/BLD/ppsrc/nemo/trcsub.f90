MODULE trcsub
   !!======================================================================
   !!                       ***  MODULE trcsubstp  ***
   !!TOP :   Averages physics variables for TOP substepping. 
   !!======================================================================
   !! History :  1.0  !  2011-10  (K. Edwards)  Original
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sub_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_sub_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sub_stp
   SUBROUTINE trc_sub_ini        ! Empty routine
      WRITE(*,*) 'trc_sub_ini: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sub_ini


   !!======================================================================
END MODULE trcsub
