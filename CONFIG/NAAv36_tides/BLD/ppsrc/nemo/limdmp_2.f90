MODULE limdmp_2
   !!======================================================================
   !!                       ***  MODULE limdmp_2   ***
   !!  LIM-2 ice model : restoring Ice thickness and Fraction leads
   !!======================================================================
   !! History :   2.0  !  2004-04 (S. Theetten) Original code
   !!             3.3  !  2010-06 (J.-M. Molines) use of fldread
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module                  No ice damping
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dmp_2( kt )        ! Dummy routine
      WRITE(*,*) 'lim_dmp_2: You should not see this print! error? ', kt
   END SUBROUTINE lim_dmp_2

   !!======================================================================
END MODULE limdmp_2
