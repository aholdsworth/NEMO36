MODULE sedchem

   !!======================================================================
   !! MODULE sedchem  :   Dummy module 
   !!======================================================================
   !! $Id: sedchem.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_chem( kt )         ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE sed_chem

   !!======================================================================


END MODULE sedchem
