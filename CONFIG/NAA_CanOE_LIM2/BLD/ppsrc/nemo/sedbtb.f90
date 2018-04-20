MODULE sedbtb
   !!======================================================================
   !! MODULE sedbtb  :   Dummy module 
   !!======================================================================
   !! $Id: sedbtb.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_btb( kt )         ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'sed_btb: You should not have seen this print! error?', kt
   END SUBROUTINE sed_btb

   !!======================================================================

END MODULE sedbtb
