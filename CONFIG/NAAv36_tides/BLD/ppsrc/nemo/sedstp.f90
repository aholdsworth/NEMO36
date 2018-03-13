MODULE sedstp
   !!======================================================================
   !! MODULE sedstp  :   Dummy module 
   !!======================================================================
   !! $Id: sedstp.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_stp( kt )         ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'sed_stp: You should not have seen this print! error?', kt
   END SUBROUTINE sed_stp
END MODULE sedstp
