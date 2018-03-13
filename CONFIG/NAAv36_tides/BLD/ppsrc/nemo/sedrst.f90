MODULE sedrst
   !!======================================================================
   !! MODULE sedrst :   Dummy module 
   !!======================================================================
   !! $Id: sedrst.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_rst_read                      ! Empty routines
   END SUBROUTINE sed_rst_read
   SUBROUTINE sed_rst_wri( kt )
      INTEGER, INTENT ( in ) :: kt
      WRITE(*,*) 'sed_rst_wri: You should not have seen this print! error?', kt
   END SUBROUTINE sed_rst_wri   

END MODULE sedrst
