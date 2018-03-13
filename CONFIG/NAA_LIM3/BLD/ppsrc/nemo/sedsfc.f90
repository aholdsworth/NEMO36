MODULE sedsfc
   !!======================================================================
   !!              ***  MODULE  sedsfc  ***
   !!    Sediment : Data at sediment surface
   !!=====================================================================
   !!======================================================================
   !! MODULE sedsfc  :   Dummy module 
   !!======================================================================
   !! $Id: sedsfc.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_sfc ( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'sed_sfc: You should not have seen this print! error?', kt
   END SUBROUTINE sed_sfc

END MODULE sedsfc
