MODULE seddsr
   !!======================================================================
   !! MODULE seddsr  :   Dummy module 
   !!======================================================================
   !! $Id: seddsr.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_dsr ( kt )
     INTEGER, INTENT(in) :: kt
     WRITE(*,*) 'sed_dsr: You should not have seen this print! error?', kt
  END SUBROUTINE sed_dsr
   
END MODULE seddsr
