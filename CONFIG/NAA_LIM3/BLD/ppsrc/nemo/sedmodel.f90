MODULE sedmodel
   !!======================================================================
   !! MODULE sedmodel  :   Dummy module 
   !!======================================================================
   LOGICAL, PUBLIC, PARAMETER ::   lk_sed = .FALSE.     !: sediment flag
   !! $Id: sedmodel.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   SUBROUTINE sed_model( kt )         ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'sed_stp: You should not have seen this print! error?', kt
   END SUBROUTINE sed_model

END MODULE sedmodel
