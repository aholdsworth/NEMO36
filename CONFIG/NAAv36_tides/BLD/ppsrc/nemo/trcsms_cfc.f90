MODULE trcsms_cfc
   !!======================================================================
   !!                      ***  MODULE trcsms_cfc  ***
   !! TOP : CFC main model
   !!======================================================================
   !! History :  OPA  !  1999-10  (JC. Dutay)  original code
   !!  NEMO      1.0  !  2004-03  (C. Ethe) free form + modularity
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  reorganisation
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No CFC tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_cfc( kt )       ! Empty routine
      WRITE(*,*) 'trc_sms_cfc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_cfc

   !!======================================================================
END MODULE trcsms_cfc
