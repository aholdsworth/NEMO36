MODULE cmocsink
   !!======================================================================
   !!                         ***  MODULE cmocsink  ***
   !! TOP :  CMOC  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking, Oct 28th 2015 simplified sink scheme according to Jim's own modifications for CMOC2
   !!           CMOC1  !  2016-02  (N. Swart) Adds calcite sinking.
   !!           CMOC1  !  2017-08  (A.Holdsworth) updated to NEMO 3.6.
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sink                    ! Empty routine
   END SUBROUTINE cmoc_sink

   !!======================================================================
END MODULE  cmocsink
