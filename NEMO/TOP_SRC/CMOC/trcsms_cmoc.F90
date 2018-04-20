MODULE trcsms_cmoc
   !!======================================================================
   !!                         ***  MODULE trcsms_cmoc  ***
   !! TOP :   CMOC Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_cmoc 
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_cmoc        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE par_cmoc
   USE cmocsms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_cmoc    ! called in trcsms.F90
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_cmoc.F90 4147 2013-11-04 11:51:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_cmoc ***
      !!
      !! ** Purpose :   Initialisation of the CMOC biochemical model
      !!----------------------------------------------------------------------


   SUBROUTINE trc_sms_cmoc( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_cmoc  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!                routines of CMOC or LOBSTER bio-model
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!---------------------------------------------------------------------
      !
      CALL cmoc_sms( kt )   !  CMOC

      !
   END SUBROUTINE trc_sms_cmoc

#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_cmoc( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_cmoc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_cmoc
#endif 

   !!======================================================================
END MODULE trcsms_cmoc 
