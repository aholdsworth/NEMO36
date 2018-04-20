MODULE trcsms
   !!======================================================================
   !!                         ***  MODULE trcsms  ***
   !! TOP :   Time loop of passive tracers sms
   !!======================================================================
   !! History :   1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc            !
   USE trc                !
   USE trcsms_pisces      ! PISCES biogeo-model
   USE trcsms_canoe      ! PISCES biogeo-model
   USE trcsms_cmoc      !CMOC biogeo-model
   USE trcsms_cfc         ! CFC 11 & 12
   USE trcsms_c14b        ! C14b tracer 
   USE trcsms_my_trc      ! MY_TRC  tracers
   USE prtctl_trc         ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms    ! called in trcstp.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms.F90 3680 2012-11-27 14:42:24Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms  ***
      !!
      !! ** Purpose :   Managment of the time loop of passive tracers sms 
      !!
      !! ** Method  : -  call the main routine of of each defined tracer model
      !! -------------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_sms')
      !
      IF( lk_pisces  )   CALL trc_sms_pisces ( kt )    ! main program of PISCES 
      IF( lk_canoe  )   CALL trc_sms_canoe ( kt )    ! main program of CanOE 
      IF( lk_cmoc  )     CALL trc_sms_cmoc ( kt )    ! main program of CMOC 
      IF( lk_cfc     )   CALL trc_sms_cfc    ( kt )    ! surface fluxes of CFC
      IF( lk_c14b    )   CALL trc_sms_c14b   ( kt )    ! surface fluxes of C14
      IF( lk_my_trc  )   CALL trc_sms_my_trc ( kt )    ! MY_TRC  tracers

      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sms ')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_sms')
      !
   END SUBROUTINE trc_sms


   !!======================================================================
END MODULE  trcsms
