MODULE trcini_cmoc
   !!======================================================================
   !!                         ***  MODULE trcini_cmoc  ***
   !! TOP :   initialisation of the CMOC biochemical model
   !!======================================================================
   !!----------------------------------------------------------------------

   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!           CMOC1  !  2013-15  (O. Riche) code editing for consistency with the rest of CMOC1 adaptation
   !!           CMOC1  !  2017  (A. Holdsworth) adapted CMOC to NEMO3.6
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_cmoc   : CMOC biochemical model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_cmoc   ! called by trcini.F90 module


#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_cmoc.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE trc_ini_cmoc
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE cmoc_ini ***
      !!
      !! ** Purpose :   Initialisation of the CMOC biochemical model
      !!----------------------------------------------------------------------
      !
      USE cmocsms         ! Main P4Z routine
      USE cmocche          !  Chemical model
      USE cmocsink         !  vertical flux of particulate matter due to sinking
      USE cmocsbc          !  Boundary conditions
      USE cmocrem          !  Remineralisation of organic matter
      USE cmocflx          !  Gas exchange
      USE cmocprod         !  Growth rate of the 2 phyto groups
      USE cmocmicro        !  Sources and sinks of microzooplankton
      USE cmocmort         !  Mortality terms for phytoplankton
      USE cmocsed          !  Sedimentation & burial
      !
      REAL(wp), SAVE :: sco2   =  2.312e-3_wp
      REAL(wp), SAVE :: alka0  =  2.423e-3_wp
      REAL(wp), SAVE :: oxyg0  =  177.6e-6_wp 
      REAL(wp), SAVE :: bioma0 =  1.000e-8_wp  
      REAL(wp), SAVE :: no3    =  31.04e-6_wp * 6.625_wp
      !
      INTEGER  ::  ji, jj, jk, ierr
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' cmoc_ini :   CMOC biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

                                                 ! Allocate CMOC arrays
      ierr =         sms_cmoc_alloc()          
      ierr = ierr +  cmoc_che_alloc()
      ierr = ierr +  cmoc_sink_alloc()
      ierr = ierr +  cmoc_flx_alloc()
      ierr = ierr +  cmoc_sed_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'cmoc_alloc: unable to allocate CMOC arrays' )
      !
      ryyss    = nyear_len(1) * rday    ! number of seconds per year
      r1_ryyss = 1. / ryyss
      !

      CALL cmoc_sms_init       !  Maint routine  -reads cmoc namelist
      !                                            ! Time-step

      ! Set biological ratios
      ! ---------------------
      rno3    =  16._wp / 122._wp
     ! po4r    =   1._wp / 122._wp
     ! o2nit   =  32._wp / 122._wp
     ! rdenit  = 105._wp /  16._wp
     ! rdenita =   3._wp /  5._wp
     ! o2ut    = 133._wp / 122._wp

      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         
         trn(:,:,:,jpdic) = sco2
         trn(:,:,:,jptal) = alka0
         trn(:,:,:,jpoxy) = oxyg0
         trn(:,:,:,jppoc) = bioma0
         trn(:,:,:,jpphy) = bioma0
         trn(:,:,:,jpzoo) = bioma0
         trn(:,:,:,jpnch) = bioma0 * 12. / 55.
         trn(:,:,:,jpno3) = no3

         ! ----------------------------------------------------
      END IF


      CALL cmoc_sink_init      !  vertical flux of particulate organic matter
      CALL cmoc_prod_init      !  phytoplankton growth rate over the global ocean.
      CALL cmoc_sbc_init       !  boundary conditions
      CALL cmoc_rem_init       !  remineralisation
      CALL cmoc_mort_init      !  phytoplankton mortality 
      CALL cmoc_micro_init     !  microzooplankton
      CALL cmoc_flx_init       !  gas exchange 

      ndayflxtr = 0

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'Initialization of CMOC tracers done'
      IF(lwp) WRITE(numout,*) 
      !
   END SUBROUTINE trc_ini_cmoc
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No CMOC biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_cmoc             ! Empty routine
   END SUBROUTINE trc_ini_cmoc
#endif

   !!======================================================================
END MODULE trcini_cmoc
