MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!             3.5  !  2012-05  (C. Ethe) Merge PISCES-LOBSTER
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_pisces   : PISCES biochemical model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_pisces   ! called by trcini.F90 module


   !!----------------------------------------------------------------------
   !!                    ***  top_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose : Statement function file: to be include in all passive tracer modules
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec) new architecture
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! s* or z*-coordinate (3D + time dependency) + use of additional now arrays (..._n)







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n





   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   Default option :                         eiv: dummy variables
   !!----------------------------------------------------------------------
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 6312 2016-02-15 11:43:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c3d' :                 aht: 3D coefficient
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: top_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_pisces.F90 6324 2016-02-18 08:39:38Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------

      IF( lk_p4z ) THEN  ;   CALL p4z_ini   !  PISCES
      ELSE               ;   CALL p2z_ini   !  LOBSTER
      ENDIF

   END SUBROUTINE trc_ini_pisces

   SUBROUTINE p4z_ini
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_ini ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------
      !
      USE p4zsms          ! Main P4Z routine
      USE p4zche          !  Chemical model
      USE p4zsink         !  vertical flux of particulate matter due to sinking
      USE p4zopt          !  optical model
      USE p4zsbc          !  Boundary conditions
      USE p4zfechem       !  Iron chemistry
      USE p4zrem          !  Remineralisation of organic matter
      USE p4zflx          !  Gas exchange
      USE p4zlim          !  Co-limitations of differents nutrients
      USE p4zprod         !  Growth rate of the 2 phyto groups
      USE p4zmicro        !  Sources and sinks of microzooplankton
      USE p4zmeso         !  Sources and sinks of mesozooplankton
      USE p4zmort         !  Mortality terms for phytoplankton
      USE p4zlys          !  Calcite saturation
      USE p4zsed          !  Sedimentation & burial
      !
      REAL(wp), SAVE :: sco2   =  2.312e-3_wp
      REAL(wp), SAVE :: alka0  =  2.426e-3_wp
      REAL(wp), SAVE :: oxyg0  =  177.6e-6_wp 
      REAL(wp), SAVE :: po4    =  2.165e-6_wp 
      REAL(wp), SAVE :: bioma0 =  1.000e-8_wp  
      REAL(wp), SAVE :: silic1 =  91.51e-6_wp  
      REAL(wp), SAVE :: no3    =  30.9e-6_wp * 7.625_wp
      !
      INTEGER  ::  ji, jj, jk, ierr
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' p4z_ini :   PISCES biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

                                                 ! Allocate PISCES arrays
      ierr =         sms_pisces_alloc()          
      ierr = ierr +  p4z_che_alloc()
      ierr = ierr +  p4z_sink_alloc()
      ierr = ierr +  p4z_opt_alloc()
      ierr = ierr +  p4z_prod_alloc()
      ierr = ierr +  p4z_rem_alloc()
      ierr = ierr +  p4z_flx_alloc()
      ierr = ierr +  p4z_sed_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'pisces_alloc: unable to allocate PISCES arrays' )
      !
      ryyss    = nyear_len(1) * rday    ! number of seconds per year
      r1_ryyss = 1. / ryyss
      !

      CALL p4z_sms_init       !  Maint routine
      !                                            ! Time-step

      ! Set biological ratios
      ! ---------------------
      rno3    =  16._wp / 122._wp
      po4r    =   1._wp / 122._wp
      o2nit   =  32._wp / 122._wp
      o2ut    = 133._wp / 122._wp
      rdenit  =  ( ( o2ut + o2nit ) * 0.80 - rno3 - rno3 * 0.60 ) / rno3
      rdenita =   3._wp /  5._wp


      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         
         trn(:,:,:,jpdic) = sco2
         trn(:,:,:,jpdoc) = bioma0
         trn(:,:,:,jptal) = alka0
         trn(:,:,:,jpoxy) = oxyg0
         trn(:,:,:,jpcal) = bioma0
         trn(:,:,:,jppo4) = po4 / po4r
         trn(:,:,:,jppoc) = bioma0
         trn(:,:,:,jpgoc) = bioma0
         trn(:,:,:,jpbfe) = bioma0 * 5.e-6
         trn(:,:,:,jpsil) = silic1
         trn(:,:,:,jpdsi) = bioma0 * 0.15
         trn(:,:,:,jpgsi) = bioma0 * 5.e-6
         trn(:,:,:,jpphy) = bioma0
         trn(:,:,:,jpdia) = bioma0
         trn(:,:,:,jpzoo) = bioma0
         trn(:,:,:,jpmes) = bioma0
         trn(:,:,:,jpfer) = 0.6E-9
         trn(:,:,:,jpsfe) = bioma0 * 5.e-6
         trn(:,:,:,jpdfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnch) = bioma0 * 12. / 55.
         trn(:,:,:,jpdch) = bioma0 * 12. / 55.
         trn(:,:,:,jpno3) = no3
         trn(:,:,:,jpnh4) = bioma0

         ! initialize the half saturation constant for silicate
         ! ----------------------------------------------------
         xksi(:,:)    = 2.e-6
         xksimax(:,:) = xksi(:,:)
      END IF


      CALL p4z_sink_init      !  vertical flux of particulate organic matter
      CALL p4z_opt_init       !  Optic: PAR in the water column
      CALL p4z_lim_init       !  co-limitations by the various nutrients
      CALL p4z_prod_init      !  phytoplankton growth rate over the global ocean.
      CALL p4z_sbc_init       !  boundary conditions
      CALL p4z_fechem_init    !  Iron chemistry
      CALL p4z_rem_init       !  remineralisation
      CALL p4z_mort_init      !  phytoplankton mortality 
      CALL p4z_micro_init     !  microzooplankton
      CALL p4z_meso_init      !  mesozooplankton
      CALL p4z_lys_init       !  calcite saturation
      CALL p4z_flx_init       !  gas exchange 

      ndayflxtr = 0

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'Initialization of PISCES tracers done'
      IF(lwp) WRITE(numout,*) 
      !
   END SUBROUTINE p4z_ini

   SUBROUTINE p2z_ini
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE p2z_ini ***
      !!
      !! ** Purpose :   Initialisation of the LOBSTER biochemical model
      !!----------------------------------------------------------------------
      !
   END SUBROUTINE p2z_ini

   !!======================================================================
END MODULE trcini_pisces
