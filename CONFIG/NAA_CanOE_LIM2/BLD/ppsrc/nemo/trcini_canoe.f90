MODULE trcini_canoe
    !!======================================================================
    !!                         ***  MODULE trcini_canoe  ***
    !! TOP :   initialisation of the CanOE biochemical model
    !!======================================================================
    !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
    !!              -   !  1999-10  (O. Aumont, C. Le Quere)
    !!              -   !  2002     (O. Aumont) PISCES
    !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
    !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.canoe.h90
    !!             3.5  !  2012-05  (C. Ethe) Merge PISCESLOBSTER
    !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
    !!----------------------------------------------------------------------

    !!----------------------------------------------------------------------
    !!   'key_canoe'                                       CanOE bio-model
    !!----------------------------------------------------------------------
    !! trc_ini_canoe   : CanOE biochemical model initialisation
    !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE trcsms_canoe    ! Main P4Z routine
   USE canche          !  Chemical model
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canopt          !  optical model
   USE cansbc          !  Boundary conditions
   USE canrem          !  Remineralisation of organic matter
   USE canflx          !  Gas exchange
   USE canprod         !  Growth rate of the 2 phyto groups
   USE canmicro        !  Sources and sinks of microzooplankton
   USE canmeso         !  Sources and sinks of mesozooplankton
   USE canmort         !  Mortality terms for phytoplankton
   USE canlys          !  Calcite saturation
   USE cansed          !  Sedimentation & burial
   USE canint          !  Temperature-dependence
   !
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_canoe   ! called by trcini.F90 module
   
   REAL(wp) :: sco2   =  2.312e-3_wp
   REAL(wp) :: alka0  =  2.423e-3_wp !mol/L
   REAL(wp) :: oxyg0  =  177.6_wp 
   REAL(wp) :: bioma0 =  1.000e-2_wp !nmol/L 
   REAL(wp) :: no3    =  31.04_wp
 
   !

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
   !! $Id: trcini_canoe.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_canoe
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_canoe ***
      !!
      !! ** Purpose :   Initialisation of the CanOE biochemical model
      !!----------------------------------------------------------------------

      !!
      !! ** Purpose :   Initialisation of the CanOE biochemical model
      !!----------------------------------------------------------------------
      !
      
      INTEGER  ::  ji, jj, jk, ierr
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_canoe :   CanOE biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

                                                 ! Allocate CanOE arrays
      ierr =         sms_canoe_alloc()          
      ierr = ierr +  can_che_alloc()
      ierr = ierr +  can_sink_alloc()
      ierr = ierr +  can_opt_alloc()
      ierr = ierr +  can_prod_alloc()
      ierr = ierr +  can_flx_alloc()
      ierr = ierr +  can_sed_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'canoe_alloc: unable to allocate CanOE arrays' )
      !
      ryyss    = nyear_len(1) * rday    ! number of seconds per year
      r1_ryyss = 1. / ryyss
      !

      CALL can_sms_init       !  Maint routine
      !                                            ! Time-step

      ! Set biological ratios
      ! ---------------------
      ! Set biological ratios
      ! ---------------------
      rr_c2n  = 6.625_wp
      rr_fe2n = 33._wp
      rr_n2c=1./rr_c2n
      rr_n2fe=1./rr_fe2n
      rr_fe2c=rr_fe2n/rr_c2n
      rr_c2fe=1./rr_fe2c

!      rno3    =  16._wp / 122._wp
!      po4r    =   1._wp / 122._wp
!      o2nit   =  32._wp / 122._wp
!      rdenit  = 105._wp /  16._wp
!      rdenita =   3._wp /  5._wp
!      o2ut    = 133._wp / 122._wp
!
      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         !all quantities are in micromolar except for iron which is in picomolar 
         trn(:,:,:,jpdic) = sco2
         trn(:,:,:,jptal) = alka0
         trn(:,:,:,jpoxy) = oxyg0
         trn(:,:,:,jpcal) = bioma0
         trn(:,:,:,jppoc) = bioma0
         trn(:,:,:,jpphy) = bioma0
         trn(:,:,:,jpnn)  = bioma0 * 16./106.
         trn(:,:,:,jpnfe) = bioma0 * 5.
         trn(:,:,:,jpnch) = bioma0 * 12. / 55.
         trn(:,:,:,jpdia) = bioma0
         trn(:,:,:,jpdn)  = bioma0 * 16./106.
         trn(:,:,:,jpdfe) = bioma0 * 5.
         trn(:,:,:,jpdch) = bioma0 * 12. / 55.
         trn(:,:,:,jpzoo) = bioma0
         trn(:,:,:,jpmes) = bioma0
         trn(:,:,:,jpfer) = 6.0E+2
         trn(:,:,:,jpgoc) = bioma0
         trn(:,:,:,jpno3) = no3
         trn(:,:,:,jpnh4) = bioma0 * 16./106.

      END IF


      CALL can_int_init       !  temperature dependence of biological rates
      CALL can_sink_init      !  vertical flux of particulate organic matter
      CALL can_opt_init       !  Optic: PAR in the water column
      CALL can_prod_init      !  phytoplankton growth rate over the global ocean.
      CALL can_sbc_init       !  boundary conditions
      CALL can_rem_init       !  remineralisation
      CALL can_mort_init      !  phytoplankton mortality 
      CALL can_micro_init     !  microzooplankton
      CALL can_meso_init      !  mesozooplankton
      CALL can_lys_init       !  calcite saturation
      CALL can_flx_init       !  gas exchange 

      ndayflxtr = 0

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'Initialization of CanOE tracers done'
      IF(lwp) WRITE(numout,*) 
      !
   END SUBROUTINE trc_ini_canoe


   !!======================================================================
END MODULE trcini_canoe
