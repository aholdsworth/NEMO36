MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficient 
   !!=====================================================================
   !! History :        ! 1997-07  (G. Madec)  from inimix.F split in 2 routines
   !!   NEMO      1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!             2.0  ! 2005-11  (G. Madec)  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE ldfslp          ! ???
   USE in_out_manager  ! I/O manager
   USE ioipsl
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_init   ! called by opa.F90

   !! * Substitutions
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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra.F90 4624 2014-04-28 12:09:03Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !! 
      !! ** Purpose :   initializations of the tracer lateral mixing coeff.
      !!
      !! ** Method  :   the Eddy diffusivity and eddy induced velocity ceoff.
      !!      are defined as follows:
      !!         default option   : constant coef. aht0, aeiv0 (namelist)
      !!        'key_traldf_c1d': depth dependent coef. defined in 
      !!                            in ldf_tra_c1d routine
      !!        'key_traldf_c2d': latitude and longitude dependent coef.
      !!                            defined in ldf_tra_c2d routine
      !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
      !!                            defined in ldf_tra_c3d routine
      !!
      !!      N.B. User defined include files.  By default, 3d and 2d coef.
      !!      are set to a constant value given in the namelist and the 1d
      !!      coefficients are initialized to a hyperbolic tangent vertical
      !!      profile.
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio               ! temporary integer
      INTEGER ::   ios                  ! temporary integer
      LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
      !! 
      NAMELIST/namtra_ldf/ ln_traldf_lap  , ln_traldf_bilap,                  &
         &                 ln_traldf_level, ln_traldf_hor  , ln_traldf_iso,   &
         &                 ln_traldf_grif , ln_traldf_gdia ,                  &
         &                 ln_triad_iso   , ln_botmix_grif ,                  &
         &                 rn_aht_0       , rn_ahtb_0      , rn_aeiv_0,       &
         &                 rn_slpmax      , rn_chsmag      ,    rn_smsh,      &
         &                 rn_aht_m
      !!----------------------------------------------------------------------

      !  Define the lateral tracer physics parameters
      ! =============================================
    

      REWIND( numnam_ref )              ! Namelist namtra_ldf in reference namelist : Lateral physics on tracers
      READ  ( numnam_ref, namtra_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namtra_ldf in configuration namelist : Lateral physics on tracers
      READ  ( numnam_cfg, namtra_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_ldf )

      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_init : lateral tracer physics'
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtra_ldf : lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      laplacian operator            ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '      bilaplacian operator          ln_traldf_bilap = ', ln_traldf_bilap
         WRITE(numout,*) '      iso-level                     ln_traldf_level = ', ln_traldf_level
         WRITE(numout,*) '      horizontal (geopotential)     ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '      iso-neutral                   ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '      iso-neutral (Griffies)        ln_traldf_grif  = ', ln_traldf_grif
         WRITE(numout,*) '      Griffies strmfn diagnostics   ln_traldf_gdia  = ', ln_traldf_gdia
         WRITE(numout,*) '      lateral eddy diffusivity      rn_aht_0        = ', rn_aht_0
         WRITE(numout,*) '      background hor. diffusivity   rn_ahtb_0       = ', rn_ahtb_0
         WRITE(numout,*) '      eddy induced velocity coef.   rn_aeiv_0       = ', rn_aeiv_0
         WRITE(numout,*) '      maximum isoppycnal slope      rn_slpmax       = ', rn_slpmax
         WRITE(numout,*) '      pure lateral mixing in ML     ln_triad_iso    = ', ln_triad_iso
         WRITE(numout,*) '      lateral mixing on bottom      ln_botmix_grif  = ', ln_botmix_grif
         WRITE(numout,*)
      ENDIF

      !                                ! convert DOCTOR namelist names into OLD names
      aht0  = rn_aht_0
      ahtb0 = rn_ahtb_0
      aeiv0 = rn_aeiv_0

      !                                ! Parameter control

      ! ... Check consistency for type and direction :
      !           ==> will be done in traldf module

      ! ... Space variation of eddy coefficients
      ioptio = 0
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio + 1
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'
        ELSEIF( ioptio > 1 ) THEN
           CALL ctl_stop('          use only one of the following keys:',   &
             &           ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d' )
      ENDIF

      IF( ln_traldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
         IF( aht0 > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal diffusivity coef. aht0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
         IF( aht0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop('The horizontal diffusivity coef. aht0 must be positive' )
      ENDIF


      !  Lateral eddy diffusivity and eddy induced velocity coefficients
      ! ================================================================
      CALL ldf_tra_c3d( ll_print )      ! aht = 3D coef. = F( longitude, latitude, depth )

        IF(lwp) WRITE(numout,*)' SMAGORINSKY DIFFUSION'
        IF(lwp .AND. rn_smsh < 1)  WRITE(numout,*)' only  shear is used '
        IF(lwp.and.ln_traldf_bilap) CALL ctl_stop(' SMAGORINSKY + BILAPLACIAN - UNSTABLE OR NON_CONSERVATIVE' )

      !
   END SUBROUTINE ldf_tra_init

   !!----------------------------------------------------------------------
   !!                        ***  ldftra_c3d.h90  ***
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_c3d.h90 3634 2012-11-22 17:32:14Z acc $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_tra_c3d( ld_print )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldftra_c3d  ***
      !!                
      !! ** Purpose :
      !!     initializations of non uniform horizontal eddy diffusivity coef.
      !!
      !! ** Method :
      !!      3D eddy diffusivity coefficients ( longitude, latitude, depth )
      !!
      !!       biharmonic operator    : ahtt = defined at T-level
      !!                                ahtu,ahtv,ahtw never used
      !!
      !!       harmonic operator (ahtt never used)
      !!         iso-model level   : ahtu, ahtv defined at u-, v-points
      !!         isopycnal         : ahtu, ahtv, ahtw defined at u-, v-, w-pts
      !!         or geopotential
      !!
      !!       eddy induced velocity
      !!         always harmonic   : aeiu, aeiv, aeiw defined at u-, v-, w-pts
      !!----------------------------------------------------------------------
      USE ioipsl
      !
      LOGICAL, INTENT (in) ::   ld_print   ! If true, output arrays on numout
      !!----------------------------------------------------------------------

      IF( lk_traldf_eiv ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' ldf_tra_c3d : 3D eddy diffusivity and eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~   --  '
         IF(lwp) WRITE(numout,*) '               Coefficients set to constant'
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' ldf_tra_c3d : 3D eddy diffusivity coefficient'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~   -- '
         IF(lwp) WRITE(numout,*) '               Coefficients set to constant'
      ENDIF


      ! biharmonic operator   (T-point)
      ! -------------------

      ahtt(:,:,:) = aht0              ! set ahtt at T-point (here no space variation)

      CALL lbc_lnk( ahtt, 'T', 1. )   ! Lateral boundary conditions on ( ahtt )

      IF(lwp .AND. ld_print ) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtt at k = 1'
         CALL prihre( ahtt(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-3, numout )
      ENDIF


      ! harmonic operator   (U-, V-, W-points)
      ! ----------------- 

      ahtu(:,:,:) = aht0              ! set ahtu = ahtv at u- and v-points,
      ahtv(:,:,:) = aht0              ! and ahtw at w-point
      ahtw(:,:,:) = aht0              ! (here example: no space variation)
      
      CALL lbc_lnk( ahtu, 'U', 1. )   ! Lateral boundary conditions
      CALL lbc_lnk( ahtv, 'V', 1. )
      CALL lbc_lnk( ahtw, 'W', 1. )

      IF(lwp .AND. ld_print ) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtu at k = 1'
         CALL prihre( ahtu(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-3, numout )
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtv at k = 1'
         CALL prihre( ahtv(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-3, numout )
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: ahtw at k = 1'
         CALL prihre( ahtw(:,:,1), jpi, jpj, 1, jpi, 1,   &
            &                                1, jpj, 1, 1.e-3, numout )
      ENDIF

   END SUBROUTINE ldf_tra_c3d

   !!======================================================================
END MODULE ldftra
