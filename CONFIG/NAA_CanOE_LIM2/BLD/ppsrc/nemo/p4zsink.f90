MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes


   INTEGER  :: ik100


   !!* Substitution
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
   !! $Id: p4zsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, jit
      INTEGER  ::   iiter1, iiter2
      REAL(wp) ::   zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zagg , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zfact, zwsmax, zmax, zstep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / 5000._wp
               wsbio4(ji,jj,jk) = wsbio2 + ( 200.- wsbio2 ) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio
      wscal (:,:,:) = wsbio4(:,:,:)
      !
      ! OA This is (I hope) a temporary solution for the problem that may 
      ! OA arise in specific situation where the CFL criterion is broken 
      ! OA for vertical sedimentation of particles. To avoid this, a time
      ! OA splitting algorithm has been coded. A specific maximum
      ! OA iteration number is provided and may be specified in the namelist 
      ! OA This is to avoid very large iteration number when explicit free
      ! OA surface is used (for instance). When niter?max is set to 1, 
      ! OA this computation is skipped. The crude old threshold method is 
      ! OA then applied. This also happens when niter exceeds nitermax.
      IF( MAX( niter1max, niter2max ) == 1 ) THEN
        iiter1 = 1
        iiter2 = 1
      ELSE
        iiter1 = 1
        iiter2 = 1
        DO jk = 1, jpkm1
          DO jj = 1, jpj
             DO ji = 1, jpi
                IF( tmask(ji,jj,jk) == 1) THEN
                   zwsmax =  0.5 * e3t_n(ji,jj,jk) / xstep
                   iiter1 =  MAX( iiter1, INT( wsbio3(ji,jj,jk) / zwsmax ) )
                   iiter2 =  MAX( iiter2, INT( wsbio4(ji,jj,jk) / zwsmax ) )
                ENDIF
             END DO
          END DO
        END DO
        IF( lk_mpp ) THEN
           CALL mpp_max( iiter1 )
           CALL mpp_max( iiter2 )
        ENDIF
        iiter1 = MIN( iiter1, niter1max )
        iiter2 = MIN( iiter2, niter2max )
      ENDIF

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 1 ) THEN
                 zwsmax = 0.5 * e3t_n(ji,jj,jk) / xstep
                 wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax * FLOAT( iiter1 ) )
                 wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax * FLOAT( iiter2 ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------
      DO jit = 1, iiter1
        CALL p4z_sink2( wsbio3, sinking , jppoc, iiter1 )
        CALL p4z_sink2( wsbio3, sinkfer , jpsfe, iiter1 )
      END DO

      DO jit = 1, iiter2
        CALL p4z_sink2( wsbio4, sinking2, jpgoc, iiter2 )
        CALL p4z_sink2( wsbio4, sinkfer2, jpbfe, iiter2 )
        CALL p4z_sink2( wsbio4, sinksil , jpgsi, iiter2 )
        CALL p4z_sink2( wscal , sinkcal , jpcal, iiter2 )
      END DO

      !  Exchange between organic matter compartments due to coagulation/disaggregation
      !  ---------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zstep = xstep 
               zfact = zstep * xdiss(ji,jj,jk)
               !  Part I : Coagulation dependent on turbulence
               zagg1 = 25.9  * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)
               zagg2 = 4452. * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)

               ! Part II : Differential settling

               !  Aggregation of small into large particles
               zagg3 =  47.1 * zstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)
               zagg4 =  3.3  * zstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)

               zagg   = zagg1 + zagg2 + zagg3 + zagg4
               zaggfe = zagg * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn )

               ! Aggregation of DOC to POC : 
               ! 1st term is shear aggregation of DOC-DOC
               ! 2nd term is shear aggregation of DOC-POC
               ! 3rd term is differential settling of DOC-POC
               zaggdoc  = ( ( 0.369 * 0.3 * trb(ji,jj,jk,jpdoc) + 102.4 * trb(ji,jj,jk,jppoc) ) * zfact       &
               &            + 2.4 * zstep * trb(ji,jj,jk,jppoc) ) * 0.3 * trb(ji,jj,jk,jpdoc)
               ! transfer of DOC to GOC : 
               ! 1st term is shear aggregation
               ! 2nd term is differential settling 
               zaggdoc2 = ( 3.53E3 * zfact + 0.1 * zstep ) * trb(ji,jj,jk,jpgoc) * 0.3 * trb(ji,jj,jk,jpdoc)
               ! tranfer of DOC to POC due to brownian motion
               zaggdoc3 =  ( 5095. * trb(ji,jj,jk,jppoc) + 114. * 0.3 * trb(ji,jj,jk,jpdoc) ) *zstep * 0.3 * trb(ji,jj,jk,jpdoc)

               !  Update the trends
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zagg + zaggdoc + zaggdoc3
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zagg + zaggdoc2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
               !
            END DO
         END DO
      END DO


     ! Total carbon export per year
     IF( iom_use( "tcexp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
        &   t_oce_co2_exp = glob_sum( ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj,      zw2d )
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "tcexp" ) )  CALL iom_put( "tcexp" , t_oce_co2_exp * zfact )   ! molC/s
          ! 
          CALL wrk_dealloc( jpi, jpj,      zw2d )
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
        ENDIF
      ELSE
         IF( ln_diatrc ) THEN
            zfact = 1.e3 * rfact2r
            trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 5) = sinking2(:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 6) = sinkfer (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 7) = sinkfer2(:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 8) = sinksil (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_pcs0_2d + 9) = sinkcal (:,:,ik100) * zfact * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk

      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
   END SUBROUTINE p4z_sink_init


   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra, kiter )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      INTEGER , INTENT(in   )                         ::   kiter     ! number of iterations for time-splitting 
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )

      zstep = rfact2 / FLOAT( kiter ) / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trb(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0
      IF( lk_degrad ) THEN
         zwsink2(:,:,:) = zwsink2(:,:,:) * facvol(:,:,:)
      ENDIF


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trb(:,:,jk-1,jp_tra) - trb(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / e3w_n(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trb(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
                  trb(ji,jj,jk,jp_tra) = trb(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      trb(:,:,:,jp_tra) = ztrb(:,:,:)
      psinkflx(:,:,:)   = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4  (jpi,jpj,jpk) , wscal(jpi,jpj,jpk) ,     &
         &      sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                      ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                      ,     &                
         &      sinkfer2(jpi,jpj,jpk)                                             ,     &                
         &      sinkfer(jpi,jpj,jpk)                                              , STAT=p4z_sink_alloc )                
         !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   

   !!======================================================================
END MODULE p4zsink
