MODULE cansink
   !!======================================================================
   !!                         ***  MODULE cansink  ***
   !! TOP :  CanOE  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   can_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   can_sink_init  :  Unitialisation of sinking speed parameters
   !!   can_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_sink         ! called in canbio.F90
   PUBLIC   can_sink_init    ! called in trcsms_canoe.F90
   PUBLIC   can_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal    !: CaCO3 and BSi sinking fluxes

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
   !! $Id: cansink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE can_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zfact,  zstep
      REAL(wp) ::   zrfact2
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      !DO jk = 1, jpkm1
      !   DO jj = 1, jpj
      !      DO ji = 1,jpi
      !         zmax  = MAX( heup(ji,jj), hmld(ji,jj) )
      !         zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / 5000._wp
      !         wsbio4(ji,jj,jk) = wsbio2 + ( 200.- wsbio2 ) * zfact
      !      END DO
      !   END DO
      !END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      !wsbio3(:,:,:) = wsbio
      !wscal (:,:,:) = wsbio4(:,:,:)
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

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0

      !   Compute the sedimentation term using cansink2 for all the sinking particles
  !   -----------------------------------------------------
      CALL can_sink2( wsbio3, sinking , jppoc )
      CALL can_sink2( wsbio4, sinking2, jpgoc )
      CALL can_sink2( wscal , sinkcal , jpcal )

    
     !
     ! IF( ln_diatrc ) THEN
     !    zrfact2 = 1.e-3 * rfact2r
     !    ik1  = iksed + 1
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj,      zw2d )
          zfact = 1.e-3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          !
          CALL wrk_dealloc( jpi, jpj,      zw2d )
        ENDIF
      ELSE
         IF( ln_diatrc ) THEN
            zfact = 1.e-3 * rfact2r
            trc2d(:,:,jp_can0_2d + 4) = sinking (:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_can0_2d + 5) = sinking2(:,:,ik100) * zfact * tmask(:,:,1)
            trc2d(:,:,jp_can0_2d + 9) = sinkcal (:,:,ik100) * zfact * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_sink')
      !
   END SUBROUTINE can_sink

   SUBROUTINE can_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk

      ik100 = 10        !  last level where depth less than 100 m
      !DO jk = jpkm1, 1, -1
      !   IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      !END DO
      !IF (lwp) WRITE(numout,*)
      !IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      !IF (lwp) WRITE(numout,*)
      !
      wsbio3(:,:,:) = wsbio
      wsbio4(:,:,:) = wsbio2
      wscal(:,:,:) = wsbioc
      !
   END SUBROUTINE can_sink_init


   SUBROUTINE can_sink2( pwsink, psinkflx, jp_tra )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_sink2  ***
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
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk,   zwsink2, ztrb )

      zstep = rfact2 /  2.

      ztrb (:,:,:) = trb(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0


      ! Vertical advective flux
      DO jk = 1, jpkm1
         DO jj = 1, jpj      
            DO ji = 1, jpi    
               zigma = zwsink2(ji,jj,jk+1) * zstep / e3w_n(ji,jj,jk+1)
               zew   = zwsink2(ji,jj,jk+1)
               psinkflx(ji,jj,jk+1) = -zew * trb(ji,jj,jk,jp_tra) * zstep
            END DO
         END DO
      END DO
      ! Boundary conditions
      psinkflx(:,:,1  ) = 0.e0
      psinkflx(:,:,jpk) = 0.e0

         

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
      CALL wrk_dealloc( jpi, jpj, jpk, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_sink2')
      !
   END SUBROUTINE can_sink2


   INTEGER FUNCTION can_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4  (jpi,jpj,jpk) , wscal(jpi,jpj,jpk) ,     &
         &      sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                      ,     &                
         &      sinkcal(jpi,jpj,jpk)   , STAT=can_sink_alloc )                
         !
      IF( can_sink_alloc /= 0 ) CALL ctl_warn('can_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION can_sink_alloc
   

   !!======================================================================
END MODULE  cansink
