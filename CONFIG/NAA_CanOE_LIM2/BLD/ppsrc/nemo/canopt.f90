MODULE canopt
   !!======================================================================
   !!                         ***  MODULE canopt  ***
   !! TOP - CanOE : Compute the light availability in the water column
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_canoe     ! Source Minus Sink of CanOE
   USE iom            ! I/O manager
   USE fldread         !  time interpolation
   USE prtctl_trc      !  print control for debugging


   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_opt        ! called in canbio.F90 module
   PUBLIC   can_opt_init   ! called in trcsms_canoe.F90 module
   PUBLIC   can_opt_alloc

   !! * Shared module variables



   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: etot      !: PAR over 24h in case of diurnal cycle
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: emoy           !: averaged PAR in the mixed layer

   INTEGER  ::   nksrp   ! levels below which the light cannot penetrate ( depth larger than 391 m)
   REAL(wp) ::   parlux = 0.43_wp / 3._wp

   REAL(wp), DIMENSION(3,61), PUBLIC ::   xkrgb   !: tabulated attenuation coefficients for RGB absorption
   
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
   !! $Id: canopt.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_opt( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl,zxsi0r
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zdepmoy, zetmp, zetmp2, zetmp1
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zekg, zekr, zekb, ze0, ze1, ze2, ze3
      !REAL(wp), POINTER, DIMENSION(:,:,:) :: zpar, ze0, ze1, ze2, ze3
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_opt')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      zdepmoy, zetmp, zetmp1, zetmp2       )
      CALL wrk_alloc( jpi, jpj, jpk, zekg, zekr, zekb, ze0, ze1, ze2, ze3 )


      !     Initialisation of variables used to compute PAR
      !     -----------------------------------------------
      ze1(:,:,:) = 0._wp
      ze2(:,:,:) = 0._wp
      ze3(:,:,:) = 0._wp
      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      DO jk = 1, jpkm1                         !  --------------------------------------------------------
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zchl = ( trb(ji,jj,jk,jpnch) + trb(ji,jj,jk,jpdch) + rtrn ) 
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !                                                         
               zekb(ji,jj,jk) = xkrgb(1,irgb) * e3t_n(ji,jj,jk)
               zekg(ji,jj,jk) = xkrgb(2,irgb) * e3t_n(ji,jj,jk)
               zekr(ji,jj,jk) = xkrgb(3,irgb) * e3t_n(ji,jj,jk)
            END DO
         END DO
      END DO

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            zc1 = parlux * qsr(ji,jj) * EXP( -0.5 * zekb(ji,jj,1) )
            zc2 = parlux * qsr(ji,jj) * EXP( -0.5 * zekg(ji,jj,1) )
            zc3 = parlux * qsr(ji,jj) * EXP( -0.5 * zekr(ji,jj,1) )
            ze1  (ji,jj,1) = zc1
            ze2  (ji,jj,1) = zc2
            ze3  (ji,jj,1) = zc3
            etot (ji,jj,1) = (       zc1 +        zc2 +       zc3 )
            enano(ji,jj,1) = ( 2.1 * zc1 + 0.42 * zc2 + 0.4 * zc3 )
            ediat(ji,jj,1) = ( 1.6 * zc1 + 0.69 * zc2 + 0.7 * zc3 )
         END DO
      END DO

    
      DO jk = 2, nksrp      
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zc1 = ze1(ji,jj,jk-1) * EXP( -0.5 * ( zekb(ji,jj,jk-1) + zekb(ji,jj,jk) ) )
               zc2 = ze2(ji,jj,jk-1) * EXP( -0.5 * ( zekg(ji,jj,jk-1) + zekg(ji,jj,jk) ) )
               zc3 = ze3(ji,jj,jk-1) * EXP( -0.5 * ( zekr(ji,jj,jk-1) + zekr(ji,jj,jk) ) )
               ze1  (ji,jj,jk) = zc1
               ze2  (ji,jj,jk) = zc2
               ze3  (ji,jj,jk) = zc3
               etot (ji,jj,jk) = (       zc1 +        zc2 +       zc3 )
               enano(ji,jj,jk) = ( 2.1 * zc1 + 0.42 * zc2 + 0.4 * zc3 )
               ediat(ji,jj,jk) = ( 1.6 * zc1 + 0.69 * zc2 + 0.7 * zc3 )
            END DO
         END DO
      END DO

      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         zxsi0r = 1.e0 / rn_si0
         !
         ze0  (:,:,1) = rn_abs * qsr(:,:)
         ze1  (:,:,1) = parlux * qsr(:,:)             ! surface value : separation in R-G-B + near surface 
         ze2  (:,:,1) = parlux * qsr(:,:)
         ze3  (:,:,1) = parlux * qsr(:,:)
         etot3(:,:,1) =          qsr(:,:) * tmask(:,:,1)
         !
         DO jk = 2, nksrp + 1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zc0 = ze0(ji,jj,jk-1) * EXP( -e3t_n(ji,jj,jk-1) * zxsi0r )
                  zc1 = ze1(ji,jj,jk-1) * EXP( -zekb(ji,jj,jk-1 ) )
                  zc2 = ze2(ji,jj,jk-1) * EXP( -zekg(ji,jj,jk-1 ) )
                  zc3 = ze3(ji,jj,jk-1) * EXP( -zekr(ji,jj,jk-1 ) )
                  ze0(ji,jj,jk) = zc0
                  ze1(ji,jj,jk) = zc1
                  ze2(ji,jj,jk) = zc2
                  ze3(ji,jj,jk) = zc3
                  etot3(ji,jj,jk) = ( zc0 + zc1 + zc2 + zc3 ) * tmask(ji,jj,jk)
              END DO
              !
            END DO
            !
        END DO
        !
      ENDIF


      !                                        !* Euphotic depth and level
      neln(:,:) = 1                            !  ------------------------
      heup(:,:) = 300.

      DO jk = 2, nksrp
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( etot(ji,jj,jk) >= 0.0043 * qsr(ji,jj) )  THEN
                 neln(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                     ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heup(ji,jj) = gdepw_n(ji,jj,jk+1)      ! Euphotic layer depth
              ENDIF
           END DO
        END DO
      END DO
 
      heup(:,:) = MIN( 300., heup(:,:) )

      !                                        !* mean light over the mixed layer
      zdepmoy(:,:)   = 0.e0                    !  -------------------------------
      zetmp  (:,:)   = 0.e0
      zetmp1 (:,:)   = 0.e0
      zetmp2 (:,:)   = 0.e0

      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zetmp  (ji,jj) = zetmp  (ji,jj) + etot (ji,jj,jk) * e3t_n(ji,jj,jk)
                  zetmp1 (ji,jj) = zetmp1 (ji,jj) + enano(ji,jj,jk) * e3t_n(ji,jj,jk)
                  zetmp2 (ji,jj) = zetmp2 (ji,jj) + ediat(ji,jj,jk) * e3t_n(ji,jj,jk)
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) + e3t_n(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO
      !
      emoy(:,:,:) = etot(:,:,:)
      !
      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  emoy (ji,jj,jk) = zetmp (ji,jj) * z1_dep
                  enano(ji,jj,jk) = zetmp1(ji,jj) * z1_dep
                  ediat(ji,jj,jk) = zetmp2(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO

      !
      !

      !
      IF( lk_iomput ) THEN
        IF( knt == nrdttrc ) THEN
           IF( iom_use( "Heup"  ) ) CALL iom_put( "Heup" , heup(:,:  ) * tmask(:,:,1) )  ! euphotic layer deptht
           IF( iom_use( "PAR"   ) ) CALL iom_put( "PAR"  , emoy(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
        ENDIF
      ELSE
         IF( ln_diatrc ) THEN        ! save output diagnostics
            trc2d(:,:,  jp_can0_2d + 10) = heup(:,:  ) * tmask(:,:,1)
            trc3d(:,:,:,jp_can0_3d + 3)  = etot(:,:,:) * tmask(:,:,:)
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zdepmoy, zetmp, zetmp1, zetmp2 )
      CALL wrk_dealloc( jpi, jpj, jpk, zekg, zekr, zekb, ze0, ze1, ze2, ze3 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_opt')
      !
   END SUBROUTINE can_opt



   SUBROUTINE can_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      !
      INTEGER :: ierr
      INTEGER :: ios                 ! Local integer output status for namelist read
      !
      !
      CALL trc_oce_rgb( xkrgb )                  ! tabulated attenuation coefficients
      nksrp = trc_oce_ext_lev( r_si2, 0.33e2 )     ! max level of light extinction (Blue Chl=0.01)
      !
      IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksrp, ' ref depth = ', gdepw_1d(nksrp+1), ' m'
      !
                         etot     (:,:,:) = 0._wp
                         enano    (:,:,:) = 0._wp
                         ediat    (:,:,:) = 0._wp
      IF( ln_qsr_bio )   etot3    (:,:,:) = 0._wp
      ! 
      IF( nn_timing == 1 )  CALL timing_stop('can_opt_init')
      !
   END SUBROUTINE can_opt_init


   INTEGER FUNCTION can_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_opt_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE(  enano(jpi,jpj,jpk)    , ediat(jpi,jpj,jpk), &
        &       etot(jpi,jpj,jpk), emoy (jpi,jpj,jpk), STAT=can_opt_alloc ) 
         !
      IF( can_opt_alloc /= 0 ) CALL ctl_warn('can_opt_alloc : failed to allocate arrays.')
      !
   END FUNCTION can_opt_alloc


   !!======================================================================
END MODULE  canopt
