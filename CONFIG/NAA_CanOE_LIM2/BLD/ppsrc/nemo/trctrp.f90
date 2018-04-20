MODULE trctrp
   !!======================================================================
   !!                       ***  MODULE trctrp  ***
   !! Ocean Physics    : manage the passive tracer transport
   !!======================================================================
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             3.3  !  2010-07 (C. Ethe) Merge TRA-TRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_trp        : passive tracer transport
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables 
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE trabbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE zdfkpp          ! KPP non-local tracer fluxes         (trc_kpp routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)
   USE trcldf          ! lateral mixing                      (trc_ldf routine)
   USE trcadv          ! advection                           (trc_adv routine)
   USE trczdf          ! vertical diffusion                  (trc_zdf routine)
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)
   USE zpshde          ! partial step: hor. derivative       (zps_hde routine)
 !! written by xiaofan Luo

    USE trcbdy          ! BDY open boundaries
    USE bdy_par, only: lk_bdy

  !!  =========xiaofanLuo





   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_trp    ! called by trc_stp

   !! * Substitutions
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
   !! $Id: trctrp.F90 5120 2015-03-03 16:11:55Z acc $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_trp  ***
      !!                      
      !! ** Purpose :   Management of passive tracers transport
      !! 
      !! ** Method  : - Compute the passive tracers trends 
      !!              - Update the passive tracers
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kstp  ! ocean time-step index
      !! ---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_trp')
      !
      IF( .NOT. lk_c1d ) THEN
         !
                                CALL trc_sbc( kstp )            ! surface boundary condition
         IF( lk_trabbl )        CALL trc_bbl( kstp )            ! advective (and/or diffusive) bottom boundary layer scheme
         IF( ln_trcdmp )        CALL trc_dmp( kstp )            ! internal damping trends
      
        IF( ln_trcdmp_clo )    CALL trc_dmp_clo( kstp )        ! internal damping trends on closed seas only
  !! written by xiaofanLuo
 #if defined 1
                               CALL trc_bdy_dmp( kstp)         !BDY damping trends 
 #endif  
  !1 xiaofanLuo
                                CALL trc_adv( kstp )            ! horizontal & vertical advection 
                                CALL trc_ldf( kstp )            ! lateral mixing
         IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
            &                   CALL trc_kpp( kstp )            ! KPP non-local tracer fluxes
                                CALL trc_zdf( kstp )            ! vertical mixing and after tracer fields
                                CALL trc_nxt( kstp )            ! tracer fields at next time step     
         IF( ln_trcrad )        CALL trc_rad( kstp )            ! Correct artificial negative concentrations


         IF( ln_zps  .AND. .NOT. ln_isfcav)        &
            &            CALL zps_hde    ( kstp, jptra, trn, gtru, gtrv )   ! Partial steps: now horizontal gradient of passive
         IF( ln_zps .AND.        ln_isfcav)        &
            &            CALL zps_hde_isf( kstp, jptra, trn, pgtu=gtru, pgtv=gtrv, pgtui=gtrui, pgtvi=gtrvi )  ! Partial steps: now horizontal gradient of passive
                                                                ! tracers at the bottom ocean level
         !
      ELSE                                               ! 1D vertical configuration
                                CALL trc_sbc( kstp )            ! surface boundary condition
         IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
            &                   CALL trc_kpp( kstp )            ! KPP non-local tracer fluxes
                                CALL trc_zdf( kstp )            ! vertical mixing and after tracer fields
                                CALL trc_nxt( kstp )            ! tracer fields at next time step     
          IF( ln_trcrad )       CALL trc_rad( kstp )            ! Correct artificial negative concentrations
         !
      END IF
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_trp')
      !
   END SUBROUTINE trc_trp

   
   !!======================================================================
END MODULE trctrp
