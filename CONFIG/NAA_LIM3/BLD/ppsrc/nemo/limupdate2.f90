MODULE limupdate2
   !!======================================================================
   !!                     ***  MODULE  limupdate2  ***
   !!   LIM-3 : Update of sea-ice global variables at the end of the time step
   !!======================================================================
   !! History :  3.0  !  2006-04  (M. Vancoppenolle) Original code
   !!            3.5  !  2014-06  (C. Rousset)       Complete rewriting/cleaning
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!    lim_update2   : computes update of sea-ice global variables from trend terms
   !!----------------------------------------------------------------------
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE dom_ice
   USE dom_oce
   USE phycst          ! physical constants
   USE ice
   USE thd_ice         ! LIM thermodynamic sea-ice variables
   USE limitd_th
   USE limvar
   USE prtctl          ! Print control
   USE lbclnk          ! lateral boundary condition - MPP exchanges
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE limcons         ! conservation tests
   USE limctl
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_update2   ! routine called by ice_step

   !! * Substitutions
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
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limupdate2.F90 7814 2017-03-20 16:21:42Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_update2( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE lim_update2  ***
      !!               
      !! ** Purpose :  Computes update of sea-ice global variables at 
      !!               the end of the time step.
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zsal
      REAL(wp) ::   zvi_b, zsmv_b, zei_b, zfs_b, zfw_b, zft_b 
      !!-------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('limupdate2')

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)''
         WRITE(numout,*)' lim_update2 '
         WRITE(numout,*)' ~~~~~~~~~~~ '
      ENDIF

      ! conservation test
      IF( ln_limdiahsb ) CALL lim_cons_hsm(0, 'limupdate2', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      !----------------------------------------------------------------------
      ! Constrain the thickness of the smallest category above himin
      !----------------------------------------------------------------------
      DO jj = 1, jpj 
         DO ji = 1, jpi
            rswitch = MAX( 0._wp , SIGN( 1._wp, a_i(ji,jj,1) - epsi20 ) )   !0 if no ice and 1 if yes
            ht_i(ji,jj,1) = v_i (ji,jj,1) / MAX( a_i(ji,jj,1) , epsi20 ) * rswitch
            IF( v_i(ji,jj,1) > 0._wp .AND. ht_i(ji,jj,1) < rn_himin ) THEN
               a_i (ji,jj,1) = a_i (ji,jj,1) * ht_i(ji,jj,1) / rn_himin
               oa_i(ji,jj,1) = oa_i(ji,jj,1) * ht_i(ji,jj,1) / rn_himin
            ENDIF
         END DO
      END DO
      
      !-----------------------------------------------------
      ! ice concentration should not exceed amax 
      !-----------------------------------------------------
      at_i(:,:) = 0._wp
      DO jl = 1, jpl
         at_i(:,:) = a_i(:,:,jl) + at_i(:,:)
      END DO

      DO jl  = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( at_i(ji,jj) > rn_amax_2d(ji,jj) .AND. a_i(ji,jj,jl) > 0._wp ) THEN
                  a_i (ji,jj,jl) = a_i (ji,jj,jl) * ( 1._wp - ( 1._wp - rn_amax_2d(ji,jj) / at_i(ji,jj) ) )
                  oa_i(ji,jj,jl) = oa_i(ji,jj,jl) * ( 1._wp - ( 1._wp - rn_amax_2d(ji,jj) / at_i(ji,jj) ) )
               ENDIF
            END DO
         END DO
      END DO

      !---------------------
      ! Ice salinity
      !---------------------
      IF (  nn_icesal == 2  ) THEN 
         DO jl = 1, jpl
            DO jj = 1, jpj 
               DO ji = 1, jpi
                  zsal            = smv_i(ji,jj,jl)
                  ! salinity stays in bounds
                  rswitch         = 1._wp - MAX( 0._wp, SIGN( 1._wp, - v_i(ji,jj,jl) ) )
                  smv_i(ji,jj,jl) = rswitch * MAX( MIN( rn_simax * v_i(ji,jj,jl), smv_i(ji,jj,jl) ), rn_simin * v_i(ji,jj,jl) )
                  ! associated salt flux
                  sfx_res(ji,jj) = sfx_res(ji,jj) - ( smv_i(ji,jj,jl) - zsal ) * rhoic * r1_rdtice
               END DO
            END DO
         END DO
      ENDIF

      !----------------------------------------------------
      ! Rebin categories with thickness out of bounds
      !----------------------------------------------------
      IF ( jpl > 1 ) CALL lim_itd_th_reb( 1, jpl )

      !-----------------
      ! zap small values
      !-----------------
      CALL lim_var_zapsmall

      !------------------------------------------------------------------------------
      ! Corrections to avoid wrong values                                        |
      !------------------------------------------------------------------------------
      ! Ice drift
      !------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            IF ( at_i(ji,jj) == 0._wp ) THEN ! what to do if there is no ice
               IF ( at_i(ji+1,jj) == 0._wp ) u_ice(ji,jj)   = 0._wp ! right side
               IF ( at_i(ji-1,jj) == 0._wp ) u_ice(ji-1,jj) = 0._wp ! left side
               IF ( at_i(ji,jj+1) == 0._wp ) v_ice(ji,jj)   = 0._wp ! upper side
               IF ( at_i(ji,jj-1) == 0._wp ) v_ice(ji,jj-1) = 0._wp ! bottom side
            ENDIF
         END DO
      END DO
      !lateral boundary conditions
      CALL lbc_lnk( u_ice(:,:), 'U', -1. )
      CALL lbc_lnk( v_ice(:,:), 'V', -1. )
      !mask velocities
      u_ice(:,:) = u_ice(:,:) * umask(:,:,1)
      v_ice(:,:) = v_ice(:,:) * vmask(:,:,1)
 
      ! -------------------------------------------------
      ! Diagnostics
      ! -------------------------------------------------
      DO jl  = 1, jpl
         oa_i(:,:,jl) = oa_i(:,:,jl) + a_i(:,:,jl) * rdt_ice / rday   ! ice natural aging
         afx_thd(:,:) = afx_thd(:,:) + ( a_i(:,:,jl) - a_i_b(:,:,jl) ) * r1_rdtice
      END DO
      afx_tot = afx_thd + afx_dyn

      DO jj = 1, jpj
         DO ji = 1, jpi            
            ! heat content variation (W.m-2)
            diag_heat(ji,jj) = diag_heat(ji,jj) -  &
               &               ( SUM( e_i(ji,jj,1:nlay_i,:) - e_i_b(ji,jj,1:nlay_i,:) ) +  & 
               &                 SUM( e_s(ji,jj,1:nlay_s,:) - e_s_b(ji,jj,1:nlay_s,:) )    &
               &               ) * r1_rdtice   
            ! salt, volume
            diag_smvi(ji,jj) = diag_smvi(ji,jj) + SUM( smv_i(ji,jj,:) - smv_i_b(ji,jj,:) ) * rhoic * r1_rdtice
            diag_vice(ji,jj) = diag_vice(ji,jj) + SUM( v_i  (ji,jj,:) - v_i_b  (ji,jj,:) ) * rhoic * r1_rdtice
            diag_vsnw(ji,jj) = diag_vsnw(ji,jj) + SUM( v_s  (ji,jj,:) - v_s_b  (ji,jj,:) ) * rhosn * r1_rdtice
         END DO
      END DO

      ! conservation test
      IF( ln_limdiahsb ) CALL lim_cons_hsm(1, 'limupdate2', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      ! necessary calls (at least for coupling)
      CALL lim_var_glo2eqv
      CALL lim_var_agg(2)

      ! -------------------------------------------------
      ! control prints
      ! -------------------------------------------------
      IF( ln_icectl )   CALL lim_prt( kt, iiceprt, jiceprt, 2, ' - Final state - ' )   ! control print

      IF(ln_ctl) THEN   ! Control print
         CALL prt_ctl_info(' ')
         CALL prt_ctl_info(' - Cell values : ')
         CALL prt_ctl_info('   ~~~~~~~~~~~~~ ')
         CALL prt_ctl(tab2d_1=e12t       , clinfo1=' lim_update2  : cell area   :')
         CALL prt_ctl(tab2d_1=at_i       , clinfo1=' lim_update2  : at_i        :')
         CALL prt_ctl(tab2d_1=vt_i       , clinfo1=' lim_update2  : vt_i        :')
         CALL prt_ctl(tab2d_1=vt_s       , clinfo1=' lim_update2  : vt_s        :')
         CALL prt_ctl(tab2d_1=strength   , clinfo1=' lim_update2  : strength    :')
         CALL prt_ctl(tab2d_1=u_ice      , clinfo1=' lim_update2  : u_ice       :', tab2d_2=v_ice      , clinfo2=' v_ice       :')
         CALL prt_ctl(tab2d_1=u_ice_b    , clinfo1=' lim_update2  : u_ice_b     :', tab2d_2=v_ice_b    , clinfo2=' v_ice_b     :')

         DO jl = 1, jpl
            CALL prt_ctl_info(' ')
            CALL prt_ctl_info(' - Category : ', ivar1=jl)
            CALL prt_ctl_info('   ~~~~~~~~~~')
            CALL prt_ctl(tab2d_1=ht_i       (:,:,jl)        , clinfo1= ' lim_update2  : ht_i        : ')
            CALL prt_ctl(tab2d_1=ht_s       (:,:,jl)        , clinfo1= ' lim_update2  : ht_s        : ')
            CALL prt_ctl(tab2d_1=t_su       (:,:,jl)        , clinfo1= ' lim_update2  : t_su        : ')
            CALL prt_ctl(tab2d_1=t_s        (:,:,1,jl)      , clinfo1= ' lim_update2  : t_snow      : ')
            CALL prt_ctl(tab2d_1=sm_i       (:,:,jl)        , clinfo1= ' lim_update2  : sm_i        : ')
            CALL prt_ctl(tab2d_1=o_i        (:,:,jl)        , clinfo1= ' lim_update2  : o_i         : ')
            CALL prt_ctl(tab2d_1=a_i        (:,:,jl)        , clinfo1= ' lim_update2  : a_i         : ')
            CALL prt_ctl(tab2d_1=a_i_b      (:,:,jl)        , clinfo1= ' lim_update2  : a_i_b       : ')
            CALL prt_ctl(tab2d_1=v_i        (:,:,jl)        , clinfo1= ' lim_update2  : v_i         : ')
            CALL prt_ctl(tab2d_1=v_i_b      (:,:,jl)        , clinfo1= ' lim_update2  : v_i_b       : ')
            CALL prt_ctl(tab2d_1=v_s        (:,:,jl)        , clinfo1= ' lim_update2  : v_s         : ')
            CALL prt_ctl(tab2d_1=v_s_b      (:,:,jl)        , clinfo1= ' lim_update2  : v_s_b       : ')
            CALL prt_ctl(tab2d_1=e_i        (:,:,1,jl)      , clinfo1= ' lim_update2  : e_i1        : ')
            CALL prt_ctl(tab2d_1=e_i_b      (:,:,1,jl)      , clinfo1= ' lim_update2  : e_i1_b      : ')
            CALL prt_ctl(tab2d_1=e_i        (:,:,2,jl)      , clinfo1= ' lim_update2  : e_i2        : ')
            CALL prt_ctl(tab2d_1=e_i_b      (:,:,2,jl)      , clinfo1= ' lim_update2  : e_i2_b      : ')
            CALL prt_ctl(tab2d_1=e_s        (:,:,1,jl)      , clinfo1= ' lim_update2  : e_snow      : ')
            CALL prt_ctl(tab2d_1=e_s_b      (:,:,1,jl)      , clinfo1= ' lim_update2  : e_snow_b    : ')
            CALL prt_ctl(tab2d_1=smv_i      (:,:,jl)        , clinfo1= ' lim_update2  : smv_i       : ')
            CALL prt_ctl(tab2d_1=smv_i_b    (:,:,jl)        , clinfo1= ' lim_update2  : smv_i_b     : ')
            CALL prt_ctl(tab2d_1=oa_i       (:,:,jl)        , clinfo1= ' lim_update2  : oa_i        : ')
            CALL prt_ctl(tab2d_1=oa_i_b     (:,:,jl)        , clinfo1= ' lim_update2  : oa_i_b      : ')

            DO jk = 1, nlay_i
               CALL prt_ctl_info(' - Layer : ', ivar1=jk)
               CALL prt_ctl(tab2d_1=t_i(:,:,jk,jl) , clinfo1= ' lim_update2  : t_i       : ')
            END DO
         END DO

         CALL prt_ctl_info(' ')
         CALL prt_ctl_info(' - Heat / FW fluxes : ')
         CALL prt_ctl_info('   ~~~~~~~~~~~~~~~~~~ ')
         CALL prt_ctl(tab2d_1=sst_m  , clinfo1= ' lim_update2 : sst   : ', tab2d_2=sss_m     , clinfo2= ' sss       : ')

         CALL prt_ctl_info(' ')
         CALL prt_ctl_info(' - Stresses : ')
         CALL prt_ctl_info('   ~~~~~~~~~~ ')
         CALL prt_ctl(tab2d_1=utau       , clinfo1= ' lim_update2 : utau      : ', tab2d_2=vtau       , clinfo2= ' vtau      : ')
         CALL prt_ctl(tab2d_1=utau_ice   , clinfo1= ' lim_update2 : utau_ice  : ', tab2d_2=vtau_ice   , clinfo2= ' vtau_ice  : ')
         CALL prt_ctl(tab2d_1=u_oce      , clinfo1= ' lim_update2 : u_oce     : ', tab2d_2=v_oce      , clinfo2= ' v_oce     : ')
      ENDIF
   
      IF( nn_timing == 1 )  CALL timing_stop('limupdate2')

   END SUBROUTINE lim_update2

END MODULE limupdate2
