MODULE traadv_eiv
   !!======================================================================
   !!                    ***  MODULE  traadv_eiv  ***
   !! Ocean tracers:  advection trend - eddy induced velocity
   !!======================================================================
   !! History :  1.0  !  2005-11 (G. Madec)  Original code, from traldf and zdf _iso
   !!            3.3  !  2010-05 (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso : update the tracer trend with the horizontal component
   !!                 of iso neutral laplacian operator or horizontal 
   !!                 laplacian operator in s-coordinate
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfslp          ! iso-neutral slopes
   USE in_out_manager  ! I/O manager
   USE iom
   USE trc_oce         ! share passive tracers/Ocean variables

   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE diaptr         ! Heat/Salt transport diagnostics
   USE trddyn
   USE trd_oce

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_eiv   ! routine called by step.F90

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

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
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
!   'key_traldf_c2d' :                 aht: 2D coefficient
   !!----------------------------------------------------------------------
   !!                   ***  ldfeiv_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaei. the eddy induced velocity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'traldf_c2d' :                           eiv: 2D coefficient
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
   !! $Id: traadv_eiv.F90 8627 2017-10-16 14:19:11Z gm $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_eiv( kt, kit000, pun, pvn, pwn, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_eiv  ***
      !! 
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The eddy induced advection is computed from the slope
      !!      of iso-neutral surfaces computed in routine ldf_slp as follows:
      !!         zu_eiv =  1/(e2u e3u)   dk[ aeiu e2u mi(wslpi) ]
      !!         zv_eiv =  1/(e1v e3v)   dk[ aeiv e1v mj(wslpj)
      !!         zw_eiv = -1/(e1t e2t) { di[ aeiu e2u mi(wslpi) ]
      !!                               + dj[ aeiv e1v mj(wslpj) ] }
      !!      add the eiv component to the model velocity:
      !!         p.n = p.n + z._eiv
      !!
      !! ** Action  : - add to p.n the eiv component
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pun      ! in : 3 ocean velocity components 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pvn      ! out: 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pwn      ! increased by the eiv
      !!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zuwk, zuwk1, zuwi, zuwi1   ! local scalars
      REAL(wp) ::   zvwk, zvwk1, zvwj, zvwj1   !   -      -
      REAL(wp) ::   zztmp                      ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_eiv, zv_eiv, zw_eiv, z2d
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d, z3d_T
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_adv_eiv')
      !
      CALL wrk_alloc( jpi, jpj, zu_eiv, zv_eiv, zw_eiv, z2d )
      CALL wrk_alloc( jpi, jpj, jpk, z3d, z3d_T )

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_eiv : eddy induced advection on ', cdtype,' :'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   add to velocity fields the eiv component'
         IF( cdtype == 'TRA') THEN
            u_eiv(:,:,:) = 0.e0
            v_eiv(:,:,:) = 0.e0
            w_eiv(:,:,:) = 0.e0
         END IF
      ENDIF

      zu_eiv(:,:) = 0.e0   ;   zv_eiv(:,:) = 0.e0   ;    zw_eiv(:,:) = 0.e0  
      
                                                    ! =================
      DO jk = 1, jpkm1                              !  Horizontal slab
         !                                          ! =================
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk  ) ) * aeiu(ji,jj) * umask(ji,jj,jk  )
               zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) ) * aeiu(ji,jj) * umask(ji,jj,jk+1)
               zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk  ) ) * aeiv(ji,jj) * vmask(ji,jj,jk  )
               zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) ) * aeiv(ji,jj) * vmask(ji,jj,jk+1)

               zu_eiv(ji,jj) = 0.5 * umask(ji,jj,jk) * ( zuwk - zuwk1 ) 
               zv_eiv(ji,jj) = 0.5 * vmask(ji,jj,jk) * ( zvwk - zvwk1 ) 
   
               pun(ji,jj,jk) = pun(ji,jj,jk) + e2u(ji,jj) * zu_eiv(ji,jj)
               pvn(ji,jj,jk) = pvn(ji,jj,jk) + e1v(ji,jj) * zv_eiv(ji,jj)
            END DO
         END DO
         IF( cdtype == 'TRA') THEN
            u_eiv(:,:,jk) = zu_eiv(:,:) / e3u_0(:,:,jk)
            v_eiv(:,:,jk) = zv_eiv(:,:) / e3v_0(:,:,jk)
         END IF
         IF( jk >=2 ) THEN                             ! jk=1 zw_eiv=0, not computed
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zuwi  = ( wslpi(ji,jj,jk)+wslpi(ji-1,jj,jk) ) * aeiu(ji-1,jj) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwi1 = ( wslpi(ji,jj,jk)+wslpi(ji+1,jj,jk) ) * aeiu(ji  ,jj) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwj  = ( wslpj(ji,jj,jk)+wslpj(ji,jj-1,jk) ) * aeiv(ji,jj-1) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwj1 = ( wslpj(ji,jj,jk)+wslpj(ji,jj+1,jk) ) * aeiv(ji,jj  ) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
  
                  zw_eiv(ji,jj) = - 0.5 * tmask(ji,jj,jk) * ( zuwi1 - zuwi + zvwj1 - zvwj ) 
                  pwn(ji,jj,jk) = pwn(ji,jj,jk) + zw_eiv(ji,jj)
               END DO
            END DO
            IF( cdtype == 'TRA')  w_eiv(:,:,jk) = zw_eiv(:,:) / ( e1t(:,:) * e2t(:,:) )
         ENDIF
         !                                          ! =================
      END DO                                        !    End of slab  
      !                                             ! =================

      IF( cdtype == 'TRA') THEN
         CALL iom_put( "uoce_eiv", u_eiv )    ! i-eiv current
         CALL iom_put( "voce_eiv", v_eiv )    ! j-eiv current
         CALL iom_put( "woce_eiv", w_eiv )    ! vert. eiv current
         !
         IF( iom_use('weiv_masstr') ) THEN   ! vertical mass transport & its square value
            z2d(:,:) = rau0 * e12t(:,:)
            DO jk = 1, jpk
               z3d(:,:,jk) = w_eiv(:,:,jk) * z2d(:,:)
            END DO
            CALL iom_put( "weiv_masstr" , z3d )  
         ENDIF
         !
         IF( iom_use("ueiv_masstr") .OR. iom_use("ueiv_heattr") .OR. iom_use('ueiv_heattr3d')        &
            &                       .OR. iom_use("ueiv_salttr") .OR. iom_use('ueiv_salttr3d') ) THEN
            z3d(:,:,jpk) = 0.e0
            z2d(:,:) = 0.e0
            DO jk = 1, jpkm1
               z3d(:,:,jk) = rau0 * u_eiv(:,:,jk) * e2u(:,:) * e3u_0(:,:,jk) * umask(:,:,jk)
               z2d(:,:) = z2d(:,:) + z3d(:,:,jk)
            END DO
            CALL iom_put( "ueiv_masstr", z3d )                  ! mass transport in i-direction
         ENDIF
         !
         IF( iom_use('ueiv_heattr') .OR. iom_use('ueiv_heattr3d') ) THEN
            zztmp = 0.5 * rcp 
            z2d(:,:) = 0.e0 
            z3d_T(:,:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
                     z2d(ji,jj) = z2d(ji,jj) + z3d_T(ji,jj,jk) 
                  END DO
               END DO
            END DO
            IF (iom_use('ueiv_heattr') ) THEN
               CALL lbc_lnk( z2d, 'U', -1. )
               CALL iom_put( "ueiv_heattr", zztmp * z2d )                  ! 2D heat transport in i-direction
            ENDIF
            IF (iom_use('ueiv_heattr3d') ) THEN
               CALL lbc_lnk( z3d_T, 'U', -1. )
               CALL iom_put( "ueiv_heattr3d", zztmp * z3d_T )              ! 3D heat transport in i-direction
            ENDIF
         ENDIF
         !
         IF( iom_use('ueiv_salttr') .OR. iom_use('ueiv_salttr3d') ) THEN
            zztmp = 0.5 * 0.001
            z2d(:,:) = 0.e0 
            z3d_T(:,:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
                     z2d(ji,jj) = z2d(ji,jj) + z3d_T(ji,jj,jk) 
                  END DO
               END DO
            END DO
            IF (iom_use('ueiv_salttr') ) THEN
               CALL lbc_lnk( z2d, 'U', -1. )
               CALL iom_put( "ueiv_salttr", zztmp * z2d )                  ! 2D salt transport in i-direction
            ENDIF
            IF (iom_use('ueiv_salttr3d') ) THEN
               CALL lbc_lnk( z3d_T, 'U', -1. )
               CALL iom_put( "ueiv_salttr3d", zztmp * z3d_T )              ! 3D salt transport in i-direction
            ENDIF
         ENDIF
         !
         IF( iom_use("veiv_masstr") .OR. iom_use("veiv_heattr") .OR. iom_use('veiv_heattr3d')       &
                                    .OR. iom_use("veiv_salttr") .OR. iom_use('veiv_salttr3d') ) THEN
            z3d(:,:,jpk) = 0.e0
            DO jk = 1, jpkm1
               z3d(:,:,jk) = rau0 * v_eiv(:,:,jk) * e1v(:,:) * e3v_0(:,:,jk) * vmask(:,:,jk)
            END DO
            CALL iom_put( "veiv_masstr", z3d )                  ! mass transport in j-direction
         ENDIF
         !   
         IF( iom_use('veiv_heattr') .OR. iom_use('veiv_heattr3d') ) THEN
            zztmp = 0.5 * rcp 
            z2d(:,:) = 0.e0 
            z3d_T(:,:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
                     z2d(ji,jj) = z2d(ji,jj) + z3d_T(ji,jj,jk) 
                  END DO
               END DO
            END DO
            IF (iom_use('veiv_heattr') ) THEN
               CALL lbc_lnk( z2d, 'V', -1. )
               CALL iom_put( "veiv_heattr", zztmp * z2d )                  ! 2D heat transport in j-direction
            ENDIF
            IF (iom_use('veiv_heattr3d') ) THEN
               CALL lbc_lnk( z3d_T, 'V', -1. )
               CALL iom_put( "veiv_heattr3d", zztmp * z3d_T )              ! 3D heat transport in j-direction
            ENDIF
         ENDIF
         !
         IF( iom_use('veiv_salttr') .OR. iom_use('veiv_salttr3d') ) THEN
            zztmp = 0.5 * 0.001
            z2d(:,:) = 0.e0 
            z3d_T(:,:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
                     z2d(ji,jj) = z2d(ji,jj) + z3d_T(ji,jj,jk)
                  END DO
               END DO
            END DO
            IF (iom_use('veiv_salttr') ) THEN
               CALL lbc_lnk( z2d, 'V', -1. )
               CALL iom_put( "veiv_salttr", zztmp * z2d )                  ! 2D salt transport in i-direction
            ENDIF
            IF (iom_use('veiv_salttr3d') ) THEN
               CALL lbc_lnk( z3d_T, 'V', -1. )
               CALL iom_put( "veiv_salttr3d", zztmp * z3d_T )              ! 3D salt transport in i-direction
            ENDIF
         ENDIF
         !
         IF( iom_use('weiv_masstr') .OR. iom_use('weiv_heattr3d') .OR. iom_use('weiv_salttr3d')) THEN   ! vertical mass transport & its square value
           z2d(:,:) = rau0 * e12t(:,:)
           DO jk = 1, jpk
              z3d(:,:,jk) = w_eiv(:,:,jk) * z2d(:,:)
           END DO
           CALL iom_put( "weiv_masstr" , z3d )                  ! mass transport in k-direction
         ENDIF
         !
         IF( iom_use('weiv_heattr3d') ) THEN
            zztmp = 0.5 * rcp 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj,jk+1,jp_tem) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z3d_T, 'T', 1. )
            CALL iom_put( "weiv_heattr3d", zztmp * z3d_T )                 ! 3D heat transport in k-direction
         ENDIF
         !
         IF( iom_use('weiv_salttr3d') ) THEN
            zztmp = 0.5 * 0.001 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d_T(ji,jj,jk) = z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj,jk+1,jp_sal) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z3d_T, 'T', 1. )
            CALL iom_put( "weiv_salttr3d", zztmp * z3d_T )                 ! 3D salt transport in k-direction
         ENDIF
         !
         IF( ln_diaptr ) THEN
            z3d(:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d(ji,jj,jk) = v_eiv(ji,jj,jk) * 0.5 * (tsn(ji,jj,jk,jp_tem)+tsn(ji,jj+1,jk,jp_tem)) &
                        &             * e1v(ji,jj) * e3v_0(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL dia_ptr_ohst_components( jp_tem, 'eiv', z3d )
            z3d(:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z3d(ji,jj,jk) = v_eiv(ji,jj,jk) * 0.5 * (tsn(ji,jj,jk,jp_sal)+tsn(ji,jj+1,jk,jp_sal)) &
                     &             * e1v(ji,jj) * e3v_0(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL dia_ptr_ohst_components( jp_sal, 'eiv', z3d )
         ENDIF
         !
!!gm add CMIP6 diag here instead of been done in trdken.F90
         !
         IF( iom_use('eketrd_eiv') ) THEN     ! tendency of EKE from parameterized eddy advection
            ! CMIP6 diagnostic tknebto = tendency of EKE from parameterized mesoscale eddy advection
            ! = vertical_integral( k (N S)^2 ) rho dz   where rho = rau0 and S = isoneutral slope.
            z2d(:,:) = 0._wp
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  DO jj = 1,jpj
                     z2d(ji,jj) = z2d(ji,jj) + rau0 * aeiw(ji,jj)                 &
                        &                    * rn2b(ji,jj,jk) * e3w_0(ji,jj,jk)        &
                        &                    * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)    &
                        &                       + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL iom_put( "eketrd_eiv", z2d )
         ENDIF
         !
!!gm  removed from trdken.F90    IF( ln_KE_trd )   CALL trd_dyn(u_eiv, v_eiv, jpdyn_eivke, kt )
         !
      ENDIF

      CALL wrk_dealloc( jpi, jpj, zu_eiv, zv_eiv, zw_eiv, z2d )
      CALL wrk_dealloc( jpi, jpj, jpk, z3d, z3d_T )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_adv_eiv')
      !
    END SUBROUTINE tra_adv_eiv


   !!==============================================================================
END MODULE traadv_eiv
