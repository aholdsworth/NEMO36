MODULE dynzdf_imp
   !!======================================================================
   !!                    ***  MODULE  dynzdf_imp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.4  !  2012-01  (H. Liu) Semi-implicit bottom friction
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_imp  : update the momentum trend with the vertical diffusion using a implicit time-stepping
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! variable volume
   USE sbc_oce         ! surface boundary condition: ocean
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE dynadv          ! dynamics: vector invariant versus flux form
   USE dynspg_oce, ONLY: lk_dynspg_ts
   USE zdfbfr          ! Bottom friction setup
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_imp   ! called by step.F90

   REAL(wp) ::  r_vvl     ! variable volume indicator, =1 if lk_vvl=T, =0 otherwise 

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
   !! $Id: dynzdf_imp.F90 8627 2017-10-16 14:19:11Z gm $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_imp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_imp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!      and the surface forcing, and add it to the general trend of 
      !!      the momentum equations.
      !!
      !! ** Method  :   The vertical momentum mixing trend is given by :
      !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
      !!      backward time stepping
      !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) arrays with the after vertical diffusive mixing trend.
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::  kt     ! ocean time-step index
      REAL(wp), INTENT(in) ::  p2dt   ! vertical profile of tracer time-step
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   z1_p2dt, zcoef, zzwi, zzws, zrhs   ! local scalars
      REAL(wp) ::   ze3ua, ze3va, zzz
      REAL(wp), POINTER, DIMENSION(:,:)   ::  z2d
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwi, zwd, zws
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_imp')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwi, zwd, zws ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
         IF( lk_vvl ) THEN   ;    r_vvl = 1._wp       ! Variable volume indicator
         ELSE                ;    r_vvl = 0._wp       
         ENDIF
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      z1_p2dt = 1._wp / p2dt      ! inverse of the timestep

      ! 1. Apply semi-implicit bottom friction
      ! --------------------------------------
      ! Only needed for semi-implicit bottom friction setup. The explicit
      ! bottom friction has been included in "u(v)a" which act as the R.H.S
      ! column vector of the tri-diagonal matrix equation
      !

      IF( ln_bfrimp ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)       ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
               avmu(ji,jj,ikbu+1) = -bfrua(ji,jj) * e3uw_n(ji,jj,ikbu+1)
               avmv(ji,jj,ikbv+1) = -bfrva(ji,jj) * e3vw_n(ji,jj,ikbv+1)
            END DO
         END DO
         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ikbu = miku(ji,jj)       ! ocean top level at u- and v-points 
                  ikbv = mikv(ji,jj)       ! (first wet ocean u- and v-points)
                  IF (ikbu .GE. 2) avmu(ji,jj,ikbu) = -tfrua(ji,jj) * e3uw_n(ji,jj,ikbu)
                  IF (ikbv .GE. 2) avmv(ji,jj,ikbv) = -tfrva(ji,jj) * e3vw_n(ji,jj,ikbv)
               END DO
            END DO
         END IF
      ENDIF

      IF( ln_dynadv_vec .OR. .NOT. lk_vvl ) THEN      ! applied on velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ub(:,:,jk) + p2dt * ua(:,:,jk) ) * umask(:,:,jk)
            va(:,:,jk) = ( vb(:,:,jk) + p2dt * va(:,:,jk) ) * vmask(:,:,jk)
         END DO
      ELSE                                            ! applied on thickness weighted velocity
         DO jk = 1, jpkm1
            ua(:,:,jk) = (          ub(:,:,jk) * e3u_b(:,:,jk)      &
               &           + p2dt * ua(:,:,jk) * e3u_n(:,:,jk)  )   &
               &                               / e3u_a(:,:,jk) * umask(:,:,jk)
            va(:,:,jk) = (          vb(:,:,jk) * e3v_b(:,:,jk)      &
               &           + p2dt * va(:,:,jk) * e3v_n(:,:,jk)  )   &
               &                               / e3v_a(:,:,jk) * vmask(:,:,jk)
         END DO
      ENDIF

      IF ( ln_bfrimp .AND.lk_dynspg_ts ) THEN
         ! remove barotropic velocities:
         DO jk = 1, jpkm1
            ua(:,:,jk) = (ua(:,:,jk) - ua_b(:,:)) * umask(:,:,jk)
            va(:,:,jk) = (va(:,:,jk) - va_b(:,:)) * vmask(:,:,jk)
         END DO
         ! Add bottom/top stress due to barotropic component only:
         DO jj = 2, jpjm1        
            DO ji = 1, jpi   ! vector opt.
               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,ikbu) + r_vvl   * e3u_a(ji,jj,ikbu)
               ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,ikbv) + r_vvl   * e3v_a(ji,jj,ikbv)
               ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * bfrua(ji,jj) * ua_b(ji,jj) / ze3ua
               va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * bfrva(ji,jj) * va_b(ji,jj) / ze3va
            END DO
         END DO
         IF ( ln_isfcav ) THEN
            DO jj = 2, jpjm1        
               DO ji = 1, jpi   ! vector opt.
                  ikbu = miku(ji,jj)         ! top ocean level at u- and v-points 
                  ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                  ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,ikbu) + r_vvl   * e3u_a(ji,jj,ikbu)
                  ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,ikbv) + r_vvl   * e3v_a(ji,jj,ikbv)
                  ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * tfrua(ji,jj) * ua_b(ji,jj) / ze3ua
                  va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * tfrva(ji,jj) * va_b(ji,jj) / ze3va
               END DO
            END DO
         END IF
      ENDIF

      ! 2. Vertical diffusion on u
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmu can take
      ! non zero value at the ocean bottom depending on the bottom friction used.
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1 
            DO ji = 1, jpi   ! vector opt.
               ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,jk) + r_vvl   * e3u_a(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3ua      
               zzwi          = zcoef * avmu  (ji,jj,jk  ) / e3uw_n(ji,jj,jk  )
               zwi(ji,jj,jk) = zzwi  * wumask(ji,jj,jk  )
               zzws          = zcoef * avmu  (ji,jj,jk+1) / e3uw_n(ji,jj,jk+1) 
               zws(ji,jj,jk) = zzws  * wumask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = 1, jpi   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in ua
      !-----------------------------------------------------------------------
      !
      !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = 1, jpi   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = 1, jpi   ! vector opt.
            ze3ua =  ( 1._wp - r_vvl ) * e3u_n(ji,jj,1) + r_vvl   * e3u_a(ji,jj,1) 
            ua(ji,jj,1) = ua(ji,jj,1) + p2dt * 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( ze3ua * rau0 ) * umask(ji,jj,1) 
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 1, jpi
               zrhs = ua(ji,jj,jk)   ! zrhs=right hand side
               ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==
         DO ji = 1, jpi   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = 1, jpi
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      ! 3. Vertical diffusion on v
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmv can take
      ! non zero value at the ocean bottom depending on the bottom friction used
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1   
            DO ji = 1, jpi   ! vector opt.
               ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,jk)  + r_vvl * e3v_a(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3va
               zzwi          = zcoef * avmv (ji,jj,jk  ) / e3vw_n(ji,jj,jk  )
               zwi(ji,jj,jk) =  zzwi * wvmask(ji,jj,jk)
               zzws          = zcoef * avmv (ji,jj,jk+1) / e3vw_n(ji,jj,jk+1)
               zws(ji,jj,jk) =  zzws * wvmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = 1, jpi   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------
      !
      !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1        
         DO jj = 2, jpjm1   
            DO ji = 1, jpi   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = 1, jpi   ! vector opt.
            ze3va =  ( 1._wp - r_vvl ) * e3v_n(ji,jj,1) + r_vvl   * e3v_a(ji,jj,1) 
            va(ji,jj,1) = va(ji,jj,1) + p2dt * 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( ze3va * rau0 ) * vmask(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 1, jpi   ! vector opt.
               zrhs = va(ji,jj,jk)   ! zrhs=right hand side
               va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==
         DO ji = 1, jpi   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = 1, jpi
               va(ji,jj,jk) = ( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( iom_use( 'dispkevfo' ) ) THEN   ! ocean kinetic energy dissipation per unit area
         !                                ! due to v friction (v=vertical) 
         !                                ! see NEMO_book appendix C, §C.8 (N.B. here averaged at t-points)
         !                                ! Note that formally, in a Leap-Frog environment, the shear**2 should be the product of 
         !                                ! now by before shears, i.e. the source term of TKE (local positivity is not ensured).
         CALL wrk_alloc(jpi,jpj,   z2d )
         z2d(:,:) = 0._wp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  z2d(ji,jj) = z2d(ji,jj)  +  (                                                                                  &
                     &   avmu(ji  ,jj,jk) * ( ua(ji  ,jj,jk-1) - ua(ji  ,jj,jk) )**2 / e3uw_n(ji  ,jj,jk) * wumask(ji  ,jj,jk)   &
                     & + avmu(ji-1,jj,jk) * ( ua(ji-1,jj,jk-1) - ua(ji-1,jj,jk) )**2 / e3uw_n(ji-1,jj,jk) * wumask(ji-1,jj,jk)   &
                     & + avmv(ji,jj  ,jk) * ( va(ji,jj  ,jk-1) - va(ji,jj  ,jk) )**2 / e3vw_n(ji,jj  ,jk) * wvmask(ji,jj  ,jk)   &
                     & + avmv(ji,jj-1,jk) * ( va(ji,jj-1,jk-1) - va(ji,jj-1,jk) )**2 / e3vw_n(ji,jj-1,jk) * wvmask(ji,jj-1,jk)   &
                     &                        )
               END DO
            END DO
         END DO
         zzz= - 0.5_wp* rau0           ! caution sign minus here
         z2d(:,:) = zzz * z2d(:,:) 
         CALL lbc_lnk( z2d,'T', 1. )
         CALL iom_put( 'dispkevfo', z2d )
         CALL wrk_dealloc(jpi,jpj,   z2d )
      ENDIF


      ! J. Chanut: Lines below are useless ?
      !! restore bottom layer avmu(v) 
      IF( ln_bfrimp ) THEN
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
              ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
              avmu(ji,jj,ikbu+1) = 0.e0
              avmv(ji,jj,ikbv+1) = 0.e0
           END DO
        END DO
        IF (ln_isfcav) THEN
           DO jj = 2, jpjm1
              DO ji = 2, jpim1
                 ikbu = miku(ji,jj)         ! ocean top level at u- and v-points 
                 ikbv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                 IF (ikbu > 1) avmu(ji,jj,ikbu) = 0.e0
                 IF (ikbv > 1) avmv(ji,jj,ikbv) = 0.e0
              END DO
           END DO
        END IF
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwi, zwd, zws) 
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_imp')
      !
   END SUBROUTINE dyn_zdf_imp

   !!==============================================================================
END MODULE dynzdf_imp
