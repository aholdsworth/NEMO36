MODULE bdydyn3d
   !!======================================================================
   !!                       ***  MODULE  bdydyn3d  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on baroclinic velocities
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite 
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   !!   bdy_dyn3d        : apply open boundary conditions to baroclinic velocities
   !!   bdy_dyn3d_frs    : apply Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdylib          ! for orlanski library routines
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  !
   Use phycst
   USE ldfdyn_oce

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn3d     ! routine called by bdy_dyn
   PUBLIC   bdy_dyn3d_dmp ! routine called by step

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdydyn3d.F90 6808 2016-07-19 08:38:35Z jamesharle $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn3d( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for baroclinic velocities
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      !
      INTEGER ::   ib_bdy   ! loop index
      !!----------------------------------------------------------------------
      !
      DO ib_bdy=1, nb_bdy
         !
         SELECT CASE( cn_dyn3d(ib_bdy) )
         CASE('none')        ;   CYCLE
         CASE('frs' )        ;   CALL bdy_dyn3d_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('specified')   ;   CALL bdy_dyn3d_spe( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('zero')        ;   CALL bdy_dyn3d_zro( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('orlanski' )   ;   CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.false. )
         CASE('orlanski_npo');   CALL bdy_dyn3d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, ll_npo=.true. )
         CASE('zerograd')    ;   CALL bdy_dyn3d_zgrad( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE('neumann')     ;   CALL bdy_dyn3d_nmn( idx_bdy(ib_bdy), ib_bdy )
         CASE DEFAULT        ;   CALL ctl_stop( 'bdy_dyn3d : unrecognised option for open boundaries for baroclinic velocities' )
         END SELECT
! ML
! apply sponge layer in the boundary for viscosity
         IF(ln_sponge_dyn(ib_bdy)) CALL bdy_dyn_sponge ( idx_bdy(ib_bdy), kt, ib_bdy  )

      END DO

      !
   END SUBROUTINE bdy_dyn3d


   SUBROUTINE bdy_dyn3d_spe( idx, dta, kt , ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_spe  ***
      !!
      !! ** Purpose : - Apply a specified value for baroclinic velocities
      !!                at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER        , INTENT(in) ::   kt      ! time step index
      TYPE(OBC_INDEX), INTENT(in) ::   idx     ! OBC indices
      TYPE(OBC_DATA) , INTENT(in) ::   dta     ! OBC external data
      INTEGER        , INTENT(in) ::   ib_bdy  ! BDY set index
      !
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_spe')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            ua(ii,ij,jk) = dta%u3d(jb,jk) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            va(ii,ij,jk) = dta%v3d(jb,jk) * vmask(ii,ij,jk)
         END DO
      END DO
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ! Boundary points should be updated  
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( kt == nit000 )   CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_spe')
      !
   END SUBROUTINE bdy_dyn3d_spe

   SUBROUTINE bdy_dyn3d_zgrad( idx, dta, kt , ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_zgrad  ***
      !!
      !! ** Purpose : - Enforce a zero gradient of normal velocity
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   fu, fv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_zgrad')
      !
      igrd = 2                      ! Copying tangential velocity into bdy points
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            fu   = ABS( ABS (NINT( idx%flagu(jb,igrd) ) ) - 1 )
            ua(ii,ij,jk) = ua(ii,ij,jk) * REAL( 1 - fu) + ( ua(ii,ij+fu,jk) * umask(ii,ij+fu,jk) &
                        &+ ua(ii,ij-fu,jk) * umask(ii,ij-fu,jk) ) * umask(ii,ij,jk) * REAL( fu )
         END DO
      END DO
      !
      igrd = 3                      ! Copying tangential velocity into bdy points
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            fv   = ABS( ABS (NINT( idx%flagv(jb,igrd) ) ) - 1 )
            va(ii,ij,jk) = va(ii,ij,jk) * REAL( 1 - fv ) + ( va(ii+fv,ij,jk) * vmask(ii+fv,ij,jk) &
                        &+ va(ii-fv,ij,jk) * vmask(ii-fv,ij,jk) ) * vmask(ii,ij,jk) * REAL( fv )
         END DO
      END DO
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ! Boundary points should be updated  
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_zgrad')

   END SUBROUTINE bdy_dyn3d_zgrad

   SUBROUTINE bdy_dyn3d_zro( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_zro  ***
      !!
      !! ** Purpose : - baroclinic velocities = 0. at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER        , INTENT(in) ::   kt      ! time step index
      TYPE(OBC_INDEX), INTENT(in) ::   idx     ! OBC indices
      TYPE(OBC_DATA) , INTENT(in) ::   dta     ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !
      INTEGER  ::   ib, ik         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_zro')
      !
      igrd = 2                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ua(ii,ij,ik) = 0._wp
         END DO
      END DO

      igrd = 3                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            va(ii,ij,ik) = 0._wp
         END DO
      END DO
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ;   CALL lbc_bdy_lnk( va, 'V', -1.,ib_bdy )   ! Boundary points should be updated
      !
      IF( kt == nit000 )   CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 )   CALL timing_stop('bdy_dyn3d_zro')
      !
   END SUBROUTINE bdy_dyn3d_zro


   SUBROUTINE bdy_dyn3d_frs( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for baroclinic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER        , INTENT(in) ::   kt      ! time step index
      TYPE(OBC_INDEX), INTENT(in) ::   idx     ! OBC indices
      TYPE(OBC_DATA) , INTENT(in) ::   dta     ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_frs')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta%u3d(jb,jk) - ua(ii,ij,jk) ) ) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta%v3d(jb,jk) - va(ii,ij,jk) ) ) * vmask(ii,ij,jk)
         END DO
      END DO 
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( kt == nit000 )   CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 )   CALL timing_stop('bdy_dyn3d_frs')
      !
   END SUBROUTINE bdy_dyn3d_frs


   SUBROUTINE bdy_dyn3d_orlanski( idx, dta, ib_bdy, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn3d_orlanski  ***
      !!             
      !!              - Apply Orlanski radiation to baroclinic velocities. 
      !!              - Wrapper routine for bdy_orlanski_3d
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! BDY set index
      LOGICAL,                      INTENT(in) ::   ll_npo  ! switch for NPO version

      INTEGER  ::   jb, igrd                               ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_orlanski')
      !
      !! Note that at this stage the ub and ua arrays contain the baroclinic velocities. 
      !
      igrd = 2      ! Orlanski bc on u-velocity; 
      !            
      CALL bdy_orlanski_3d( idx, igrd, ub, ua, dta%u3d, ll_npo )

      igrd = 3      ! Orlanski bc on v-velocity
      !  
      CALL bdy_orlanski_3d( idx, igrd, vb, va, dta%v3d, ll_npo )
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )   
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_orlanski')
      !
   END SUBROUTINE bdy_dyn3d_orlanski


   SUBROUTINE bdy_dyn3d_dmp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_dmp  ***
      !!
      !! ** Purpose : Apply damping for baroclinic velocities at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ib_bdy         ! loop index
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('bdy_dyn3d_dmp')
      !
      DO ib_bdy=1, nb_bdy
         IF ( ln_dyn3d_dmp(ib_bdy) .and. cn_dyn3d(ib_bdy) /= 'none' ) THEN
            igrd = 2                      ! Relaxation of zonal velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%u3d(jb,jk) - &
                                   ub(ii,ij,jk) + ub_b(ii,ij)) ) * umask(ii,ij,jk)
               END DO
            END DO
            !
            igrd = 3                      ! Relaxation of meridional velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%v3d(jb,jk) -  &
                                   vb(ii,ij,jk) + vb_b(ii,ij)) ) * vmask(ii,ij,jk)
               END DO
            END DO
         ENDIF
      END DO
      !
      CALL lbc_lnk( ua, 'U', -1. )   ;   CALL lbc_lnk( va, 'V', -1. )   ! Boundary points should be updated
      !
      IF( nn_timing == 1 )   CALL timing_stop('bdy_dyn3d_dmp')
      !
   END SUBROUTINE bdy_dyn3d_dmp

   SUBROUTINE bdy_dyn3d_nmn( idx, ib_bdy )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn3d_nmn  ***
      !!             
      !!              - Apply Neumann condition to baroclinic velocities. 
      !!              - Wrapper routine for bdy_nmn
      !! 
      !!
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      INTEGER,                      INTENT(in) ::   ib_bdy  ! BDY set index

      INTEGER  ::   jb, igrd                               ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_nmn')
      !
      !! Note that at this stage the ub and ua arrays contain the baroclinic velocities. 
      !
      igrd = 2      ! Neumann bc on u-velocity; 
      !            
      CALL bdy_nmn( idx, igrd, ua )

      igrd = 3      ! Neumann bc on v-velocity
      !  
      CALL bdy_nmn( idx, igrd, va )
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( va, 'V', -1., ib_bdy )
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_nmn')
      !
   END SUBROUTINE bdy_dyn3d_nmn
   SUBROUTINE bdy_dyn_sponge( idx, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for baroclinic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER        , INTENT(in) ::   kt      ! time step index
      TYPE(OBC_INDEX), INTENT(in) ::   idx     ! OBC indices
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !
      INTEGER  ::   jb, jk ,ji,jj  ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn_sponge')
      !
      IF(ln_dynldf_lap)THEN
      igrd = 1                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = 1. + (rn_sponge_dyn(igrd)-1.) * idx%nbw(jb,igrd)

            ahm1(ii,ij,:) = ahm1(ii,ij,:)  * zwgt 



      END DO
! in point F no boundary points defined
      DO jj = 1,jpjm1
       DO ji = 1, jpim1

            ahm2(ii,ij,:) = 0.25 * (ahm1(ji+1,jj,:) + ahm1(ji,jj+1,:)  &
                                 +  ahm1(ji,jj,:) + ahm1(ji+1,jj+1,:) ) 





         
       ENDDO
      ENDDO  
     ENDIF
         
      IF(ln_dynldf_bilap)THEN
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = 1. + (rn_sponge_dynb(igrd)-1.) * idx%nbw(jb,igrd)

            ahm3(ii,ij,:) = ahm3(ii,ij,:)  * zwgt



      END DO
      igrd = 3

      DO jb = 1, idx%nblen(igrd)
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = 1. + (rn_sponge_dynb(igrd)-1.) * idx%nbw(jb,igrd)

            ahm4(ii,ij,:) = ahm4(ii,ij,:)  * zwgt



      END DO


     ENDIF

      CALL lbc_lnk( ahm2, 'F', 1. )    ! Boundary points should be updated
      !

      CALL lbc_bdy_lnk( ahm3, 'U', -1., ib_bdy )    ! Boundary points should be updated
      CALL lbc_bdy_lnk( ahm4, 'V', -1., ib_bdy )
      !
      IF( kt == nit000 )   CLOSE( unit = 102 )
      !
      IF( nn_timing == 1 )   CALL timing_stop('bdy_dyn_sponge')
      !
   END SUBROUTINE bdy_dyn_sponge


   !!======================================================================
END MODULE bdydyn3d