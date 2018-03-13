MODULE bdyice_lim
   !!======================================================================
   !!                       ***  MODULE  bdyice_lim  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (LIM2 and LIM3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2011    (D. Storkey)  rewrite in preparation for OBC-BDY merge
   !!              -   !  2012-01 (C. Rousset)  add lim3 and remove useless jk loop 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_bdy'            and                 Unstructured Open Boundary Conditions
   !!   'key_lim2'                                                 LIM-2 sea ice model
   !!   'key_lim3'                                                 LIM-3 sea ice model
   !!----------------------------------------------------------------------
   !!   bdy_ice_lim        : Application of open boundaries to ice
   !!   bdy_ice_frs        : Application of Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE phycst          ! physical constant
   USE eosbn2          ! equation of state
   USE oce             ! ocean dynamics and tracers variables

   USE par_ice_2
   USE ice_2           ! LIM_2 ice variables
   USE dom_ice_2       ! sea-ice domain





   USE par_oce         ! ocean parameters
   USE dom_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE bdy_oce         ! ocean open boundary conditions
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! write to numout file
   USE lib_mpp         ! distributed memory computing
   USE lib_fortran     ! to use key_nosignedzero

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_ice_lim     ! routine called in sbcmod
   PUBLIC   bdy_ice_lim_dyn ! routine called in limrhg

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdyice_lim.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_ice_lim( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ice_lim  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for ice (LIM2 and LIM3)
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt     ! Main time step counter
      INTEGER               :: ib_bdy ! Loop index





      DO ib_bdy=1, nb_bdy

         SELECT CASE( cn_ice_lim(ib_bdy) )
         CASE('none')
            CYCLE
         CASE('frs')
            CALL bdy_ice_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE DEFAULT
            CALL ctl_stop( 'bdy_ice_lim : unrecognised option for open boundaries for ice fields' )
         END SELECT

      END DO






   END SUBROUTINE bdy_ice_lim

   SUBROUTINE bdy_ice_frs( idx, dta, kt, ib_bdy )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_frs  ***
      !!                    
      !! ** Purpose : Apply the Flow Relaxation Scheme for sea-ice fields in the case 
      !!              of unstructured open boundaries.
      !! 
      !! Reference : Engedahl H., 1995: Use of the flow relaxation scheme in a three-
      !!             dimensional baroclinic ocean model with realistic topography. Tellus, 365-382.
      !!------------------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   kt   ! main time-step counter
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index

      INTEGER  ::   jpbound            ! 0 = incoming ice
                                       ! 1 = outgoing ice
      INTEGER  ::   jb, jk, jgrd, jl   ! dummy loop indices
      INTEGER  ::   ji, jj, ii, ij     ! local scalar
      REAL(wp) ::   zwgt, zwgt1        ! local scalar
      REAL(wp) ::   ztmelts, zdh





      !!------------------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_ice_frs')
      !
      jgrd = 1      ! Everything is at T-points here
      !

      DO jb = 1, idx%nblenrim(jgrd)
         ji    = idx%nbi(jb,jgrd)
         jj    = idx%nbj(jb,jgrd)
         zwgt  = idx%nbw(jb,jgrd)
         zwgt1 = 1.e0 - idx%nbw(jb,jgrd)
         frld (ji,jj) = ( frld (ji,jj) * zwgt1 + dta%frld (jb) * zwgt ) * tmask(ji,jj,1)     ! Leads fraction 
         hicif(ji,jj) = ( hicif(ji,jj) * zwgt1 + dta%hicif(jb) * zwgt ) * tmask(ji,jj,1)     ! Ice depth 
         hsnif(ji,jj) = ( hsnif(ji,jj) * zwgt1 + dta%hsnif(jb) * zwgt ) * tmask(ji,jj,1)     ! Snow depth
      END DO 

      CALL lbc_bdy_lnk( frld,  'T', 1., ib_bdy )                                         ! lateral boundary conditions
      CALL lbc_bdy_lnk( hicif, 'T', 1., ib_bdy )
      CALL lbc_bdy_lnk( hsnif, 'T', 1., ib_bdy )

      vt_i(:,:) = hicif(:,:) * frld(:,:)
      vt_s(:,:) = hsnif(:,:) * frld(:,:)
      !
      !      
      IF( nn_timing == 1 ) CALL timing_stop('bdy_ice_frs')
      !
   END SUBROUTINE bdy_ice_frs


   SUBROUTINE bdy_ice_lim_dyn( cd_type )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_lim_dyn  ***
      !!                    
      !! ** Purpose : Apply dynamics boundary conditions for sea-ice in the cas of unstructured open boundaries.
      !!              u_ice and v_ice are equal to the value of the adjacent grid point if this latter is not ice free
      !!              if adjacent grid point is ice free, then u_ice and v_ice are equal to ocean velocities
      !!
      !! 2013-06 : C. Rousset
      !!------------------------------------------------------------------------------
      !!
      CHARACTER(len=1), INTENT(in)  ::   cd_type   ! nature of velocity grid-points
      INTEGER  ::   jb, jgrd           ! dummy loop indices
      INTEGER  ::   ji, jj             ! local scalar
      INTEGER  ::   ib_bdy             ! Loop index
      REAL(wp) ::   zmsk1, zmsk2, zflag
     !!------------------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_ice_lim_dyn')
      !
      DO ib_bdy=1, nb_bdy
         !
         SELECT CASE( cn_ice_lim(ib_bdy) )

         CASE('none')

            CYCLE
            
         CASE('frs')
            
            IF( nn_ice_lim_dta(ib_bdy) == 0 ) CYCLE            ! case ice boundaries = initial conditions 
                                                               !      do not change ice velocity (it is only computed by rheology)
 
            SELECT CASE ( cd_type )
               
            CASE ( 'U' )
               
               jgrd = 2      ! u velocity
               DO jb = 1, idx_bdy(ib_bdy)%nblenrim(jgrd)
                  ji    = idx_bdy(ib_bdy)%nbi(jb,jgrd)
                  jj    = idx_bdy(ib_bdy)%nbj(jb,jgrd)
                  zflag = idx_bdy(ib_bdy)%flagu(jb,jgrd)
                  
                  IF ( ABS( zflag ) == 1. ) THEN  ! eastern and western boundaries
                     ! one of the two zmsk is always 0 (because of zflag)
                     zmsk1 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji+1,jj) ) ) ! 0 if no ice
                     zmsk2 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji-1,jj) ) ) ! 0 if no ice
                     
                     ! u_ice = u_ice of the adjacent grid point except if this grid point is ice-free (then u_ice = u_oce)
                     u_ice (ji,jj) = u_ice(ji+1,jj) * 0.5_wp * ABS( zflag + 1._wp ) * zmsk1 + &
                        &            u_ice(ji-1,jj) * 0.5_wp * ABS( zflag - 1._wp ) * zmsk2 + &
                        &            u_oce(ji  ,jj) * ( 1._wp - MIN( 1._wp, zmsk1 + zmsk2 ) )
                  ELSE                             ! everywhere else
                     !u_ice(ji,jj) = u_oce(ji,jj)
                     u_ice(ji,jj) = 0._wp
                  ENDIF
                  ! mask ice velocities
                  rswitch = MAX( 0.0_wp , SIGN ( 1.0_wp , at_i(ji,jj) - 0.01_wp ) ) ! 0 if no ice
                  u_ice(ji,jj) = rswitch * u_ice(ji,jj)
                  
               ENDDO
               
               CALL lbc_bdy_lnk( u_ice(:,:), 'U', -1., ib_bdy )
               
            CASE ( 'V' )
               
               jgrd = 3      ! v velocity
               DO jb = 1, idx_bdy(ib_bdy)%nblenrim(jgrd)
                  ji    = idx_bdy(ib_bdy)%nbi(jb,jgrd)
                  jj    = idx_bdy(ib_bdy)%nbj(jb,jgrd)
                  zflag = idx_bdy(ib_bdy)%flagv(jb,jgrd)
                  
                  IF ( ABS( zflag ) == 1. ) THEN  ! northern and southern boundaries
                     ! one of the two zmsk is always 0 (because of zflag)
                     zmsk1 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji,jj+1) ) ) ! 0 if no ice
                     zmsk2 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji,jj-1) ) ) ! 0 if no ice
                     
                     ! u_ice = u_ice of the adjacent grid point except if this grid point is ice-free (then u_ice = u_oce)
                     v_ice (ji,jj) = v_ice(ji,jj+1) * 0.5_wp * ABS( zflag + 1._wp ) * zmsk1 + &
                        &            v_ice(ji,jj-1) * 0.5_wp * ABS( zflag - 1._wp ) * zmsk2 + &
                        &            v_oce(ji,jj  ) * ( 1._wp - MIN( 1._wp, zmsk1 + zmsk2 ) )
                  ELSE                             ! everywhere else
                     !v_ice(ji,jj) = v_oce(ji,jj)
                     v_ice(ji,jj) = 0._wp
                  ENDIF
                  ! mask ice velocities
                  rswitch = MAX( 0.0_wp , SIGN ( 1.0_wp , at_i(ji,jj) - 0.01 ) ) ! 0 if no ice
                  v_ice(ji,jj) = rswitch * v_ice(ji,jj)
                  
               ENDDO
               
               CALL lbc_bdy_lnk( v_ice(:,:), 'V', -1., ib_bdy )
                  
            END SELECT
            
         CASE DEFAULT
            CALL ctl_stop( 'bdy_ice_lim_dyn : unrecognised option for open boundaries for ice fields' )
         END SELECT
         
      ENDDO

      IF( nn_timing == 1 ) CALL timing_stop('bdy_ice_lim_dyn')
      
    END SUBROUTINE bdy_ice_lim_dyn


   !!=================================================================================
END MODULE bdyice_lim
