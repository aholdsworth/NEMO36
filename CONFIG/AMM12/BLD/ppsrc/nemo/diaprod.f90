MODULE diaprod
! Requires key_iom_put

   !!======================================================================
   !!                     ***  MODULE  diaprod  ***
   !! Ocean diagnostics :  write ocean product diagnostics
   !!=====================================================================
   !! History :  3.4  ! 2012  (D. Storkey)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_prod      : calculate and write out product diagnostics
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! for thickness weighted diagnostics if key_vvl
   USE eosbn2          ! equation of state  (eos call)
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE diadimg         ! dimg direct access file format output
   USE iom
   USE ioipsl
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE wrk_nemo        ! working array
   USE diaptr

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_prod                 ! routines called by step.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   'key_zdfddm' :                      avs: 3D array defined in zdfddm module
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 8026 2017-05-15 15:54:57Z lovato $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
   !! NEMO/OPA 3.4 , NEMO Consortium (2012)
   !! $Id $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_prod( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_prod  ***
      !!                   
      !! ** Purpose :   Write out product diagnostics (uT, vS etc.)
      !!
      !! ** Method  :  use iom_put
      !!               Product diagnostics are not thickness-weighted in this routine.
      !!               They should be thickness-weighted using XIOS if key_vvl is set. 
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zrhop    ! potential density
      !!----------------------------------------------------------------------
      ! 
      IF( nn_timing == 1 )   CALL timing_start('dia_prod')
      ! 
      CALL wrk_alloc( jpi , jpj      , z2d )
      CALL wrk_alloc( jpi , jpj, jpk , z3d )
      CALL wrk_alloc( jpi , jpj, jpk , zrhop )
      !

      IF( iom_use("urhop") .OR. iom_use("vrhop") .OR. iom_use("wrhop") &
     &  .OR. iom_use("rhop") &
     & ) THEN
         CALL eos( tsn, z3d, zrhop )                 ! now in situ and potential density
         zrhop(:,:,:) = zrhop(:,:,:)-1000.e0         ! reference potential density to 1000 to avoid precision issues in rhop2 calculation
         zrhop(:,:,jpk) = 0._wp
         CALL iom_put( 'rhop', zrhop )
      ENDIF

      IF( iom_use("ut") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "ut", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vt") .OR. iom_use("sopht_vt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "vt", z3d )                  ! product of temperature and meridional velocity at V points
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = z3d(ji,jj,jk) * e3v_0(ji,jj,jk) * e1v(ji,jj)
               END DO
            END DO
         END DO
         IF(ln_diaptr) CALL dia_ptr_ohst_components( jp_tem, 'vts', z3d)
      ENDIF

      IF( iom_use("wt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * tsn(ji,jj,1,jp_tem)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk-1,jp_tem) + tsn(ji,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "wt", z3d )                  ! product of temperature and vertical velocity at W points
      ENDIF

      IF( iom_use("us") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "us", z3d )                  ! product of salinity and zonal velocity at U points
      ENDIF

      IF( iom_use("vs") .OR. iom_use("sopst_vs") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "vs", z3d )                  ! product of salinity and meridional velocity at V points
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = z3d(ji,jj,jk) * e3v_0(ji,jj,jk) * e1v(ji,jj)
               END DO
            END DO
         END DO
         IF(ln_diaptr) CALL dia_ptr_ohst_components( jp_sal, 'vts', z3d)
      ENDIF

      IF( iom_use("ws") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * tsn(ji,jj,1,jp_sal)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk-1,jp_sal) + tsn(ji,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "ws", z3d )                  ! product of salinity and vertical velocity at W points
      ENDIF

      IF( iom_use("urhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji+1,jj,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "urhop", z3d )                  ! product of density and zonal velocity at U points
      ENDIF

      IF( iom_use("vrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji,jj+1,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "vrhop", z3d )                  ! product of density and meridional velocity at V points
      ENDIF

      IF( iom_use("wrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * zrhop(ji,jj,1)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk-1) + zrhop(ji,jj,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "wrhop", z3d )                  ! product of density and vertical velocity at W points
      ENDIF

      !
      CALL wrk_dealloc( jpi , jpj      , z2d )
      CALL wrk_dealloc( jpi , jpj, jpk , z3d )
      CALL wrk_dealloc( jpi , jpj, jpk , zrhop )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_prod')
      !
   END SUBROUTINE dia_prod
   !!======================================================================
END MODULE diaprod
