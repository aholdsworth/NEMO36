MODULE limadv 
   !!======================================================================
   !!                       ***  MODULE limadv   ***
   !! LIM sea-ice model : sea-ice advection
   !!======================================================================
   !! History :  LIM  ! 2008-03 (M. Vancoppenolle)  LIM-3 from LIM-2 code
   !!            3.2  ! 2009-06 (F. Dupont)  correct a error in the North fold b. c.
   !!            4.0  ! 2011-02 (G. Madec) dynamical allocation
   !!--------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_lim3'                                     LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_adv_x  : advection of sea ice on x axis
   !!   lim_adv_y  : advection of sea ice on y axis
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean domain
   USE dom_ice          ! LIM-3 domain
   USE ice              ! LIM-3 variables
   USE lbclnk           ! lateral boundary condition - MPP exchanges
   USE in_out_manager   ! I/O manager
   USE prtctl           ! Print control
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE lib_fortran      ! to use key_nosignedzero

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_adv_x   ! called by lim_trp
   PUBLIC   lim_adv_y   ! called by lim_trp

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
   !! $Id: limadv.F90 5429 2015-06-16 09:57:07Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_adv_x( pdf, put , pcrh, psm , ps0 ,   &
      &                  psx, psxx, psy , psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine lim_adv_x  **
      !!  
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!              variable on i-axis
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!              but also theirquadratic forms. The method preserves
      !!              tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!--------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdf                ! reduction factor for the time step
      REAL(wp)                    , INTENT(in   ) ::   pcrh               ! call lim_adv_x then lim_adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   put                ! i-direction ice velocity at U-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   ps0                ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psx , psy          ! 1st moments 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !! 
      INTEGER  ::   ji, jj                               ! dummy loop indices
      REAL(wp) ::   zs1max, zrdt, zslpmax, ztemp         ! local scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !   -      -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zf0 , zfx  , zfy   , zbet   ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zfm , zfxx , zfyy  , zfxy   !  -      -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zalg, zalg1, zalg1q         !  -      -
      !---------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zf0 , zfx , zfy , zbet, zfm )
      CALL wrk_alloc( jpi, jpj, zfxx, zfyy, zfxy, zalg, zalg1, zalg1q )

      ! Limitation of moments.                                           

      zrdt = rdt_ice * pdf      ! If ice drift field is too fast, use an appropriate time step for advection.

      DO jj = 1, jpj
         DO ji = 1, jpi
            zslpmax = MAX( 0._wp, ps0(ji,jj) )
            zs1max  = 1.5 * zslpmax
            zs1new  = MIN( zs1max, MAX( -zs1max, psx(ji,jj) ) )
            zs2new  = MIN(  2.0 * zslpmax - 0.3334 * ABS( zs1new ),      &
               &            MAX( ABS( zs1new ) - zslpmax, psxx(ji,jj) )  )
            rswitch = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zslpmax) ) ) * tmask(ji,jj,1)   ! Case of empty boxes & Apply mask

            ps0 (ji,jj) = zslpmax  
            psx (ji,jj) = zs1new      * rswitch
            psxx(ji,jj) = zs2new      * rswitch
            psy (ji,jj) = psy (ji,jj) * rswitch
            psyy(ji,jj) = psyy(ji,jj) * rswitch
            psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * rswitch
         END DO
      END DO

      !  Initialize volumes of boxes  (=area if adv_x first called, =psm otherwise)                                     
      psm (:,:)  = MAX( pcrh * e12t(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20 )

      !  Calculate fluxes and moments between boxes i<-->i+1              
      DO jj = 1, jpj                      !  Flux from i to i+1 WHEN u GT 0 
         DO ji = 1, jpi
            zbet(ji,jj)  =  MAX( 0._wp, SIGN( 1._wp, put(ji,jj) ) )
            zalf         =  MAX( 0._wp, put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji,jj)
            zalfq        =  zalf * zalf
            zalf1        =  1.0 - zalf
            zalf1q       =  zalf1 * zalf1
            !
            zfm (ji,jj)  =  zalf  *   psm (ji,jj)
            zf0 (ji,jj)  =  zalf  * ( ps0 (ji,jj) + zalf1 * ( psx(ji,jj) + (zalf1 - zalf) * psxx(ji,jj) )  )
            zfx (ji,jj)  =  zalfq * ( psx (ji,jj) + 3.0 * zalf1 * psxx(ji,jj) )
            zfxx(ji,jj)  =  zalf  *   psxx(ji,jj) * zalfq
            zfy (ji,jj)  =  zalf  * ( psy (ji,jj) + zalf1 * psxy(ji,jj) )
            zfxy(ji,jj)  =  zalfq *   psxy(ji,jj)
            zfyy(ji,jj)  =  zalf  *   psyy(ji,jj)

            !  Readjust moments remaining in the box.
            psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
            ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
            psx (ji,jj)  =  zalf1q * ( psx(ji,jj) - 3.0 * zalf * psxx(ji,jj) )
            psxx(ji,jj)  =  zalf1  * zalf1q * psxx(ji,jj)
            psy (ji,jj)  =  psy (ji,jj) - zfy(ji,jj)
            psyy(ji,jj)  =  psyy(ji,jj) - zfyy(ji,jj)
            psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
         END DO
      END DO

      DO jj = 1, jpjm1                      !  Flux from i+1 to i when u LT 0.
         DO ji = 1, jpi
            zalf          = MAX( 0._wp, -put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji+1,jj) 
            zalg  (ji,jj) = zalf
            zalfq         = zalf * zalf
            zalf1         = 1.0 - zalf
            zalg1 (ji,jj) = zalf1
            zalf1q        = zalf1 * zalf1
            zalg1q(ji,jj) = zalf1q
            !
            zfm   (ji,jj) = zfm (ji,jj) + zalf  *   psm (ji+1,jj)
            zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0 (ji+1,jj) - zalf1 * ( psx(ji+1,jj) - (zalf1 - zalf ) * psxx(ji+1,jj) ) )
            zfx   (ji,jj) = zfx (ji,jj) + zalfq * ( psx (ji+1,jj) - 3.0 * zalf1 * psxx(ji+1,jj) )
            zfxx  (ji,jj) = zfxx(ji,jj) + zalf  *   psxx(ji+1,jj) * zalfq
            zfy   (ji,jj) = zfy (ji,jj) + zalf  * ( psy (ji+1,jj) - zalf1 * psxy(ji+1,jj) )
            zfxy  (ji,jj) = zfxy(ji,jj) + zalfq *   psxy(ji+1,jj)
            zfyy  (ji,jj) = zfyy(ji,jj) + zalf  *   psyy(ji+1,jj)
         END DO
      END DO

      DO jj = 2, jpjm1                     !  Readjust moments remaining in the box. 
         DO ji = 1, jpi
            zbt  =       zbet(ji-1,jj)
            zbt1 = 1.0 - zbet(ji-1,jj)
            !
            psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji-1,jj) )
            ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji-1,jj) )
            psx (ji,jj) = zalg1q(ji-1,jj) * ( psx(ji,jj) + 3.0 * zalg(ji-1,jj) * psxx(ji,jj) )
            psxx(ji,jj) = zalg1 (ji-1,jj) * zalg1q(ji-1,jj) * psxx(ji,jj)
            psy (ji,jj) = zbt * psy (ji,jj) + zbt1 * ( psy (ji,jj) - zfy (ji-1,jj) )
            psyy(ji,jj) = zbt * psyy(ji,jj) + zbt1 * ( psyy(ji,jj) - zfyy(ji-1,jj) )
            psxy(ji,jj) = zalg1q(ji-1,jj) * psxy(ji,jj)
         END DO
      END DO

      !   Put the temporary moments into appropriate neighboring boxes.    
      DO jj = 2, jpjm1                     !   Flux from i to i+1 IF u GT 0.
         DO ji = 1, jpi
            zbt  =       zbet(ji-1,jj)
            zbt1 = 1.0 - zbet(ji-1,jj)
            psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji-1,jj) ) + zbt1 * psm(ji,jj)
            zalf        = zbt * zfm(ji-1,jj) / psm(ji,jj)
            zalf1       = 1.0 - zalf
            ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji-1,jj)
            !
            ps0 (ji,jj) = zbt * ( ps0(ji,jj) + zf0(ji-1,jj) ) + zbt1 * ps0(ji,jj)
            psx (ji,jj) = zbt * ( zalf * zfx(ji-1,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp ) + zbt1 * psx(ji,jj)
            psxx(ji,jj) = zbt * ( zalf * zalf * zfxx(ji-1,jj) + zalf1 * zalf1 * psxx(ji,jj)                               &
               &                + 5.0 * ( zalf * zalf1 * ( psx (ji,jj) - zfx(ji-1,jj) ) - ( zalf1 - zalf ) * ztemp )  )   &
               &                                                + zbt1 * psxx(ji,jj)
            psxy(ji,jj) = zbt * ( zalf * zfxy(ji-1,jj) + zalf1 * psxy(ji,jj)             &
               &                + 3.0 * (- zalf1*zfy(ji-1,jj)  + zalf * psy(ji,jj) ) )   &
               &                                                + zbt1 * psxy(ji,jj)
            psy (ji,jj) = zbt * ( psy (ji,jj) + zfy (ji-1,jj) ) + zbt1 * psy (ji,jj)
            psyy(ji,jj) = zbt * ( psyy(ji,jj) + zfyy(ji-1,jj) ) + zbt1 * psyy(ji,jj)
         END DO
      END DO

      DO jj = 2, jpjm1                     !  Flux from i+1 to i IF u LT 0.
         DO ji = 1, jpi
            zbt  =       zbet(ji,jj)
            zbt1 = 1.0 - zbet(ji,jj)
            psm(ji,jj)  = zbt * psm(ji,jj)  + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
            zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
            zalf1       = 1.0 - zalf
            ztemp       = - zalf * ps0(ji,jj) + zalf1 * zf0(ji,jj)
            !
            ps0(ji,jj)  = zbt * ps0 (ji,jj) + zbt1 * ( ps0(ji,jj) + zf0(ji,jj) )
            psx(ji,jj)  = zbt * psx (ji,jj) + zbt1 * ( zalf * zfx(ji,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp )
            psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( zalf * zalf * zfxx(ji,jj)  + zalf1 * zalf1 * psxx(ji,jj)  &
               &                                      + 5.0 *( zalf * zalf1 * ( - psx(ji,jj) + zfx(ji,jj) )      &
               &                                      + ( zalf1 - zalf ) * ztemp ) )
            psxy(ji,jj) = zbt * psxy(ji,jj) + zbt1 * (  zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)  &
               &                                      + 3.0 * ( zalf1 * zfy(ji,jj) - zalf * psy(ji,jj) )  )
            psy(ji,jj)  = zbt * psy (ji,jj)  + zbt1 * ( psy (ji,jj) + zfy (ji,jj) )
            psyy(ji,jj) = zbt * psyy(ji,jj)  + zbt1 * ( psyy(ji,jj) + zfyy(ji,jj) )
         END DO
      END DO

      !-- Lateral boundary conditions
      CALL lbc_lnk_multi( psm , 'T',  1., ps0 , 'T',  1.   &
         &              , psx , 'T', -1., psy , 'T', -1.   &   ! caution gradient ==> the sign changes
         &              , psxx, 'T',  1., psyy, 'T',  1.   &
         &              , psxy, 'T',  1. )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_x: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
         CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_x: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
         CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_x: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
         CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_x: psxy :')
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zf0 , zfx , zfy , zbet, zfm )
      CALL wrk_dealloc( jpi, jpj, zfxx, zfyy, zfxy, zalg, zalg1, zalg1q )
      !
   END SUBROUTINE lim_adv_x


   SUBROUTINE lim_adv_y( pdf, pvt , pcrh, psm , ps0 ,   &
      &                  psx, psxx, psy , psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine lim_adv_y  **
      !!            
      !! ** purpose :   Computes and adds the advection trend to sea-ice 
      !!              variable on y axis
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!              but also their quadratic forms. The method preserves 
      !!              tracer structures by conserving second order moments.
      !! 
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!---------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdf                ! reduction factor for the time step
      REAL(wp)                    , INTENT(in   ) ::   pcrh               ! call lim_adv_x then lim_adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pvt                ! j-direction ice velocity at V-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   ps0                ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psx , psy          ! 1st moments 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!
      INTEGER  ::   ji, jj                               ! dummy loop indices
      REAL(wp) ::   zs1max, zrdt, zslpmax, ztemp         ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zf0, zfx , zfy , zbet   ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zfm, zfxx, zfyy, zfxy   !  -      -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zalg, zalg1, zalg1q     !  -      -
      !---------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zf0 , zfx , zfy , zbet, zfm )
      CALL wrk_alloc( jpi, jpj, zfxx, zfyy, zfxy, zalg, zalg1, zalg1q )

      ! Limitation of moments.

      zrdt = rdt_ice * pdf ! If ice drift field is too fast, use an appropriate time step for advection.

      DO jj = 1, jpj
         DO ji = 1, jpi
            zslpmax = MAX( 0._wp, ps0(ji,jj) )
            zs1max  = 1.5 * zslpmax
            zs1new  = MIN( zs1max, MAX( -zs1max, psy(ji,jj) ) )
            zs2new  = MIN(  ( 2.0 * zslpmax - 0.3334 * ABS( zs1new ) ),   &
               &             MAX( ABS( zs1new )-zslpmax, psyy(ji,jj) )  )
            rswitch = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zslpmax) ) ) * tmask(ji,jj,1)   ! Case of empty boxes & Apply mask
            !
            ps0 (ji,jj) = zslpmax  
            psx (ji,jj) = psx (ji,jj) * rswitch
            psxx(ji,jj) = psxx(ji,jj) * rswitch
            psy (ji,jj) = zs1new * rswitch
            psyy(ji,jj) = zs2new * rswitch
            psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * rswitch
         END DO
      END DO

      !  Initialize volumes of boxes (=area if adv_x first called, =psm otherwise)
      psm(:,:)  = MAX(  pcrh * e12t(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20  )

      !  Calculate fluxes and moments between boxes j<-->j+1              
      DO jj = 1, jpj                     !  Flux from j to j+1 WHEN v GT 0   
         DO ji = 1, jpi
            zbet(ji,jj)  =  MAX( 0._wp, SIGN( 1._wp, pvt(ji,jj) ) )
            zalf         =  MAX( 0._wp, pvt(ji,jj) ) * zrdt * e1v(ji,jj) / psm(ji,jj)
            zalfq        =  zalf * zalf
            zalf1        =  1.0 - zalf
            zalf1q       =  zalf1 * zalf1
            !
            zfm (ji,jj)  =  zalf  * psm(ji,jj)
            zf0 (ji,jj)  =  zalf  * ( ps0(ji,jj) + zalf1 * ( psy(ji,jj)  + (zalf1-zalf) * psyy(ji,jj)  ) ) 
            zfy (ji,jj)  =  zalfq *( psy(ji,jj) + 3.0*zalf1*psyy(ji,jj) )
            zfyy(ji,jj)  =  zalf  * zalfq * psyy(ji,jj)
            zfx (ji,jj)  =  zalf  * ( psx(ji,jj) + zalf1 * psxy(ji,jj) )
            zfxy(ji,jj)  =  zalfq * psxy(ji,jj)
            zfxx(ji,jj)  =  zalf  * psxx(ji,jj)
            !
            !  Readjust moments remaining in the box.
            psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
            ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
            psy (ji,jj)  =  zalf1q * ( psy(ji,jj) -3.0 * zalf * psyy(ji,jj) )
            psyy(ji,jj)  =  zalf1 * zalf1q * psyy(ji,jj)
            psx (ji,jj)  =  psx (ji,jj) - zfx(ji,jj)
            psxx(ji,jj)  =  psxx(ji,jj) - zfxx(ji,jj)
            psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
         END DO
      END DO
      !
      DO jj = 1, jpjm1                   !  Flux from j+1 to j when v LT 0.
         DO ji = 1, jpi
            zalf          = ( MAX(0._wp, -pvt(ji,jj) ) * zrdt * e1v(ji,jj) ) / psm(ji,jj+1) 
            zalg  (ji,jj) = zalf
            zalfq         = zalf * zalf
            zalf1         = 1.0 - zalf
            zalg1 (ji,jj) = zalf1
            zalf1q        = zalf1 * zalf1
            zalg1q(ji,jj) = zalf1q
            !
            zfm   (ji,jj) = zfm (ji,jj) + zalf  *   psm (ji,jj+1)
            zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0 (ji,jj+1) - zalf1 * (psy(ji,jj+1) - (zalf1 - zalf ) * psyy(ji,jj+1) ) )
            zfy   (ji,jj) = zfy (ji,jj) + zalfq * ( psy (ji,jj+1) - 3.0 * zalf1 * psyy(ji,jj+1) )
            zfyy  (ji,jj) = zfyy(ji,jj) + zalf  *   psyy(ji,jj+1) * zalfq
            zfx   (ji,jj) = zfx (ji,jj) + zalf  * ( psx (ji,jj+1) - zalf1 * psxy(ji,jj+1) )
            zfxy  (ji,jj) = zfxy(ji,jj) + zalfq *   psxy(ji,jj+1)
            zfxx  (ji,jj) = zfxx(ji,jj) + zalf  *   psxx(ji,jj+1)
         END DO
      END DO

      !  Readjust moments remaining in the box. 
      DO jj = 2, jpj
         DO ji = 1, jpi
            zbt  =         zbet(ji,jj-1)
            zbt1 = ( 1.0 - zbet(ji,jj-1) )
            !
            psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji,jj-1) )
            ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji,jj-1) )
            psy (ji,jj) = zalg1q(ji,jj-1) * ( psy(ji,jj) + 3.0 * zalg(ji,jj-1) * psyy(ji,jj) )
            psyy(ji,jj) = zalg1 (ji,jj-1) * zalg1q(ji,jj-1) * psyy(ji,jj)
            psx (ji,jj) = zbt * psx (ji,jj) + zbt1 * ( psx (ji,jj) - zfx (ji,jj-1) )
            psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) - zfxx(ji,jj-1) )
            psxy(ji,jj) = zalg1q(ji,jj-1) * psxy(ji,jj)
         END DO
      END DO

      !   Put the temporary moments into appropriate neighboring boxes.    
      DO jj = 2, jpjm1                    !   Flux from j to j+1 IF v GT 0.
         DO ji = 1, jpi
            zbt  =         zbet(ji,jj-1)
            zbt1 = ( 1.0 - zbet(ji,jj-1) )
            psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji,jj-1) ) + zbt1 * psm(ji,jj) 
            zalf        = zbt * zfm(ji,jj-1) / psm(ji,jj) 
            zalf1       = 1.0 - zalf
            ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji,jj-1)
            !
            ps0(ji,jj)  = zbt  * ( ps0(ji,jj) + zf0(ji,jj-1) ) + zbt1 * ps0(ji,jj)
            psy(ji,jj)  = zbt  * ( zalf * zfy(ji,jj-1) + zalf1 * psy(ji,jj) + 3.0 * ztemp )   &
               &                                               + zbt1 * psy(ji,jj)  
            psyy(ji,jj) = zbt  * ( zalf * zalf * zfyy(ji,jj-1) + zalf1 * zalf1 * psyy(ji,jj)                             &
               &                 + 5.0 * ( zalf * zalf1 * ( psy(ji,jj) - zfy(ji,jj-1) ) - ( zalf1 - zalf ) * ztemp ) )   & 
               &                                               + zbt1 * psyy(ji,jj)
            psxy(ji,jj) = zbt  * (  zalf * zfxy(ji,jj-1) + zalf1 * psxy(ji,jj)               &
               &                  + 3.0 * (- zalf1 * zfx(ji,jj-1) + zalf * psx(ji,jj) )  )   &
               &                                                + zbt1 * psxy(ji,jj)
            psx (ji,jj) = zbt * ( psx (ji,jj) + zfx (ji,jj-1) ) + zbt1 * psx (ji,jj)
            psxx(ji,jj) = zbt * ( psxx(ji,jj) + zfxx(ji,jj-1) ) + zbt1 * psxx(ji,jj)
         END DO
      END DO

      DO jj = 2, jpjm1                   !  Flux from j+1 to j IF v LT 0.
         DO ji = 1, jpi
            zbt  =         zbet(ji,jj)
            zbt1 = ( 1.0 - zbet(ji,jj) )
            psm(ji,jj)  = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
            zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
            zalf1       = 1.0 - zalf
            ztemp       = - zalf * ps0 (ji,jj) + zalf1 * zf0(ji,jj)
            ps0 (ji,jj) =   zbt  * ps0 (ji,jj) + zbt1  * ( ps0(ji,jj) + zf0(ji,jj) )
            psy (ji,jj) =   zbt  * psy (ji,jj) + zbt1  * ( zalf * zfy(ji,jj) + zalf1 * psy(ji,jj) + 3.0 * ztemp )
            psyy(ji,jj) =   zbt  * psyy(ji,jj) + zbt1 * (  zalf * zalf * zfyy(ji,jj) + zalf1 * zalf1 * psyy(ji,jj)   &
               &                                         + 5.0 *( zalf *zalf1 *( -psy(ji,jj) + zfy(ji,jj) )          &
               &                                         + ( zalf1 - zalf ) * ztemp )                                )
            psxy(ji,jj) =   zbt  * psxy(ji,jj) + zbt1 * (  zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)       &
               &                                         + 3.0 * ( zalf1 * zfx(ji,jj) - zalf * psx(ji,jj) )  )
            psx (ji,jj) =   zbt  * psx (ji,jj) + zbt1 * ( psx (ji,jj) + zfx (ji,jj) )
            psxx(ji,jj) =   zbt  * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) + zfxx(ji,jj) )
         END DO
      END DO

      !-- Lateral boundary conditions
      CALL lbc_lnk_multi( psm , 'T',  1.,  ps0 , 'T',  1.   &
         &              , psx , 'T', -1.,  psy , 'T', -1.   &   ! caution gradient ==> the sign changes
         &              , psxx, 'T',  1.,  psyy, 'T',  1.   &
         &              , psxy, 'T',  1. )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_y: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
         CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_y: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
         CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_y: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
         CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_y: psxy :')
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zf0 , zfx , zfy , zbet, zfm )
      CALL wrk_dealloc( jpi, jpj, zfxx, zfyy, zfxy, zalg, zalg1, zalg1q )
      !
   END SUBROUTINE lim_adv_y


   !!======================================================================
END MODULE limadv
