MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces or key_pisces_reduced'                     PISCES model
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE sms_pisces  ! PISCES variables
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

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
CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)           :: cltra
      REAL(wp)                     :: zfact
      INTEGER                      :: ji, jj, jk, jn
      REAL(wp), DIMENSION(jpi,jpj) :: zdic, zo2min, zdepo2min
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn = jp_pcs0, jp_pcs1
         zfact = 1.0e+6 
         IF( jn == jpno3 .OR. jn == jpnh4 ) zfact = rno3 * 1.0e+6 
         IF( jn == jppo4  )                 zfact = po4r * 1.0e+6
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         IF( iom_use( cltra ) )  CALL iom_put( cltra, trn(:,:,:,jn) * zfact )
      END DO

      IF( iom_use( "INTDIC" ) ) THEN                     !   DIC content in kg/m2
         zdic(:,:) = 0.
         DO jk = 1, jpkm1
            zdic(:,:) = zdic(:,:) + trn(:,:,jk,jpdic) * e3t_n(:,:,jk) * tmask(:,:,jk) * 12.
         ENDDO
         CALL iom_put( 'INTDIC', zdic )     
      ENDIF
      !
      IF( iom_use( "O2MIN" ) .OR. iom_use ( "ZO2MIN" ) ) THEN  ! Oxygen minimum concentration and depth 
         zo2min   (:,:) = trn(:,:,1,jpoxy) * tmask(:,:,1)
         zdepo2min(:,:) = gdepw_n(:,:,1)    * tmask(:,:,1)
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( tmask(ji,jj,jk) == 1 ) then
                     IF( trn(ji,jj,jk,jpoxy) < zo2min(ji,jj) ) then
                        zo2min   (ji,jj) = trn(ji,jj,jk,jpoxy)
                        zdepo2min(ji,jj) = gdepw_n(ji,jj,jk)
                     ENDIF
                  ENDIF
               END DO
            END DO
         END DO
         !
         CALL iom_put('O2MIN' , zo2min     )                              ! oxygen minimum concentration
         CALL iom_put('ZO2MIN', zdepo2min  )                              ! depth of oxygen minimum concentration
          !
      ENDIF
      !
   END SUBROUTINE trc_wri_pisces


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_pisces.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
