MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zfechem
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

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
   !! $Id: p4zbio.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_bio')
      !
      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

          
      CALL p4z_opt  ( kt, knt )     ! Optic: PAR in the water column
      CALL p4z_sink ( kt, knt )     ! vertical flux of particulate organic matter
      CALL p4z_fechem(kt, knt )     ! Iron chemistry/scavenging
      CALL p4z_lim  ( kt, knt )     ! co-limitations by the various nutrients
      CALL p4z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
      !                             ! (for each element : C, Si, Fe, Chl )
      CALL p4z_mort ( kt      )     ! phytoplankton mortality
     !                             ! zooplankton sources/sinks routines 
      CALL p4z_micro( kt, knt )           ! microzooplankton
      CALL p4z_meso ( kt, knt )           ! mesozooplankton
      CALL p4z_rem  ( kt, knt )     ! remineralization terms of organic matter+scavenging of Fe
      !                             ! test if tracers concentrations fall below 0.
      !                                                             !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_bio')
      !
   END SUBROUTINE p4z_bio


   !!======================================================================
END MODULE p4zbio
