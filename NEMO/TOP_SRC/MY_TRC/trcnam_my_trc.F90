MODULE trcnam_my_trc
   !!======================================================================
   !!                      ***  MODULE trcnam_my_trc  ***
   !! TOP :   initialisation of some run parameters for MY_TRC bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'   :                                       MY_TRC model
   !!----------------------------------------------------------------------
   !! trc_nam_my_trc      : MY_TRC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_my_trc   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_my_trc.F90 8353 2017-07-19 14:41:00Z lovato $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_my_trc
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_my_trc  ***  
      !!
      !! ** Purpose :   read MY_TRC namelist
      !!
      !!----------------------------------------------------------------------
      INTEGER :: jn
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_my_trc : read MY_TRC namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      ! 
      do jn = jp_myt0 , jp_myt1
         ctrcnm    (jn) = 'NONAME'
         ctrcln    (jn) = 'NO Long Name'
         ctrcun    (jn) = 'NOUNIT'
         ln_trc_ini(jn) = .false.
         ln_trc_wri(jn) = .true.
      enddo
      !
      END SUBROUTINE trc_nam_my_trc
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No MY_TRC
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_my_trc                      ! Empty routine
   END  SUBROUTINE  trc_nam_my_trc
#endif  

   !!======================================================================
END MODULE trcnam_my_trc
