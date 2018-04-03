MODULE trcnam_pisces
   !!======================================================================
   !!                      ***  MODULE trcnam_pisces  ***
   !! TOP :   initialisation of some run parameters for PISCES bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.pisces.h90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'   :                                   PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_pisces       : PISCES model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_pisces      ! sms trends
   USE trdtrc_oce
   USE iom             ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_pisces   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_pisces.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_pisces  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !!
      INTEGER :: jl, jn
      INTEGER :: ios                 ! Local integer output status for namelist read
      TYPE(DIAG), DIMENSION(jp_pisces_2d)  :: pisdia2d
      TYPE(DIAG), DIMENSION(jp_pisces_3d)  :: pisdia3d
      TYPE(DIAG), DIMENSION(jp_pisces_trd) :: pisdiabio
      CHARACTER(LEN=20)   ::   clname
      !!
      NAMELIST/nampisdia/ pisdia3d, pisdia2d     ! additional diagnostics




      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_pisces'

      IF(lwp) WRITE(numout,*) ' trc_nam_pisces : read PISCES namelist'



      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.pis' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist nampisdia
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist nampisdia in reference namelist : Pisces diagnostics
         READ  ( numnatp_ref, nampisdia, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist nampisdia in configuration namelist : Pisces diagnostics
         READ  ( numnatp_cfg, nampisdia, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisdia )

         DO jl = 1, jp_pisces_2d
            jn = jp_pcs0_2d + jl - 1
            ctrc2d(jn) = pisdia2d(jl)%sname
            ctrc2l(jn) = pisdia2d(jl)%lname
            ctrc2u(jn) = pisdia2d(jl)%units
         END DO

         DO jl = 1, jp_pisces_3d
            jn = jp_pcs0_3d + jl - 1
            ctrc3d(jn) = pisdia3d(jl)%sname
            ctrc3l(jn) = pisdia3d(jl)%lname
            ctrc3u(jn) = pisdia3d(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : natadd'
            DO jl = 1, jp_pisces_3d
               jn = jp_pcs0_3d + jl - 1
               WRITE(numout,*) '  3d diag nb : ', jn, '    short name : ', ctrc3d(jn), &
                 &             '  long name  : ', ctrc3l(jn), '   unit : ', ctrc3u(jn)
            END DO
            WRITE(numout,*) ' '

            DO jl = 1, jp_pisces_2d
               jn = jp_pcs0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 &             '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF


   END SUBROUTINE trc_nam_pisces


   !!======================================================================
END MODULE trcnam_pisces
