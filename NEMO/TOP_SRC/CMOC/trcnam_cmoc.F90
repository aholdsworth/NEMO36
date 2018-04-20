MODULE trcnam_cmoc
   !!======================================================================
   !!                      ***  MODULE trcnam_cmoc  ***
   !! TOP :   initialisation of some run parameters for CMOC bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.cmoc.h90
   !!           CMOC1  !  2013-15 (O. Riche) include CMOC variables and CMOC namelist + NEMO list bug fix
   !!           CMOC1  !  2017-08 (A. Holdsworth) adapted for NEMO 3.6
   !!----------------------------------------------------------------------
#if defined key_cmoc 
   !!----------------------------------------------------------------------
   !!   'key_cmoc'   :                                   CMOC bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_cmoc       : CMOC model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_cmoc      ! sms trends
   USE trdtrc_oce
   USE iom             ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_cmoc   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_cmoc.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_cmoc
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_cmoc  ***  
      !!
      !! ** Purpose :   read CMOC namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !!
      INTEGER :: jl, jn
      INTEGER :: ios                 ! Local integer output status for namelist read
      TYPE(DIAG), DIMENSION(jp_cmoc_2d)  :: pisdia2d
      TYPE(DIAG), DIMENSION(jp_cmoc_3d)  :: pisdia3d
      CHARACTER(LEN=20)   ::   clname
      !!
      NAMELIST/nampisbio/ nrdttrc, xkmort, wsbio,wsbio2
      NAMELIST/nampisdia/ pisdia3d, pisdia2d     ! additional diagnostics

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_cmoc'
      IF(lwp) WRITE(numout,*) ' trc_nam_cmoc : read CMOC (CMOC) namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.pis' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist nampisbio
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : CMOC diagnostics
         READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )

     

         REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : CMOC diagnostics
         READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisbio )
     
         IF(lwp) THEN                         ! control print
            WRITE(numout,*) ' Namelist : nampisbio'
            WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
            WRITE(numout,*) '    half saturation constant for mortality    xkmort    =', xkmort
            WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
            WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
         ENDIF
!      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist nampisdia
         ! -------------------
!         DO jl = 1, jp_pisces_2d
!            WRITE(pisdia2d(jl)%sname,'("2D_",I2.2)') jl                    ! short name ! <NEMO BUG OR 02/18/2014>  [nemo_ticket] [nemo] [4362] v3.4beta:bugfix to correct I/O format of passive tracer name , see ticket #1144
!            WRITE(pisdia2d(jl)%lname,'("2D DIAGNOSTIC NUMBER ",I2.2)') jl  ! long name  ! <NEMO BUG OR 02/18/2014>  [nemo_ticket] [nemo] [4362] v3.4beta:bugfix to correct I/O format of passive tracer na    me , see ticket #1144
!            pisdia2d(jl)%units = ' '                                       ! units
!         END DO
!         !                                 ! 3D output arrays
!         DO jl = 1, jp_pisces_3d
!            WRITE(pisdia3d(jl)%sname,'("3D_",I2.2)') jl                    ! short name ! <NEMO BUG OR 02/18/2014>  [nemo_ticket] [nemo] [4362] v3.4beta:bugfix to correct I/O format of passive tracer na    me , see ticket #1144
!            WRITE(pisdia3d(jl)%lname,'("3D DIAGNOSTIC NUMBER ",I2.2)') jl  ! long name  ! <NEMO BUG OR 02/18/2014>  [nemo_ticket] [nemo] [4362] v3.4beta:bugfix to correct I/O format of passive tracer na        me , see ticket #1144
!            pisdia3d(jl)%units = ' '                                       ! units
!         END DO
!
         REWIND( numnatp_ref )              ! Namelist nampisdia in reference namelist : Pisces diagnostics
         READ  ( numnatp_ref, nampisdia, IOSTAT = ios, ERR = 903)
903      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in reference namelist', lwp )

     

         REWIND( numnatp_cfg )              ! Namelist nampisdia in configuration namelist : Pisces diagnostics
         READ  ( numnatp_cfg, nampisdia, IOSTAT = ios, ERR = 904 )
904      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in configuration namelist', lwp )
     

         DO jl = 1, jp_pisces_2d
            jn = jp_cm0_2d + jl - 1
            ctrc2d(jn) = pisdia2d(jl)%sname
            ctrc2l(jn) = pisdia2d(jl)%lname
            ctrc2u(jn) = pisdia2d(jl)%units
         END DO

         DO jl = 1, jp_pisces_3d
            jn = jp_cm0_3d + jl - 1
            ctrc3d(jn) = pisdia3d(jl)%sname
            ctrc3l(jn) = pisdia3d(jl)%lname
            ctrc3u(jn) = pisdia3d(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : natadd'
            DO jl = 1, jp_pisces_3d
               jn = jp_cm0_3d + jl - 1
               WRITE(numout,*) '  3d diag nb : ', jn, '    short name : ', ctrc3d(jn), &
                 &             '  long name  : ', ctrc3l(jn), '   unit : ', ctrc3u(jn)
            END DO
            WRITE(numout,*) ' '

            DO jl = 1, jp_pisces_2d
               jn = jp_cm0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 &             '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF


   END SUBROUTINE trc_nam_cmoc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No CMOC bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_cmoc                      ! Empty routine
   END  SUBROUTINE  trc_nam_cmoc
#endif  

   !!======================================================================
END MODULE trcnam_cmoc
