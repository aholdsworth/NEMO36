MODULE trcnam_canoe
   !!======================================================================
   !!                      ***  MODULE trcnam_canoe  ***
   !! TOP :   initialisation of some run parameters for CanOE bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.canoe.h90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'   :                                   CanOE bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_canoe       : CanOE model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_canoe      ! sms trends
   USE trdtrc_oce
   USE iom             ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_canoe   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_canoe.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_canoe
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_canoe  ***  
      !!
      !! ** Purpose :   read CanOE namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !!
      INTEGER :: jl, jn
      INTEGER :: ios                 ! Local integer output status for namelist read
      TYPE(DIAG), DIMENSION(jp_canoe_2d)  :: pisdia2d
      TYPE(DIAG), DIMENSION(jp_canoe_3d)  :: pisdia3d
      TYPE(DIAG), DIMENSION(jp_canoe_trd) :: pisdiabio
      CHARACTER(LEN=20)   ::   clname
      !!
      NAMELIST/nampisbio/ nrdttrc, wsbio, wsbio2, wsbioc
      NAMELIST/nampisdia/ pisdia3d, pisdia2d     ! additional diagnostics
      NAMELIST/nampisdmp/ ln_pisdmp, nn_pisdmp!, ln_pisclo

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_canoe'
!#if defined 1      
      IF(lwp) WRITE(numout,*) ' trc_nam_canoe : read CanOE namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
!#endif
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonp     , 'output.namelist.can' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
         !
         ! Namelist nampisbio
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : Pisces biognostics
         READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : Pisces biognostics
         READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisbio )
         
         IF(lwp) THEN                         ! control print
            WRITE(numout,*) ' Namelist : nampisbio'
            WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
            WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
            WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
            WRITE(numout,*) '    CaCO3 sinking speed                       wsbioc    =', wsbioc
         ENDIF

       ! Namelist nampisdmp
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist nampisdmp in reference namelist : Pisces dmpgnostics
         READ  ( numnatp_ref, nampisdmp, IOSTAT = ios, ERR = 905)
905      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist nampisdmp in configuration namelist : Pisces dmpgnostics
         READ  ( numnatp_cfg, nampisdmp, IOSTAT = ios, ERR = 906 )
906      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisdmp )
         IF(lwp) THEN                         ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : nampisdmp'
            WRITE(numout,*) '    Relaxation of tracer to glodap mean value             ln_pisdmp      =', ln_pisdmp
            WRITE(numout,*) '    Frequency of Relaxation                               nn_pisdmp      =', nn_pisdmp
!            WRITE(numout,*) '    Restoring of tracer to initial value  on closed seas  ln_pisclo      =', ln_pisclo
            WRITE(numout,*) ' '
        ENDIF


      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         ! Namelist nampisdia
         ! -------------------
         REWIND( numnatp_ref )              ! Namelist nampisdia in reference namelist : CanOE diagnostics
         READ  ( numnatp_ref, nampisdia, IOSTAT = ios, ERR = 903)
903      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist nampisdia in configuration namelist : CanOE diagnostics
         READ  ( numnatp_cfg, nampisdia, IOSTAT = ios, ERR = 904 )
904      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdia in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisdia )

         
         DO jl = 1, jp_canoe_2d
            jn = jp_can0_2d + jl - 1
            ctrc2d(jn) = pisdia2d(jl)%sname
            ctrc2l(jn) = pisdia2d(jl)%lname
            ctrc2u(jn) = pisdia2d(jl)%units
         END DO

         DO jl = 1, jp_canoe_3d
            jn = jp_can0_3d + jl - 1
            ctrc3d(jn) = pisdia3d(jl)%sname
            ctrc3l(jn) = pisdia3d(jl)%lname
            ctrc3u(jn) = pisdia3d(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : natadd'
            DO jl = 1, jp_canoe_3d
               jn = jp_can0_3d + jl - 1
               WRITE(numout,*) '  3d diag nb : ', jn, '    short name : ', ctrc3d(jn), &
                 &             '  long name  : ', ctrc3l(jn), '   unit : ', ctrc3u(jn)
            END DO
            WRITE(numout,*) ' '

            DO jl = 1, jp_canoe_2d
               jn = jp_can0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 &             '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
     ENDIF
  !
      

   END SUBROUTINE trc_nam_canoe


   !!======================================================================
END MODULE trcnam_canoe
