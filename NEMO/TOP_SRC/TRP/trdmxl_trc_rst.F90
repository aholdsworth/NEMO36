MODULE trdmxl_trc_rst
   !!======================================================================
   !!                       ***  MODULE  trdmxl_rst  ***
   !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
   !!======================================================================
   !! History :  9.0  ! 07-03 (C. Deltel) Original code
   !!----------------------------------------------------------------------
  
#if defined key_top && defined key_trdmxl_trc
   !!----------------------------------------------------------------------
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE trc             ! for nn_dttrc ctrcnm
   USE trdmxl_trc_oce  ! for lk_trdmxl_trc

   IMPLICIT NONE
   PRIVATE
  
   PUBLIC   trd_mxl_trc_rst_read    ! routine called by trd_mxl_init
   PUBLIC   trd_mxl_trc_rst_write   ! routine called by step.F90
  
   INTEGER ::   nummldw_trc               ! logical unit for mld restart
   !!---------------------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trdmxl_trc_rst.F90 5341 2015-06-03 14:59:46Z davestorkey $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------
  
CONTAINS
  

    SUBROUTINE trd_mxl_trc_rst_write( kt )
      !!--------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE trd_mxl_rst_wri  ***
      !!                
      !! ** Purpose :   Write mixed-layer diagnostics restart fields.
      !!--------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
      !
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step deine as a character
      CHARACTER(LEN=50)   ::   clname   ! output restart file name
      CHARACTER(LEN=256)  ::   clpath   ! full path to restart file
      CHARACTER (len=35) :: charout
      INTEGER :: jl,  jk, jn               ! loop indice
      !!--------------------------------------------------------------------------------

      IF( kt == nitrst - nn_dttrc .OR. nitend - nit000 + 1 < 2 * nn_dttrc ) THEN ! idem trcrst.F90
         IF( nitrst > 1.0e9 ) THEN
            WRITE(clkt,*) nitrst
         ELSE
           WRITE(clkt,'(i8.8)') nitrst
         ENDIF
         clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_trdrst_trc_out)
         clpath = TRIM(cn_trcrst_outdir)
         IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
         IF(lwp) WRITE(numout,*) '             open ocean restart_mld_trc NetCDF  'TRIM(clpath)//TRIM(clname)
         CALL iom_open( TRIM(clpath)//TRIM(clname), nummldw_trc, ldwrt = .TRUE., kiolib = jprstlib )
      ENDIF

      IF( kt == nitend .AND. lk_trdmxl_trc ) THEN

         IF( kt == nitend .AND. lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'trdmxl_trc_rst: output for ML diags. restart, with trd_mxl_trc_rst_write routine'
            WRITE(numout,*) '~~~~~~~~~~~~~~'
            WRITE(numout,*)
         ENDIF

         IF( ln_trdmxl_trc_instant ) THEN 
            !
            DO jn = 1, jptra
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlbb_trc_'  //ctrcnm(jn), tmlbb_trc  (:,:,jn) )
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlbn_trc_'  //ctrcnm(jn), tmlbn_trc  (:,:,jn) )
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlatfb_trc_'//ctrcnm(jn), tmlatfb_trc(:,:,jn) )
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlradb_trc_'//ctrcnm(jn), tmlradb_trc(:,:,jn) )
            END DO
            !
         ELSE
            !
            CALL iom_rstput( kt, nitrst, nummldw_trc, 'rmldbn_trc', rmldbn_trc )  ! 2D x 1
            
            !                                                          ! ===========
            DO jn = 1, jptra                                           ! tracer loop
               !                                                       ! ===========

               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlatfb_trc_' //ctrcnm(jn), tmlatfb_trc (:,:,jn) ) ! 2D x jptra
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlbb_trc_'   //ctrcnm(jn), tmlbb_trc   (:,:,jn) ) ! 2D x jptra
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlradb_trc_' //ctrcnm(jn), tmlradb_trc (:,:,jn) ) ! 2D x jptra

               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmlbn_trc_'   //ctrcnm(jn), tmlbn_trc   (:,:,jn) ) ! 2D x jptra
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tml_sumb_trc_'//ctrcnm(jn), tml_sumb_trc(:,:,jn) ) ! 2D x jptra
               
               DO jk = 1, jpltrd_trc
                  IF( jk < 10 )   THEN
                     WRITE(charout,FMT="('tmltrd_csum_ub_trc_', A3, '_', I1)") ctrcnm(jn), jk
                  ELSE
                     WRITE(charout,FMT="('tmltrd_csum_ub_trc_', A3, '_', I2)") ctrcnm(jn), jk
                  ENDIF
                  CALL iom_rstput( kt, nitrst, nummldw_trc, charout, tmltrd_csum_ub_trc(:,:,jk,jn) )
               END DO
               
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmltrd_atf_sumb_trc_'//ctrcnm(jn) , &
                    &           tmltrd_atf_sumb_trc(:,:,jn) ) ! 2D x jptra

               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmltrd_rad_sumb_trc_'//ctrcnm(jn) , &
                    &           tmltrd_rad_sumb_trc(:,:,jn) ) ! 2D x jptra
               !                                                       ! ===========
            END DO                                                     ! tracer loop
            !                                                          ! ===========
#if defined key_pisces_reduced
            DO jl = 1, jp_pisces_trd
               CALL iom_rstput( kt, nitrst, nummldw_trc, 'tmltrd_csum_ub_bio'//ctrd_bio(jl,2), tmltrd_csum_ub_bio(:,:,jl) )
            ENDDO
#endif

         ENDIF
         
         CALL iom_close( nummldw_trc )
         lrst_trc = .TRUE.

      ENDIF

    END SUBROUTINE trd_mxl_trc_rst_write


    SUBROUTINE trd_mxl_trc_rst_read
      !!----------------------------------------------------------------------------
      !!                   ***  SUBROUTINE trd_mxl_rst_lec  ***
      !!                   
      !! ** Purpose :   Read file for mixed-layer diagnostics restart.
      !!----------------------------------------------------------------------------
      INTEGER  ::  inum       ! temporary logical unit
      !
      CHARACTER (len=35) :: charout
      INTEGER ::  jk, jn, jl     ! loop indice
      INTEGER ::  jlibalt = jprstlib
      LOGICAL ::  llok
      CHARACTER(LEN=256)  ::   clpath   ! full path to restart file
      !!-----------------------------------------------------------------------------
      
      IF(lwp)  THEN
         WRITE(numout,*)
         WRITE(numout,*) ' trd_mxl_trc_rst_read : read the NetCDF MLD restart file'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~'
      ENDIF
      
      clpath = TRIM(cn_trcrst_indir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'

      IF ( jprstlib == jprstdimg ) THEN
        ! eventually read netcdf file (monobloc)  for restarting on different number of processors
        ! if {cn_trdrst_trc_in}.nc exists, then set jlibalt to jpnf90
        INQUIRE( FILE = TRIM(clpath)//TRIM(cn_trdrst_trc_in)//'.nc', EXIST = llok )
        IF ( llok ) THEN ; jlibalt = jpnf90  ; ELSE ; jlibalt = jprstlib ; ENDIF
      ENDIF

      CALL iom_open( TRIM(clpath)//TRIM(cn_trdrst_trc_in), inum, kiolib = jlibalt ) 
      
      IF( ln_trdmxl_trc_instant ) THEN 
         
         DO jn = 1, jptra
            CALL iom_get( inum, jpdom_autoglo, 'tmlbb_trc_'  //ctrcnm(jn), tmlbb_trc  (:,:,jn) )
            CALL iom_get( inum, jpdom_autoglo, 'tmlbn_trc_'  //ctrcnm(jn), tmlbn_trc  (:,:,jn) )
            CALL iom_get( inum, jpdom_autoglo, 'tmlatfb_trc_'//ctrcnm(jn), tmlatfb_trc(:,:,jn) )
            CALL iom_get( inum, jpdom_autoglo, 'tmlradb_trc_'//ctrcnm(jn), tmlradb_trc(:,:,jn) )
         END DO
         
      ELSE
         CALL iom_get( inum, jpdom_autoglo, 'rmldbn_trc', rmldbn_trc ) ! needed for rmld_sum
         
         !                                                          ! ===========
         DO jn = 1, jptra                                           ! tracer loop
            !                                                       ! ===========
            CALL iom_get( inum, jpdom_autoglo, 'tmlatfb_trc_' //ctrcnm(jn), tmlatfb_trc(:,:,jn) )
            CALL iom_get( inum, jpdom_autoglo, 'tmlbb_trc_'   //ctrcnm(jn), tmlbb_trc  (:,:,jn) )
            CALL iom_get( inum, jpdom_autoglo, 'tmlradb_trc_' //ctrcnm(jn), tmlradb_trc(:,:,jn) )

            CALL iom_get( inum, jpdom_autoglo, 'tmlbn_trc_'   //ctrcnm(jn), tmlbn_trc   (:,:,jn) ) ! needed for tml_sum
            CALL iom_get( inum, jpdom_autoglo, 'tml_sumb_trc_'//ctrcnm(jn), tml_sumb_trc(:,:,jn) )
            
            DO jk = 1, jpltrd_trc
               IF( jk < 10 )   THEN
                  WRITE(charout,FMT="('tmltrd_csum_ub_trc_', A3, '_', I1)") ctrcnm(jn), jk
               ELSE
                  WRITE(charout,FMT="('tmltrd_csum_ub_trc_', A3, '_', I2)") ctrcnm(jn), jk
               ENDIF
               CALL iom_get( inum, jpdom_autoglo, charout, tmltrd_csum_ub_trc(:,:,jk,jn) )
            END DO
            
            CALL iom_get( inum, jpdom_autoglo, 'tmltrd_atf_sumb_trc_'//ctrcnm(jn) , &
                 &        tmltrd_atf_sumb_trc(:,:,jn) )

            CALL iom_get( inum, jpdom_autoglo, 'tmltrd_rad_sumb_trc_'//ctrcnm(jn) , &
                 &        tmltrd_rad_sumb_trc(:,:,jn) )
            !                                                       ! ===========
         END DO                                                     ! tracer loop
         !                                                          ! ===========

#if defined key_pisces_reduced
         DO jl = 1, jp_pisces_trd
            CALL iom_get( inum, jpdom_autoglo, 'tmltrd_csum_ub_bio'//ctrd_bio(jl,2), tmltrd_csum_ub_bio(:,:,jl) )
         ENDDO
#endif
         
         CALL iom_close( inum )
      ENDIF
      
    END SUBROUTINE trd_mxl_trc_rst_read
  
#else
  !!=================================================================================
  !!                       ***  MODULE  trdmxl_rst  ***
  !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
  !!=================================================================================
CONTAINS
  SUBROUTINE trd_mxl_trc_rst_opn( kt )
    WRITE(*,*) 'trd_mxl_trc_rst_opn: You should not have seen this print! error?', kt
  END SUBROUTINE trd_mxl_trc_rst_opn
  SUBROUTINE trd_mxl_trc_rst_write( kt )           !  No ML diags ==> empty routine
    WRITE(*,*) 'trd_mxl_trc_rst_wri: You should not have seen this print! error?', kt
  END SUBROUTINE trd_mxl_trc_rst_write
  SUBROUTINE trd_mxl_trc_rst_read                  !  No ML Diags ==> empty routine
  END SUBROUTINE trd_mxl_trc_rst_read
#endif

  !!=================================================================================
END MODULE trdmxl_trc_rst
