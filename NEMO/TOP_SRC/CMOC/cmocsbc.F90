MODULE cmocsbc
   !!======================================================================
   !!                         ***  MODULE p4sbc  ***
   !! TOP :   CMOC surface boundary conditions of external inputs of nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   cmoc_sbc_init   :  Initialization of p4z_sbc
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE iom             !  I/O manager
   USE fldread         !  time interpolation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_sbc
   PUBLIC   cmoc_sbc_init   

   !! * Shared module variables
   LOGICAL , PUBLIC  :: ln_river    !: boolean for river input of nutrients

   LOGICAL , PUBLIC  :: ll_sbc

   !! * Module variables
   LOGICAL  ::  ll_solub
   INTEGER , PARAMETER  :: jpriv  = 2   !: Maximum number of river input fields
   INTEGER , PARAMETER  :: jr_dic = 1   !: index of dissolved inorganic carbon
   INTEGER , PARAMETER  :: jr_doc = 2   !: index of dissolved organic carbon

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_river  ! structure of input riverdic


   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file
   INTEGER  :: ntimes_riv


   !! * File structures
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_fmsk      ! structure of input iron limitation mask ! <CMOC code OR 10/22/2015>
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivdic, rivalk    !: river input fields
   REAL(wp), PUBLIC ::  rivalkinput, rivdicinput

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocsbc.F90 5507 2015-06-29 15:19:49Z aumont $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE cmoc_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine cmoc_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER  :: ji,jj 
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sbc')
      !
      ! N release due to coastal rivers
      ! Compute river at nit000 or only if there is more than 1 time record in river file
      ! -----------------------------------------
      IF( ln_river ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_riv > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_river )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zcoef = ryyss * e1e2t(ji,jj) * h_rnf(ji,jj) 
                  rivalk(ji,jj) =   sf_river(jr_dic)%fnow(ji,jj,1)                                    &
                     &              * 1.E3        / ( 12. * zcoef + rtrn )
                  rivdic(ji,jj) = ( sf_river(jr_dic)%fnow(ji,jj,1) + sf_river(jr_doc)%fnow(ji,jj,1) ) &
                     &              * 1.E3         / ( 12. * zcoef + rtrn )
               END DO
            END DO
         ENDIF
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sbc')
!      !
   END SUBROUTINE cmoc_sbc



   SUBROUTINE cmoc_sbc_init

      !!----------------------------------------------------------------------
      !!                  ***  routine cmoc_sbc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!Adapted from CMOC NEMO3.4 AMH 08-2017
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ji, jj, jk, jm,ifpr
      INTEGER  :: numfmsk, numriv                            ! <CMOC code OR 10/22/2015> add numfmsk (iron limitation mask)
      INTEGER  :: ios                 ! Local integer output status for namelist read
      INTEGER  :: ierr, ierr1, ierr2
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      REAL(wp) ::  ztimes_riv 
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zdust, zndepo, zriver, zcmask
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: zfmask      ! <CMOC code OR 10/22/2015> add zfmask for iron limitation mask
      REAL(wp), DIMENSION(:), ALLOCATABLE :: rivinput
      !
      CHARACTER(len=100) ::  cn_dir                          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(jpriv) ::  slf_river    ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_fmsk, sn_riverdoc, sn_riverdic     ! <CMOC code OR 10/22/2015> add sn_fmsk (iron limitation mask)        ! informations about the fields to be read

      NAMELIST/nampissbc/cn_dir, sn_riverdic, sn_riverdoc,                         &
        &                ln_river,                                                 &
        &                sn_fmsk
      NAMELIST/namcmocnfx/ phinf_cmoc, phi0_cmoc, anf_cmoc, pnf_cmoc, inf_cmoc, tnfMa_cmoc, tnfmi_cmoc
  !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sbc_init')
      !
             
      !                            !* set file information
      REWIND( numnatp_ref )              ! Namelist nampissbc in reference namelist : Pisces external sources of nutrients
      READ  ( numnatp_ref, nampissbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampissbc in configuration namelist : Pisces external sources of nutrients
      READ  ( numnatp_cfg, nampissbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissbc )

      REWIND( numnatp_ref )
      READ  ( numnatp_ref,  namcmocnfx, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocnfx in reference namelist', lwp )
      
      REWIND( numnatp_cfg )
      READ  ( numnatp_cfg,  namcmocnfx, IOSTAT = ios, ERR = 904)
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocnfx in reference namelist', lwp )

      IF(lwm) WRITE ( numonp, namcmocnfx )

      IF(lwp) THEN
         WRITE(numout,*) ' namelist : nampissbc '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    river input of nutrients          ln_river    = ', ln_river
      END IF

      IF(  ln_river  ) THEN  ;  ll_sbc = .TRUE.
      ELSE                                            ;  ll_sbc = .FALSE.
      ENDIF



     ! nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         !
         slf_river(jr_dic) = sn_riverdic  ;  slf_river(jr_doc) = sn_riverdoc   
         !
         ALLOCATE(rivdic(jpi,jpj), rivalk(jpi,jpj)) 
         !
         ALLOCATE( sf_river(jpriv), rivinput(jpriv), STAT=ierr1 )           !* allocate and fill sf_river (forcing structure) with sn_river_
         rivinput(:) = 0.0

         IF( ierr1 > 0 )   CALL ctl_stop( 'STOP', 'cmoc_sbc_init: unable to allocate sf_irver structure' )
         !
         CALL fld_fill( sf_river, slf_river, cn_dir, 'cmoc_sbc_init', 'Input from river ', 'nampissbc' )
         DO ifpr = 1, jpriv
                                          ALLOCATE( sf_river(ifpr)%fnow(jpi,jpj,1  ) )
            IF( slf_river(ifpr)%ln_tint ) ALLOCATE( sf_river(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         IF( Agrif_Root() ) THEN   !  Only on the master grid
            ! Get total input rivers ; need to compute total river supply in a year
            DO ifpr = 1, jpriv
               CALL iom_open ( TRIM( slf_river(ifpr)%clname ), numriv )
               CALL iom_gettime( numriv, zsteps, kntime=ntimes_riv)
               ALLOCATE( zriver(jpi,jpj,ntimes_riv) )
               DO jm = 1, ntimes_riv
                  CALL iom_get( numriv, jpdom_data, TRIM( slf_river(ifpr)%clvar ), zriver(:,:,jm), jm )
               END DO
               CALL iom_close( numriv )
               ztimes_riv = 1._wp / FLOAT(ntimes_riv) 
               DO jm = 1, ntimes_riv
                  rivinput(ifpr) = rivinput(ifpr) + glob_sum( zriver(:,:,jm) * tmask(:,:,1) * ztimes_riv ) 
               END DO
               DEALLOCATE( zriver)
            END DO
            ! N/P and Si releases due to coastal rivers
            ! -----------------------------------------
            rivdicinput = (rivinput(jr_dic) + rivinput(jr_doc) ) * 1E3 / 12._wp
            rivalkinput = rivinput(jr_dic) * 1E3 / 12._wp
            !
         ENDIF
      ELSE
         rivdicinput = 0._wp
         rivalkinput = 0._wp
      END IF 
      
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampissbc '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    river input of nutrients                 ln_river    = ', ln_river
         WRITE(numout,*) '    Alk Supply : ', rivalkinput*1E3/1E12,' Teq/yr'
         WRITE(numout,*) '    DIC Supply : ', rivdicinput*1E3*12./1E12,'TgC/yr'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for dinitrogen fix. , namcmocnfx'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum ref. diazotroph concentration    phinf_cmoc =',  phinf_cmoc
         WRITE(numout,*) '    Surface ref. diazotroph concentration     phi0_cmoc =',   phi0_cmoc
         WRITE(numout,*) '    Inverse depth of diazotroph conc.max.      anf_cmoc =',    anf_cmoc
         WRITE(numout,*) '    Maximum ref. rate of dinitrogen fix.       pnf_cmoc =',    pnf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix surf. irr.     inf_cmoc =',    inf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix SST          tnfMa_cmoc =',  tnfMa_cmoc
         WRITE(numout,*) '    Minimum ref. dinitrogen fix SST          tnfmi_cmoc =',  tnfmi_cmoc
         WRITE(numout,*) ' '
       END IF

      IF( ln_river ) THEN
           ll_sbc = .TRUE.
      ELSE
           ll_sbc = .FALSE.
      ENDIF

      ! iron mask/iron limitation for the ocean ! <CMOC code OR 10/22/2015>
      ! ---------------------------------------
         IF(lwp) WRITE(numout,*) '    initialize CMOC iron limitation'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ierr = 0 
         ALLOCATE( sf_fmsk(1), STAT=ierr )           !* allocate and fill sf_fmsk (forcing structure) with sn_fmsk
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'cmoc_sbc_init: unable to allocate sf_fmsk structure' )
         ierr = 0
         
         ! I think this section copy the file structure in memory
         CALL fld_fill( sf_fmsk, (/ sn_fmsk /), cn_dir, 'cmoc_sbc_init', 'Iron limitation mask', 'nampissbc' )
                                   ALLOCATE( sf_fmsk(1)%fnow(jpi,jpj,1), STAT=ierr )  ! fnow current values based on interpolation (OR)?
                                   IF( ierr > 0 ) THEN
                                            CALL ctl_stop('cmocsbc: iron limitation mask,                  & 
                                                          unable to allocate iron limitation array' )     &
                                            ;    RETURN 
                                   ENDIF 
         !
         ! Open the the channel numfmsk associated with  file 'sn_fmsk%clname'
         CALL iom_open (  TRIM( sn_fmsk%clname ) , numfmsk )
         IF(lwp) WRITE(numout,*) '    iron mask file opened' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         ! allocate dynamically mem space for the mask and get the data with iom_get
         ! jpdom_data must hold some information on the grid tile
         ALLOCATE(zfmask(jpi,jpj))
         CALL iom_get  ( numfmsk, jpdom_data, TRIM( sn_fmsk%clvar ), zfmask(:,:), 1 ) ! <CMOC code OR 10/22/2015> get the data from the iron mask file without worrying about the time step
         ! because data time scale is yearly
         IF(lwp) WRITE(numout,*) '    iron mask data read' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout) 
         CALL iom_close( numfmsk )
           DO ji = 1, jpi
              DO jj = 1, jpj
                    xlimnfecmoc(ji,jj) = zfmask(ji,jj)
              END DO
           END DO
         IF(lwp) WRITE(numout,*) '    iron mask data save in memory' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         DEALLOCATE(zfmask) 
         IF(lwp) WRITE(numout,*) '    iron limitation initialization: done'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)  

!      !
      IF( ll_sbc ) CALL cmoc_sbc( nit000 ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sbc_init')
      !
   END SUBROUTINE cmoc_sbc_init
#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sbc                         ! Empty routine
   END SUBROUTINE cmoc_sbc
#endif 

   !!======================================================================
END MODULE  cmocsbc
