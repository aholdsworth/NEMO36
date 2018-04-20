MODULE canflx
   !!======================================================================
   !!                         ***  MODULE canflx  ***
   !! TOP :   CanOE CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   can_flx_init  :   Read the namelist
   !!   can_patm      :   Read sfc atm pressure [atm] for each grid cell
   !!----------------------------------------------------------------------
   USE oce_trc                      !  shared variables between ocean and passive tracers 
   USE trc                          !  passive tracers common variables
   USE sms_canoe                   !  CanOE Source Minus Sink variables
   USE canche                       !  Chemical model
   USE prtctl_trc                   !  print control for debugging
   USE iom                          !  I/O manager
   USE fldread                      !  read input fields




   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_flx  
   PUBLIC   can_flx_init  
   PUBLIC   can_flx_alloc  

   !                               !!** Namelist  nampisext  **
   REAL(wp)          ::  atcco2     !: pre-industrial atmospheric [co2] (ppm) 	
   LOGICAL           ::  ln_co2int  !: flag to read in a file and interpolate atmospheric pco2 or not
   CHARACTER(len=34) ::  clname     !: filename of pco2 values
   INTEGER           ::  nn_offset  !: Offset model-data start year (default = 0) 

   !!  Variables related to reading atmospheric CO2 time history    
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) :: atcco2h, years
   INTEGER  :: nmaxrec, numco2

   !                               !!* nampisatm namelist (Atmospheric PRessure) *
   LOGICAL, PUBLIC ::   ln_presatm  !: ref. pressure: global mean Patm (F) or a constant (F)

   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:)  ::  patm      ! atmospheric pressure at kt                 [atm]
   TYPE(FLD), ALLOCATABLE,       DIMENSION(:)    ::  sf_patm   ! structure of input fields (file informations, fields read)


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_co2   !: ocean carbon flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmco2   !: atmospheric pco2 

   REAL(wp) ::  t_oce_co2_flx               !: Total ocean carbon flux 
   REAL(wp) ::  t_atm_co2_flx               !: global mean of atmospheric pco2
   REAL(wp) ::  area                        !: ocean surface
   REAL(wp) ::  xconv  = 0.01_wp / 3600._wp !: coefficients for conversion 

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
   !! $Id: canflx.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_flx ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : 
      !!              - Include total atm P correction via Esbensen & Kushnir (1981) 
      !!              - Pressure correction NOT done for key_cpl_carbon_cycle
      !!              - Remove Wanninkhof chemical enhancement;
      !!              - Add option for time-interpolation of atcco2.txt  
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zph, zah2, zbot, zdic, zalk, zsch_o2, zalka, zsch_co2
      REAL(wp) ::   zyr_dec, zdco2dt
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:) :: zkgco2, zkgo2, zh2co3, zoflx, zw2d 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_flx')
      !
      CALL wrk_alloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx )
      !

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      IF( kt /= nit000 .AND. knt == 1 ) CALL can_patm( kt )    ! Get sea-level pressure (E&K [1981] climatology) for use in flux calcs

      IF( ln_co2int ) THEN 
         ! Linear temporal interpolation  of atmospheric pco2.  atcco2.txt has annual values.
         ! Caveats: First column of .txt must be in years, decimal  years preferably. 
         ! For nn_offset, if your model year is iyy, nn_offset=(years(1)-iyy) 
         ! then the first atmospheric CO2 record read is at years(1)
         zyr_dec = REAL( nyear + nn_offset, wp ) + REAL( nday_year, wp ) / REAL( nyear_len(1), wp )
         jm = 2
         DO WHILE( jm <= nmaxrec .AND. years(jm-1) < zyr_dec .AND. years(jm) >= zyr_dec ) ;  jm = jm + 1 ;  END DO
         iind = jm  ;   iindm1 = jm - 1
         zdco2dt = ( atcco2h(iind) - atcco2h(iindm1) ) / ( years(iind) - years(iindm1) + rtrn )
         atcco2  = zdco2dt * ( zyr_dec - years(iindm1) ) + atcco2h(iindm1)
         satmco2(:,:) = atcco2 
      ENDIF


      DO jm = 1, 10
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
               zbot  = borat(ji,jj,1)
               zfact = rhop(ji,jj,1) / 1000. + rtrn
               zdic  = trb(ji,jj,1,jpdic) / zfact
               zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact
               zalka = trb(ji,jj,1,jptal) / zfact

               ! CALCULATE [ALK]([CO3--], [HCO3-])
               zalk  = zalka - (  akw3(ji,jj,1) / zph - zph + zbot / ( 1.+ zph / akb3(ji,jj,1) )  )

               ! CALCULATE [H+] AND [H2CO3]
               zah2   = SQRT(  (zdic-zalk)**2 + 4.* ( zalk * ak23(ji,jj,1)   &
                  &                                        / ak13(ji,jj,1) ) * ( 2.* zdic - zalk )  )
               zah2   = 0.5 * ak13(ji,jj,1) / zalk * ( ( zdic - zalk ) + zah2 )
               zh2co3(ji,jj) = ( 2.* zdic - zalk ) / ( 2.+ ak13(ji,jj,1) / zah2 ) * zfact
               hi(ji,jj,1)   = zah2 * zfact
            END DO
         END DO
      END DO


      ! --------------
      ! COMPUTE FLUXES
      ! --------------

      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! -------------------------------------------

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztc  = MIN( 35., tsn(ji,jj,1,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2073.1 - 125.62 * ztc + 3.6276 * ztc2 - 0.043126 * ztc3
            zsch_o2  = 1953.4 - 128.0  * ztc + 3.9918 * ztc2 - 0.050091 * ztc3
            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
            zkgwan = 0.3 * zws  + 2.5 * ( 0.5246 + 0.016256 * ztc + 0.00049946  * ztc2 )
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1)
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Compute CO2 flux for the sea and air
            zfld = satmco2(ji,jj) * patm(ji,jj) * tmask(ji,jj,1) * chemc(ji,jj,1) * zkgco2(ji,jj)   ! (mol/L) * (m/s)
            zflu = zh2co3(ji,jj) * tmask(ji,jj,1) * zkgco2(ji,jj)                                   ! (mol/L) (m/s) ?
            oce_co2(ji,jj) = ( zfld - zflu ) * rfact2 * e1e2t(ji,jj) * tmask(ji,jj,1) * 1000.
            ! compute the trend
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + ( zfld - zflu )*rfact2 / e3t_n(ji,jj,1)

            ! Compute O2 flux 
            zfld16 = atcox * patm(ji,jj) * chemc(ji,jj,2) * tmask(ji,jj,1) * zkgo2(ji,jj)          ! (mol/L) * (m/s)
            zflu16 = trb(ji,jj,1,jpoxy) * tmask(ji,jj,1) * zkgo2(ji,jj)*1.e-6
!convert to mol/L
            zoflx(ji,jj) = zfld16 - zflu16
            tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy) + zoflx(ji,jj) *rfact2 / e3t_n(ji,jj,1)*1.e6
         END DO
      END DO

      t_oce_co2_flx = t_oce_co2_flx + glob_sum( oce_co2(:,:) )            ! Cumulative Total Flux of Carbon
 
      IF( kt == nitend ) THEN
         t_atm_co2_flx = glob_sum( satmco2(:,:) * patm(:,:) * e1e2t(:,:) )            ! Total atmospheric pCO2
         !
         t_oce_co2_flx = (-1.) * t_oce_co2_flx  * 12. / 1.e15             ! Conversion in PgC ; negative for out of the ocean
         t_atm_co2_flx = t_atm_co2_flx  / area                            ! global mean of atmospheric pCO2
         !
         IF( lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' Global mean of atmospheric pCO2 (ppm) at it= ', kt, ' date= ', ndastp
            WRITE(numout,*) '------------------------------------------------------- :  ',t_atm_co2_flx
            WRITE(numout,*)
            WRITE(numout,*) ' Cumulative total Flux of Carbon out of the ocean (PgC) :'
            WRITE(numout,*) '-------------------------------------------------------  ',t_oce_co2_flx
         ENDIF
     ENDIF
         !
         !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         CALL wrk_alloc( jpi, jpj, zw2d )  
         IF( iom_use( "Cflx"  ) )  THEN
           !zw2d(:,:) =  satmco2(:,:) * tmask(:,:,1) 
            zw2d(:,:) = oce_co2(:,:) / e1e2t(:,:) / rfact2r
            CALL iom_put( "Cflx"     , zw2d ) 
         ENDIF
         IF( iom_use( "Oflx"  ) )  THEN
            zw2d(:,:) =  zoflx(:,:) * 1000 * tmask(:,:,1)
           !zw2d(:,:) =  chemc(:,:,1)  * tmask(:,:,1) 
            CALL iom_put( "Oflx" , zw2d )
         ENDIF
         IF( iom_use( "Kg"    ) )  THEN
            zw2d(:,:) =  zkgco2(:,:) * tmask(:,:,1)
            CALL iom_put( "Kg"   , zw2d )
         ENDIF
         IF( iom_use( "Dpo2" ) ) THEN
           zw2d(:,:) = (atcox * patm(:,:) - trb(:,:,1,jpoxy) / ( chemc(:,:,2) + rtrn ) ) * tmask(:,:,1)
           !zw2d(:,:) =  patm(:,:)  * tmask(:,:,1) 
           CALL iom_put( "Dpo2" ,  zw2d )
         ENDIF
         IF( iom_use( "Dpco2" ) )  THEN
           zw2d(:,:) = ( zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) - satmco2(:,:) * patm(:,:) ) * tmask(:,:,1) 
           CALL iom_put( "Dpco2"  , zw2d )
         ENDIF
         IF( iom_use( "pco2" ) )  THEN
            !zw2d(:,:) =   zh2co3(:,:)*tmask(:,:,1)   ! molC/s
            zw2d(:,:) =   (zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) )*tmask(:,:,1)   ! molC/s
         CALL iom_put("pco2" , zw2d )      ! molC
         ENDIF
         !
         CALL wrk_dealloc( jpi, jpj, zw2d )
      ELSE
         IF( ln_diatrc ) THEN
            trc2d(:,:,jp_can0_2d    ) = oce_co2(:,:) / e1e2t(:,:) * rfact2r 
            trc2d(:,:,jp_can0_2d + 1) = zoflx(:,:) * 1000 * tmask(:,:,1) 
            trc2d(:,:,jp_can0_2d + 2) = zkgco2(:,:) * tmask(:,:,1) 
            trc2d(:,:,jp_can0_2d + 3) = ( satmco2(:,:) * patm(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1) 
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx )
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_flx')
      !
   END SUBROUTINE can_flx


   SUBROUTINE can_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !! ** input   :   Namelist nampisext
      !!----------------------------------------------------------------------
      NAMELIST/nampisext/ln_co2int, atcco2, clname, nn_offset
      INTEGER :: jm
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !

      REWIND( numnatp_ref )              ! Namelist nampisext in reference namelist : Pisces atm. conditions
      READ  ( numnatp_ref, nampisext, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisext in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisext in configuration namelist : Pisces atm. conditions
      READ  ( numnatp_cfg, nampisext, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisext in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for air-sea exchange, nampisext'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Choice for reading in the atm pCO2 file or constant value, ln_co2int =', ln_co2int
         WRITE(numout,*) ' '
      ENDIF
      IF( .NOT.ln_co2int ) THEN
         IF(lwp) THEN                         ! control print
            WRITE(numout,*) '    Constant Atmospheric pCO2 value  atcco2    =', atcco2
            WRITE(numout,*) ' '
         ENDIF
         satmco2(:,:)  = atcco2      ! Initialisation of atmospheric pco2
      ELSE
         IF(lwp)  THEN
            WRITE(numout,*) '    Atmospheric pCO2 value  from file clname      =', TRIM( clname )
            WRITE(numout,*) '    Offset model-data start year      nn_offset   =', nn_offset
            WRITE(numout,*) ' '
         ENDIF
         CALL ctl_opn( numco2, TRIM( clname) , 'OLD', 'FORMATTED', 'SEQUENTIAL', -1 , numout, lwp )
         jm = 0                      ! Count the number of record in co2 file
         DO
           READ(numco2,*,END=100) 
           jm = jm + 1
         END DO
 100     nmaxrec = jm - 1 
         ALLOCATE( years  (nmaxrec) )     ;      years  (:) = 0._wp
         ALLOCATE( atcco2h(nmaxrec) )     ;      atcco2h(:) = 0._wp

         REWIND(numco2)
         DO jm = 1, nmaxrec          ! get  xCO2 data
            READ(numco2, *)  years(jm), atcco2h(jm)
            IF(lwp) WRITE(numout, '(f6.0,f7.2)')  years(jm), atcco2h(jm)
         END DO
         CLOSE(numco2)
      ENDIF
      !
      area = glob_sum( e1e2t(:,:) )        ! interior global domain surface
      oce_co2(:,:)  = 0._wp                ! Initialization of Flux of Carbon
      t_oce_co2_flx = 0._wp
      t_atm_co2_flx = 0._wp
      !
      CALL can_patm( nit000 )
      !
   END SUBROUTINE can_flx_init

   SUBROUTINE can_patm( kt )

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_atm  ***
      !!
      !! ** Purpose :   Read and interpolate the external atmospheric sea-levl pressure
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !
      INTEGER            ::  ierr
      INTEGER            ::  ios      ! Local integer output status for namelist read
      CHARACTER(len=100) ::  cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::  sn_patm  ! informations about the fields to be read
      !!
      NAMELIST/nampisatm/ ln_presatm, sn_patm, cn_dir

      !                                         ! ----------------------- !
      IF( kt == nit000 ) THEN                   ! First call kt=nittrc000 !

         REWIND( numnatp_ref )              ! Namelist nampisatm in reference namelist : Pisces atm. sea level pressure file
         READ  ( numnatp_ref, nampisatm, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisatm in reference namelist', lwp )

         REWIND( numnatp_cfg )              ! Namelist nampisatm in configuration namelist : Pisces atm. sea level pressure file 
         READ  ( numnatp_cfg, nampisatm, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisatm in configuration namelist', lwp )
         IF(lwm) WRITE ( numonp, nampisatm )
         !
         !
         IF(lwp) THEN                                 !* control print
            WRITE(numout,*)
            WRITE(numout,*) '   Namelist nampisatm : Atmospheric Pressure as external forcing'
            WRITE(numout,*) '      constant atmopsheric pressure (F) or from a file (T)  ln_presatm = ', ln_presatm
            WRITE(numout,*)
         ENDIF
         !
         IF( ln_presatm ) THEN
            ALLOCATE( sf_patm(1), STAT=ierr )           !* allocate and fill sf_patm (forcing structure) with sn_patm
            IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'can_flx: unable to allocate sf_patm structure' )
            !
            CALL fld_fill( sf_patm, (/ sn_patm /), cn_dir, 'can_flx', 'Atmospheric pressure ', 'nampisatm' )
                                   ALLOCATE( sf_patm(1)%fnow(jpi,jpj,1)   )
            IF( sn_patm%ln_tint )  ALLOCATE( sf_patm(1)%fdta(jpi,jpj,1,2) )
         ENDIF
         !                                         
         IF( .NOT.ln_presatm )   patm(:,:) = 1.e0    ! Initialize patm if no reading from a file
         !
      ENDIF
      !
      IF( ln_presatm ) THEN
            
         CALL fld_read( kt, 1, sf_patm )               !* input Patm provided at kt + 1/2
         patm(:,:) = sf_patm(1)%fnow(:,:,1) /101325._wp                       ! atmospheric pressure converted from Pa to atmospheres
      ENDIF
      !
   END SUBROUTINE can_patm

   INTEGER FUNCTION can_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( oce_co2(jpi,jpj), satmco2(jpi,jpj), patm(jpi,jpj), STAT=can_flx_alloc )
      !
      IF( can_flx_alloc /= 0 )   CALL ctl_warn('can_flx_alloc : failed to allocate arrays')
      !
   END FUNCTION can_flx_alloc


   !!======================================================================
END MODULE  canflx
