MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :   PISCES CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   p4z_flx_init  :   Read the namelist
   !!   p4z_patm      :   Read sfc atm pressure [atm] for each grid cell
   !!----------------------------------------------------------------------
   USE oce_trc                      !  shared variables between ocean and passive tracers 
   USE trc                          !  passive tracers common variables
   USE sms_pisces                   !  PISCES Source Minus Sink variables
   USE p4zche                       !  Chemical model
   USE prtctl_trc                   !  print control for debugging
   USE iom                          !  I/O manager
   USE fldread                      !  read input fields




   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_flx  
   PUBLIC   p4z_flx_init  
   PUBLIC   p4z_flx_alloc  

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

   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:)  ::  patm      ! atmospheric pressure at kt                 [N/m2]
   TYPE(FLD), ALLOCATABLE,       DIMENSION(:)    ::  sf_patm   ! structure of input fields (file informations, fields read)


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_co2   !: ocean carbon flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmco2   !: atmospheric pco2 

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
   !! $Id: p4zflx.F90 7607 2017-01-25 15:37:31Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_flx ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
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
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zvapsw, zsal, zfco2, zxc2, xCO2approx, ztkel, zfugcoeff
      REAL(wp) ::   zph, zdic, zsch_o2, zsch_co2
      REAL(wp) ::   zyr_dec, zdco2dt
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:) :: zkgco2, zkgo2, zh2co3, zoflx, zw2d, zpco2atm 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_flx')
      !
      CALL wrk_alloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx, zpco2atm )
      !

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      IF( kt /= nit000 .AND. knt == 1 ) CALL p4z_patm( kt )    ! Get sea-level pressure (E&K [1981] climatology) for use in flux calcs

      IF( ln_co2int ) THEN 
         ! Linear temporal interpolation  of atmospheric pco2.  atcco2.txt has annual values.
         ! Caveats: First column of .txt must be in years, decimal  years preferably. 
         ! For nn_offset, if your model year is iyy, nn_offset=(years(1)-iyy) 
         ! then the first atmospheric CO2 record read is at years(1)
         zyr_dec = REAL( nyear + nn_offset, wp ) + REAL( nday_year, wp ) / REAL( nyear_len(1), wp )
         jm = 1
         DO WHILE( jm <= nmaxrec .AND. years(jm) < zyr_dec ) ;  jm = jm + 1 ;  END DO
         iind = jm  ;   iindm1 = jm - 1
         zdco2dt = ( atcco2h(iind) - atcco2h(iindm1) ) / ( years(iind) - years(iindm1) + rtrn )
         atcco2  = zdco2dt * ( zyr_dec - years(iindm1) ) + atcco2h(iindm1)
         satmco2(:,:) = atcco2 
      ENDIF


      DO jj = 1, jpj
         DO ji = 1, jpi
            ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
            zfact = rhop(ji,jj,1) / 1000. + rtrn
            zdic  = trb(ji,jj,1,jpdic)
            zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact
            ! CALCULATE [H2CO3]
            zh2co3(ji,jj) = zdic/(1. + ak13(ji,jj,1)/zph + ak13(ji,jj,1)*ak23(ji,jj,1)/zph**2)
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
            ztc4 = ztc2 * ztc2 
            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2116.8 - 136.25 * ztc + 4.7353 * ztc2 - 0.092307 * ztc3 + 0.0007555 * ztc4
            zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.109390 * ztc3 + 0.0009377 * ztc4
            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
            zkgwan = 0.251 * zws
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1)
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ztkel  = tsn(ji,jj,1,jp_tem) + 273.15
            zsal   = tsn(ji,jj,1,jp_sal) + ( 1.- tmask(ji,jj,1) ) * 35.
            zvapsw = EXP(24.4543 - 67.4509*(100.0/ztkel) - 4.8489*LOG(ztkel/100) - 0.000544*zsal)
            zpco2atm(ji,jj) = satmco2(ji,jj) * ( patm(ji,jj) - zvapsw )
            zxc2 = (1.0 - zpco2atm(ji,jj) * 1E-6 )**2
            zfugcoeff = EXP(patm(ji,jj) * (chemc(ji,jj,2) + 2.0 * zxc2 * chemc(ji,jj,3) )   &
            &           / (82.05736 * ztkel))
            zfco2 = zpco2atm(ji,jj) * zfugcoeff

            ! Compute CO2 flux for the sea and air
            zfld = zfco2 * chemc(ji,jj,1) * zkgco2(ji,jj)  ! (mol/L) * (m/s)
            zflu = zh2co3(ji,jj) * zkgco2(ji,jj)                                   ! (mol/L) (m/s) ?
            oce_co2(ji,jj) = ( zfld - zflu ) * rfact2 * e1e2t(ji,jj) * tmask(ji,jj,1) * 1000.
            ! compute the trend
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + ( zfld - zflu ) * rfact2 / e3t_n(ji,jj,1) * tmask(ji,jj,1)

            ! Compute O2 flux 
            zfld16 = patm(ji,jj) * chemo2(ji,jj,1) * zkgo2(ji,jj)          ! (mol/L) * (m/s)
            zflu16 = trb(ji,jj,1,jpoxy) * zkgo2(ji,jj)
            zoflx(ji,jj) = ( zfld16 - zflu16 ) * tmask(ji,jj,1)
            tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy) + zoflx(ji,jj) * rfact2 / e3t_n(ji,jj,1)
         END DO
      END DO

      t_oce_co2_flx     = glob_sum( oce_co2(:,:) )                    !  Total Flux of Carbon
      t_oce_co2_flx_cum = t_oce_co2_flx_cum + t_oce_co2_flx       !  Cumulative Total Flux of Carbon
!      t_atm_co2_flx     = glob_sum( satmco2(:,:) * e1e2t(:,:) )       ! Total atmospheric pCO2
      t_atm_co2_flx     =  atcco2      ! Total atmospheric pCO2
 
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         CALL wrk_alloc( jpi, jpj, zw2d )  
         IF( iom_use( "Cflx"  ) )  THEN
            zw2d(:,:) = oce_co2(:,:) / e1e2t(:,:) * rfact2r
            CALL iom_put( "Cflx"     , zw2d ) 
         ENDIF
         IF( iom_use( "Oflx"  ) )  THEN
            zw2d(:,:) =  zoflx(:,:) * 1000 * tmask(:,:,1)
            CALL iom_put( "Oflx" , zw2d )
         ENDIF
         IF( iom_use( "Kg"    ) )  THEN
            zw2d(:,:) =  zkgco2(:,:) * tmask(:,:,1)
            CALL iom_put( "Kg"   , zw2d )
         ENDIF
         IF( iom_use( "Dpco2" ) ) THEN
           zw2d(:,:) = ( zpco2atm(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1)
           CALL iom_put( "Dpco2" ,  zw2d )
         ENDIF
         IF( iom_use( "Dpo2" ) )  THEN
           zw2d(:,:) = ( atcox * patm(:,:) - atcox * trn(:,:,1,jpoxy) / ( chemo2(:,:,1) + rtrn ) ) * tmask(:,:,1)
           CALL iom_put( "Dpo2"  , zw2d )
         ENDIF
         IF( iom_use( "tcflx" ) )  CALL iom_put( "tcflx"    , t_oce_co2_flx * rfact2r )   ! molC/s
         CALL iom_put( "tcflxcum" , t_oce_co2_flx_cum )      ! molC
         !
         CALL wrk_dealloc( jpi, jpj, zw2d )
      ELSE
         IF( ln_diatrc ) THEN
            trc2d(:,:,jp_pcs0_2d    ) = oce_co2(:,:) / e1e2t(:,:) * rfact2r 
            trc2d(:,:,jp_pcs0_2d + 1) = zoflx(:,:) * 1000 * tmask(:,:,1) 
            trc2d(:,:,jp_pcs0_2d + 2) = zkgco2(:,:) * tmask(:,:,1) 
            trc2d(:,:,jp_pcs0_2d + 3) = ( zpco2atm(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      !
      CALL wrk_dealloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx, zpco2atm )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_flx')
      !
   END SUBROUTINE p4z_flx


   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
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
      oce_co2(:,:)  = 0._wp                ! Initialization of Flux of Carbon
      t_oce_co2_flx = 0._wp
      t_atm_co2_flx = 0._wp
      !
      CALL p4z_patm( nit000 )
      !
   END SUBROUTINE p4z_flx_init

   SUBROUTINE p4z_patm( kt )

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_atm  ***
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
            IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_patm structure' )
            !
            CALL fld_fill( sf_patm, (/ sn_patm /), cn_dir, 'p4z_flx', 'Atmospheric pressure ', 'nampisatm' )
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
         patm(:,:) = sf_patm(1)%fnow(:,:,1)                        ! atmospheric pressure
      ENDIF
      !
   END SUBROUTINE p4z_patm

   INTEGER FUNCTION p4z_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( oce_co2(jpi,jpj), satmco2(jpi,jpj), patm(jpi,jpj), STAT=p4z_flx_alloc )
      !
      IF( p4z_flx_alloc /= 0 )   CALL ctl_warn('p4z_flx_alloc : failed to allocate arrays')
      !
   END FUNCTION p4z_flx_alloc


   !!======================================================================
END MODULE p4zflx
