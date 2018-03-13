MODULE sbccpl
   !!======================================================================
   !!                       ***  MODULE  sbccpl  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in coupled mode
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!            3.4  ! 2011_11  (C. Harris) more flexibility + multi-category fields
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   namsbc_cpl      : coupled formulation namlist
   !!   sbc_cpl_init    : initialisation of the coupled exchanges
   !!   sbc_cpl_rcv     : receive fields from the atmosphere over the ocean (ocean only)
   !!                     receive stress from the atmosphere over the ocean (ocean-ice case)
   !!   sbc_cpl_ice_tau : receive stress from the atmosphere over ice
   !!   sbc_cpl_ice_flx : receive fluxes from the atmosphere over ice
   !!   sbc_cpl_snd     : send     fields to the atmosphere
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE sbcapr
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE phycst          ! physical constants







   USE cpl_oasis3      ! OASIS3 coupling
   USE geo2ocean       ! 
   USE oce   , ONLY : tsn, un, vn, sshn, ub, vb, sshb, fraqsr_1lev
   USE albedo          !
   USE in_out_manager  ! I/O manager
   USE iom             ! NetCDF library
   USE lib_mpp         ! distribued memory computing library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE eosbn2
   USE sbcrnf   , ONLY : l_rnfcpl
   USE sbcisf   , ONLY : l_isfcpl

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_cpl_init       ! routine called by sbcmod.F90
   PUBLIC   sbc_cpl_rcv        ! routine called by sbc_ice_lim(_2).F90
   PUBLIC   sbc_cpl_snd        ! routine called by step.F90
   PUBLIC   sbc_cpl_ice_tau    ! routine called by sbc_ice_lim(_2).F90
   PUBLIC   sbc_cpl_ice_flx    ! routine called by sbc_ice_lim(_2).F90
   PUBLIC   sbc_cpl_alloc      ! routine called in sbcice_cice.F90

   INTEGER, PARAMETER ::   jpr_otx1   =  1            ! 3 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jpr_oty1   =  2            ! 
   INTEGER, PARAMETER ::   jpr_otz1   =  3            ! 
   INTEGER, PARAMETER ::   jpr_otx2   =  4            ! 3 atmosphere-ocean stress components on grid 2
   INTEGER, PARAMETER ::   jpr_oty2   =  5            ! 
   INTEGER, PARAMETER ::   jpr_otz2   =  6            ! 
   INTEGER, PARAMETER ::   jpr_itx1   =  7            ! 3 atmosphere-ice   stress components on grid 1
   INTEGER, PARAMETER ::   jpr_ity1   =  8            ! 
   INTEGER, PARAMETER ::   jpr_itz1   =  9            ! 
   INTEGER, PARAMETER ::   jpr_itx2   = 10            ! 3 atmosphere-ice   stress components on grid 2
   INTEGER, PARAMETER ::   jpr_ity2   = 11            ! 
   INTEGER, PARAMETER ::   jpr_itz2   = 12            ! 
   INTEGER, PARAMETER ::   jpr_qsroce = 13            ! Qsr above the ocean
   INTEGER, PARAMETER ::   jpr_qsrice = 14            ! Qsr above the ice
   INTEGER, PARAMETER ::   jpr_qsrmix = 15 
   INTEGER, PARAMETER ::   jpr_qnsoce = 16            ! Qns above the ocean
   INTEGER, PARAMETER ::   jpr_qnsice = 17            ! Qns above the ice
   INTEGER, PARAMETER ::   jpr_qnsmix = 18
   INTEGER, PARAMETER ::   jpr_rain   = 19            ! total liquid precipitation (rain)
   INTEGER, PARAMETER ::   jpr_snow   = 20            ! solid precipitation over the ocean (snow)
   INTEGER, PARAMETER ::   jpr_tevp   = 21            ! total evaporation
   INTEGER, PARAMETER ::   jpr_ievp   = 22            ! solid evaporation (sublimation)
   INTEGER, PARAMETER ::   jpr_sbpr   = 23            ! sublimation - liquid precipitation - solid precipitation
   INTEGER, PARAMETER ::   jpr_semp   = 24            ! solid freshwater budget (sublimation - snow)
   INTEGER, PARAMETER ::   jpr_oemp   = 25            ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jpr_w10m   = 26            ! 10m wind
   INTEGER, PARAMETER ::   jpr_dqnsdt = 27            ! d(Q non solar)/d(temperature)
   INTEGER, PARAMETER ::   jpr_rnf    = 28            ! runoffs
   INTEGER, PARAMETER ::   jpr_cal    = 29            ! calving
   INTEGER, PARAMETER ::   jpr_taum   = 30            ! wind stress module
   INTEGER, PARAMETER ::   jpr_co2    = 31
   INTEGER, PARAMETER ::   jpr_topm   = 32            ! topmeltn
   INTEGER, PARAMETER ::   jpr_botm   = 33            ! botmeltn
   INTEGER, PARAMETER ::   jpr_sflx   = 34            ! salt flux
   INTEGER, PARAMETER ::   jpr_toce   = 35            ! ocean temperature
   INTEGER, PARAMETER ::   jpr_soce   = 36            ! ocean salinity
   INTEGER, PARAMETER ::   jpr_ocx1   = 37            ! ocean current on grid 1
   INTEGER, PARAMETER ::   jpr_ocy1   = 38            !
   INTEGER, PARAMETER ::   jpr_ssh    = 39            ! sea surface height
   INTEGER, PARAMETER ::   jpr_fice   = 40            ! ice fraction          
   INTEGER, PARAMETER ::   jpr_e3t1st = 41            ! first T level thickness 
   INTEGER, PARAMETER ::   jpr_fraqsr = 42            ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jpr_isf    = 43
   INTEGER, PARAMETER ::   jpr_icb    = 44
   INTEGER, PARAMETER ::   jprcv      = 44            ! total number of fields received

   INTEGER, PARAMETER ::   jps_fice   =  1            ! ice fraction sent to the atmosphere
   INTEGER, PARAMETER ::   jps_toce   =  2            ! ocean temperature
   INTEGER, PARAMETER ::   jps_tice   =  3            ! ice   temperature
   INTEGER, PARAMETER ::   jps_tmix   =  4            ! mixed temperature (ocean+ice)
   INTEGER, PARAMETER ::   jps_albice =  5            ! ice   albedo
   INTEGER, PARAMETER ::   jps_albmix =  6            ! mixed albedo
   INTEGER, PARAMETER ::   jps_hice   =  7            ! ice  thickness
   INTEGER, PARAMETER ::   jps_hsnw   =  8            ! snow thickness
   INTEGER, PARAMETER ::   jps_ocx1   =  9            ! ocean current on grid 1
   INTEGER, PARAMETER ::   jps_ocy1   = 10            !
   INTEGER, PARAMETER ::   jps_ocz1   = 11            !
   INTEGER, PARAMETER ::   jps_ivx1   = 12            ! ice   current on grid 1
   INTEGER, PARAMETER ::   jps_ivy1   = 13            !
   INTEGER, PARAMETER ::   jps_ivz1   = 14            !
   INTEGER, PARAMETER ::   jps_co2    = 15
   INTEGER, PARAMETER ::   jps_soce   = 16            ! ocean salinity
   INTEGER, PARAMETER ::   jps_ssh    = 17            ! sea surface height
   INTEGER, PARAMETER ::   jps_qsroce = 18            ! Qsr above the ocean
   INTEGER, PARAMETER ::   jps_qnsoce = 19            ! Qns above the ocean
   INTEGER, PARAMETER ::   jps_oemp   = 20            ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jps_sflx   = 21            ! salt flux
   INTEGER, PARAMETER ::   jps_otx1   = 22            ! 2 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jps_oty1   = 23            ! 
   INTEGER, PARAMETER ::   jps_rnf    = 24            ! runoffs
   INTEGER, PARAMETER ::   jps_taum   = 25            ! wind stress module
   INTEGER, PARAMETER ::   jps_fice2  = 26            ! ice fraction sent to OPA (by SAS when doing SAS-OPA coupling)
   INTEGER, PARAMETER ::   jps_e3t1st = 27            ! first level depth (vvl)
   INTEGER, PARAMETER ::   jps_fraqsr = 28            ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jpsnd      = 28            ! total number of fields sended

   !                                                         !!** namelist namsbc_cpl **
   TYPE ::   FLD_C
      CHARACTER(len = 32) ::   cldes                  ! desciption of the coupling strategy
      CHARACTER(len = 32) ::   clcat                  ! multiple ice categories strategy
      CHARACTER(len = 32) ::   clvref                 ! reference of vector ('spherical' or 'cartesian')
      CHARACTER(len = 32) ::   clvor                  ! orientation of vector fields ('eastward-northward' or 'local grid')
      CHARACTER(len = 32) ::   clvgrd                 ! grids on which is located the vector fields
   END TYPE FLD_C
   ! Send to the atmosphere                           !
   TYPE(FLD_C) ::   sn_snd_temp, sn_snd_alb, sn_snd_thick, sn_snd_crt, sn_snd_co2                        
   ! Received from the atmosphere                     !
   TYPE(FLD_C) ::   sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_dqnsdt, sn_rcv_qsr, sn_rcv_qns, sn_rcv_emp, sn_rcv_rnf
   TYPE(FLD_C) ::   sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2, sn_rcv_icb, sn_rcv_isf                               
   ! Other namelist parameters                        !
   INTEGER     ::   nn_cplmodel            ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
   LOGICAL     ::   ln_usecplmask          !  use a coupling mask file to merge data received from several models
                                           !   -> file cplmask.nc with the float variable called cplmask (jpi,jpj,nn_cplmodel)
   TYPE ::   DYNARR     
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   z3   
   END TYPE DYNARR

   TYPE( DYNARR ), SAVE, DIMENSION(jprcv) ::   frcv                      ! all fields recieved from the atmosphere

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   albedo_oce_mix     ! ocean albedo sent to atmosphere (mix clear/overcast sky)

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(    :) ::   nrcvinfo           ! OASIS info argument

   !! Substitution
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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: sbccpl.F90 7607 2017-01-25 15:37:31Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
  
   INTEGER FUNCTION sbc_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION sbc_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( albedo_oce_mix(jpi,jpj), nrcvinfo(jprcv),  STAT=ierr(1) )
      
      ALLOCATE( a_i(jpi,jpj,1) , STAT=ierr(2) )  ! used in sbcice_if.F90 (done here as there is no sbc_ice_if_init)
      ALLOCATE( xcplmask(jpi,jpj,0:nn_cplmodel) , STAT=ierr(3) )
      !
      sbc_cpl_alloc = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_cpl_alloc )
      IF( sbc_cpl_alloc > 0 )   CALL ctl_warn('sbc_cpl_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_cpl_alloc


   SUBROUTINE sbc_cpl_init( k_ice )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_init  ***
      !!
      !! ** Purpose :   Initialisation of send and received information from
      !!                the atmospheric component
      !!
      !! ** Method  : * Read namsbc_cpl namelist 
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_ice       ! ice management in the sbc (=0/1/2/3)
      !!
      INTEGER ::   jn   ! dummy loop index
      INTEGER ::   ios  ! Local integer output status for namelist read
      INTEGER ::   inum 
      REAL(wp), POINTER, DIMENSION(:,:) ::   zacs, zaos
      !!
      NAMELIST/namsbc_cpl/  sn_snd_temp, sn_snd_alb   , sn_snd_thick, sn_snd_crt   , sn_snd_co2,      &
         &                  sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau  , sn_rcv_dqnsdt, sn_rcv_qsr,      &
         &                  sn_rcv_qns , sn_rcv_emp   , sn_rcv_rnf  , sn_rcv_cal   , sn_rcv_iceflx,   &
         &                  sn_rcv_co2 , sn_rcv_icb , sn_rcv_isf, nn_cplmodel  , ln_usecplmask
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_cpl_init')
      !
      CALL wrk_alloc( jpi,jpj, zacs, zaos )

      ! ================================ !
      !      Namelist informations       !
      ! ================================ !

      REWIND( numnam_ref )              ! Namelist namsbc_cpl in reference namelist : Variables for OASIS coupling
      READ  ( numnam_ref, namsbc_cpl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_cpl in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc_cpl in configuration namelist : Variables for OASIS coupling
      READ  ( numnam_cfg, namsbc_cpl, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_cpl in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_cpl )

      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'sbc_cpl_init : namsbc_cpl namelist '
         WRITE(numout,*)'~~~~~~~~~~~~'
      ENDIF
      IF( lwp .AND. ln_cpl ) THEN                        ! control print
         WRITE(numout,*)'  received fields (mutiple ice categogies)'
         WRITE(numout,*)'      10m wind module                 = ', TRIM(sn_rcv_w10m%cldes  ), ' (', TRIM(sn_rcv_w10m%clcat  ), ')'
         WRITE(numout,*)'      stress module                   = ', TRIM(sn_rcv_taumod%cldes), ' (', TRIM(sn_rcv_taumod%clcat), ')'
         WRITE(numout,*)'      surface stress                  = ', TRIM(sn_rcv_tau%cldes   ), ' (', TRIM(sn_rcv_tau%clcat   ), ')'
         WRITE(numout,*)'                     - referential    = ', sn_rcv_tau%clvref
         WRITE(numout,*)'                     - orientation    = ', sn_rcv_tau%clvor
         WRITE(numout,*)'                     - mesh           = ', sn_rcv_tau%clvgrd
         WRITE(numout,*)'      non-solar heat flux sensitivity = ', TRIM(sn_rcv_dqnsdt%cldes), ' (', TRIM(sn_rcv_dqnsdt%clcat), ')'
         WRITE(numout,*)'      solar heat flux                 = ', TRIM(sn_rcv_qsr%cldes   ), ' (', TRIM(sn_rcv_qsr%clcat   ), ')'
         WRITE(numout,*)'      non-solar heat flux             = ', TRIM(sn_rcv_qns%cldes   ), ' (', TRIM(sn_rcv_qns%clcat   ), ')'
         WRITE(numout,*)'      freshwater budget               = ', TRIM(sn_rcv_emp%cldes   ), ' (', TRIM(sn_rcv_emp%clcat   ), ')'
         WRITE(numout,*)'      runoffs                         = ', TRIM(sn_rcv_rnf%cldes   ), ' (', TRIM(sn_rcv_rnf%clcat   ), ')'
         WRITE(numout,*)'      calving                         = ', TRIM(sn_rcv_cal%cldes   ), ' (', TRIM(sn_rcv_cal%clcat   ), ')'
         WRITE(numout,*)'      iceberg                         = ', TRIM(sn_rcv_icb%cldes   ), ' (', TRIM(sn_rcv_icb%clcat   ), ')'
         WRITE(numout,*)'      ice shelf                       = ', TRIM(sn_rcv_isf%cldes   ), ' (', TRIM(sn_rcv_isf%clcat   ), ')'
         WRITE(numout,*)'      sea ice heat fluxes             = ', TRIM(sn_rcv_iceflx%cldes), ' (', TRIM(sn_rcv_iceflx%clcat), ')'
         WRITE(numout,*)'      atm co2                         = ', TRIM(sn_rcv_co2%cldes   ), ' (', TRIM(sn_rcv_co2%clcat   ), ')'
         WRITE(numout,*)'  sent fields (multiple ice categories)'
         WRITE(numout,*)'      surface temperature             = ', TRIM(sn_snd_temp%cldes  ), ' (', TRIM(sn_snd_temp%clcat  ), ')'
         WRITE(numout,*)'      albedo                          = ', TRIM(sn_snd_alb%cldes   ), ' (', TRIM(sn_snd_alb%clcat   ), ')'
         WRITE(numout,*)'      ice/snow thickness              = ', TRIM(sn_snd_thick%cldes ), ' (', TRIM(sn_snd_thick%clcat ), ')'
         WRITE(numout,*)'      surface current                 = ', TRIM(sn_snd_crt%cldes   ), ' (', TRIM(sn_snd_crt%clcat   ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crt%clvref 
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crt%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crt%clvgrd
         WRITE(numout,*)'      oce co2 flux                    = ', TRIM(sn_snd_co2%cldes   ), ' (', TRIM(sn_snd_co2%clcat   ), ')'
         WRITE(numout,*)'  nn_cplmodel                         = ', nn_cplmodel
         WRITE(numout,*)'  ln_usecplmask                       = ', ln_usecplmask
      ENDIF

      !                                   ! allocate sbccpl arrays
      IF( sbc_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_cpl_alloc : unable to allocate arrays' )
     
      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress 

      ! for each field: define the OASIS name                              (srcv(:)%clname)
      !                 define receive or not from the namelist parameters (srcv(:)%laction)
      !                 define the north fold type of lbc                  (srcv(:)%nsgn)

      ! default definitions of srcv
      srcv(:)%laction = .FALSE.   ;   srcv(:)%clgrid = 'T'   ;   srcv(:)%nsgn = 1.   ;   srcv(:)%nct = 1

      !                                                      ! ------------------------- !
      !                                                      ! ice and ocean wind stress !   
      !                                                      ! ------------------------- !
      !                                                           ! Name 
      srcv(jpr_otx1)%clname = 'O_OTaux1'      ! 1st ocean component on grid ONE (T or U)
      srcv(jpr_oty1)%clname = 'O_OTauy1'      ! 2nd   -      -         -     - 
      srcv(jpr_otz1)%clname = 'O_OTauz1'      ! 3rd   -      -         -     - 
      srcv(jpr_otx2)%clname = 'O_OTaux2'      ! 1st ocean component on grid TWO (V)
      srcv(jpr_oty2)%clname = 'O_OTauy2'      ! 2nd   -      -         -     - 
      srcv(jpr_otz2)%clname = 'O_OTauz2'      ! 3rd   -      -         -     - 
      !
      srcv(jpr_itx1)%clname = 'O_ITaux1'      ! 1st  ice  component on grid ONE (T, F, I or U)
      srcv(jpr_ity1)%clname = 'O_ITauy1'      ! 2nd   -      -         -     - 
      srcv(jpr_itz1)%clname = 'O_ITauz1'      ! 3rd   -      -         -     - 
      srcv(jpr_itx2)%clname = 'O_ITaux2'      ! 1st  ice  component on grid TWO (V)
      srcv(jpr_ity2)%clname = 'O_ITauy2'      ! 2nd   -      -         -     - 
      srcv(jpr_itz2)%clname = 'O_ITauz2'      ! 3rd   -      -         -     - 
      ! 
      ! Vectors: change of sign at north fold ONLY if on the local grid
      IF( TRIM( sn_rcv_tau%clvor ) == 'local grid' )   srcv(jpr_otx1:jpr_itz2)%nsgn = -1.
      
      !                                                           ! Set grid and action
      SELECT CASE( TRIM( sn_rcv_tau%clvgrd ) )      !  'T', 'U,V', 'U,V,I', 'U,V,F', 'T,I', 'T,F', or 'T,U,V'
      CASE( 'T' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'U,V' ) 
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_itz2)%laction = .TRUE.     ! receive oce and ice components on both grid 1 & 2
      CASE( 'U,V,T' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'T'        ! ice components given at T-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,I' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,F' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'T,I' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'T,F' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'T,U,V' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'T'        ! oce components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 only
         srcv(jpr_itx1:jpr_itz2)%laction = .TRUE.     ! receive ice components on grid 1 & 2
      CASE default   
         CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_tau%clvgrd' )
      END SELECT
      !
      IF( TRIM( sn_rcv_tau%clvref ) == 'spherical' )   &           ! spherical: 3rd component not received
         &     srcv( (/jpr_otz1, jpr_otz2, jpr_itz1, jpr_itz2/) )%laction = .FALSE. 
      !
      IF( TRIM( sn_rcv_tau%clvor  ) == 'local grid' ) THEN        ! already on local grid -> no need of the second grid
            srcv(jpr_otx2:jpr_otz2)%laction = .FALSE. 
            srcv(jpr_itx2:jpr_itz2)%laction = .FALSE. 
            srcv(jpr_oty1)%clgrid = srcv(jpr_oty2)%clgrid   ! not needed but cleaner...
            srcv(jpr_ity1)%clgrid = srcv(jpr_ity2)%clgrid   ! not needed but cleaner...
      ENDIF
      !
      IF( TRIM( sn_rcv_tau%cldes ) /= 'oce and ice' ) THEN        ! 'oce and ice' case ocean stress on ocean mesh used
         srcv(jpr_itx1:jpr_itz2)%laction = .FALSE.    ! ice components not received
         srcv(jpr_itx1)%clgrid = 'U'                  ! ocean stress used after its transformation
         srcv(jpr_ity1)%clgrid = 'V'                  ! i.e. it is always at U- & V-points for i- & j-comp. resp.
      ENDIF
       
      !                                                      ! ------------------------- !
      !                                                      !    freshwater budget      !   E-P
      !                                                      ! ------------------------- !
      ! we suppose that atmosphere modele do not make the difference between precipiration (liquide or solid)
      ! over ice of free ocean within the same atmospheric cell.cd 
      srcv(jpr_rain)%clname = 'OTotRain'      ! Rain = liquid precipitation
      srcv(jpr_snow)%clname = 'OTotSnow'      ! Snow = solid precipitation
      srcv(jpr_tevp)%clname = 'OTotEvap'      ! total evaporation (over oce + ice sublimation)
      srcv(jpr_ievp)%clname = 'OIceEvap'      ! evaporation over ice = sublimation
      srcv(jpr_sbpr)%clname = 'OSubMPre'      ! sublimation - liquid precipitation - solid precipitation 
      srcv(jpr_semp)%clname = 'OISubMSn'      ! ice solid water budget = sublimation - solid precipitation
      srcv(jpr_oemp)%clname = 'OOEvaMPr'      ! ocean water budget = ocean Evap - ocean precip
      SELECT CASE( TRIM( sn_rcv_emp%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(                                 jpr_oemp   )%laction = .TRUE. 
      CASE( 'conservative'  )
         srcv( (/jpr_rain, jpr_snow, jpr_ievp, jpr_tevp/) )%laction = .TRUE.
         IF ( k_ice <= 1 )  srcv(jpr_ievp)%laction = .FALSE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_ievp, jpr_sbpr, jpr_semp, jpr_oemp/) )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_emp%cldes' )
      END SELECT


      !                                                      ! ---------------------------------------------------- !
      !                                                      !     Runoffs, Calving, Iceberg, Iceshelf cavities     !   
      !                                                      ! ---------------------------------------------------- !
      srcv(jpr_rnf   )%clname = 'O_Runoff'
      IF( TRIM( sn_rcv_rnf%cldes ) == 'coupled' ) THEN
         srcv(jpr_rnf)%laction = .TRUE.
         l_rnfcpl              = .TRUE.                      ! -> no need to read runoffs in sbcrnf
         ln_rnf                = nn_components /= jp_iam_sas ! -> force to go through sbcrnf if not sas
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   runoffs received from oasis -> force ln_rnf = ', ln_rnf
      ENDIF
      !
      srcv(jpr_cal)%clname = 'OCalving'   ;  IF( TRIM( sn_rcv_cal%cldes) == 'coupled' )   srcv(jpr_cal)%laction = .TRUE.
      srcv(jpr_isf)%clname = 'OIcshelf'   ;  IF( TRIM( sn_rcv_isf%cldes) == 'coupled' )   srcv(jpr_isf)%laction = .TRUE.
      srcv(jpr_icb)%clname = 'OIceberg'   ;  IF( TRIM( sn_rcv_icb%cldes) == 'coupled' )   srcv(jpr_icb)%laction = .TRUE.

      IF( srcv(jpr_isf)%laction .AND. nn_isf > 0 ) THEN
         l_isfcpl             = .TRUE.                      ! -> no need to read isf in sbcisf
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   iceshelf received from oasis '
      ENDIF

      !                                                      ! ------------------------- !
      !                                                      !    non solar radiation    !   Qns
      !                                                      ! ------------------------- !
      srcv(jpr_qnsoce)%clname = 'O_QnsOce'
      srcv(jpr_qnsice)%clname = 'O_QnsIce'
      srcv(jpr_qnsmix)%clname = 'O_QnsMix'
      SELECT CASE( TRIM( sn_rcv_qns%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(               jpr_qnsoce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qnsice, jpr_qnsmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qnsice, jpr_qnsoce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qnsmix   )%laction = .TRUE. 
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qns%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' .AND. jpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qns%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !                                                      ! ------------------------- !
      !                                                      !    solar radiation        !   Qsr
      !                                                      ! ------------------------- !
      srcv(jpr_qsroce)%clname = 'O_QsrOce'
      srcv(jpr_qsrice)%clname = 'O_QsrIce'
      srcv(jpr_qsrmix)%clname = 'O_QsrMix'
      SELECT CASE( TRIM( sn_rcv_qsr%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(               jpr_qsroce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qsrice, jpr_qsrmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qsrice, jpr_qsroce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qsrmix   )%laction = .TRUE. 
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qsr%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qsr%cldes ) == 'mixed oce-ice' .AND. jpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qsr%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !                                                      ! ------------------------- !
      !                                                      !   non solar sensitivity   !   d(Qns)/d(T)
      !                                                      ! ------------------------- !
      srcv(jpr_dqnsdt)%clname = 'O_dQnsdT'   
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'coupled' )   srcv(jpr_dqnsdt)%laction = .TRUE.
      !
      ! non solar sensitivity mandatory for LIM ice model
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'none' .AND. k_ice /= 0 .AND. k_ice /= 4 .AND. nn_components /= jp_iam_sas ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_dqnsdt%cldes must be coupled in namsbc_cpl namelist' )
      ! non solar sensitivity mandatory for mixed oce-ice solar radiation coupling technique
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'none' .AND. TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' ) &
         CALL ctl_stop( 'sbc_cpl_init: namsbc_cpl namelist mismatch between sn_rcv_qns%cldes and sn_rcv_dqnsdt%cldes' )
      !                                                      ! ------------------------- !
      !                                                      !      10m wind module      !   
      !                                                      ! ------------------------- !
      srcv(jpr_w10m)%clname = 'O_Wind10'   ;   IF( TRIM(sn_rcv_w10m%cldes  ) == 'coupled' )   srcv(jpr_w10m)%laction = .TRUE. 
      !
      !                                                      ! ------------------------- !
      !                                                      !   wind stress module      !   
      !                                                      ! ------------------------- !
      srcv(jpr_taum)%clname = 'O_TauMod'   ;   IF( TRIM(sn_rcv_taumod%cldes) == 'coupled' )   srcv(jpr_taum)%laction = .TRUE.
      lhftau = srcv(jpr_taum)%laction

      !                                                      ! ------------------------- !
      !                                                      !      Atmospheric CO2      !
      !                                                      ! ------------------------- !
      srcv(jpr_co2 )%clname = 'O_AtmCO2'   ;   IF( TRIM(sn_rcv_co2%cldes   ) == 'coupled' )    srcv(jpr_co2 )%laction = .TRUE.
      !                                                      ! ------------------------- !
      !                                                      !   topmelt and botmelt     !   
      !                                                      ! ------------------------- !
      srcv(jpr_topm )%clname = 'OTopMlt'
      srcv(jpr_botm )%clname = 'OBotMlt'
      IF( TRIM(sn_rcv_iceflx%cldes) == 'coupled' ) THEN
         IF ( TRIM( sn_rcv_iceflx%clcat ) == 'yes' ) THEN
            srcv(jpr_topm:jpr_botm)%nct = jpl
         ELSE
            CALL ctl_stop( 'sbc_cpl_init: sn_rcv_iceflx%clcat should always be set to yes currently' )
         ENDIF
         srcv(jpr_topm:jpr_botm)%laction = .TRUE.
      ENDIF
      !                                                      ! ------------------------------- !
      !                                                      !   OPA-SAS coupling - rcv by opa !   
      !                                                      ! ------------------------------- !
      srcv(jpr_sflx)%clname = 'O_SFLX'
      srcv(jpr_fice)%clname = 'RIceFrc'
      !
      IF( nn_components == jp_iam_opa ) THEN    ! OPA coupled to SAS via OASIS: force received field by OPA (sent by SAS)
         srcv(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         srcv(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         srcv(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv( (/jpr_qsroce, jpr_qnsoce, jpr_oemp, jpr_sflx, jpr_fice, jpr_otx1, jpr_oty1, jpr_taum/) )%laction = .TRUE.
         srcv(jpr_otx1)%clgrid = 'U'        ! oce components given at U-point
         srcv(jpr_oty1)%clgrid = 'V'        !           and           V-point
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv( (/jpr_otx1,jpr_oty1/) )%nsgn = -1.
         sn_rcv_tau%clvgrd = 'U,V'
         sn_rcv_tau%clvor = 'local grid'
         sn_rcv_tau%clvref = 'spherical'
         sn_rcv_emp%cldes = 'oce only'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OPA coupling  '
            WRITE(numout,*)'               OPA component  '
            WRITE(numout,*)
            WRITE(numout,*)'  received fields from SAS component '
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V on local grid and sperical coordinates '
            WRITE(numout,*)'                  wind stress module'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! -------------------------------- !
      !                                                      !   OPA-SAS coupling - rcv by sas  !   
      !                                                      ! -------------------------------- !
      srcv(jpr_toce  )%clname = 'I_SSTSST'
      srcv(jpr_soce  )%clname = 'I_SSSal'
      srcv(jpr_ocx1  )%clname = 'I_OCurx1'
      srcv(jpr_ocy1  )%clname = 'I_OCury1'
      srcv(jpr_ssh   )%clname = 'I_SSHght'
      srcv(jpr_e3t1st)%clname = 'I_E3T1st'   
      srcv(jpr_fraqsr)%clname = 'I_FraQsr'   
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) srcv(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv( (/jpr_toce, jpr_soce, jpr_ssh, jpr_fraqsr, jpr_ocx1, jpr_ocy1/) )%laction = .TRUE.
         srcv( jpr_e3t1st )%laction = lk_vvl
         srcv(jpr_ocx1)%clgrid = 'U'        ! oce components given at U-point
         srcv(jpr_ocy1)%clgrid = 'V'        !           and           V-point
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv(jpr_ocx1:jpr_ocy1)%nsgn = -1.
         ! Change first letter to couple with atmosphere if already coupled OPA
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_Runoff received by OPA from SAS and therefore O_Runoff received by SAS from the Atmosphere
         DO jn = 1, jprcv
            IF ( srcv(jn)%clname(1:1) == "O" ) srcv(jn)%clname = "S"//srcv(jn)%clname(2:LEN(srcv(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OPA coupling  '
            WRITE(numout,*)'               SAS component  '
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  received fields from OPA component '
            ELSE
               WRITE(numout,*)'  Additional received fields from OPA component : '
            ENDIF
            WRITE(numout,*)'               sea surface temperature (Celcius) '
            WRITE(numout,*)'               sea surface salinity ' 
            WRITE(numout,*)'               surface currents ' 
            WRITE(numout,*)'               sea surface height ' 
            WRITE(numout,*)'               thickness of first ocean T level '        
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF
      
      ! =================================================== !
      ! Allocate all parts of frcv used for received fields !
      ! =================================================== !
      DO jn = 1, jprcv
         IF ( srcv(jn)%laction ) ALLOCATE( frcv(jn)%z3(jpi,jpj,srcv(jn)%nct) )
      END DO
      ! Allocate taum part of frcv which is used even when not received as coupling field
      IF ( .NOT. srcv(jpr_taum)%laction ) ALLOCATE( frcv(jpr_taum)%z3(jpi,jpj,srcv(jpr_taum)%nct) )
      ! Allocate w10m part of frcv which is used even when not received as coupling field
      IF ( .NOT. srcv(jpr_w10m)%laction ) ALLOCATE( frcv(jpr_w10m)%z3(jpi,jpj,srcv(jpr_w10m)%nct) )
      ! Allocate jpr_otx1 part of frcv which is used even when not received as coupling field
      IF ( .NOT. srcv(jpr_otx1)%laction ) ALLOCATE( frcv(jpr_otx1)%z3(jpi,jpj,srcv(jpr_otx1)%nct) )
      IF ( .NOT. srcv(jpr_oty1)%laction ) ALLOCATE( frcv(jpr_oty1)%z3(jpi,jpj,srcv(jpr_oty1)%nct) )
      ! Allocate itx1 and ity1 as they are used in sbc_cpl_ice_tau even if srcv(jpr_itx1)%laction = .FALSE.
      IF( k_ice /= 0 ) THEN
         IF ( .NOT. srcv(jpr_itx1)%laction ) ALLOCATE( frcv(jpr_itx1)%z3(jpi,jpj,srcv(jpr_itx1)%nct) )
         IF ( .NOT. srcv(jpr_ity1)%laction ) ALLOCATE( frcv(jpr_ity1)%z3(jpi,jpj,srcv(jpr_ity1)%nct) )
      END IF

      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (ssnd(:)%clname)
      !                 define send or not from the namelist parameters (ssnd(:)%laction)
      !                 define the north fold type of lbc               (ssnd(:)%nsgn)
      
      ! default definitions of nsnd
      ssnd(:)%laction = .FALSE.   ;   ssnd(:)%clgrid = 'T'   ;   ssnd(:)%nsgn = 1.  ; ssnd(:)%nct = 1
         
      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !
      !                                                      ! ------------------------- !
      ssnd(jps_toce)%clname = 'O_SSTSST'
      ssnd(jps_tice)%clname = 'O_TepIce'
      ssnd(jps_tmix)%clname = 'O_TepMix'
      SELECT CASE( TRIM( sn_snd_temp%cldes ) )
      CASE( 'none'                                 )       ! nothing to do
      CASE( 'oce only'                             )   ;   ssnd( jps_toce )%laction = .TRUE.
      CASE( 'oce and ice' , 'weighted oce and ice' )
         ssnd( (/jps_toce, jps_tice/) )%laction = .TRUE.
         IF ( TRIM( sn_snd_temp%clcat ) == 'yes' )  ssnd(jps_tice)%nct = jpl
      CASE( 'mixed oce-ice'                        )   ;   ssnd( jps_tmix )%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_temp%cldes' )
      END SELECT
           
      !                                                      ! ------------------------- !
      !                                                      !          Albedo           !
      !                                                      ! ------------------------- !
      ssnd(jps_albice)%clname = 'O_AlbIce' 
      ssnd(jps_albmix)%clname = 'O_AlbMix'
      SELECT CASE( TRIM( sn_snd_alb%cldes ) )
      CASE( 'none'                 )     ! nothing to do
      CASE( 'ice' , 'weighted ice' )   ; ssnd(jps_albice)%laction = .TRUE.
      CASE( 'mixed oce-ice'        )   ; ssnd(jps_albmix)%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_alb%cldes' )
      END SELECT
      !
      ! Need to calculate oceanic albedo if
      !     1. sending mixed oce-ice albedo or
      !     2. receiving mixed oce-ice solar radiation 
      IF ( TRIM ( sn_snd_alb%cldes ) == 'mixed oce-ice' .OR. TRIM ( sn_rcv_qsr%cldes ) == 'mixed oce-ice' ) THEN
         CALL albedo_oce( zaos, zacs )
         ! Due to lack of information on nebulosity : mean clear/overcast sky
         albedo_oce_mix(:,:) = ( zacs(:,:) + zaos(:,:) ) * 0.5
      ENDIF

      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness ! 
      !                                                      ! ------------------------- !
      ssnd(jps_fice)%clname = 'OIceFrc'
      ssnd(jps_hice)%clname = 'OIceTck'
      ssnd(jps_hsnw)%clname = 'OSnwTck'
      IF( k_ice /= 0 ) THEN
         ssnd(jps_fice)%laction = .TRUE.                  ! if ice treated in the ocean (even in climato case)
! Currently no namelist entry to determine sending of multi-category ice fraction so use the thickness entry for now
         IF ( TRIM( sn_snd_thick%clcat ) == 'yes' ) ssnd(jps_fice)%nct = jpl
      ENDIF
      
      SELECT CASE ( TRIM( sn_snd_thick%cldes ) )
      CASE( 'none'         )       ! nothing to do
      CASE( 'ice and snow' ) 
         ssnd(jps_hice:jps_hsnw)%laction = .TRUE.
         IF ( TRIM( sn_snd_thick%clcat ) == 'yes' ) THEN
            ssnd(jps_hice:jps_hsnw)%nct = jpl
         ENDIF
      CASE ( 'weighted ice and snow' ) 
         ssnd(jps_hice:jps_hsnw)%laction = .TRUE.
         IF ( TRIM( sn_snd_thick%clcat ) == 'yes' ) ssnd(jps_hice:jps_hsnw)%nct = jpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_thick%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !      Surface current      !
      !                                                      ! ------------------------- !
      !        ocean currents              !            ice velocities
      ssnd(jps_ocx1)%clname = 'O_OCurx1'   ;   ssnd(jps_ivx1)%clname = 'O_IVelx1'
      ssnd(jps_ocy1)%clname = 'O_OCury1'   ;   ssnd(jps_ivy1)%clname = 'O_IVely1'
      ssnd(jps_ocz1)%clname = 'O_OCurz1'   ;   ssnd(jps_ivz1)%clname = 'O_IVelz1'
      !
      ssnd(jps_ocx1:jps_ivz1)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( sn_snd_crt%clvgrd == 'U,V' ) THEN
         ssnd(jps_ocx1)%clgrid = 'U' ; ssnd(jps_ocy1)%clgrid = 'V'
      ELSE IF( sn_snd_crt%clvgrd /= 'T' ) THEN  
         CALL ctl_stop( 'sn_snd_crt%clvgrd must be equal to T' )
         ssnd(jps_ocx1:jps_ivz1)%clgrid  = 'T'      ! all oce and ice components on the same unique grid
      ENDIF
      ssnd(jps_ocx1:jps_ivz1)%laction = .TRUE.   ! default: all are send
      IF( TRIM( sn_snd_crt%clvref ) == 'spherical' )   ssnd( (/jps_ocz1, jps_ivz1/) )%laction = .FALSE. 
      IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) ssnd(jps_ocx1:jps_ivz1)%nsgn = 1.
      SELECT CASE( TRIM( sn_snd_crt%cldes ) )
      CASE( 'none'                 )   ;   ssnd(jps_ocx1:jps_ivz1)%laction = .FALSE.
      CASE( 'oce only'             )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE( 'weighted oce and ice' )   !   nothing to do
      CASE( 'mixed oce-ice'        )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_crt%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !          CO2 flux         !
      !                                                      ! ------------------------- !
      ssnd(jps_co2)%clname = 'O_CO2FLX' ;  IF( TRIM(sn_snd_co2%cldes) == 'coupled' )    ssnd(jps_co2 )%laction = .TRUE.

      !                                                      ! ------------------------------- !
      !                                                      !   OPA-SAS coupling - snd by opa !   
      !                                                      ! ------------------------------- !
      ssnd(jps_ssh   )%clname = 'O_SSHght' 
      ssnd(jps_soce  )%clname = 'O_SSSal' 
      ssnd(jps_e3t1st)%clname = 'O_E3T1st'   
      ssnd(jps_fraqsr)%clname = 'O_FraQsr'
      !
      IF( nn_components == jp_iam_opa ) THEN
         ssnd(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd( (/jps_toce, jps_soce, jps_ssh, jps_fraqsr, jps_ocx1, jps_ocy1/) )%laction = .TRUE.
         ssnd( jps_e3t1st )%laction = lk_vvl
         ! vector definition: not used but cleaner...
         ssnd(jps_ocx1)%clgrid  = 'U'        ! oce components given at U-point
         ssnd(jps_ocy1)%clgrid  = 'V'        !           and           V-point
         sn_snd_crt%clvgrd = 'U,V'
         sn_snd_crt%clvor = 'local grid'
         sn_snd_crt%clvref = 'spherical'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'  sent fields to SAS component '
            WRITE(numout,*)'               sea surface temperature (T before, Celcius) '
            WRITE(numout,*)'               sea surface salinity ' 
            WRITE(numout,*)'               surface currents U,V on local grid and spherical coordinates' 
            WRITE(numout,*)'               sea surface height ' 
            WRITE(numout,*)'               thickness of first ocean T level '        
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! ------------------------------- !
      !                                                      !   OPA-SAS coupling - snd by sas !   
      !                                                      ! ------------------------------- !
      ssnd(jps_sflx  )%clname = 'I_SFLX'     
      ssnd(jps_fice2 )%clname = 'IIceFrc'
      ssnd(jps_qsroce)%clname = 'I_QsrOce'   
      ssnd(jps_qnsoce)%clname = 'I_QnsOce'   
      ssnd(jps_oemp  )%clname = 'IOEvaMPr' 
      ssnd(jps_otx1  )%clname = 'I_OTaux1'   
      ssnd(jps_oty1  )%clname = 'I_OTauy1'   
      ssnd(jps_rnf   )%clname = 'I_Runoff'   
      ssnd(jps_taum  )%clname = 'I_TauMod'   
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) ssnd(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd( (/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/) )%laction = .TRUE.
         !
         ! Change first letter to couple with atmosphere if already coupled with sea_ice
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_SSTSST sent by OPA to SAS and therefore S_SSTSST sent by SAS to the Atmosphere
         DO jn = 1, jpsnd
            IF ( ssnd(jn)%clname(1:1) == "O" ) ssnd(jn)%clname = "S"//ssnd(jn)%clname(2:LEN(ssnd(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  sent fields to OPA component '
            ELSE
               WRITE(numout,*)'  Additional sent fields to OPA component : '
            ENDIF
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V components'
            WRITE(numout,*)'                  wind stress module'
         ENDIF
      ENDIF

      !
      ! ================================ !
      !   initialisation of the coupler  !
      ! ================================ !

      CALL cpl_define(jprcv, jpsnd, nn_cplmodel)
      
      IF (ln_usecplmask) THEN 
         xcplmask(:,:,:) = 0.
         CALL iom_open( 'cplmask', inum )
         CALL iom_get( inum, jpdom_unknown, 'cplmask', xcplmask(1:nlci,1:nlcj,1:nn_cplmodel),   &
            &          kstart = (/ mig(1),mjg(1),1 /), kcount = (/ nlci,nlcj,nn_cplmodel /) )
         CALL iom_close( inum )
      ELSE
         xcplmask(:,:,:) = 1.
      ENDIF
      xcplmask(:,:,0) = 1. - SUM( xcplmask(:,:,1:nn_cplmodel), dim = 3 )
      !
      ncpl_qsr_freq = cpl_freq( 'O_QsrOce' ) + cpl_freq( 'O_QsrMix' ) + cpl_freq( 'I_QsrOce' ) + cpl_freq( 'I_QsrMix' )
      IF( ln_dm2dc .AND. ln_cpl .AND. ncpl_qsr_freq /= 86400 )   &
         &   CALL ctl_stop( 'sbc_cpl_init: diurnal cycle reconstruction (ln_dm2dc) needs daily couping for solar radiation' )
      ncpl_qsr_freq = 86400 / ncpl_qsr_freq

      CALL wrk_dealloc( jpi,jpj, zacs, zaos )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_cpl_init')
      !
   END SUBROUTINE sbc_cpl_init


   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in frcv array). called at each time step.
      !!                OASIS controls if there is something do receive or not. nrcvinfo contains the info
      !!                to know if the field was really received or not
      !!
      !!              --> If ocean stress was really received:
      !!
      !!                  - transform the received ocean stress vector from the received
      !!                 referential and grid into an atmosphere-ocean stress in 
      !!                 the (i,j) ocean referencial and at the ocean velocity point. 
      !!                    The received stress are :
      !!                     - defined by 3 components (if cartesian coordinate)
      !!                            or by 2 components (if spherical)
      !!                     - oriented along geographical   coordinate (if eastward-northward)
      !!                            or  along the local grid coordinate (if local grid)
      !!                     - given at U- and V-point, resp.   if received on 2 grids
      !!                            or at T-point               if received on 1 grid
      !!                    Therefore and if necessary, they are successively 
      !!                  processed in order to obtain them 
      !!                     first  as  2 components on the sphere 
      !!                     second as  2 components oriented along the local grid
      !!                     third  as  2 components on the U,V grid 
      !!
      !!              --> 
      !!
      !!              - In 'ocean only' case, non solar and solar ocean heat fluxes 
      !!             and total ocean freshwater fluxes  
      !!
      !! ** Method  :   receive all fields from the atmosphere and transform 
      !!              them into ocean surface boundary condition fields 
      !!
      !! ** Action  :   update  utau, vtau   ocean stress at U,V grid 
      !!                        taum         wind stress module at T-point
      !!                        wndm         wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!                        qns          non solar heat fluxes including emp heat content    (ocean only case)
      !!                                     and the latent heat flux of solid precip. melting
      !!                        qsr          solar ocean heat fluxes   (ocean only case)
      !!                        emp          upward mass flux [evap. - precip. (- runoffs) (- calving)] (ocean only case)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           ::   kt          ! ocean model time step index
      INTEGER, INTENT(in)           ::   k_fsbc      ! frequency of sbc (-> ice model) computation 
      INTEGER, INTENT(in)           ::   k_ice       ! ice management in the sbc (=0/1/2/3)

      !!
      LOGICAL  ::   llnewtx, llnewtau      ! update wind stress components and module??
      INTEGER  ::   ji, jj, jn             ! dummy loop indices
      INTEGER  ::   isec                   ! number of seconds since nit000 (assuming rdttra did not change since nit000)
      REAL(wp) ::   zcumulneg, zcumulpos   ! temporary scalars     
      REAL(wp) ::   zcoef                  ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3        ! drag coefficient
      REAL(wp) ::   zzx, zzy               ! temporary variables
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztx, zty, zmsk, zemp, zqns, zqsr
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_cpl_rcv')
      !
      CALL wrk_alloc( jpi,jpj, ztx, zty, zmsk, zemp, zqns, zqsr )
      !
      IF( ln_mixcpl )   zmsk(:,:) = 1. - xcplmask(:,:,0)
      !
      !                                                      ! ======================================================= !
      !                                                      ! Receive all the atmos. fields (including ice information)
      !                                                      ! ======================================================= !
      isec = ( kt - nit000 ) * NINT( rdttra(1) )                ! date of exchanges
      DO jn = 1, jprcv                                          ! received fields sent by the atmosphere
         IF( srcv(jn)%laction )   CALL cpl_rcv( jn, isec, frcv(jn)%z3, xcplmask(:,:,1:nn_cplmodel), nrcvinfo(jn) )
      END DO

      !                                                      ! ========================= !
      IF( srcv(jpr_otx1)%laction ) THEN                      !  ocean stress components  !
         !                                                   ! ========================= !
         ! define frcv(jpr_otx1)%z3(:,:,1) and frcv(jpr_oty1)%z3(:,:,1): stress at U/V point along model grid
         ! => need to be done only when we receive the field
         IF(  nrcvinfo(jpr_otx1) == OASIS_Rcv ) THEN
            !
            IF( TRIM( sn_rcv_tau%clvref ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               !
               CALL geo2oce( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), frcv(jpr_otz1)%z3(:,:,1),   &
                  &          srcv(jpr_otx1)%clgrid, ztx, zty )
               frcv(jpr_otx1)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 1st grid
               frcv(jpr_oty1)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 1st grid
               !
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL geo2oce( frcv(jpr_otx2)%z3(:,:,1), frcv(jpr_oty2)%z3(:,:,1), frcv(jpr_otz2)%z3(:,:,1),   &
                     &          srcv(jpr_otx2)%clgrid, ztx, zty )
                  frcv(jpr_otx2)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 2nd grid
                  frcv(jpr_oty2)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 2nd grid
               ENDIF
               !
            ENDIF
            !
            IF( TRIM( sn_rcv_tau%clvor ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), srcv(jpr_otx1)%clgrid, 'en->i', ztx )   
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL rot_rep( frcv(jpr_otx2)%z3(:,:,1), frcv(jpr_oty2)%z3(:,:,1), srcv(jpr_otx2)%clgrid, 'en->j', zty )   
               ELSE	
                  CALL rot_rep( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), srcv(jpr_otx1)%clgrid, 'en->j', zty )  
               ENDIF
               frcv(jpr_otx1)%z3(:,:,1) = ztx(:,:)      ! overwrite 1st component on the 1st grid
               frcv(jpr_oty1)%z3(:,:,1) = zty(:,:)      ! overwrite 2nd component on the 2nd grid
            ENDIF
            !                              
            IF( srcv(jpr_otx1)%clgrid == 'T' ) THEN
               DO jj = 2, jpjm1                                          ! T ==> (U,V)
                  DO ji = 2, jpim1   ! vector opt.
                     frcv(jpr_otx1)%z3(ji,jj,1) = 0.5 * ( frcv(jpr_otx1)%z3(ji+1,jj  ,1) + frcv(jpr_otx1)%z3(ji,jj,1) )
                     frcv(jpr_oty1)%z3(ji,jj,1) = 0.5 * ( frcv(jpr_oty1)%z3(ji  ,jj+1,1) + frcv(jpr_oty1)%z3(ji,jj,1) )
                  END DO
               END DO
               CALL lbc_lnk( frcv(jpr_otx1)%z3(:,:,1), 'U',  -1. )   ;   CALL lbc_lnk( frcv(jpr_oty1)%z3(:,:,1), 'V',  -1. )
            ENDIF
            llnewtx = .TRUE.
         ELSE
            llnewtx = .FALSE.
         ENDIF
         !                                                   ! ========================= !
      ELSE                                                   !   No dynamical coupling   !
         !                                                   ! ========================= !
         frcv(jpr_otx1)%z3(:,:,1) = 0.e0                               ! here simply set to zero 
         frcv(jpr_oty1)%z3(:,:,1) = 0.e0                               ! an external read in a file can be added instead
         llnewtx = .TRUE.
         !
      ENDIF
      !                                                      ! ========================= !
      !                                                      !    wind stress module     !   (taum)
      !                                                      ! ========================= !
      !
      IF( .NOT. srcv(jpr_taum)%laction ) THEN                    ! compute wind stress module from its components if not received 
         ! => need to be done only when otx1 was changed
         IF( llnewtx ) THEN
!CDIR NOVERRCHK
            DO jj = 2, jpjm1
!CDIR NOVERRCHK
               DO ji = 2, jpim1   ! vect. opt.
                  zzx = frcv(jpr_otx1)%z3(ji-1,jj  ,1) + frcv(jpr_otx1)%z3(ji,jj,1)
                  zzy = frcv(jpr_oty1)%z3(ji  ,jj-1,1) + frcv(jpr_oty1)%z3(ji,jj,1)
                  frcv(jpr_taum)%z3(ji,jj,1) = 0.5 * SQRT( zzx * zzx + zzy * zzy )
               END DO
            END DO
            CALL lbc_lnk( frcv(jpr_taum)%z3(:,:,1), 'T', 1. )
            llnewtau = .TRUE.
         ELSE
            llnewtau = .FALSE.
         ENDIF
      ELSE
         llnewtau = nrcvinfo(jpr_taum) == OASIS_Rcv
         ! Stress module can be negative when received (interpolation problem)
         IF( llnewtau ) THEN 
            frcv(jpr_taum)%z3(:,:,1) = MAX( 0._wp, frcv(jpr_taum)%z3(:,:,1) )
         ENDIF
      ENDIF
      !
      !                                                      ! ========================= !
      !                                                      !      10 m wind speed      !   (wndm)
      !                                                      ! ========================= !
      !
      IF( .NOT. srcv(jpr_w10m)%laction ) THEN                    ! compute wind spreed from wind stress module if not received  
         ! => need to be done only when taumod was changed
         IF( llnewtau ) THEN 
            zcoef = 1. / ( zrhoa * zcdrag ) 
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi 
                  frcv(jpr_w10m)%z3(ji,jj,1) = SQRT( frcv(jpr_taum)%z3(ji,jj,1) * zcoef )
               END DO
            END DO
         ENDIF
      ENDIF

      ! u(v)tau and taum will be modified by ice model
      ! -> need to be reset before each call of the ice/fsbc      
      IF( MOD( kt-1, k_fsbc ) == 0 ) THEN
         !
         IF( ln_mixcpl ) THEN
            utau(:,:) = utau(:,:) * xcplmask(:,:,0) + frcv(jpr_otx1)%z3(:,:,1) * zmsk(:,:)
            vtau(:,:) = vtau(:,:) * xcplmask(:,:,0) + frcv(jpr_oty1)%z3(:,:,1) * zmsk(:,:)
            taum(:,:) = taum(:,:) * xcplmask(:,:,0) + frcv(jpr_taum)%z3(:,:,1) * zmsk(:,:)
            wndm(:,:) = wndm(:,:) * xcplmask(:,:,0) + frcv(jpr_w10m)%z3(:,:,1) * zmsk(:,:)
         ELSE
            utau(:,:) = frcv(jpr_otx1)%z3(:,:,1)
            vtau(:,:) = frcv(jpr_oty1)%z3(:,:,1)
            taum(:,:) = frcv(jpr_taum)%z3(:,:,1)
            wndm(:,:) = frcv(jpr_w10m)%z3(:,:,1)
         ENDIF
         CALL iom_put( "taum_oce", taum )   ! output wind stress module
         !  
      ENDIF


      !  Fields received by SAS when OASIS coupling
      !  (arrays no more filled at sbcssm stage)
      !                                                      ! ================== !
      !                                                      !        SSS         !
      !                                                      ! ================== !
      IF( srcv(jpr_soce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         sss_m(:,:) = frcv(jpr_soce)%z3(:,:,1)
         CALL iom_put( 'sss_m', sss_m )
      ENDIF
      !                                               
      !                                                      ! ================== !
      !                                                      !        SST         !
      !                                                      ! ================== !
      IF( srcv(jpr_toce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         sst_m(:,:) = frcv(jpr_toce)%z3(:,:,1)
         IF( srcv(jpr_soce)%laction .AND. ln_useCT ) THEN    ! make sure that sst_m is the potential temperature
            sst_m(:,:) = eos_pt_from_ct( sst_m(:,:), sss_m(:,:) )
         ENDIF
      ENDIF
      !                                                      ! ================== !
      !                                                      !        SSH         !
      !                                                      ! ================== !
      IF( srcv(jpr_ssh )%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssh_m(:,:) = frcv(jpr_ssh )%z3(:,:,1)
         CALL iom_put( 'ssh_m', ssh_m )
      ENDIF
      !                                                      ! ================== !
      !                                                      !  surface currents  !
      !                                                      ! ================== !
      IF( srcv(jpr_ocx1)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssu_m(:,:) = frcv(jpr_ocx1)%z3(:,:,1)
         ub (:,:,1) = ssu_m(:,:)                             ! will be used in sbcice_lim in the call of lim_sbc_tau
         un (:,:,1) = ssu_m(:,:)                             ! will be used in sbc_cpl_snd if atmosphere coupling
         CALL iom_put( 'ssu_m', ssu_m )
      ENDIF
      IF( srcv(jpr_ocy1)%laction ) THEN
         ssv_m(:,:) = frcv(jpr_ocy1)%z3(:,:,1)
         vb (:,:,1) = ssv_m(:,:)                             ! will be used in sbcice_lim in the call of lim_sbc_tau
         vn (:,:,1) = ssv_m(:,:)                             ! will be used in sbc_cpl_snd if atmosphere coupling
         CALL iom_put( 'ssv_m', ssv_m )
      ENDIF
      !                                                      ! ======================== !
      !                                                      !  first T level thickness !
      !                                                      ! ======================== !
      IF( srcv(jpr_e3t1st )%laction ) THEN                   ! received by sas in case of opa <-> sas coupling
         e3t_m(:,:) = frcv(jpr_e3t1st )%z3(:,:,1)
         CALL iom_put( 'e3t_m', e3t_m(:,:) )
      ENDIF
      !                                                      ! ================================ !
      !                                                      !  fraction of solar net radiation !
      !                                                      ! ================================ !
      IF( srcv(jpr_fraqsr)%laction ) THEN                    ! received by sas in case of opa <-> sas coupling
         frq_m(:,:) = frcv(jpr_fraqsr)%z3(:,:,1)
         CALL iom_put( 'frq_m', frq_m )
      ENDIF
      
      !                                                      ! ========================= !
      IF( k_ice <= 1 .AND. MOD( kt-1, k_fsbc ) == 0 ) THEN   !  heat & freshwater fluxes ! (Ocean only case)
         !                                                   ! ========================= !
         !
         !                                                       ! total freshwater fluxes over the ocean (emp)
         IF( srcv(jpr_oemp)%laction .OR. srcv(jpr_rain)%laction ) THEN
            SELECT CASE( TRIM( sn_rcv_emp%cldes ) )                                    ! evaporation - precipitation
            CASE( 'conservative' )
               zemp(:,:) = frcv(jpr_tevp)%z3(:,:,1) - ( frcv(jpr_rain)%z3(:,:,1) + frcv(jpr_snow)%z3(:,:,1) )
            CASE( 'oce only', 'oce and ice' )
               zemp(:,:) = frcv(jpr_oemp)%z3(:,:,1)
            CASE default
               CALL ctl_stop( 'sbc_cpl_rcv: wrong definition of sn_rcv_emp%cldes' )
            END SELECT
         ELSE
            zemp(:,:) = 0._wp
         ENDIF
         !
         !   
         !                                                        ! runoffs and calving (added in emp)
         IF( srcv(jpr_rnf)%laction )     rnf(:,:)  = frcv(jpr_rnf)%z3(:,:,1)
         IF( srcv(jpr_cal)%laction )     zemp(:,:) = zemp(:,:) - frcv(jpr_cal)%z3(:,:,1)

         IF( srcv(jpr_icb)%laction )  THEN 
             fwficb(:,:) = frcv(jpr_icb)%z3(:,:,1)
             rnf(:,:)    = rnf(:,:) + fwficb(:,:)   ! iceberg added to runfofs
         ENDIF
         IF( srcv(jpr_isf)%laction )  fwfisf(:,:) = - frcv(jpr_isf)%z3(:,:,1)  ! fresh water flux from the isf (fwfisf <0 mean melting)  
         
         IF( ln_mixcpl ) THEN   ;   emp(:,:) = emp(:,:) * xcplmask(:,:,0) + zemp(:,:) * zmsk(:,:)
         ELSE                   ;   emp(:,:) =                              zemp(:,:)
         ENDIF
         !
         !                                                       ! non solar heat flux over the ocean (qns)
         IF(      srcv(jpr_qnsoce)%laction ) THEN   ;   zqns(:,:) = frcv(jpr_qnsoce)%z3(:,:,1)
         ELSE IF( srcv(jpr_qnsmix)%laction ) THEN   ;   zqns(:,:) = frcv(jpr_qnsmix)%z3(:,:,1)
         ELSE                                       ;   zqns(:,:) = 0._wp
         END IF
         ! update qns over the free ocean with:
         IF( nn_components /= jp_iam_opa ) THEN
            zqns(:,:) =  zqns(:,:) - zemp(:,:) * sst_m(:,:) * rcp         ! remove heat content due to mass flux (assumed to be at SST)
            IF( srcv(jpr_snow  )%laction ) THEN
               zqns(:,:) = zqns(:,:) - frcv(jpr_snow)%z3(:,:,1) * lfus    ! energy for melting solid precipitation over the free ocean
            ENDIF
         ENDIF
         !
         IF( srcv(jpr_icb)%laction )  zqns(:,:) = zqns(:,:) - frcv(jpr_icb)%z3(:,:,1) * lfus ! remove heat content associated to iceberg melting
         !
         IF( ln_mixcpl ) THEN   ;   qns(:,:) = qns(:,:) * xcplmask(:,:,0) + zqns(:,:) * zmsk(:,:)
         ELSE                   ;   qns(:,:) =                              zqns(:,:)
         ENDIF

         !                                                       ! solar flux over the ocean          (qsr)
         IF     ( srcv(jpr_qsroce)%laction ) THEN   ;   zqsr(:,:) = frcv(jpr_qsroce)%z3(:,:,1)
         ELSE IF( srcv(jpr_qsrmix)%laction ) then   ;   zqsr(:,:) = frcv(jpr_qsrmix)%z3(:,:,1)
         ELSE                                       ;   zqsr(:,:) = 0._wp
         ENDIF
         IF( ln_dm2dc .AND. ln_cpl )   zqsr(:,:) = sbc_dcy( zqsr )   ! modify qsr to include the diurnal cycle
         IF( ln_mixcpl ) THEN   ;   qsr(:,:) = qsr(:,:) * xcplmask(:,:,0) + zqsr(:,:) * zmsk(:,:)
         ELSE                   ;   qsr(:,:) =                              zqsr(:,:)
         ENDIF
         !
         ! salt flux over the ocean (received by opa in case of opa <-> sas coupling)
         IF( srcv(jpr_sflx )%laction )   sfx(:,:) = frcv(jpr_sflx  )%z3(:,:,1)
         ! Ice cover  (received by opa in case of opa <-> sas coupling)
         IF( srcv(jpr_fice )%laction )   fr_i(:,:) = frcv(jpr_fice )%z3(:,:,1)
         !

      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, ztx, zty, zmsk, zemp, zqns, zqsr )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_cpl_rcv')
      !
   END SUBROUTINE sbc_cpl_rcv
   

   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_tau  ***
      !!
      !! ** Purpose :   provide the stress over sea-ice in coupled mode 
      !!
      !! ** Method  :   transform the received stress from the atmosphere into
      !!             an atmosphere-ice stress in the (i,j) ocean referencial
      !!             and at the velocity point of the sea-ice model (cp_ice_msh):
      !!                'C'-grid : i- (j-) components given at U- (V-) point 
      !!                'I'-grid : B-grid lower-left corner: both components given at I-point 
      !!
      !!                The received stress are :
      !!                 - defined by 3 components (if cartesian coordinate)
      !!                        or by 2 components (if spherical)
      !!                 - oriented along geographical   coordinate (if eastward-northward)
      !!                        or  along the local grid coordinate (if local grid)
      !!                 - given at U- and V-point, resp.   if received on 2 grids
      !!                        or at a same point (T or I) if received on 1 grid
      !!                Therefore and if necessary, they are successively 
      !!             processed in order to obtain them 
      !!                 first  as  2 components on the sphere 
      !!                 second as  2 components oriented along the local grid
      !!                 third  as  2 components on the cp_ice_msh point 
      !!
      !!                Except in 'oce and ice' case, only one vector stress field 
      !!             is received. It has already been processed in sbc_cpl_rcv
      !!             so that it is now defined as (i,j) components given at U-
      !!             and V-points, respectively. Therefore, only the third
      !!             transformation is done and only if the ice-grid is a 'I'-grid. 
      !!
      !! ** Action  :   return ptau_i, ptau_j, the stress over the ice at cp_ice_msh point
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_tauj   ! at I-point (B-grid) or U & V-point (C-grid)
      !!
      INTEGER ::   ji, jj                          ! dummy loop indices
      INTEGER ::   itx                             ! index of taux over ice
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztx, zty 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_cpl_ice_tau')
      !
      CALL wrk_alloc( jpi,jpj, ztx, zty )

      IF( srcv(jpr_itx1)%laction ) THEN   ;   itx =  jpr_itx1   
      ELSE                                ;   itx =  jpr_otx1
      ENDIF

      ! do something only if we just received the stress from atmosphere
      IF(  nrcvinfo(itx) == OASIS_Rcv ) THEN

         !                                                      ! ======================= !
         IF( srcv(jpr_itx1)%laction ) THEN                      !   ice stress received   !
            !                                                   ! ======================= !
            !  
            IF( TRIM( sn_rcv_tau%clvref ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               CALL geo2oce(  frcv(jpr_itx1)%z3(:,:,1), frcv(jpr_ity1)%z3(:,:,1), frcv(jpr_itz1)%z3(:,:,1),   &
                  &          srcv(jpr_itx1)%clgrid, ztx, zty )
               frcv(jpr_itx1)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 1st grid
               frcv(jpr_ity1)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 1st grid
               !
               IF( srcv(jpr_itx2)%laction ) THEN
                  CALL geo2oce( frcv(jpr_itx2)%z3(:,:,1), frcv(jpr_ity2)%z3(:,:,1), frcv(jpr_itz2)%z3(:,:,1),   &
                     &          srcv(jpr_itx2)%clgrid, ztx, zty )
                  frcv(jpr_itx2)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 2nd grid
                  frcv(jpr_ity2)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 2nd grid
               ENDIF
               !
            ENDIF
            !
            IF( TRIM( sn_rcv_tau%clvor ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(jpr_itx1)%z3(:,:,1), frcv(jpr_ity1)%z3(:,:,1), srcv(jpr_itx1)%clgrid, 'en->i', ztx )   
               IF( srcv(jpr_itx2)%laction ) THEN
                  CALL rot_rep( frcv(jpr_itx2)%z3(:,:,1), frcv(jpr_ity2)%z3(:,:,1), srcv(jpr_itx2)%clgrid, 'en->j', zty )   
               ELSE
                  CALL rot_rep( frcv(jpr_itx1)%z3(:,:,1), frcv(jpr_ity1)%z3(:,:,1), srcv(jpr_itx1)%clgrid, 'en->j', zty )  
               ENDIF
               frcv(jpr_itx1)%z3(:,:,1) = ztx(:,:)      ! overwrite 1st component on the 1st grid
               frcv(jpr_ity1)%z3(:,:,1) = zty(:,:)      ! overwrite 2nd component on the 1st grid
            ENDIF
            !                                                   ! ======================= !
         ELSE                                                   !     use ocean stress    !
            !                                                   ! ======================= !
            frcv(jpr_itx1)%z3(:,:,1) = frcv(jpr_otx1)%z3(:,:,1)
            frcv(jpr_ity1)%z3(:,:,1) = frcv(jpr_oty1)%z3(:,:,1)
            !
         ENDIF
         !                                                      ! ======================= !
         !                                                      !     put on ice grid     !
         !                                                      ! ======================= !
         !    
         !                                                  j+1   j     -----V---F
         ! ice stress on ice velocity point (cp_ice_msh)                 !       |
         ! (C-grid ==>(U,V) or B-grid ==> I or F)                 j      |   T   U
         !                                                               |       |
         !                                                   j    j-1   -I-------|
         !                                               (for I)         |       |
         !                                                              i-1  i   i
         !                                                               i      i+1 (for I)
         SELECT CASE ( cp_ice_msh )
            !
         CASE( 'I' )                                         ! B-grid ==> I
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               DO jj = 2, jpjm1                                   ! (U,V) ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(jpr_itx1)%z3(ji-1,jj  ,1) + frcv(jpr_itx1)%z3(ji-1,jj-1,1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(jpr_ity1)%z3(ji  ,jj-1,1) + frcv(jpr_ity1)%z3(ji-1,jj-1,1) )
                  END DO
               END DO
            CASE( 'F' )
               DO jj = 2, jpjm1                                   ! F ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = frcv(jpr_itx1)%z3(ji-1,jj-1,1)
                     p_tauj(ji,jj) = frcv(jpr_ity1)%z3(ji-1,jj-1,1)
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.25 * ( frcv(jpr_itx1)%z3(ji,jj  ,1) + frcv(jpr_itx1)%z3(ji-1,jj  ,1)   &
                        &                   + frcv(jpr_itx1)%z3(ji,jj-1,1) + frcv(jpr_itx1)%z3(ji-1,jj-1,1) ) 
                     p_tauj(ji,jj) = 0.25 * ( frcv(jpr_ity1)%z3(ji,jj  ,1) + frcv(jpr_ity1)%z3(ji-1,jj  ,1)   &
                        &                   + frcv(jpr_oty1)%z3(ji,jj-1,1) + frcv(jpr_ity1)%z3(ji-1,jj-1,1) )
                  END DO
               END DO
            CASE( 'I' )
               p_taui(:,:) = frcv(jpr_itx1)%z3(:,:,1)                   ! I ==> I
               p_tauj(:,:) = frcv(jpr_ity1)%z3(:,:,1)
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'I' ) THEN 
               CALL lbc_lnk( p_taui, 'I',  -1. )   ;   CALL lbc_lnk( p_tauj, 'I',  -1. )
            ENDIF
            !
         CASE( 'F' )                                         ! B-grid ==> F
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               DO jj = 2, jpjm1                                   ! (U,V) ==> F
                  DO ji = 2, jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(jpr_itx1)%z3(ji,jj,1) + frcv(jpr_itx1)%z3(ji  ,jj+1,1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(jpr_ity1)%z3(ji,jj,1) + frcv(jpr_ity1)%z3(ji+1,jj  ,1) )
                  END DO
               END DO
            CASE( 'I' )
               DO jj = 2, jpjm1                                   ! I ==> F
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = frcv(jpr_itx1)%z3(ji+1,jj+1,1)
                     p_tauj(ji,jj) = frcv(jpr_ity1)%z3(ji+1,jj+1,1)
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> F
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.25 * ( frcv(jpr_itx1)%z3(ji,jj  ,1) + frcv(jpr_itx1)%z3(ji+1,jj  ,1)   &
                        &                   + frcv(jpr_itx1)%z3(ji,jj+1,1) + frcv(jpr_itx1)%z3(ji+1,jj+1,1) ) 
                     p_tauj(ji,jj) = 0.25 * ( frcv(jpr_ity1)%z3(ji,jj  ,1) + frcv(jpr_ity1)%z3(ji+1,jj  ,1)   &
                        &                   + frcv(jpr_ity1)%z3(ji,jj+1,1) + frcv(jpr_ity1)%z3(ji+1,jj+1,1) )
                  END DO
               END DO
            CASE( 'F' )
               p_taui(:,:) = frcv(jpr_itx1)%z3(:,:,1)                   ! F ==> F
               p_tauj(:,:) = frcv(jpr_ity1)%z3(:,:,1)
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'F' ) THEN 
               CALL lbc_lnk( p_taui, 'F',  -1. )   ;   CALL lbc_lnk( p_tauj, 'F',  -1. )
            ENDIF
            !
         CASE( 'C' )                                         ! C-grid ==> U,V
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               p_taui(:,:) = frcv(jpr_itx1)%z3(:,:,1)                   ! (U,V) ==> (U,V)
               p_tauj(:,:) = frcv(jpr_ity1)%z3(:,:,1)
            CASE( 'F' )
               DO jj = 2, jpjm1                                   ! F ==> (U,V)
                  DO ji = 2, jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(jpr_itx1)%z3(ji,jj,1) + frcv(jpr_itx1)%z3(ji  ,jj-1,1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(jpr_ity1)%z3(jj,jj,1) + frcv(jpr_ity1)%z3(ji-1,jj  ,1) )
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> (U,V)
                  DO ji = 2, jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(jpr_itx1)%z3(ji+1,jj  ,1) + frcv(jpr_itx1)%z3(ji,jj,1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(jpr_ity1)%z3(ji  ,jj+1,1) + frcv(jpr_ity1)%z3(ji,jj,1) )
                  END DO
               END DO
            CASE( 'I' )
               DO jj = 2, jpjm1                                   ! I ==> (U,V)
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(jpr_itx1)%z3(ji+1,jj+1,1) + frcv(jpr_itx1)%z3(ji+1,jj  ,1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(jpr_ity1)%z3(ji+1,jj+1,1) + frcv(jpr_ity1)%z3(ji  ,jj+1,1) )
                  END DO
               END DO
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'U' ) THEN 
               CALL lbc_lnk( p_taui, 'U',  -1. )   ;   CALL lbc_lnk( p_tauj, 'V',  -1. )
            ENDIF
         END SELECT

      ENDIF
      !   
      CALL wrk_dealloc( jpi,jpj, ztx, zty )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_cpl_ice_tau')
      !
   END SUBROUTINE sbc_cpl_ice_tau
   

   SUBROUTINE sbc_cpl_ice_flx( p_frld, palbi, psst, pist )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_flx  ***
      !!
      !! ** Purpose :   provide the heat and freshwater fluxes of the ocean-ice system
      !!
      !! ** Method  :   transform the fields received from the atmosphere into
      !!             surface heat and fresh water boundary condition for the 
      !!             ice-ocean system. The following fields are provided:
      !!               * total non solar, solar and freshwater fluxes (qns_tot, 
      !!             qsr_tot and emp_tot) (total means weighted ice-ocean flux)
      !!             NB: emp_tot include runoffs and calving.
      !!               * fluxes over ice (qns_ice, qsr_ice, emp_ice) where
      !!             emp_ice = sublimation - solid precipitation as liquid
      !!             precipitation are re-routed directly to the ocean and 
      !!             calving directly enter the ocean (runoffs are read but included in trasbc.F90)
      !!               * solid precipitation (sprecip), used to add to qns_tot 
      !!             the heat lost associated to melting solid precipitation
      !!             over the ocean fraction.
      !!               * heat content of rain, snow and evap can also be provided,
      !!             otherwise heat flux associated with these mass flux are
      !!             guessed (qemp_oce, qemp_ice)
      !!
      !!             - the fluxes have been separated from the stress as
      !!               (a) they are updated at each ice time step compare to
      !!               an update at each coupled time step for the stress, and
      !!               (b) the conservative computation of the fluxes over the
      !!               sea-ice area requires the knowledge of the ice fraction
      !!               after the ice advection and before the ice thermodynamics,
      !!               so that the stress is updated before the ice dynamics
      !!               while the fluxes are updated after it.
      !!
      !! ** Details
      !!             qns_tot = pfrld * qns_oce + ( 1 - pfrld ) * qns_ice   => provided
      !!                     + qemp_oce + qemp_ice                         => recalculated and added up to qns
      !!
      !!             qsr_tot = pfrld * qsr_oce + ( 1 - pfrld ) * qsr_ice   => provided
      !!
      !!             emp_tot = emp_oce + emp_ice                           => calving is provided and added to emp_tot (and emp_oce)
      !!                                                                      river runoff (rnf) is provided but not included here
      !!
      !! ** Action  :   update at each nf_ice time step:
      !!                   qns_tot, qsr_tot  non-solar and solar total heat fluxes
      !!                   qns_ice, qsr_ice  non-solar and solar heat fluxes over the ice
      !!                   emp_tot           total evaporation - precipitation(liquid and solid) (-calving)
      !!                   emp_ice           ice sublimation - solid precipitation over the ice
      !!                   dqns_ice          d(non-solar heat flux)/d(Temperature) over the ice
      !!                   sprecip           solid precipitation over the ocean  
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ), DIMENSION(:,:)   ::   p_frld     ! lead fraction                [0 to 1]
      ! optional arguments, used only in 'mixed oce-ice' case
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   palbi      ! all skies ice albedo 
      REAL(wp), INTENT(in   ), DIMENSION(:,:  ), OPTIONAL ::   psst       ! sea surface temperature     [Celsius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   pist       ! ice surface temperature     [Kelvin]
      !
      INTEGER ::   jl         ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zcptn, ztmp, zicefr, zmsk, zsnw
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip, zevap_oce, zevap_ice, zdevap_ice
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_cpl_ice_flx')
      !
      CALL wrk_alloc( jpi,jpj,     zcptn, ztmp, zicefr, zmsk, zsnw )
      CALL wrk_alloc( jpi,jpj,     zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip, zevap_oce, zevap_ice, zdevap_ice )
      CALL wrk_alloc( jpi,jpj,     zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice )
      CALL wrk_alloc( jpi,jpj,jpl, zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice )

      IF( ln_mixcpl )   zmsk(:,:) = 1. - xcplmask(:,:,0)
      zicefr(:,:) = 1.- p_frld(:,:)
      zcptn(:,:) = rcp * sst_m(:,:)
      !
      !                                                      ! ========================= !
      !                                                      !    freshwater budget      !   (emp_tot)
      !                                                      ! ========================= !
      !
      !                                                           ! solid Precipitation                                (sprecip)
      !                                                           ! liquid + solid Precipitation                       (tprecip)
      !                                                           ! total Evaporation - total Precipitation            (emp_tot)
      !                                                           ! sublimation - solid precipitation (cell average)   (emp_ice)
      SELECT CASE( TRIM( sn_rcv_emp%cldes ) )
      CASE( 'conservative' )   ! received fields: jpr_rain, jpr_snow, jpr_ievp, jpr_tevp
         zsprecip(:,:) =   frcv(jpr_snow)%z3(:,:,1)                  ! May need to ensure positive here
         ztprecip(:,:) =   frcv(jpr_rain)%z3(:,:,1) + zsprecip(:,:)  ! May need to ensure positive here
         zemp_tot(:,:) =   frcv(jpr_tevp)%z3(:,:,1) - ztprecip(:,:)
         zemp_ice(:,:) = ( frcv(jpr_ievp)%z3(:,:,1) - frcv(jpr_snow)%z3(:,:,1) ) * zicefr(:,:)
         IF( iom_use('precip') )          &
            &  CALL iom_put( 'precip'       ,   frcv(jpr_rain)%z3(:,:,1) + frcv(jpr_snow)%z3(:,:,1)                              )  ! total  precipitation
         IF( iom_use('rain') )            &
            &  CALL iom_put( 'rain'         ,   frcv(jpr_rain)%z3(:,:,1)                                                         )  ! liquid precipitation 
         IF( iom_use('rain_ao_cea') )   &
            &  CALL iom_put( 'rain_ao_cea'  , frcv(jpr_rain)%z3(:,:,1)* p_frld(:,:) * tmask(:,:,1)      )   ! liquid precipitation 
         IF( iom_use('hflx_rain_cea') )   &
            CALL iom_put( 'hflx_rain_cea', frcv(jpr_rain)%z3(:,:,1) * zcptn(:,:) * tmask(:,:,1))   ! heat flux from liq. precip. 
         IF( iom_use('hflx_prec_cea') )   &
            CALL iom_put( 'hflx_prec_cea', ztprecip * zcptn(:,:) * tmask(:,:,1) * p_frld(:,:) )   ! heat content flux from all precip  (cell avg)
         IF( iom_use('evap_ao_cea') .OR. iom_use('hflx_evap_cea') )   &
            ztmp(:,:) = frcv(jpr_tevp)%z3(:,:,1) - frcv(jpr_ievp)%z3(:,:,1) * zicefr(:,:)
         IF( iom_use('evap_ao_cea'  ) )   &
            CALL iom_put( 'evap_ao_cea'  , ztmp * tmask(:,:,1)                  )   ! ice-free oce evap (cell average)
         IF( iom_use('hflx_evap_cea') )   &
            CALL iom_put( 'hflx_evap_cea', ztmp(:,:) * zcptn(:,:) * tmask(:,:,1) )   ! heat flux from from evap (cell average)
      CASE( 'oce and ice'   )   ! received fields: jpr_sbpr, jpr_semp, jpr_oemp, jpr_ievp
         zemp_tot(:,:) = p_frld(:,:) * frcv(jpr_oemp)%z3(:,:,1) + zicefr(:,:) * frcv(jpr_sbpr)%z3(:,:,1)
         zemp_ice(:,:) = frcv(jpr_semp)%z3(:,:,1) * zicefr(:,:)
         zsprecip(:,:) = frcv(jpr_ievp)%z3(:,:,1) - frcv(jpr_semp)%z3(:,:,1)
         ztprecip(:,:) = frcv(jpr_semp)%z3(:,:,1) - frcv(jpr_sbpr)%z3(:,:,1) + zsprecip(:,:)
      END SELECT

      ! runoffs and calving (put in emp_tot)
      IF( srcv(jpr_rnf)%laction )   rnf(:,:) = frcv(jpr_rnf)%z3(:,:,1)
      IF( iom_use('hflx_rnf_cea') )   &
         CALL iom_put( 'hflx_rnf_cea' , rnf(:,:) * zcptn(:,:) )
      IF( srcv(jpr_cal)%laction ) THEN 
         zemp_tot(:,:) = zemp_tot(:,:) - frcv(jpr_cal)%z3(:,:,1)
         CALL iom_put( 'calving_cea', frcv(jpr_cal)%z3(:,:,1) )
      ENDIF


      IF( srcv(jpr_icb)%laction )  THEN 
         fwficb(:,:) = frcv(jpr_icb)%z3(:,:,1)
         rnf(:,:)    = rnf(:,:) + fwficb(:,:)   ! iceberg added to runoffs
         CALL iom_put( 'iceberg_cea', frcv(jpr_icb)%z3(:,:,1) )
      ENDIF
      IF( srcv(jpr_isf)%laction )  THEN
        fwfisf(:,:) = - frcv(jpr_isf)%z3(:,:,1)  ! fresh water flux from the isf (fwfisf <0 mean melting)  
        CALL iom_put( 'iceshelf_cea', frcv(jpr_isf)%z3(:,:,1) )
      ENDIF


      IF( ln_mixcpl ) THEN
         emp_tot(:,:) = emp_tot(:,:) * xcplmask(:,:,0) + zemp_tot(:,:) * zmsk(:,:)
         emp_ice(:,:) = emp_ice(:,:) * xcplmask(:,:,0) + zemp_ice(:,:) * zmsk(:,:)
         sprecip(:,:) = sprecip(:,:) * xcplmask(:,:,0) + zsprecip(:,:) * zmsk(:,:)
         tprecip(:,:) = tprecip(:,:) * xcplmask(:,:,0) + ztprecip(:,:) * zmsk(:,:)
      ELSE
         emp_tot(:,:) =                                  zemp_tot(:,:)
         emp_ice(:,:) =                                  zemp_ice(:,:)
         sprecip(:,:) =                                  zsprecip(:,:)
         tprecip(:,:) =                                  ztprecip(:,:)
      ENDIF

      IF( iom_use('subl_ai_cea') )  CALL iom_put( 'subl_ai_cea', frcv(jpr_ievp)%z3(:,:,1) * zicefr(:,:) )  ! Sublimation over sea-ice (cell average)
                                    CALL iom_put( 'snowpre'    , sprecip(:,:)               )   ! Snow
      IF( iom_use('snow_ao_cea') )  CALL iom_put( 'snow_ao_cea', sprecip(:,:) * p_frld(:,:) )   ! Snow over ice-free ocean  (cell average)
      IF( iom_use('snow_ai_cea') )  CALL iom_put( 'snow_ai_cea', sprecip(:,:) * zicefr(:,:) )   ! Snow over sea-ice         (cell average)

      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_qns%cldes ) )                !   non solar heat fluxes   !   (qns)
      !                                                      ! ========================= !
      CASE( 'oce only' )         ! the required field is directly provided
         zqns_tot(:,:) = frcv(jpr_qnsoce)%z3(:,:,1)
      CASE( 'conservative' )     ! the required fields are directly provided
         zqns_tot(:,:) = frcv(jpr_qnsmix)%z3(:,:,1)
         IF ( TRIM(sn_rcv_qns%clcat) == 'yes' ) THEN
            zqns_ice(:,:,1:jpl) = frcv(jpr_qnsice)%z3(:,:,1:jpl)
         ELSE
            DO jl=1,jpl
               zqns_ice(:,:,jl) = frcv(jpr_qnsice)%z3(:,:,1) ! Set all category values equal
            ENDDO
         ENDIF
      CASE( 'oce and ice' )      ! the total flux is computed from ocean and ice fluxes
         zqns_tot(:,:) =  p_frld(:,:) * frcv(jpr_qnsoce)%z3(:,:,1)
         IF ( TRIM(sn_rcv_qns%clcat) == 'yes' ) THEN
            DO jl=1,jpl
               zqns_tot(:,:   ) = zqns_tot(:,:) + a_i(:,:,jl) * frcv(jpr_qnsice)%z3(:,:,jl)   
               zqns_ice(:,:,jl) = frcv(jpr_qnsice)%z3(:,:,jl)
            ENDDO
         ELSE
            qns_tot(:,:) = qns_tot(:,:) + zicefr(:,:) * frcv(jpr_qnsice)%z3(:,:,1)
            DO jl=1,jpl
               zqns_tot(:,:   ) = zqns_tot(:,:) + zicefr(:,:) * frcv(jpr_qnsice)%z3(:,:,1)
               zqns_ice(:,:,jl) = frcv(jpr_qnsice)%z3(:,:,1)
            ENDDO
         ENDIF
      CASE( 'mixed oce-ice' )    ! the ice flux is cumputed from the total flux, the SST and ice informations
! ** NEED TO SORT OUT HOW THIS SHOULD WORK IN THE MULTI-CATEGORY CASE - CURRENTLY NOT ALLOWED WHEN INTERFACE INITIALISED **
         zqns_tot(:,:  ) = frcv(jpr_qnsmix)%z3(:,:,1)
         zqns_ice(:,:,1) = frcv(jpr_qnsmix)%z3(:,:,1)    &
            &            + frcv(jpr_dqnsdt)%z3(:,:,1) * ( pist(:,:,1) - ( (rt0 + psst(:,:  ) ) * p_frld(:,:)   &
            &                                           + pist(:,:,1) * zicefr(:,:) ) )
      END SELECT
!!gm
!!    currently it is taken into account in leads budget but not in the zqns_tot, and thus not in 
!!    the flux that enter the ocean....
!!    moreover 1 - it is not diagnose anywhere.... 
!!             2 - it is unclear for me whether this heat lost is taken into account in the atmosphere or not...
!!
!! similar job should be done for snow and precipitation temperature
      !                                     
      IF( srcv(jpr_cal)%laction ) THEN   ! Iceberg melting 
         zqns_tot(:,:) = zqns_tot(:,:) - frcv(jpr_cal)%z3(:,:,1) * lfus  ! add the latent heat of iceberg melting
                                                                         ! we suppose it melts at 0deg, though it should be temp. of surrounding ocean
         IF( iom_use('hflx_cal_cea') )   CALL iom_put( 'hflx_cal_cea', - frcv(jpr_cal)%z3(:,:,1) * lfus )   ! heat flux from calving
      ENDIF

!!chris     
!!    The heat content associated to the ice shelf in removed in the routine sbcisf.F90
      !
      IF( srcv(jpr_icb)%laction )  zqns_tot(:,:) = zqns_tot(:,:) - frcv(jpr_icb)%z3(:,:,1) * lfus ! remove heat content associated to iceberg melting
      !
!!      !

      ! clem: this formulation is certainly wrong... but better than it was...
      zqns_tot(:,:) = zqns_tot(:,:)                       &            ! zqns_tot update over free ocean with:
         &          - ztmp(:,:)                           &            ! remove the latent heat flux of solid precip. melting
         &          - (  zemp_tot(:,:)                    &            ! remove the heat content of mass flux (assumed to be at SST)
         &             - zemp_ice(:,:) ) * zcptn(:,:) 

     IF( ln_mixcpl ) THEN
         qns_tot(:,:) = qns(:,:) * p_frld(:,:) + SUM( qns_ice(:,:,:) * a_i(:,:,:), dim=3 )   ! total flux from blk
         qns_tot(:,:) = qns_tot(:,:) * xcplmask(:,:,0) +  zqns_tot(:,:)* zmsk(:,:)
         DO jl=1,jpl
            qns_ice(:,:,jl) = qns_ice(:,:,jl) * xcplmask(:,:,0) +  zqns_ice(:,:,jl)* zmsk(:,:)
         ENDDO
      ELSE
         qns_tot(:,:  ) = zqns_tot(:,:  )
         qns_ice(:,:,:) = zqns_ice(:,:,:)
      ENDIF

      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_qsr%cldes ) )                !      solar heat fluxes    !   (qsr)
      !                                                      ! ========================= !
      CASE( 'oce only' )
         zqsr_tot(:,:  ) = MAX( 0._wp , frcv(jpr_qsroce)%z3(:,:,1) )
      CASE( 'conservative' )
         zqsr_tot(:,:  ) = frcv(jpr_qsrmix)%z3(:,:,1)
         IF ( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            zqsr_ice(:,:,1:jpl) = frcv(jpr_qsrice)%z3(:,:,1:jpl)
         ELSE
            ! Set all category values equal for the moment
            DO jl=1,jpl
               zqsr_ice(:,:,jl) = frcv(jpr_qsrice)%z3(:,:,1)
            ENDDO
         ENDIF
         zqsr_tot(:,:  ) = frcv(jpr_qsrmix)%z3(:,:,1)
         zqsr_ice(:,:,1) = frcv(jpr_qsrice)%z3(:,:,1)
      CASE( 'oce and ice' )
         zqsr_tot(:,:  ) =  p_frld(:,:) * frcv(jpr_qsroce)%z3(:,:,1)
         IF ( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            DO jl=1,jpl
               zqsr_tot(:,:   ) = zqsr_tot(:,:) + a_i(:,:,jl) * frcv(jpr_qsrice)%z3(:,:,jl)   
               zqsr_ice(:,:,jl) = frcv(jpr_qsrice)%z3(:,:,jl)
            ENDDO
         ELSE
            qsr_tot(:,:   ) = qsr_tot(:,:) + zicefr(:,:) * frcv(jpr_qsrice)%z3(:,:,1)
            DO jl=1,jpl
               zqsr_tot(:,:   ) = zqsr_tot(:,:) + zicefr(:,:) * frcv(jpr_qsrice)%z3(:,:,1)
               zqsr_ice(:,:,jl) = frcv(jpr_qsrice)%z3(:,:,1)
            ENDDO
         ENDIF
      CASE( 'mixed oce-ice' )
         zqsr_tot(:,:  ) = frcv(jpr_qsrmix)%z3(:,:,1)
! ** NEED TO SORT OUT HOW THIS SHOULD WORK IN THE MULTI-CATEGORY CASE - CURRENTLY NOT ALLOWED WHEN INTERFACE INITIALISED **
!       Create solar heat flux over ice using incoming solar heat flux and albedos
!       ( see OASIS3 user guide, 5th edition, p39 )
         zqsr_ice(:,:,1) = frcv(jpr_qsrmix)%z3(:,:,1) * ( 1.- palbi(:,:,1) )   &
            &            / (  1.- ( albedo_oce_mix(:,:  ) * p_frld(:,:)       &
            &                     + palbi         (:,:,1) * zicefr(:,:) ) )
      END SELECT
      IF( ln_dm2dc .AND. ln_cpl ) THEN   ! modify qsr to include the diurnal cycle
         zqsr_tot(:,:  ) = sbc_dcy( zqsr_tot(:,:  ) )
         DO jl=1,jpl
            zqsr_ice(:,:,jl) = sbc_dcy( zqsr_ice(:,:,jl) )
         ENDDO
      ENDIF


      IF( ln_mixcpl ) THEN
         qsr_tot(:,:) = qsr(:,:) * p_frld(:,:) + SUM( qsr_ice(:,:,:) * a_i(:,:,:), dim=3 )   ! total flux from blk
         qsr_tot(:,:) = qsr_tot(:,:) * xcplmask(:,:,0) +  zqsr_tot(:,:)* zmsk(:,:)
         DO jl=1,jpl
            qsr_ice(:,:,jl) = qsr_ice(:,:,jl) * xcplmask(:,:,0) +  zqsr_ice(:,:,jl)* zmsk(:,:)
         ENDDO
      ELSE
         qsr_tot(:,:  ) = zqsr_tot(:,:  )
         qsr_ice(:,:,:) = zqsr_ice(:,:,:)
      ENDIF

      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_dqnsdt%cldes ) )             !          d(qns)/dt        !
      !                                                      ! ========================= !
      CASE ('coupled')
         IF ( TRIM(sn_rcv_dqnsdt%clcat) == 'yes' ) THEN
            zdqns_ice(:,:,1:jpl) = frcv(jpr_dqnsdt)%z3(:,:,1:jpl)
         ELSE
            ! Set all category values equal for the moment
            DO jl=1,jpl
               zdqns_ice(:,:,jl) = frcv(jpr_dqnsdt)%z3(:,:,1)
            ENDDO
         ENDIF
      END SELECT
      
      IF( ln_mixcpl ) THEN
         DO jl=1,jpl
            dqns_ice(:,:,jl) = dqns_ice(:,:,jl) * xcplmask(:,:,0) + zdqns_ice(:,:,jl) * zmsk(:,:)
         ENDDO
      ELSE
         dqns_ice(:,:,:) = zdqns_ice(:,:,:)
      ENDIF
      
      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_iceflx%cldes ) )             !    topmelt and botmelt    !
      !                                                      ! ========================= !
      CASE ('coupled')
         topmelt(:,:,:)=frcv(jpr_topm)%z3(:,:,:)
         botmelt(:,:,:)=frcv(jpr_botm)%z3(:,:,:)
      END SELECT

      ! Surface transimission parameter io (Maykut Untersteiner , 1971 ; Ebert and Curry, 1993 )
      ! Used for LIM2 and LIM3
      ! Coupled case: since cloud cover is not received from atmosphere 
      !               ===> used prescribed cloud fraction representative for polar oceans in summer (0.81)
      fr1_i0(:,:) = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )
      fr2_i0(:,:) = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )

      CALL wrk_dealloc( jpi,jpj,     zcptn, ztmp, zicefr, zmsk, zsnw )
      CALL wrk_dealloc( jpi,jpj,     zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip, zevap_oce, zevap_ice, zdevap_ice )
      CALL wrk_dealloc( jpi,jpj,     zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice )
      CALL wrk_dealloc( jpi,jpj,jpl, zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_cpl_ice_flx')
      !
   END SUBROUTINE sbc_cpl_ice_flx
   
   
   SUBROUTINE sbc_cpl_snd( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_snd  ***
      !!
      !! ** Purpose :   provide the ocean-ice informations to the atmosphere
      !!
      !! ** Method  :   send to the atmosphere through a call to cpl_snd
      !!              all the needed fields (as defined in sbc_cpl_init)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      !
      INTEGER ::   ji, jj, jl   ! dummy loop indices
      INTEGER ::   isec, info   ! local integer
      REAL(wp) ::   zumax, zvmax
      REAL(wp), POINTER, DIMENSION(:,:)   ::   zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztmp3, ztmp4   
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_cpl_snd')
      !
      CALL wrk_alloc( jpi,jpj, zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1 )
      CALL wrk_alloc( jpi,jpj,jpl, ztmp3, ztmp4 )

      isec = ( kt - nit000 ) * NINT(rdttra(1))        ! date of exchanges

      zfr_l(:,:) = 1.- fr_i(:,:)
      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !   in Kelvin
      !                                                      ! ------------------------- !
      IF( ssnd(jps_toce)%laction .OR. ssnd(jps_tice)%laction .OR. ssnd(jps_tmix)%laction ) THEN
         
         IF ( nn_components == jp_iam_opa ) THEN
            ztmp1(:,:) = tsn(:,:,1,jp_tem)   ! send temperature as it is (potential or conservative) -> use of ln_useCT on the received part
         ELSE
            ! we must send the surface potential temperature 
            IF( ln_useCT )  THEN    ;   ztmp1(:,:) = eos_pt_from_ct( tsn(:,:,1,jp_tem), tsn(:,:,1,jp_sal) )
            ELSE                    ;   ztmp1(:,:) = tsn(:,:,1,jp_tem)
            ENDIF
            !
            SELECT CASE( sn_snd_temp%cldes)
            CASE( 'oce only'             )   ;   ztmp1(:,:) =   ztmp1(:,:) + rt0
            CASE( 'oce and ice'          )   ;   ztmp1(:,:) =   ztmp1(:,:) + rt0
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )   
                  ztmp3(:,:,1:jpl) = tn_ice(:,:,1:jpl)
               CASE( 'no' )
                  WHERE( SUM( a_i, dim=3 ) /= 0. )
                     ztmp3(:,:,1) = SUM( tn_ice * a_i, dim=3 ) / SUM( a_i, dim=3 )
                  ELSEWHERE
                     ztmp3(:,:,1) = rt0
                  END WHERE
               CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'weighted oce and ice' )   ;   ztmp1(:,:) = ( ztmp1(:,:) + rt0 ) * zfr_l(:,:)   
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )   
                  ztmp3(:,:,1:jpl) = tn_ice(:,:,1:jpl) * a_i(:,:,1:jpl)
               CASE( 'no' )
                  ztmp3(:,:,:) = 0.0
                  DO jl=1,jpl
                     ztmp3(:,:,1) = ztmp3(:,:,1) + tn_ice(:,:,jl) * a_i(:,:,jl)
                  ENDDO
               CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'mixed oce-ice'        )   
               ztmp1(:,:) = ( ztmp1(:,:) + rt0 ) * zfr_l(:,:) 
               DO jl=1,jpl
                  ztmp1(:,:) = ztmp1(:,:) + tn_ice(:,:,jl) * a_i(:,:,jl)
               ENDDO
            CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%cldes' )
            END SELECT
         ENDIF
         IF( ssnd(jps_toce)%laction )   CALL cpl_snd( jps_toce, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
         IF( ssnd(jps_tice)%laction )   CALL cpl_snd( jps_tice, isec, ztmp3, info )
         IF( ssnd(jps_tmix)%laction )   CALL cpl_snd( jps_tmix, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !           Albedo          !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_albice)%laction ) THEN                         ! ice 
          SELECT CASE( sn_snd_alb%cldes )
          CASE( 'ice' )
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )   
                ztmp3(:,:,1:jpl) = alb_ice(:,:,1:jpl)
             CASE( 'no' )
                WHERE( SUM( a_i, dim=3 ) /= 0. )
                   ztmp1(:,:) = SUM( alb_ice (:,:,1:jpl) * a_i(:,:,1:jpl), dim=3 ) / SUM( a_i(:,:,1:jpl), dim=3 )
                ELSEWHERE
                   ztmp1(:,:) = albedo_oce_mix(:,:)
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%clcat' )
             END SELECT
          CASE( 'weighted ice' )   ;
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )   
                ztmp3(:,:,1:jpl) =  alb_ice(:,:,1:jpl) * a_i(:,:,1:jpl)
             CASE( 'no' )
                WHERE( fr_i (:,:) > 0. )
                   ztmp1(:,:) = SUM (  alb_ice(:,:,1:jpl) * a_i(:,:,1:jpl), dim=3 )
                ELSEWHERE
                   ztmp1(:,:) = 0.
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_ice%clcat' )
             END SELECT
          CASE default      ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%cldes' )
         END SELECT

         SELECT CASE( sn_snd_alb%clcat )
            CASE( 'yes' )   
               CALL cpl_snd( jps_albice, isec, ztmp3, info )      !-> MV this has never been checked in coupled mode
            CASE( 'no'  )   
               CALL cpl_snd( jps_albice, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info ) 
         END SELECT
      ENDIF

      IF( ssnd(jps_albmix)%laction ) THEN                         ! mixed ice-ocean
         ztmp1(:,:) = albedo_oce_mix(:,:) * zfr_l(:,:)
         DO jl=1,jpl
            ztmp1(:,:) = ztmp1(:,:) + alb_ice(:,:,jl) * a_i(:,:,jl)
         ENDDO
         CALL cpl_snd( jps_albmix, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness ! 
      !                                                      ! ------------------------- !
      ! Send ice fraction field to atmosphere
      IF( ssnd(jps_fice)%laction ) THEN
         SELECT CASE( sn_snd_thick%clcat )
         CASE( 'yes' )   ;   ztmp3(:,:,1:jpl) =  a_i(:,:,1:jpl)
         CASE( 'no'  )   ;   ztmp3(:,:,1    ) = fr_i(:,:      )
         CASE default    ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
         END SELECT
         IF( ssnd(jps_fice)%laction )   CALL cpl_snd( jps_fice, isec, ztmp3, info )
      ENDIF
      
      ! Send ice fraction field to OPA (sent by SAS in SAS-OPA coupling)
      IF( ssnd(jps_fice2)%laction ) THEN
         ztmp3(:,:,1) = fr_i(:,:)
         IF( ssnd(jps_fice2)%laction )   CALL cpl_snd( jps_fice2, isec, ztmp3, info )
      ENDIF

      ! Send ice and snow thickness field 
      IF( ssnd(jps_hice)%laction .OR. ssnd(jps_hsnw)%laction ) THEN 
         SELECT CASE( sn_snd_thick%cldes)
         CASE( 'none'                  )       ! nothing to do
         CASE( 'weighted ice and snow' )   
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )   
               ztmp3(:,:,1:jpl) =  ht_i(:,:,1:jpl) * a_i(:,:,1:jpl)
               ztmp4(:,:,1:jpl) =  ht_s(:,:,1:jpl) * a_i(:,:,1:jpl)
            CASE( 'no' )
               ztmp3(:,:,:) = 0.0   ;  ztmp4(:,:,:) = 0.0
               DO jl=1,jpl
                  ztmp3(:,:,1) = ztmp3(:,:,1) + ht_i(:,:,jl) * a_i(:,:,jl)
                  ztmp4(:,:,1) = ztmp4(:,:,1) + ht_s(:,:,jl) * a_i(:,:,jl)
               ENDDO
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE( 'ice and snow'         )   
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )
               ztmp3(:,:,1:jpl) = ht_i(:,:,1:jpl)
               ztmp4(:,:,1:jpl) = ht_s(:,:,1:jpl)
            CASE( 'no' )
               WHERE( SUM( a_i, dim=3 ) /= 0. )
                  ztmp3(:,:,1) = SUM( ht_i * a_i, dim=3 ) / SUM( a_i, dim=3 )
                  ztmp4(:,:,1) = SUM( ht_s * a_i, dim=3 ) / SUM( a_i, dim=3 )
               ELSEWHERE
                 ztmp3(:,:,1) = 0.
                 ztmp4(:,:,1) = 0.
               END WHERE
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%cldes' )
         END SELECT
         IF( ssnd(jps_hice)%laction )   CALL cpl_snd( jps_hice, isec, ztmp3, info )
         IF( ssnd(jps_hsnw)%laction )   CALL cpl_snd( jps_hsnw, isec, ztmp4, info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_ocx1)%laction ) THEN                      !      Surface current      !
         !                                                   ! ------------------------- !
         !    
         !                                                  j+1   j     -----V---F
         ! surface velocity always sent from T point                     !       |
         !                                                        j      |   T   U
         !                                                               |       |
         !                                                   j    j-1   -I-------|
         !                                               (for I)         |       |
         !                                                              i-1  i   i
         !                                                               i      i+1 (for I)
         IF( nn_components == jp_iam_opa ) THEN
            zotx1(:,:) = un(:,:,1)  
            zoty1(:,:) = vn(:,:,1)  
         ELSE        
            SELECT CASE( TRIM( sn_snd_crt%cldes ) )
            CASE( 'oce only'             )      ! C-grid ==> T
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     zotx1(ji,jj) = 0.5 * ( un(ji,jj,1) + un(ji-1,jj  ,1) )
                     zoty1(ji,jj) = 0.5 * ( vn(ji,jj,1) + vn(ji  ,jj-1,1) ) 
                  END DO
               END DO
            CASE( 'weighted oce and ice' )   
               SELECT CASE ( cp_ice_msh )
               CASE( 'C' )                      ! Ocean and Ice on C-grid ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        zotx1(ji,jj) = 0.5 * ( un   (ji,jj,1) + un   (ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                        zoty1(ji,jj) = 0.5 * ( vn   (ji,jj,1) + vn   (ji  ,jj-1,1) ) * zfr_l(ji,jj)
                        zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                        zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
                     END DO
                  END DO
               CASE( 'I' )                      ! Ocean on C grid, Ice on I-point (B-grid) ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! NO vector opt.
                        zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                        zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)  
                        zitx1(ji,jj) = 0.25 * ( u_ice(ji+1,jj+1) + u_ice(ji,jj+1)                     &
                           &                  + u_ice(ji+1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                        zity1(ji,jj) = 0.25 * ( v_ice(ji+1,jj+1) + v_ice(ji,jj+1)                     &
                           &                  + v_ice(ji+1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     END DO
                  END DO
               CASE( 'F' )                      ! Ocean on C grid, Ice on F-point (B-grid) ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! NO vector opt.
                        zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                        zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)  
                        zitx1(ji,jj) = 0.25 * ( u_ice(ji-1,jj-1) + u_ice(ji,jj-1)                     &
                           &                  + u_ice(ji-1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                        zity1(ji,jj) = 0.25 * ( v_ice(ji-1,jj-1) + v_ice(ji,jj-1)                     &
                           &                  + v_ice(ji-1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     END DO
                  END DO
               END SELECT
               CALL lbc_lnk( zitx1, 'T', -1. )   ;   CALL lbc_lnk( zity1, 'T', -1. )
            CASE( 'mixed oce-ice'        )
               SELECT CASE ( cp_ice_msh )
               CASE( 'C' )                      ! Ocean and Ice on C-grid ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        zotx1(ji,jj) = 0.5 * ( un   (ji,jj,1) + un   (ji-1,jj  ,1) ) * zfr_l(ji,jj)   &
                           &         + 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                        zoty1(ji,jj) = 0.5 * ( vn   (ji,jj,1) + vn   (ji  ,jj-1,1) ) * zfr_l(ji,jj)   &
                           &         + 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
                     END DO
                  END DO
               CASE( 'I' )                      ! Ocean on C grid, Ice on I-point (B-grid) ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! NO vector opt.
                        zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)   &   
                           &         + 0.25 * ( u_ice(ji+1,jj+1) + u_ice(ji,jj+1)                     &
                           &                  + u_ice(ji+1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                        zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)   & 
                           &         + 0.25 * ( v_ice(ji+1,jj+1) + v_ice(ji,jj+1)                     &
                           &                  + v_ice(ji+1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     END DO
                  END DO
               CASE( 'F' )                      ! Ocean on C grid, Ice on F-point (B-grid) ==> T
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! NO vector opt.
                        zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)   &   
                           &         + 0.25 * ( u_ice(ji-1,jj-1) + u_ice(ji,jj-1)                     &
                           &                  + u_ice(ji-1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                        zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)   & 
                           &         + 0.25 * ( v_ice(ji-1,jj-1) + v_ice(ji,jj-1)                     &
                           &                  + v_ice(ji-1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     END DO
                  END DO
               END SELECT
            END SELECT
            CALL lbc_lnk( zotx1, ssnd(jps_ocx1)%clgrid, -1. )   ;   CALL lbc_lnk( zoty1, ssnd(jps_ocy1)%clgrid, -1. )
            !
         ENDIF
         !
         !
         IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) THEN             ! Rotation of the components
            !                                                                     ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->e', ztmp1 )       ! 1st component 
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->n', ztmp2 )       ! 2nd component 
            zotx1(:,:) = ztmp1(:,:)                                                   ! overwrite the components 
            zoty1(:,:) = ztmp2(:,:)
            IF( ssnd(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component 
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component 
               zitx1(:,:) = ztmp1(:,:)                                                ! overwrite the components 
               zity1(:,:) = ztmp2(:,:)
            ENDIF
         ENDIF
         !
         ! spherical coordinates to cartesian -> 2 components to 3 components
         IF( TRIM( sn_snd_crt%clvref ) == 'cartesian' ) THEN
            ztmp1(:,:) = zotx1(:,:)                     ! ocean currents
            ztmp2(:,:) = zoty1(:,:)
            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
            !
            IF( ssnd(jps_ivx1)%laction ) THEN           ! ice velocities
               ztmp1(:,:) = zitx1(:,:)
               ztmp1(:,:) = zity1(:,:)
               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
            ENDIF
         ENDIF
         !
         IF( ssnd(jps_ocx1)%laction )   CALL cpl_snd( jps_ocx1, isec, RESHAPE ( zotx1, (/jpi,jpj,1/) ), info )   ! ocean x current 1st grid
         IF( ssnd(jps_ocy1)%laction )   CALL cpl_snd( jps_ocy1, isec, RESHAPE ( zoty1, (/jpi,jpj,1/) ), info )   ! ocean y current 1st grid
         IF( ssnd(jps_ocz1)%laction )   CALL cpl_snd( jps_ocz1, isec, RESHAPE ( zotz1, (/jpi,jpj,1/) ), info )   ! ocean z current 1st grid
         !
         IF( ssnd(jps_ivx1)%laction )   CALL cpl_snd( jps_ivx1, isec, RESHAPE ( zitx1, (/jpi,jpj,1/) ), info )   ! ice   x current 1st grid
         IF( ssnd(jps_ivy1)%laction )   CALL cpl_snd( jps_ivy1, isec, RESHAPE ( zity1, (/jpi,jpj,1/) ), info )   ! ice   y current 1st grid
         IF( ssnd(jps_ivz1)%laction )   CALL cpl_snd( jps_ivz1, isec, RESHAPE ( zitz1, (/jpi,jpj,1/) ), info )   ! ice   z current 1st grid
         ! 
      ENDIF
      !
      !
      !  Fields sent by OPA to SAS when doing OPA<->SAS coupling
      !                                                        ! SSH
      IF( ssnd(jps_ssh )%laction )  THEN
         !                          ! removed inverse barometer ssh when Patm
         !                          forcing is used (for sea-ice dynamics)
         IF( ln_apr_dyn ) THEN   ;   ztmp1(:,:) = sshb(:,:) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
         ELSE                    ;   ztmp1(:,:) = sshn(:,:)
         ENDIF
         CALL cpl_snd( jps_ssh   , isec, RESHAPE ( ztmp1            , (/jpi,jpj,1/) ), info )

      ENDIF
      !                                                        ! SSS
      IF( ssnd(jps_soce  )%laction )  THEN
         CALL cpl_snd( jps_soce  , isec, RESHAPE ( tsn(:,:,1,jp_sal), (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                        ! first T level thickness 
      IF( ssnd(jps_e3t1st )%laction )  THEN
         CALL cpl_snd( jps_e3t1st, isec, RESHAPE ( e3t_n(:,:,1)   , (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                        ! Qsr fraction
      IF( ssnd(jps_fraqsr)%laction )  THEN
         CALL cpl_snd( jps_fraqsr, isec, RESHAPE ( fraqsr_1lev(:,:) , (/jpi,jpj,1/) ), info )
      ENDIF
      !
      !  Fields sent by SAS to OPA when OASIS coupling
      !                                                        ! Solar heat flux
      IF( ssnd(jps_qsroce)%laction )  CALL cpl_snd( jps_qsroce, isec, RESHAPE ( qsr , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_qnsoce)%laction )  CALL cpl_snd( jps_qnsoce, isec, RESHAPE ( qns , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_oemp  )%laction )  CALL cpl_snd( jps_oemp  , isec, RESHAPE ( emp , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_sflx  )%laction )  CALL cpl_snd( jps_sflx  , isec, RESHAPE ( sfx , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_otx1  )%laction )  CALL cpl_snd( jps_otx1  , isec, RESHAPE ( utau, (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_oty1  )%laction )  CALL cpl_snd( jps_oty1  , isec, RESHAPE ( vtau, (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_rnf   )%laction )  CALL cpl_snd( jps_rnf   , isec, RESHAPE ( rnf , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_taum  )%laction )  CALL cpl_snd( jps_taum  , isec, RESHAPE ( taum, (/jpi,jpj,1/) ), info )

      CALL wrk_dealloc( jpi,jpj, zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1 )
      CALL wrk_dealloc( jpi,jpj,jpl, ztmp3, ztmp4 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_cpl_snd')
      !
   END SUBROUTINE sbc_cpl_snd
   
   !!======================================================================
END MODULE sbccpl
