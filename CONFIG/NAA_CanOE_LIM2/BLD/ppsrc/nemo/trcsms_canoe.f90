MODULE trcsms_canoe
   !!======================================================================
   !!                         ***  MODULE trcsms_canoe  ***
   !! TOP :   CanOE Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_canoe        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE trcdta
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE iom             !  I/O manager
   USE trd_oce         !  Ocean trends variables
   USE trdtrc          !  TOP trends variables
   USE sedmodel        !  Sediment model
   USE prtctl_trc      !  print control for debugging

   USE canbio          !  Biological model
   USE canche          !  Chemical model
   USE canlys          !  Calcite saturation
   USE canflx          !  Gas exchange
   USE cansed          !  Sedimentation
   USE cansbc          !  Sedimentation
   USE canint          !  time interpolation
   USE par_canoe

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_canoe    ! called in trcsms.F90
   PUBLIC   can_sms_init    ! called in trcsms.F90
    
   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: ???
   INTEGER ::  numno3  !: logical unit for NO3 budget
   INTEGER ::  numalk  !: logical unit for talk budget
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_canoe.F90 4147 2013-11-04 11:51:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_canoe ***
      !!
      !! ** Purpose :   Initialisation of the CanOE biochemical model
      !!----------------------------------------------------------------------


   SUBROUTINE trc_sms_canoe( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_canoe  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of CanOE bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   ji, jj, jk, knt, jn, jl
      REAL(wp) ::  ztra
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_canoe')
      !
      IF( kt == nittrc000 ) THEN
        !
      !* Array used to indicate negative tracer values  
        ALLOCATE( xnegtr(jpi,jpj,jpk) )
        !
        CALL can_che                              ! initialize the chemical constants
        !
        IF( .NOT. ln_rsttr ) THEN  ;   CALL can_ph_ini   !  set PH at kt=nit000 
        ELSE                       ;   CALL can_rst( nittrc000, 'READ' )  !* read or initialize all required fields 
        ENDIF
        !
      ENDIF
      !
      IF( ln_pisdmp .AND. MOD( kt - nn_dttrc, nn_pisdmp ) == 0 )   CALL can_dmp( kt )      ! Relaxation of some tracers
      !
      !                                                                    !   set time step size (Euler/Leapfrog)
      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN   ;    rfact = rdttrc(1)     !  at nittrc000
      ELSEIF( kt <= nittrc000 + nn_dttrc )                          THEN   ;    rfact = 2. * rdttrc(1)   ! at nittrc000 or nittrc000+nn_dttrc (Leapfrog)
      ENDIF
      !
      IF( ( ln_top_euler .AND. kt == nittrc000 )  .OR. ( .NOT.ln_top_euler .AND. kt <= nittrc000 + nn_dttrc ) ) THEN
         rfactr  = 1. / rfact
         rfact2  = rfact / FLOAT( nrdttrc )
         rfact2r = 1. / rfact2
         xstep = rfact2 / rday         ! Time step duration for biology
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rfact  = ', rfact, ' rdt = ', rdttra(1)
         IF(lwp) write(numout,*) '    CanOE  Biology time step    rfact2 = ', rfact2
         IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN
         DO jn = jp_can0, jp_can1              !   SMS on tracer without Asselin time-filter
            trb(:,:,:,jn) = trn(:,:,:,jn)
         END DO
      ENDIF
      !
      IF( ndayflxtr /= nday_year ) THEN      ! New days
         !
         ndayflxtr = nday_year

         IF(lwp) write(numout,*)
         IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
         IF(lwp) write(numout,*) '~~~~~~'

         CALL can_che              ! computation of chemical constants
         CALL can_int( kt )        ! computation of various rates for biogeochemistry
         !
      ENDIF

      IF( ll_sbc ) CALL can_sbc( kt )   ! external sources of nutrients 

      DO knt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL can_bio( kt, knt )   ! Biology
         CALL can_sed( kt, knt )   ! Sedimentation
         CALL can_lys( kt, knt )   ! Compute CaCO3 saturation
         CALL can_flx( kt, knt )   ! Compute surface fluxes
         !
         xnegtr(:,:,:) = 1.e0
         DO jn = jp_can0, jp_can1
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( ( trb(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN
                        ztra             = ABS( trb(ji,jj,jk,jn) ) / ( ABS( tra(ji,jj,jk,jn) ) + rtrn )
                        xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                     ENDIF
                 END DO
               END DO
            END DO
         END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
         DO jn = jp_can0, jp_can1
           trb(:,:,:,jn) = trb(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
         END DO
        !
         DO jn = jp_can0, jp_can1
            tra(:,:,:,jn) = 0._wp
         END DO
         !
         IF( ln_top_euler ) THEN
            DO jn = jp_can0, jp_can1
               trn(:,:,:,jn) = trb(:,:,:,jn)
            END DO
         ENDIF
      END DO

      !
      !
      IF( l_trdtrc ) THEN
         DO jn = jp_can0, jp_can1
           CALL trd_trc( tra(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
         END DO
      END IF
      !
      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
         DO jn = jp_can0, jp_can1
           CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF
      !
      IF( lrst_trc )  CALL can_rst( kt, 'WRITE' )  !* Write CanOE informations in restart file 
      !

      IF( lk_iomput .OR. ln_check_mass )  CALL can_chk_mass( kt ) ! Mass conservation checking

      IF ( lwm .AND. kt == nittrc000 ) CALL FLUSH    ( numonp )     ! flush output namelist CanOE
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_canoe')
      !
      !
   END SUBROUTINE trc_sms_canoe

   SUBROUTINE can_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  initan_sms_init  ***  
      !!
      !! ** Purpose :   read CanOE namelist
      !!
      !! ** input   :   file 'namelist.trc.s' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !NAMELIST/nampisbio/ nrdttrc, wsbio,  wsbio2,wsbioc
!      NAMELIST/nampisdmp/ ln_pisdmp, nn_pisdmp
      NAMELIST/nampismass/ ln_check_mass
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
!
!      REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : Pisces variables
!      READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
!901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )
!
!      REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : Pisces variables
!      READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
!902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
!      IF(lwm) WRITE ( numonp, nampisbio )
!
!      IF(lwp) THEN                         ! control print
!         WRITE(numout,*) ' Namelist : nampisbio'
!         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
!         WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
!         WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
!      ENDIF
!
!
!      REWIND( numnatp_ref )              ! Namelist nampisdmp in reference namelist : Pisces damping
!      READ  ( numnatp_ref, nampisdmp, IOSTAT = ios, ERR = 905)
!905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in reference namelist', lwp )
!
!      REWIND( numnatp_cfg )              ! Namelist nampisdmp in configuration namelist : Pisces damping
!      READ  ( numnatp_cfg, nampisdmp, IOSTAT = ios, ERR = 906 )
!906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in configuration namelist', lwp )
!      IF(lwm) WRITE ( numonp, nampisdmp )
!
!      IF(lwp) THEN                         ! control print
!         WRITE(numout,*)
!         WRITE(numout,*) ' Namelist : nampisdmp'
!         WRITE(numout,*) '    Relaxation of tracer to glodap mean value             ln_pisdmp      =', ln_pisdmp
!         WRITE(numout,*) '    Frequency of Relaxation                               nn_pisdmp      =', nn_pisdmp
!         WRITE(numout,*) ' '
!      ENDIF
!
      REWIND( numnatp_ref )              ! Namelist nampismass in reference namelist : Pisces mass conservation check
      READ  ( numnatp_ref, nampismass, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismass in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismass in configuration namelist : Pisces mass conservation check 
      READ  ( numnatp_cfg, nampismass, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismass in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismass )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameter for mass conservation checking'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Flag to check mass conservation of NO3/Si/TALK ln_check_mass = ', ln_check_mass
      ENDIF

   END SUBROUTINE can_sms_init

   SUBROUTINE can_ph_ini
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE can_ini_ph  ***
      !!
      !!  ** Purpose : Initialization of chemical variables of the carbon cycle
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      ! Set PH from  total alkalinity, borat (???), akb3 (???) and ak23 (???)
      ! --------------------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmas   = tmask(ji,jj,jk)
               ztmas1  = 1. - tmask(ji,jj,jk)
               zcaralk = trb(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
               zco3    = ( zcaralk - trb(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
               zbicarb = ( 2. * trb(ji,jj,jk,jpdic) - zcaralk )
               hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
            END DO
         END DO
     END DO
     !
   END SUBROUTINE can_ph_ini

   SUBROUTINE can_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE can_rst  ***
      !!
      !!  ** Purpose : Read or write variables in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!---------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      IF( TRIM(cdrw) == 'READ' ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' can_rst : Read specific variables from canoe model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         ! 
         IF( iom_varid( numrtr, 'PH', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'PH' , hi(:,:,:)  )
         ELSE
!            hi(:,:,:) = 1.e-9 
            CALL can_ph_ini
         ENDIF
         !
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         IF( kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'can_rst : write canoe restart file  kt =', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'PH', hi(:,:,:) )
      ENDIF
      !
   END SUBROUTINE can_rst

   SUBROUTINE can_dmp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  can_dmp  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt ! time step
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      !
      REAL(wp) :: zdntrsum, zdnfsum, ztau, nsum
      !!---------------------------------------------------------------------


      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' can_dmp : Restoring of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA configuration (not 1D) !
         !                                                    ! --------------------------- !
         ! adjust NO3 according to difference between global total rates of denitrification and N2 fixation
         ! this adjustment must be multiplicative rather than additive to prevent negative concentrations

         nsum = 1. / glob_sum( trn(:,:,:,jpno3)*cvol(:,:,:)*0.0010008 )              ! inverse global total N in mol^-1 (1.0008 is an approximate correction for non-NO3 N)
         zdnfsum = glob_sum( zdnf(:,:,:) * cvol(:,:,:)  )                           ! global total in molN s^-1
         zdntrsum = glob_sum( denitr(:,:,:) * cvol(:,:,:)  )         
         ztau = 1.+(zdntrsum-zdnfsum)*nsum*FLOAT(nn_pisdmp)*rfact       
         !IF(lwp) WRITE(numout,*) '       Totals  : ', zdnfsum, zdntrsum, ztau, nsum

         !trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksum
         !IF(lwp) WRITE(numout,*) '       NO3N  mean : ', zno3sumn
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * ztau

        !
      ENDIF
        !
   END SUBROUTINE can_dmp


   SUBROUTINE can_chk_mass( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_chk_mass  ***
      !!
      !! ** Purpose :  Mass conservation check 
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      REAL(wp) :: zalkbudget, zno3budget
      CHARACTER(LEN=100)   ::   cltxt
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zvol
      INTEGER :: jk
      !!----------------------------------------------------------------------

      !
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN 
         IF( ln_check_mass .AND. lwp) THEN      !   Open budget file of NO3, ALK, Si, Fer
            CALL ctl_opn( numno3, 'no3.budget'  , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numalk, 'talk.budget', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         ENDIF
      ENDIF

      !
      IF(  ln_check_mass .AND. kt == nitend )   THEN
         !   Compute the budget of NO3, ALK, Si, Fer
         zno3budget = glob_sum( (   trn(:,:,:,jpno3) + trn(:,:,:,jpnh4)  &
            &                    + trn(:,:,:,jpphy) + trn(:,:,:,jpdia)  &
            &                    + trn(:,:,:,jpzoo) + trn(:,:,:,jpmes)  &
            &                    + trn(:,:,:,jppoc)                     &
            &                    + trn(:,:,:,jpgoc))  * cvol(:,:,:)  )
         !
         zalkbudget = glob_sum( (   trn(:,:,:,jpno3)                     &
            &                     + trn(:,:,:,jptal)                     &
            &                     + trn(:,:,:,jpcal) * 2.                ) * cvol(:,:,:)  )

         IF( lwp ) THEN
            WRITE(numno3,9000) kt,  zno3budget / areatot
            WRITE(numalk,9100) kt,  zalkbudget / areatot
         ENDIF
      ENDIF
      !


 9100  FORMAT(i10,e18.10)
 9000  FORMAT(i10,e18.10)

       !
   END SUBROUTINE can_chk_mass



   !!======================================================================
END MODULE trcsms_canoe 
