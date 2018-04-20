MODULE cmocsms
   !!======================================================================
   !!                         ***  MODULE cmocsms  ***
   !! TOP :   CMOC Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!           CMOC1  !  2017-08 (A. Holdsworth) adapted CMOC for NEMO3.6
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_cmoc         :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE par_cmoc
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE trcdta
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE cmocbio          !  Biological model
   USE cmocche          !  Chemical model
   USE cmocflx          !  Gas exchange
   USE cmocsbc          !  External source of nutrients
   USE cmocsed          !  Sedimentation
   USE cmocrem          !  remineralisation
   USE iom             !  I/O manager
   USE trd_oce         !  Ocean trends variables
   USE trdtrc          !  TOP trends variables
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_sms        ! called in trcsms_cmoc
   PUBLIC   cmoc_sms_init   ! called in p4zsms.F90

   REAL(wp) :: alkbudget, no3budget 
   REAL(wp) :: xfact1, xfact2, xfact3
   INTEGER ::  numco2, numnut, numnit  !: logical unit for co2 budget

   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: ???


   !! * Substitutions
#  include "top_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_cmoc.F90 4147 2013-11-04 11:51:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_cmoc ***
      !!
      !! ** Purpose :   Initialisation of the CMOC biochemical model
      !!----------------------------------------------------------------------


   SUBROUTINE cmoc_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sms  ***
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of CMOC bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   ji, jj, jk, jnt, jn, jl
      REAL(wp) ::  ztra
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sms')
      !The first passive tracer time step
      IF( kt == nittrc000 ) THEN
        !
        ALLOCATE( xnegtr(jpi,jpj,jpk) )
        !
        CALL cmoc_che                              ! initialize the chemical constants
        !
        IF( .NOT. ln_rsttr ) THEN  ;   CALL cmoc_ph_ini   !  set PH at kt=nit000 
        ELSE                       ;   CALL cmoc_rst( nittrc000, 'READ' )  !* read or initialize all required fields 
        ENDIF
        !
      ENDIF
      !
      IF( ln_pisdmp .AND. MOD( kt - nn_dttrc, nn_pisdmp ) == 0 )   CALL cmoc_dmp( kt )      ! Relaxation of some tracers
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
         IF(lwp) write(numout,*) '    CMOC  Biology time step    rfact2 = ', rfact2
         IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN
         DO jn = jp_cm0, jp_cm1              !   SMS on tracer without Asselin time-filter
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

         CALL cmoc_che              ! computation of chemical constants
         !
      ENDIF

      IF( ll_sbc ) CALL cmoc_sbc( kt )   ! external sources of nutrients 

      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL cmoc_bio( kt, jnt )   ! Biology soft tissue production
         CALL cmoc_sed( kt, jnt )   ! Sedimentation - soft tissue remineralization
         CALL cmoc_flx( kt, jnt )   ! Compute surface fluxes
         !
         xnegtr(:,:,:) = 1.e0
         DO jn = jp_cm0, jp_cm1
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
         DO jn = jp_cm0, jp_cm1
           trb(:,:,:,jn) = trb(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
         END DO
        !
         DO jn = jp_cm0, jp_cm1
            tra(:,:,:,jn) = 0._wp
         END DO
         !
         IF( ln_top_euler ) THEN
            DO jn = jp_cm0, jp_cm1
               trn(:,:,:,jn) = trb(:,:,:,jn)
            END DO
         ENDIF
      END DO

      !
      !
      IF( l_trdtrc ) THEN
         DO jn = jp_cm0, jp_cm1
           CALL trd_trc( tra(:,:,:,jn), jn, jptra_sms, kt )   ! save trends
         END DO
      END IF
      !
      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !set lateral boundary conditions
         DO jn = jp_cm0, jp_cm1
           CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF
      !
      IF( lrst_trc )  CALL cmoc_rst( kt, 'WRITE' )  !* Write CMOC informations in restart file 
      !


      IF ( lwm .AND. kt == nittrc000 ) CALL FLUSH    ( numonp )     ! flush output namelist CMOC
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sms')
      !
      !

      !
   END SUBROUTINE cmoc_sms
   
   SUBROUTINE cmoc_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  cmoc_sms_init  ***  
      !!
      !! ** Purpose :   read CMOC namelist
      !!
      !! ** input   :   file 'namelist.trc.s' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      NAMELIST/nampisbio/ nrdttrc,  xkmort,  wsbio2, wsbio
      NAMELIST/nampisdmp/ ln_pisdmp, nn_pisdmp
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisbio in reference namelist : Pisces variables
      READ  ( numnatp_ref, nampisbio, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisbio in configuration namelist : Pisces variables
      READ  ( numnatp_cfg, nampisbio, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisbio in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisbio )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' Namelist : nampisbio'
         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
         WRITE(numout,*) '    half saturation constant for mortality    xkmort    =', xkmort
         WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
         WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
      ENDIF
      



      REWIND( numnatp_ref )              ! Namelist nampisdmp in reference namelist : Pisces damping
      READ  ( numnatp_ref, nampisdmp, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisdmp in configuration namelist : Pisces damping
      READ  ( numnatp_cfg, nampisdmp, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisdmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisdmp )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampisdmp'
         WRITE(numout,*) '    Relaxation of tracer to glodap mean value             ln_pisdmp      =', ln_pisdmp
         WRITE(numout,*) '    Frequency of Relaxation                               nn_pisdmp      =', nn_pisdmp
         WRITE(numout,*) ' '
      ENDIF


   END SUBROUTINE cmoc_sms_init

   SUBROUTINE cmoc_ph_ini
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cmoc_ini_ph  ***
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
   END SUBROUTINE cmoc_ph_ini

   SUBROUTINE cmoc_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cmoc_rst  ***
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
      !!---------------------------------------------------------------------

      IF( TRIM(cdrw) == 'READ' ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' cmoc_rst : Read specific variables from cmoc model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         ! 
         IF( iom_varid( numrtr, 'PH', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'PH' , hi(:,:,:)  )
         ELSE
!            hi(:,:,:) = 1.e-9 
            CALL cmoc_ph_ini
         ENDIF
         !
         IF( iom_varid( numrtr, 'tcflxcum', ldstop = .FALSE. ) > 0 ) THEN  ! cumulative total flux of carbon
            CALL iom_get( numrtr, 'tcflxcum' , t_oce_co2_flx_cum  )
         ELSE
            t_oce_co2_flx_cum = 0._wp
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         IF( kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'cmoc_rst : write cmoc restart file  kt =', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'PH', hi(:,:,:) )
      !   CALL iom_rstput( kt, nitrst, numrtw, 'tcflxcum', t_oce_co2_flx_cum )
      ENDIF
      !
   END SUBROUTINE cmoc_rst

   SUBROUTINE cmoc_dmp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  cmoc_dmp  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt ! time step
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.165     ! mean value of phosphates
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      REAL(wp) ::  silmean = 91.51     ! mean value of silicate
      !
      REAL(wp) :: zarea, zalksumn, zno3sumn 
      REAL(wp) :: zalksumb,  zno3sumb 
      !!---------------------------------------------------------------------


      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' cmoc_dmp : Restoring of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA configuration (not 1D) !
         !                                                    ! --------------------------- !
         ! set total alkalinity, phosphate, nitrate & silicate
         zarea          = 1._wp / glob_sum( cvol(:,:,:) ) * 1e6              

         zalksumn = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
         zno3sumn = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea * rno3
 
         IF(lwp) WRITE(numout,*) '       TALKN mean : ', zalksumn
         trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksumn


         IF(lwp) WRITE(numout,*) '       NO3N  mean : ', zno3sumn
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * no3mean / zno3sumn

         !
         !
         IF( .NOT. ln_top_euler ) THEN
            zalksumb = glob_sum( trb(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
            zno3sumb = glob_sum( trb(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea * rno3
 
            IF(lwp) WRITE(numout,*) ' '
            IF(lwp) WRITE(numout,*) '       TALKB mean : ', zalksumb
            trb(:,:,:,jptal) = trb(:,:,:,jptal) * alkmean / zalksumb


            IF(lwp) WRITE(numout,*) '       NO3B  mean : ', zno3sumb
            trb(:,:,:,jpno3) = trb(:,:,:,jpno3) * no3mean / zno3sumb

        ENDIF
        !
      ENDIF
        !
   END SUBROUTINE cmoc_dmp


#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'cmoc_sms: You should not have seen this print! error?', kt
   END SUBROUTINE cmoc_sms
#endif 

   !!======================================================================
END MODULE cmocsms 
