MODULE cmocsink
   !!======================================================================
   !!                         ***  MODULE cmocsink  ***
   !! TOP :  CMOC  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking, Oct 28th 2015 simplified sink scheme according to Jim's own modifications for CMOC2
   !!           CMOC1  !  2016-02  (N. Swart) Adds calcite sinking.
   !!           CMOC1  !  2017-08  (A.Holdsworth) updated to NEMO 3.6.
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   cmoc_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   cmoc_sink_init  :  Unitialisation of sinking speed parameters
   !!   cmoc_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp
   USE cmocsbc
   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_sink         ! called in p4zbio.F90
   PUBLIC   cmoc_sink_init    ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   INTEGER  :: iksed  = 10



   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE cmoc_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk, jit
      INTEGER  ::   iiter1
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zfpon         ! Calcite export flx at the bottom of the euphotic zone
      REAL(wp) ::   zwsmax,zrfact2
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zdeup, zideup ! Euphotic zone depth, inverse depth
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zcalbotflx    ! Calcite flux to sediments
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zcalflxexp    ! exponential decay of calcite flux with depth
      REAL(wp) ::   zcaldiv                                  ! divergence of the calcite flux
      REAL(wp), POINTER, DIMENSION(:) :: zdepw ! computation of depths between t-grid cells.

      INTEGER  ::   ik100
      
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sink')
      !CMOC AMH2017
      CALL wrk_alloc( jpi, jpj, zfpon, zcalbotflx, zdeup, zideup)
      CALL wrk_alloc( jpi, jpj, jpk, zcalflxexp )
      CALL wrk_alloc( jpk, zdepw)
      !

      !CMOC AMH2017    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio  !set in namelist 
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8 * fse3t(ji,jj,jk) / xstep !dz/dt 
               wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO
      

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      zfpon(:,:) = 0._wp

      !   Compute the sedimentation term using cmocsink2 for  POC
      !   -----------------------------------------------------
      DO jit = 1, iiter1
        CALL cmoc_sink2( wsbio3, sinking , jppoc, iiter1 )
      END DO
      !     Calcite sinking flux
      !     --------------------------------------------------------------------
      ! Define an open ocean mask, based on where mbkt > nk_bal_cmoc == 15 (or as define in namelist)
      oomask(:,:) = 0._wp
      WHERE ( mbkt(:,:) >= nk_bal_cmoc ) oomask = 1._wp

      ! The euphotic zone depth and inverse depth, used for computing averages.
      zdeup(:,:) = 0._wp
      DO jk =1, jk_eud_cmoc
       zdeup(:,:) = zdeup(:,:) +  fse3t(:,:,jk)
      ENDDO
      zideup(:,:) = 1._wp / zdeup(:,:)
      ! Computing t-grid bounding depths more accurate then using w-grid depths, despite supposed identity.
      zdepw(:) = 0._wp
      DO jk = 2, jpk
         zdepw(jk) = zdepw(jk-1) + fse3t(ji,jj,jk-1)
      ENDDO

      !  Rain ratio at level jk_eud_cmoc - bottom of the euphotic zone:
      !  Temperature is however taken from the 1st layer (confirmed with RJC, 16/02/2016)
      xrcico(:,:) = rmcico_cmoc * exp(aci_cmoc * ( tsn(:,:,1,jp_tem)  - trcico_cmoc ) )                   &
      &                         / (1._wp + exp(aci_cmoc *( tsn(:,:,1,jp_tem) - trcico_cmoc ) + rtrn ) )

      ! PIC export at the bottom of the euphotic zone based on Zahariev et al 2008 p.59
      ! Time stepping is included with xstep, so units are in mol/m2/step
      zfpon(:,:) = xrcico(:,:) * wsbio3(:,:,jk_eud_cmoc) * xstep * trb(:,:,jk_eud_cmoc,jppoc)             &
      &                        *  tmask(:,:,jk_eud_cmoc) * oomask(:,:)

      ! Exponential decay of calcite flux with depth, at w-points.
      zcalflxexp(:,:,:) = 0._wp
      DO jk = jk_eud_cmoc+1, jpk                                                         
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcalflxexp(ji,jj,jk) = zfpon(ji,jj) * exp(-1._wp*(zdepw(jk)-zdeup(ji,jj)) / dci_cmoc)      &
               &                                   * tmask(ji,jj,jk-1)                ! Mask at jk-1 ensures bottom flux
            ENDDO                                                                     ! is included.
         ENDDO
      ENDDO
                
      ! Bottom PIC flux into sediments, which is removed from the deepest layer and
      ! added back to the surface (below).
      DO jj = 1, jpj
         DO ji = 1,jpi
             zcalbotflx(ji,jj) = zcalflxexp(ji,jj,mbkt(ji,jj)+1)

             ! Set the bottom boundary condition on the calcite flux to zero (no flux to sediment),
             ! and deal with the sediment flux separately from the sinking in the sections below.
               zcalflxexp(ji,jj,mbkt(ji,jj)+1) = 0._wp
         ENDDO
      ENDDO

      ! Over the levels of the euphotic zone, remove the euphotic-zone averaged
      ! PIC flux (mol/m3) from each level. 
      !the last terms in eq 1 and 2 on page 57
      DO jk =1, jk_eud_cmoc
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) -                                   &
         &                              zfpon(:,:) * zideup(:,:) 

         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) -                                   &
         &                      2._wp * zfpon(:,:) * zideup(:,:) 
      END DO

      ! Compute the divergence of the calcite flux and distribute it over the t-cell. 
      ! No sinking flux through the bottom here, given the masking above.
      DO jk = jk_eud_cmoc+1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcaldiv =  ( zcalflxexp(ji,jj,jk) - zcalflxexp(ji,jj,jk+1) ) / fse3t(ji,jj,jk) * tmask(ji,jj,jk)

               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +         zcaldiv 
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2._wp * zcaldiv                      
            ENDDO
         ENDDO
      ENDDO

      ! Do the bottom sedimentation of calcite. The sedimenting flux is added back
      ! to the surface layer (psuedo "river flux") for conservation.
      !AMH2017 removed the conservation equation if we treat rivers separately.
      DO jj = 1, jpj
         DO ji = 1,jpi
            tra(ji,jj,mbkt(ji,jj),jpdic) = tra(ji,jj,mbkt(ji,jj),jpdic) - zcalbotflx(ji,jj) / fse3t(ji,jj,mbkt(ji,jj))

            tra(ji,jj,mbkt(ji,jj),jptal) = tra(ji,jj,mbkt(ji,jj),jptal) - 2._wp * zcalbotflx(ji,jj) / fse3t(ji,jj,mbkt(ji,jj))
            
            !if ln_river then we don't need conservation
            IF( .NOT. ln_river ) THEN 
                tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic)  + zcalbotflx(ji,jj) / fse3t(ji,jj, 1) 
                tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal)  + 2._wp * zcalbotflx(ji,jj) / fse3t(ji,jj, 1) 
            ENDIF    
         ENDDO
      ENDDO
      !AMH 2017 update the trends

      !     --------------------------------------------------------------------

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         ik100  = iksed + 1
         IF( lk_iomput ) THEN
           IF( knt == nrdttrc ) THEN
              IF( iom_use( "oomask" ) )  THEN
                  !CALL iom_put( "EPC100"  , zw2d )
                  CALL iom_put( "oomask", oomask(:,:))
              ENDIF
              IF( iom_use( "EPC100" ) )  THEN
                  CALL iom_put( "EPC100", sinking(:,:,ik100) * zrfact2 * tmask(:,:,1) )
              ENDIF
              IF( iom_use("EPCALC100"))  THEN
                  CALL iom_put( "EPCALC100",    zfpon(:,:) * zrfact2 * tmask(:,:,1) ) !
              ENDIF
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUPOC"  , wsbio3(:,:,11) /rday * zbpoc(:,:) * 1e+3_wp  )  ! POC burial flux
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUCALC" , zfpon(:,:) * 1e+3_wp * rfact2r * zbpon(:,:)  )  ! <CMOC code OR 12/11/2015> *rfact2r replaces /rfact2 ! PIC burial flux
           ENDIF
         ELSE
           trc2d(:,:,jp_cm0_2d + 4) = sinking (:,:,ik100) * zrfact2 * tmask(:,:,1)

         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !

      CALL wrk_dealloc( jpi, jpj, zfpon, zcalbotflx, zdeup, zideup)
      CALL wrk_dealloc( jpi, jpj, jpk, zcalflxexp )
      CALL wrk_dealloc( jpk, zdepw)
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sink')

      !
   END SUBROUTINE cmoc_sink

   SUBROUTINE cmoc_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_sink_init  ***
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk
      INTEGER  ::   ios                 ! Local integer output status for namelist read

      NAMELIST/namcmocrr/  cnrr_cmoc, ncrr_cmoc
      NAMELIST/nampiscal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      NAMELIST/namcmocdeu/ jk_eud_cmoc, nk_bal_cmoc

      !
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist namcmocrr in reference namelist 
      READ  ( numnatp_ref,namcmocrr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocrr in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namcmocrr in configuration namelist 
      READ  ( numnatp_cfg, namcmocrr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocrr in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namcmocrr )
      
      REWIND( numnatp_ref )              ! Namelist nampiscal in reference namelist 
      READ  ( numnatp_ref,nampiscal, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscal in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiscal in configuration namelist 
      READ  ( numnatp_cfg, nampiscal, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscal in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiscal )

      REWIND( numnatp_ref )              ! Namelist namcmocdeu in reference namelist 
      READ  ( numnatp_ref,namcmocdeu, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocdeu in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namcmocdeu in configuration namelist 
      READ  ( numnatp_cfg, namcmocdeu, IOSTAT = ios, ERR = 906 )
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocdeu in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namcmocdeu )




      ! <CMOC code OR 10/22/2015> CMOC namelist end
      
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Redfield ratios, namcmocrr'    
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    C:N ratio                                 cnrr_cmoc =', cnrr_cmoc
         WRITE(numout,*) '    N:C ratio                                 ncrr_cmoc =', ncrr_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for calcite export  , namcmoccal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum rain ratio                      rmcico_cmoc =', rmcico_cmoc
         WRITE(numout,*) '    Rain ratio half-point temperature       trcico_cmoc =', trcico_cmoc
         WRITE(numout,*) '    Rain ratio scaling factor                  aci_cmoc =',    aci_cmoc
         WRITE(numout,*) '    CaCO3 redissolution depth scale            dci_cmoc =',    dci_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Euphotic zone depth , namcmocdeu'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Index of the bottom of the euphotic zone jk_eud_cmoc =', jk_eud_cmoc
         WRITE(numout,*) '    Open ocean, min. number of vert. layers  nk_bal_cmoc =', nk_bal_cmoc
      ENDIF

   END SUBROUTINE cmoc_sink_init


   SUBROUTINE cmoc_sink2( pwsink, psinkflx, jp_tra, kiter )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      INTEGER , INTENT(in   )                         ::   kiter     ! number of iterations for time-splitting 
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zew, zflx 
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sink2')
      !
      ! Allocate temporary workspace
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      CALL wrk_alloc( jpi, jpj, jpk,   zwsink2, ztrb )

      ztrb (:,:,:) = trb(:,:,:,jp_tra)
      !convert to m/s and apply the t land mask
      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0


      ! Vertical advective flux
      DO jk = 1, jpkm1
        DO jj = 1, jpj      
            DO ji = 1, jpi    
              ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
              zew   = zwsink2(ji,jj,jk+1)
              ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
              psinkflx(ji,jj,jk+1) = -zew * trb(ji,jj,jk,jp_tra) * rfact2
            END DO
        END DO
      END DO
      !
      ! Boundary conditions
      psinkflx(:,:,1  ) = 0.e0
      psinkflx(:,:,jpk) = 0.e0
      
      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
                ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + zflx
            END DO
         END DO
      END DO
      !update tracer concentrations based on advection
      trb(:,:,:,jp_tra) = ztrb(:,:,:)
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme      
      ! psinkflx(:,:,:)        = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwsink2, ztrb )      
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sink2')
   END SUBROUTINE cmoc_sink2


   INTEGER FUNCTION cmoc_sink_alloc()
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) ,      &
         &      sinking(jpi,jpj,jpk) ,      &                
         &       STAT=cmoc_sink_alloc )                
         !
      IF( cmoc_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION cmoc_sink_alloc
   
#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sink                    ! Empty routine
   END SUBROUTINE cmoc_sink
#endif 

   !!======================================================================
END MODULE  cmocsink
