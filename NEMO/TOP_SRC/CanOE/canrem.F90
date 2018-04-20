MODULE canrem
   !!======================================================================
   !!                         ***  MODULE canrem  ***
   !! TOP :   CanOE Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
    !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------
#if defined key_canoe
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_rem       :  Compute remineralization/dissolution of organic compounds
   !!   can_rem_init  :  Initialisation of parameters for remineralisation
   !!   can_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE canopt          !  optical model
   USE canche          !  chemical model
   USE canprod         !  Growth rate of the 2 phyto groups
   USE canmeso         !  Sources and sinks of mesozooplankton
   USE canint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_rem         ! called in canbio.F90
   PUBLIC   can_rem_init    ! called in trcsms_canoe.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  xremik    = 0.25_wp    !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xremip    = 0.025_wp   !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::  nitrif    = 0.05_wp    !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::  xlam1     = 0.0001_wp  !: scavenging rate of iron (low concentrations)
   REAL(wp), PUBLIC ::  xlam2     = 2.5_wp     !: scavenging rate of iron (high concentrations)
   REAL(wp), PUBLIC ::  ligand    = 6.0E+2_wp  !: ligand concentration
   REAL(wp), PUBLIC ::  pocfctr   = 0.65574_wp !: multiplier for POC-dependent scavenging
   REAL(wp), PUBLIC ::  o2thresh  = 6._wp      !: O2 threshold for denitrification
   REAL(wp), PUBLIC ::  nh4frx    = 02.5_wp    !: annamox fraction of denitrification
   REAL(wp), PUBLIC ::  oxymin    = 1._wp      !: half saturation constant for anoxia 
   REAL(wp), PUBLIC ::  nyld      = 0.8_wp     !: denitrification stoichiometric coefficient


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik, Tf
      REAL(wp) ::   zkeq, zfeequi
      REAL(wp) ::   zsatur, zsatur2, zdep, zfactdep
      REAL(wp) ::   zorem, zorem2, zofer, zofer2
      REAL(wp) ::   zscave, zscavex, fexs, zcoag
      REAL(wp) ::   zlamfac, zonitr, zstep
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_rem')
      !

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.  - trb(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
            END DO
         END DO
      END DO

      nh4ox(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
               !    NH4 nitrification to NO3. Ceased for oxygen concentrations
               !    below 2 umol/L. Inhibited at strong light 
               !    ----------------------------------------------------------
               zonitr  =nitrif * zstep * trb(ji,jj,jk,jpnh4) / ( 1.+ emoy(ji,jj,jk) ) * ( 1.- nitrfac(ji,jj,jk) ) 
               !denitnh4(ji,jj,jk) = nitrif * zstep * trn(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk) 
               !   Update of the tracers trends
               !   ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - 2. * zonitr
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2.e-6 * zonitr
               nh4ox(ji,jj,jk) = zonitr
            END DO
         END DO
      END DO


       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF


      denitr(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               Tf = tgfuncr(ji,jj,jk)
               zorem  = xremik * xstep * Tf * trb(ji,jj,jk,jppoc)
               zofer  = zorem * rr_fe2c
               zorem2 = xremik * xstep * Tf * trb(ji,jj,jk,jpgoc)
               zofer2 = zorem2 * rr_fe2c

! denitrification is assumed to remove NO3 as a fraction of remineralization increasing linearly from 0 to 1 with declining [O2] for [O2]<10 uM
! NO3 fraction is then divided 0.75/0.25 between NO3 and NH4 (50% classical denitrification and 50% annamox)
               zonitr=1.-MIN(trb(ji,jj,jk,jpoxy),o2thresh)/o2thresh
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zorem + zorem2)*rr_n2c - (zorem + zorem2)*nyld*zonitr*0.5*nh4frx
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - (zorem + zorem2)*nyld*zonitr*(1.-0.5*nh4frx)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - (zorem + zorem2)*(1.-zonitr)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zofer + zofer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (zorem + zorem2)*1.e-6
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zorem
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zorem2
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 1.e-6 * (zorem + zorem2)*rr_n2c                         ! 1 mol of alkalinity per mol of N
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 1.e-6 * (zorem + zorem2)*nyld*zonitr*(1.-nh4frx)        ! +1 mol if denitrification, 0 if annamox
               denitr(ji,jj,jk) = (zorem + zorem2)*zonitr*nyld

            END DO
         END DO
      END DO
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('rem3')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji = 1, jpi
               zstep   = xstep
! irreversible scavenging as in Christian et al 2002
               zcoag = MIN((trb(ji,jj,jk,jppoc)+trb(ji,jj,jk,jpgoc))*pocfctr,1.)
               fexs = MAX(trb(ji,jj,jk,jpfer)-ligand,0.)
               zscave = xlam1 * xstep * (trb(ji,jj,jk,jpfer)-fexs) * zcoag
               zscavex = xlam2 * xstep * fexs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - (zscave+zscavex)
            END DO
         END DO
      END DO
      !

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem5')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !     --------------------------------------------------------------------

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e-3 * rfact2r  ! conversion from umol/L/timestep into mol/m3/s
         denitr(:,:,:) = denitr(:,:,:) * zrfact2
         nh4ox(:,:,:) = nh4ox(:,:,:) * zrfact2
         IF( knt == nrdttrc ) THEN
          CALL iom_put( "Denitr"   , denitr(:,:,:) * tmask(:,:,:) )  ! rate of denitrification
          CALL iom_put( "Nitrif"   , nh4ox(:,:,:) * tmask(:,:,:) )  ! rate of nitrification
         ENDIF
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_rem')
      !
   END SUBROUTINE can_rem


   SUBROUTINE can_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, xremip, nitrif, xlam1, xlam2, ligand, pocfctr, o2thresh, &
                        & nh4frx, oxymin
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, nampisrem'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    remineralisation rate of POC(1)           xremik    =', xremik
         WRITE(numout,*) '    remineralisation rate of POC(2)           xremip    =', xremip
         WRITE(numout,*) '    Nitrification maximum rate                nitrif    =', nitrif
         WRITE(numout,*) '    Fe scavenging rate (low concentrations)   xlam1     =', xlam1
         WRITE(numout,*) '    Fe scavenging rate (high concentrations)  xlam2     =', xlam2
         WRITE(numout,*) '    Fe-binding ligand concentration           ligand    =', ligand
         WRITE(numout,*) '    POC-dependence of Fe scavenging           pocfctr   =', pocfctr
         WRITE(numout,*) '    O2 threshold for denitrification          o2thresh  =', o2thresh
         WRITE(numout,*) '    Annamox fraction of denitrification       nh4frx    =', nh4frx
         WRITE(numout,*) '    O2 dependence of nitrification            oxymin    =', oxymin
      ENDIF
      !
      !
      nitrfac (:,:,:) = 0._wp
      denitr  (:,:,:) = 0._wp
      nh4ox   (:,:,:) = 0._wp
      !
   END SUBROUTINE can_rem_init



#else
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_rem                    ! Empty routine
   END SUBROUTINE can_rem
#endif 

   !!======================================================================
END MODULE canrem
