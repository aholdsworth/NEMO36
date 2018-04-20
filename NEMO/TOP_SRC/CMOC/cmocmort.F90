MODULE cmocmort
   !!======================================================================
   !!                         ***  MODULE cmocmort  ***
   !! TOP :   CMOC Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!          CMOC 1  !  2013-2015(O. Riche) phytoplankton mortality, based on Zahariev et al 2008
   !!          CMOC 1  !  2013-2015(A.Holdsworth) updated for 3.6
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_mort       :   Compute the mortality terms for phytoplankton
   !!   cmoc_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE cmocprod         !  Primary productivity 
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_mort    
   PUBLIC   cmoc_mort_init    



   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocmort.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE cmoc_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------


      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zsizerat, zcompaph
      REAL(wp) :: zfactch
      REAL(wp) :: ztortp , zrespp , zmortp 
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_mort')
      !
      ! <CMOC code OR 10/19/2015> Note: replacing xstep (and removing facvol/key_grad instances) by xstep time step in days AMH2017
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Note: use the PISCES practice of minimum phyto biomass (see zcompaph below)
               zcompaph = MAX( ( trb(ji,jj,jk,jpphy) - 1e-8 ), 0.e0 )
               ! Quadratic mortality
               ! <CMOC code OR 10/19/2015> 1e+3_wp convert L^-1 to m^-3 (); use xstep the global constant to convert to d^-1
               zrespp = mpd2_cmoc * ncrr_cmoc * 1.0e3_wp * xstep * zcompaph * trb(ji,jj,jk,jpphy)

               !  Linear mortality
               ztortp = mpd_cmoc * xstep * zcompaph

               !  Total mortality
               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks
               !   Calculate the chlorophyll to phytoplankton ratio
               zfactch = trb(ji,jj,jk,jpnch)/(trb(ji,jj,jk,jpphy)+rtrn)

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('mort')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_mort')
      !
   END SUBROUTINE cmoc_mort

      SUBROUTINE cmoc_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/19/2015> CMOC namelist
      NAMELIST/nampismort/ mpd_cmoc, mpd2_cmoc
      ! <CMOC code OR 10/19/2015> CMOC namelist end 
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )
      READ  (numnatp_ref, nampismort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocmor in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : cmoc phytoplankton
      READ  ( numnatp_cfg, nampismort, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcmocmor in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namcmocmor'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Phytoplankton mortality to detritus       mpd_cmoc  =', mpd_cmoc
         WRITE(numout,*) '    Phytoplankton quadratic mortality         mpd2_cmoc =', mpd2_cmoc
      ENDIF

   END SUBROUTINE cmoc_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_mort                    ! Empty routine
   END SUBROUTINE cmoc_mort
#endif 

   !!======================================================================
END MODULE  cmocmort
