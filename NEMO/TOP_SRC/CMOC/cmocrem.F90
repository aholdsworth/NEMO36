MODULE cmocrem
   !!======================================================================
   !!                         ***  MODULE cmocrem  ***
   !! TOP :   CMOC Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_rem       :  Compute remineralization/dissolution of organic compounds
   !!   cmoc_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   !USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE cmocche          !  chemical model
   USE cmocprod         !  Growth rate of the 2 phyto groups
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE cmocsink         !


   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_rem         ! called in p4zbio.F90
   PUBLIC   cmoc_rem_init    ! called in trcsms_cmoc.F90

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cmoc_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !!CMOC AMH2017 ** Method  : Temperature dependent remineralization based on
      !!              Zahariev et al. (2008)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_rem')
      !
      ! Initialization of CMOC arrays
       redet   (:,:,:) = 0._wp
       redettot(:,:)   = 0._wp

      ! Remineralisation rate of detritus
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               redet (ji,jj,jk) = reref_cmoc * xstep &
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp *        &
               &                ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp  &
               &                 + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp )  &
               &                )      ) * trb(ji,jj,jk,jppoc) * tmask(ji,jj,jk)
            END DO
          END DO 
      END DO      

      ! Integration of remineralization below the euphotic zone (used for dentrification scaling)
      DO jk = jk_eud_cmoc+1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
                redettot(ji,jj) = redettot(ji,jj) + redet(ji,jj,jk)         &
                &                                 * fse3t(ji,jj,jk)         &
                &                                 * tmask(ji,jj,jk)
            END DO
          END DO 
      END DO      

      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      DO jk = 1, jpkm1
         tra(:,:,jk,jppoc) = tra(:,:,jk,jppoc) - redet(:,:,jk) 
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) + redet(:,:,jk) 
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - redet(:,:,jk) 
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + redet(:,:,jk) 
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) - redet(:,:,jk) * ncrr_cmoc
      END DO


       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_rem')
      !
   END SUBROUTINE cmoc_rem


   SUBROUTINE cmoc_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      ! <CMOC code OR 10/15/2015> CMOC namelist
      NAMELIST/nampisrem/ ed_cmoc, reref_cmoc

      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisrem )

      ! <CMOC code OR 10/15/2015> CMOC namelist end 
      !!----------------------------------------------------------------------
      
      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocpoc'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Remineralisation rate of POC              reref_cmoc=', reref_cmoc
         WRITE(numout,*) '    Activation energy for remineralization    ed_cmoc   =', ed_cmoc   
         WRITE(numout,*) ' '
      ENDIF
      !
      !
   END SUBROUTINE cmoc_rem_init



#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_rem                    ! Empty routine
   END SUBROUTINE cmoc_rem
#endif 

   !!======================================================================
END MODULE cmocrem
