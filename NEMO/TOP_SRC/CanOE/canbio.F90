MODULE canbio
   !!======================================================================
   !!                         ***  MODULE canbio  ***
   !! TOP :   CanOE bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_canoe
   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_bio        :   computes the interactions between the different
   !!                      compartments of CanOE
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canopt          !  optical model
   USE canprod         !  Growth rate of the 2 phyto groups
   USE canmort         !  Mortality terms for phytoplankton
   USE canmicro        !  Sources and sinks of microzooplankton
   USE canmeso         !  Sources and sinks of mesozooplankton
   USE canrem          !  Remineralisation of organic matter
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  can_bio    

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canbio.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE can_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of CanOE
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_bio')
      !
      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

          
      CALL can_sink ( kt, knt )     ! vertical flux of particulate organic matter
      CALL can_opt  ( kt, knt )     ! Optic: PAR in the water column
      CALL can_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
      !                             ! (for each element : C, Si, Fe, Chl )
      CALL can_rem  ( kt, knt )     ! remineralization terms of organic matter+scavenging of Fe
      CALL can_mort ( kt      )     ! phytoplankton mortality
     !                             ! zooplankton sources/sinks routines 
      CALL can_micro( kt, knt )           ! microzooplankton
      CALL can_meso ( kt, knt )           ! mesozooplankton
      !                             ! test if tracers concentrations fall below 0.
      !                                                             !
      !CALL total_element(totfe,totn)
      !WRITE(numout,*) totfe, totn
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_bio')
      !
   END SUBROUTINE can_bio

   SUBROUTINE total_element(totfe,totn)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE total_element  ***
      !!
      !! ** Purpose :   calculate model total N and Fe   2015/03/23 JRC
      !!---------------------------------------------------------------------
      REAL(wp) :: totfe, totn
      !!---------------------------------------------------------------------

      totfe=0.
      totn=0.

            totn = glob_sum( (   (trb(:,:,:,jpno3) + trb(:,:,:,jpnh4)  &
            &                     + trb(:,:,:,jpnn) + trb(:,:,:,jpdn))*6.625  &
            &                     + trb(:,:,:,jpzoo) + trb(:,:,:,jpmes)  &
            &                     + trb(:,:,:,jppoc) + trb(:,:,:,jpgoc)  ) * cvol(:,:,:)  )
            totn=totn/6.625
 
! zoo and poc concs are in C units, 4.981=33/6.625

            totfe = glob_sum( (   trb(:,:,:,jpfer) + trb(:,:,:,jpdfe) + trb(:,:,:,jpnfe)   &
            &                     + trb(:,:,:,jpzoo)*4.981 + trb(:,:,:,jpmes)*4.981  &
            &                     + trb(:,:,:,jppoc)*4.981 + trb(:,:,:,jpgoc)*4.981  ) * cvol(:,:,:)  )

      !
   END SUBROUTINE total_element

#else
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_bio                         ! Empty routine
   END SUBROUTINE can_bio
#endif 

   !!======================================================================
END MODULE  canbio

