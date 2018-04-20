MODULE cmocbio
   !!======================================================================
   !!                         ***  MODULE cmocbio  ***
   !! TOP :   CMOC bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!          CMOC 1  !  2013-2015(O. Riche) calls to ecosystem procedures and some extra diagnostics
   !!          CMOC 1  !  2017(A.Holdsworth) updated for 3.6
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_bio        :   computes the interactions between the different
   !!                      compartments of CMOC
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE cmocsink         !  vertical flux of particulate matter due to sinking
   USE cmocprod         !  Growth rate of the 2 phyto groups
   USE cmocmort         !  Mortality terms for phytoplankton
   USE cmocmicro        !  Sources and sinks of microzooplankton
   USE cmocrem          !  Remineralisation of organic matter
!   USE cmocfechem
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  cmoc_bio    

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocbio.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE cmoc_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of CMOC
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('cmoc_bio')
      !trim for CMOC AMH 2017

      CALL cmoc_sink ( kt, knt )     ! vertical flux of particulate organic matter
      !CALL cmoc_fechem(kt, knt )     ! Iron chemistry/scavenging
      CALL cmoc_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
      CALL cmoc_rem  ( kt, knt )     ! remineralization terms of organic matter+scavenging of Fe
      !                             ! (for each element : C, Si, Fe, Chl )
      CALL cmoc_mort ( kt      )     ! phytoplankton mortality
     !                             ! zooplankton sources/sinks routines 
      CALL cmoc_micro( kt, knt )           ! microzooplankton
      !                             ! test if tracers concentrations fall below 0.
      !                                                             !
      ! <CMOC code OR AMH 2017> Diagnostics of Normalized Alkalinity and DIC
       IF( ln_diatrc ) THEN                      
         IF( lk_iomput ) THEN
          IF( knt == nrdttrc ) THEN             

               CALL iom_put( "sDIC"   , trn(:,:,:,jpdic) / tsn(:,:,:,jp_sal) * 35 * 1e+3_wp * tmask(:,:,:) )
               CALL iom_put( "sAlk"   , trn(:,:,:,jptal) / tsn(:,:,:,jp_sal) * 35 * 1e+3_wp * tmask(:,:,:) )

          ENDIF
         ENDIF
       ENDIF 

      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_bio')
      !
   END SUBROUTINE cmoc_bio

#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_bio                         ! Empty routine
   END SUBROUTINE cmoc_bio
#endif 

   !!======================================================================
END MODULE  cmocbio

