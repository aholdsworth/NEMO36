MODULE cmocsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   CMOC Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
#if defined key_cmoc
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_cmoc      !  CMOC Source Minus Sink variables
   USE cmocsink         !  vertical flux of particulate matter due to sinking
   USE cmocsbc          !  External source of nutrients 
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_sed  
   PUBLIC   cmoc_sed_alloc
 

   !! * Module variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocsed.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cmoc_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  ::   ji, jj, jk, ikt
      REAL(wp) ::    rivconv 
          
!
      REAL(wp) ::   zwsbio3, zdep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zirondep, zsoufer
      !!---------------------------------------------------------------------
      !
      rivconv=1./2.631 
      IF( nn_timing == 1 )  CALL timing_start('cmoc_sed')
      !
      IF( kt == nittrc000 .AND. knt == 1 )   r1_rday  = 1. / rday
      !


     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
          tra(ji,jj,jk,jpno3) = tra(ji,jj,1,jpno3) +  (rivconv*rivdic(ji,jj)) * rfact2
          tra(ji,jj,jk,jpdic) = tra(ji,jj,1,jpdic) +  rivdic(ji,jj) * rfact2
          tra(ji,jj,jk,jptal) = tra(ji,jj,1,jptal) +  ( rivalk(ji,jj) - ncrr_cmoc *rivconv*rivdic(ji,jj))   * rfact2
      
      ENDIF
      !
      ! Fate of POC reaching the ocean floor: complete remineralization into DIC, DIN
      ! and sink of O2 and TALK
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwsbio3 = wsbio3(ji,jj,ikt) * zdep

            tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic)                       &
               &                             + tra(ji,jj,ikt,jppoc) * zwsbio3 
            tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal)                       &
               &                             - tra(ji,jj,ikt,jppoc) * zwsbio3 * ncrr_cmoc
            tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3)                       &
               &                             + tra(ji,jj,ikt,jppoc) * zwsbio3 
            tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy)                       &
               &                             - tra(ji,jj,ikt,jppoc) * zwsbio3 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc)                       &
                                             - tra(ji,jj,ikt,jppoc) * zwsbio3 
         END DO
      END DO

      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('cmoc_sed')
      !
 9100  FORMAT(i8,3f10.5)
      !
   END SUBROUTINE cmoc_sed


   INTEGER FUNCTION cmoc_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(jpi,jpj,jpk), sdenit(jpi,jpj), STAT=cmoc_sed_alloc )
      !
      IF( cmoc_sed_alloc /= 0 )   CALL ctl_warn('p4z_sed_alloc: failed to allocate arrays')
      !
   END FUNCTION cmoc_sed_alloc


#else
   !!======================================================================
   !!  Dummy module :                                   No CMOC bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE cmoc_sed                         ! Empty routine
   END SUBROUTINE cmoc_sed
#endif 

   !!======================================================================
END MODULE  cmocsed
