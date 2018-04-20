MODULE canmort
   !!======================================================================
   !!                         ***  MODULE canmort  ***
   !! TOP :   CanOE Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_mort       :   Compute the mortality terms for phytoplankton
   !!   can_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canprod         !  Primary productivity 
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_mort    
   PUBLIC   can_mort_init    

   !! * Shared module variables

   REAL(wp), PUBLIC :: mprat   = 5.E-2_wp   !: phytoplankton mortality rate 
   REAL(wp), PUBLIC :: mprat2  = 2.E-1_wp   !: Diatoms mortality rate
   REAL(wp), PUBLIC :: mpratm  = 5.E-2_wp   !: Phytoplankton minimum mortality rate
   REAL(wp), PUBLIC :: mpqua   = 1.E-09_wp  !: quadratic mortality of phytoplankton
   REAL(wp), PUBLIC :: mpquad  = 2.E-08_wp  !: maximum quadratic mortality of diatoms
   REAL(wp), PUBLIC :: chldegr = 2.E-2_wp   !: Chlorophyll photooxidation rate
   REAL(wp), PUBLIC :: picfrx  = 1.E-1_wp   !: CaCO3 fraction of mortality (0.1 implies 1 mol caCO3 for each 10 mol POC)
   REAL(wp), PUBLIC :: xminp   = 0.01       !: minimum phytoplankton concentration for linear mortality


   !!* Substitution
   !!----------------------------------------------------------------------
   !!                    ***  top_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose : Statement function file: to be include in all passive tracer modules
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec) new architecture
   !!----------------------------------------------------------------------
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
!   Default option :                         eiv: dummy variables
   !!----------------------------------------------------------------------
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 6312 2016-02-15 11:43:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c3d' :                 aht: 3D coefficient
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
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: top_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canmort.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE can_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      CALL can_nano            ! nanophytoplankton

      CALL can_diat            ! diatoms

   END SUBROUTINE can_mort


   SUBROUTINE can_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zmortp,zmortz
      REAL(wp) :: spc,spn,spf,szc,chl
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n,thetac
      REAL(wp) :: cxs,nxs1,nxs2,fexs1,fexs2
      REAL(wp) :: csw1,csw2
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_nano')
      !
      prodcal(:,:,:) = 0.  !: calcite production variable set to zero
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               spc = MAX(trb(ji,jj,jk,jpphy),0.)
               spn = MAX(trb(ji,jj,jk,jpnn),0.)
               spf = MAX(trb(ji,jj,jk,jpnfe),0.)
               szc = MAX(trb(ji,jj,jk,jpzoo),0.)
               chl = MAX(trb(ji,jj,jk,jpnch),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)
               thetac=chl/(spc+rtrn)

! simplified CMOC type mortality: sum of linear and quadratic terms
               zmortp = mprat * xstep * spc + mpqua * xstep * spc * spc
               if (spc.le.xminp) zmortp = mpqua * xstep * spc * spc           ! no linear mortality below biomass threshold xminp
               zmortz = mprat * xstep * szc + mpqua * xstep * szc * szc
               if (szc.le.xminp) zmortz = mpqua * xstep * szc * szc
! reduce mortality to what can support detritus production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the detritus ratio
               zmortp=zmortp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
               zmortpn(ji,jj,jk) = zmortp
! calculate "excess" relative to grazer RR
               cxs=zmortp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=zmortp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=zmortp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=zmortp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=zmortp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)

               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2

               !   Update the arrays TRA which contains the biological sources and sinks

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp - cxs
               tra(ji,jj,jk,jpnn) = tra(ji,jj,jk,jpnn) - zmortp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zmortp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - (zmortp+cxs) * thetac - chl*chldegr*xstep
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + cxs*1.E-6
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - cxs
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + nxs1 + nxs2
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + fexs1 + fexs2
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp + zmortz
! Calcification
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - picfrx*(zmortp + zmortz)*1.E-6
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2.*picfrx*(zmortp + zmortz)*1.E-6
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + picfrx*(zmortp + zmortz)
               prodcal(ji,jj,jk) = picfrx*(zmortp + zmortz)         ! diagnostic array should be in mmol/m^-3/s but conversion is in canmeso for now
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_nano')
      !
      !
   END SUBROUTINE can_nano

   SUBROUTINE can_diat
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) :: zmortp,zmortz
      REAL(wp) :: spc,spn,spf,szc,chl
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n,thetac
      REAL(wp) :: cxs,nxs1,nxs2,fexs1,fexs2
      REAL(wp) :: csw1,csw2
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_diat')
      !

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               !     Phytoplankton mortality. 
               !     ------------------------
               spc = MAX(trb(ji,jj,jk,jpdia),0.)
               spn = MAX(trb(ji,jj,jk,jpdn),0.)
               spf = MAX(trb(ji,jj,jk,jpdfe),0.)
               szc = MAX(trb(ji,jj,jk,jpmes),0.)
               chl = MAX(trb(ji,jj,jk,jpdch),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)
               thetac=chl/(spc+rtrn)

               zmortp = mpratm * xstep * spc + mpqua * xstep * spc * spc
               if (spc.le.xminp) zmortp = mpqua * xstep * spc * spc           ! no linear mortality below biomass threshold xminp
               zmortz = mprat2 * xstep * szc + mpquad * xstep * szc * szc
               if (szc.le.xminp) zmortz = mpquad * xstep * szc * szc
! reduce mortality to what can support detritus production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the detritus ratio
               zmortp=zmortp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
               zmortpd(ji,jj,jk) = zmortp
! calculate "excess" relative to grazer RR
               cxs=zmortp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=zmortp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=zmortp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=zmortp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=zmortp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)

               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2

               !   Update the arrays tra which contains the biological sources and sinks
               !   ---------------------------------------------------------------------
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zmortp - cxs 
               tra(ji,jj,jk,jpdn) = tra(ji,jj,jk,jpdn) - zmortp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zmortp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - (zmortp+cxs) * thetac - chl*chldegr*xstep
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + cxs*1.E-6
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - cxs
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + nxs1 + nxs2
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + fexs1 + fexs2
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz 
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortp + zmortz
            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_diat')
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_diat')
      !
   END SUBROUTINE can_diat

   SUBROUTINE can_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismort/ mprat, mprat2, mpratm, mpqua, mpquad, chldegr, picfrx, xminp
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismort in reference namelist : Pisces phytoplankton
      READ  ( numnatp_ref, nampismort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : Pisces phytoplankton
      READ  ( numnatp_cfg, nampismort, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, nampismort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      mpqua     =', mpqua
         WRITE(numout,*) '    maximum quadratic mortality of diatoms    mpquad    =', mpquad
         WRITE(numout,*) '    phytoplankton mortality rate              mprat     =', mprat
         WRITE(numout,*) '    Diatoms mortality rate                    mprat2    =', mprat2
         WRITE(numout,*) '    Phytoplankton minimum mortality rate      mpratm    =', mpratm
         WRITE(numout,*) '    Chlorophyll photooxidation rate           chldegr   =', chldegr
         WRITE(numout,*) '    CaCO3 production rate                     picfrx    =', picfrx
         WRITE(numout,*) '    Biomass threshold for linear mortality    xminp     =', xminp
      ENDIF


   END SUBROUTINE can_mort_init


   !!======================================================================
END MODULE  canmort
