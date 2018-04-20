MODULE canmeso
   !!======================================================================
   !!                         ***  MODULE canmeso  ***
   !! TOP :   CanOE Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_meso       :   Compute the sources/sinks for mesozooplankton
   !!   can_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canint          !  interpolation and computation of various fields
   USE canprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_meso              ! called in canbio.F90
   PUBLIC   can_meso_init         ! called in trcsms_canoe.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts

   REAL(wp), PUBLIC ::  gmax2      = 0.85_wp         !: maximum grazing rate rate
   REAL(wp), PUBLIC ::  apl        = 0.075_wp        !: large zooplankton functional response parameter
   REAL(wp), PUBLIC ::  zsr2       = 0.3_wp          !: specific respiration rate
   REAL(wp), PUBLIC ::  lambda2    = 0.8_wp          !: assimilation efficiency

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
   !! $Id: canmeso.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE can_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,lpc,lpn,lpf,chl,grazt,grazp,grazz,R,szc,itfc
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2
      !!---------------------------------------------------------------------
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_meso')
      !
      grazing2(:,:,:) = 0.  !: grazing set to zero
      grazing3(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztn = tsn(ji,jj,jk,jp_tem)
               Tf = tgfuncz2(ji,jj,jk)

               lpc = MAX(trb(ji,jj,jk,jpdia),0.)
               lpn = MAX(trb(ji,jj,jk,jpdn),0.)
               lpf = MAX(trb(ji,jj,jk,jpdfe),0.)
               chl = MAX(trb(ji,jj,jk,jpdch),0.)
               szc = MAX(trb(ji,jj,jk,jpzoo),0.)
               itfc=1./(lpc+szc+rtrn)

               c2n=lpc/(lpn+rtrn)
               n2c=lpn/(lpc+rtrn)
               c2fe=lpc/(lpf+rtrn)
               fe2c=lpf/(lpc+rtrn)
               n2fe=lpn/(lpf+rtrn)
               fe2n=lpf/(lpn+rtrn)

! assume grazing hyperbola is determined by total food concentration and the two food types are consumed in proportion to their concentrations (in C units)
               grazt=gmax2*(1.-EXP(-apl*(lpc+szc)))*trb(ji,jj,jk,jpmes)*xstep
               grazz=grazt*szc*itfc
               grazp=grazt*lpc*itfc
! reduce phytoplankton fraction to what can support grazer biomass production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the grazer ratio
               grazp=grazp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
! calculate "excess" relative to grazer RR
               cxs=grazp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=grazp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=grazp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=grazp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=grazp*rr_fe2c*(fe2n*rr_n2fe-1.)
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
! calculate zooplankton respiration (in carbon units)
               R = MAX(zsr2*Tf*trb(ji,jj,jk,jpmes)*xstep-cxs,0.)

               !   Update the arrays TRA which contain the biological sources and sinks
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + R*rr_n2c + nxs1 + nxs2
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - R - cxs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + R*rr_fe2c + fexs1 + fexs2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (R + cxs)*1.E-6
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - grazz
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + lambda2*(grazz+grazp) - R
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - grazp - cxs
               tra(ji,jj,jk,jpdn) = tra(ji,jj,jk,jpdn) - grazp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - grazp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - (grazp+cxs)*chl/(lpc+rtrn)
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + (1.-lambda2)*(grazz+grazp)

! I am excluding cxs from grazing diagnostics (represent zooplankton gains rather than phytoplankton losses)
               grazing2(ji,jj,jk) = grazp
               grazing3(ji,jj,jk) = grazz

            END DO
         END DO
      END DO
      !
      IF( ln_diatrc .AND. lk_iomput ) THEN
         zrfact2 = 1.e-3 * rfact2r
!         grazing(:,:,:) = grazing(:,:,:) * zrfact2 * tmask(:,:,:)   ! Total grazing of phyto by zoo
         prodcal(:,:,:) = prodcal(:,:,:) * zrfact2 * tmask(:,:,:)   ! Calcite production
         IF( knt == nrdttrc ) THEN
            CALL iom_put( "GRAZ2" , grazing2 * zrfact2 * tmask(:,:,:) )  ! Grazing of large phytoplankton
            CALL iom_put( "GRAZ3" , grazing3 * zrfact2 * tmask(:,:,:) )  ! Grazing of microzooplankton
            CALL iom_put( "PCAL" , prodcal  )  ! Calcite production
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_meso')
      !
   END SUBROUTINE can_meso

   SUBROUTINE can_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ part2, gmax2, apl, zsr2, lambda2
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, nampismes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, nampismes, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismes )



      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '    Mesozooplankton maximum grazing rate           gmax2        =', gmax2
         WRITE(numout,*) '    Mesozooplankton grazing initial slope          apl          =', apl
         WRITE(numout,*) '    Mesozooplankton specific respiration           zsr2         =', zsr2
         WRITE(numout,*) '    Mesozooplankton unassimilated fraction         lambda2      =', lambda2
      ENDIF


   END SUBROUTINE can_meso_init



   !!======================================================================
END MODULE  canmeso
