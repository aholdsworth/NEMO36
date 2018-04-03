MODULE canmicro
   !!======================================================================
   !!                         ***  MODULE canmicro  ***
   !! TOP :   CanOE Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_micro       :   Compute the sources/sinks for microzooplankton
   !!   can_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canint          !  interpolation and computation of various fields
   USE canprod         !  production
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_micro         ! called in canbio.F90
   PUBLIC   can_micro_init    ! called in trcsms_canoe.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part      = 0.5_wp          !: maximum grazing rate
   REAL(wp), PUBLIC ::  gmax1     = 1.7_wp          !: maximum grazing rate
   REAL(wp), PUBLIC ::  aps       = 0.075_wp        !: small zooplankton functional response parameter
   REAL(wp), PUBLIC ::  zsr1      = 0.3_wp          !: specific respiration rate
   REAL(wp), PUBLIC ::  lambda1   = 0.8_wp          !: assimilation efficiency

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
   !! $Id: canmicro.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE can_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  knt 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,spc,spn,spf,chl,grazp,R
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      REAL(wp) :: zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zw3d
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_micro')
      !
      !
      grazing1(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ztn = tsn(ji,jj,jk,jp_tem)
               Tf = tgfuncz(ji,jj,jk) 

               spc = MAX(trb(ji,jj,jk,jpphy),0.)
               spn = MAX(trb(ji,jj,jk,jpnn),0.)
               spf = MAX(trb(ji,jj,jk,jpnfe),0.)
               chl = MAX(trb(ji,jj,jk,jpnch),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)

! Micrograzer functional response is determined by phytoplankton C
               grazp=gmax1*(1.-EXP(-aps*spc))*trb(ji,jj,jk,jpzoo)*xstep
! reduce phytoplankton consumption to what can support grazer biomass production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the grazer ratio
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
               R = MAX(zsr1*Tf*trb(ji,jj,jk,jpzoo)*xstep-cxs,0.)

               ! Grazing by microzooplankton
               grazing1(ji,jj,jk) = grazp

               !  Update of the TRA arrays
               !  ------------------------
               !zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + R*rr_n2c + nxs1 + nxs2
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - R - cxs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + R*rr_fe2c + fexs1 + fexs2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (R + cxs)*1.E-6
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + lambda1*grazp - R
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - grazp - cxs
               tra(ji,jj,jk,jpnn) = tra(ji,jj,jk,jpnn) - grazp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - grazp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - (grazp+cxs)*chl/(spc+rtrn)
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + (1.-lambda1)*grazp
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, zw3d )
         IF( iom_use( "GRAZ1" ) ) THEN
            zw3d(:,:,:) =grazing1(:,:,:) * 1.e-3 * rfact2r * tmask(:,:,:)  !  Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ1", zw3d )
         ENDIF
         CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_micro')
      !
   END SUBROUTINE can_micro


   SUBROUTINE can_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiszoo/ part, gmax1, aps, zsr1, lambda1
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : Pisces microzooplankton
      READ  ( numnatp_ref, nampiszoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : Pisces microzooplankton
      READ  ( numnatp_cfg, nampiszoo, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '    Microzooplankton maximum grazing rate           gmax1       =', gmax1
         WRITE(numout,*) '    Microzooplankton grazing initial slope          aps         =', aps
         WRITE(numout,*) '    Microzooplankton specific respiration           zsr1        =', zsr1
         WRITE(numout,*) '    Microzooplankton unassimilated fraction         lambda1     =', lambda1
      ENDIF


   END SUBROUTINE can_micro_init


   !!======================================================================
END MODULE  canmicro