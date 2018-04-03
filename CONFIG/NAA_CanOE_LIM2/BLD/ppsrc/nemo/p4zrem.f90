MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zint          !  interpolation and computation of various fields
   USE p4zlim
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xremip     !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::  nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::  xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsilab     !: fraction of labile biogenic silica 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr     !: denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitnh4   !: -    -    -    -   -

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
   !! $Id: p4zrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik, zsiremin, zammonic 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zorem, zorem2, zofer, zolimit
      REAL(wp) ::   zosil, ztem
      REAL(wp) ::   zofer2
      REAL(wp) ::   zonitr, zstep, zfact
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: ztempbac
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdepbac, zolimi, zdepprod, zw3d
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zoxyrem
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      ztempbac                  )
      CALL wrk_alloc( jpi, jpj, jpk, zdepbac, zdepprod, zolimi, zoxyrem )

      ! Initialisation of temprary arrys
      zdepprod(:,:,:) = 1._wp
      ztempbac(:,:)   = 0._wp

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( gdept_n(ji,jj,jk) < zdep ) THEN
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trb(ji,jj,jk,jpzoo) + 2.* trb(ji,jj,jk,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepmin = MIN( 1., zdep / gdept_n(ji,jj,jk) )
                  zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj)
                  zdepprod(ji,jj,jk) = zdepmin**0.273
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
               ! DOC ammonification. Depends on depth, phytoplankton biomass
               ! and a limitation term which is supposed to be a parameterization
               !     of the bacterial activity. 
               zremik = xremik * zstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
               zremik = MAX( zremik, 2.74e-4 * xstep )
               ! Ammonification in oxic waters with oxygen consumption
               ! -----------------------------------------------------
               zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
               zolimi(ji,jj,jk) = MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) 
               ! Ammonification in suboxic waters with denitrification
               ! -------------------------------------------------------
               zammonic = zremik * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc)
               denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
               zoxyrem(ji,jj,jk) = zammonic *        nitrfac2(ji,jj,jk)
               !
               zolimi (ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
               denitr (ji,jj,jk) = MAX( 0.e0, denitr (ji,jj,jk) )
               zoxyrem(ji,jj,jk) = MAX( 0.e0, zoxyrem(ji,jj,jk) )
               !
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
               ! NH4 nitrification to NO3. Ceased for oxygen concentrations
               ! below 2 umol/L. Inhibited at strong light 
               ! ----------------------------------------------------------
               zonitr  =nitrif * zstep * trb(ji,jj,jk,jpnh4) / ( 1.+ emoy(ji,jj,jk) ) * ( 1.- nitrfac(ji,jj,jk) ) 
               denitnh4(ji,jj,jk) = nitrif * zstep * trb(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk) 
               ! Update of the tracers trends
               ! ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr - denitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr - rdenita * denitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2nit * zonitr
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * rno3 * zonitr &
                  &                + rno3 * ( rdenita - 1. ) * denitnh4(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Bacterial uptake of iron. No iron is available in DOC. So
               ! Bacteries are obliged to take up iron from the water. Some
               ! studies (especially at Papa) have shown this uptake to be significant
               ! ----------------------------------------------------------
               zbactfer = 10.e-6 *  rfact2 * prmax(ji,jj,jk) * xlimbacl(ji,jj,jk)             &
                  &              * trb(ji,jj,jk,jpfer) / ( 2.5E-10 + trb(ji,jj,jk,jpfer) )    &
                  &              * zdepprod(ji,jj,jk) * zdepbac(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer*0.16
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer*0.12
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer*0.04
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
               ! POC disaggregation by turbulence and bacterial activity. 
               ! --------------------------------------------------------
               zremip = xremip * zstep * tgfunc(ji,jj,jk) * ( 1.- 0.55 * nitrfac(ji,jj,jk) ) 

               ! POC disaggregation rate is reduced in anoxic zone as shown by
               ! sediment traps data. In oxic area, the exponent of the martin s
               ! law is around -0.87. In anoxic zone, it is around -0.35. This
               ! means a disaggregation constant about 0.5 the value in oxic zones
               ! -----------------------------------------------------------------
               zorem  = zremip * trb(ji,jj,jk,jppoc)
               zofer  = zremip * trb(ji,jj,jk,jpsfe)
               zorem2 = zremip * trb(ji,jj,jk,jpgoc)
               zofer2 = zremip * trb(ji,jj,jk,jpbfe)

               ! Update the appropriate tracers trends
               ! -------------------------------------

               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zorem
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zofer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zorem2 - zorem
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zorem2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zofer2 - zofer
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) - zofer2

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
               ! Remineralization rate of BSi depedant on T and saturation
               ! ---------------------------------------------------------
               zsatur   = ( sio3eq(ji,jj,jk) - trb(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn )
               zsatur   = MAX( rtrn, zsatur )
               zsatur2  = ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**37
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
               znusil2  = 0.225  * ( 1. + tsn(ji,jj,1,jp_tem) / 15.) + 0.775 * zsatur2

               ! Two classes of BSi are considered : a labile fraction and 
               ! a more refractory one. The ratio between both fractions is
               ! constant and specified in the namelist.
               ! ----------------------------------------------------------
               zdep     = MAX( hmld(ji,jj), heup(ji,jj) ) 
               zdep     = MAX( 0., gdept_n(ji,jj,jk) - zdep )
               ztem     = MAX( tsn(ji,jj,1,jp_tem), 0. )
               zfactdep = xsilab * EXP(-( xsiremlab - xsirem ) * znusil2 * zdep / wsbio2 ) * ztem / ( ztem + 10. )
               zsiremin = ( xsiremlab * zfactdep + xsirem * ( 1. - zfactdep ) ) * zstep * znusil
               zosil    = zsiremin * trb(ji,jj,jk,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
               !
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem4')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      ! Update the arrays TRA which contain the biological sources and sinks
      ! --------------------------------------------------------------------

      DO jk = 1, jpkm1
         tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zolimi (:,:,jk) + denitr(:,:,jk) + zoxyrem(:,:,jk)
         tra(:,:,jk,jpnh4) = tra(:,:,jk,jpnh4) + zolimi (:,:,jk) + denitr(:,:,jk) + zoxyrem(:,:,jk)
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) - denitr (:,:,jk) * rdenit
         tra(:,:,jk,jpdoc) = tra(:,:,jk,jpdoc) - zolimi (:,:,jk) - denitr(:,:,jk) - zoxyrem(:,:,jk)
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - zolimi (:,:,jk) * o2ut
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + zolimi (:,:,jk) + denitr(:,:,jk) + zoxyrem(:,:,jk)
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) + rno3 * ( zolimi(:,:,jk) + zoxyrem(:,:,jk) &
              &                                        + ( rdenit + 1.) * denitr(:,:,jk) )
      END DO

      IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN
              zw3d(:,:,:) = zolimi(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN"  , zw3d )
          ENDIF
          IF( iom_use( "DENIT" ) )  THEN
              zw3d(:,:,:) = denitr(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENIT"  , zw3d )
          ENDIF
          !
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
       ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      ztempbac                  )
      CALL wrk_dealloc( jpi, jpj, jpk, zdepbac, zdepprod, zolimi, zoxyrem )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, xremip, nitrif, xsirem, xsiremlab, xsilab
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
         WRITE(numout,*) '    remineralisation rate of POC              xremip    =', xremip
         WRITE(numout,*) '    remineralization rate of DOC              xremik    =', xremik
         WRITE(numout,*) '    remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '    fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '    fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '    NH4 nitrification rate                    nitrif    =', nitrif
      ENDIF
      !
      denitr  (:,:,:) = 0._wp
      denitnh4(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(jpi,jpj,jpk), denitnh4(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_warn('p4z_rem_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_rem_alloc


   !!======================================================================
END MODULE p4zrem
