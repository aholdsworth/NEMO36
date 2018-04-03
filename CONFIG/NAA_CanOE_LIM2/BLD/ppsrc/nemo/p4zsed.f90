MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed  
   PUBLIC   p4z_sed_alloc
 

   !! * Module variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday

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
   !! $Id: p4zsed.F90 8610 2017-10-11 09:50:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
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
      REAL(wp) ::   zsumsedsi, zsumsedpo4, zsumsedcal
      REAL(wp) ::   zrivalk, zrivsil, zrivno3
      REAL(wp) ::  zwflux, zfminus, zfplus
      REAL(wp) ::  zlim, zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zflx, zpdenit, z1pdenit, zolimit
      REAL(wp) ::  zsiloss, zcaloss, zws3, zws4, zwsc, zdep, zwstpoc
      REAL(wp) ::  ztrfer, ztrpo4, zwdust, zlight
      !
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zpdep, zsidep, zwork1, zwork2, zwork3
      REAL(wp), POINTER, DIMENSION(:,:)   :: zsedcal, zsedsi, zsedc
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zdenit2d, zironice, zbureff
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zwsbio3, zwsbio4, zwscal
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zirondep, zsoufer
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed')
      !
      IF( kt == nittrc000 .AND. knt == 1 )   r1_rday  = 1. / rday
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, zdenit2d, zwork1, zwork2, zwork3, zbureff )
      CALL wrk_alloc( jpi, jpj, zsedcal,  zsedsi, zsedc )
      CALL wrk_alloc( jpi, jpj, zwsbio3, zwsbio4, zwscal )
      CALL wrk_alloc( jpi, jpj, jpk, zsoufer )

      zdenit2d(:,:) = 0.e0
      zbureff (:,:) = 0.e0
      zwork1  (:,:) = 0.e0
      zwork2  (:,:) = 0.e0
      zwork3  (:,:) = 0.e0
      zsedsi   (:,:) = 0.e0
      zsedcal  (:,:) = 0.e0
      zsedc    (:,:) = 0.e0

      ! Iron input/uptake due to sea ice : Crude parameterization based on Lancelot et al.
      ! ----------------------------------------------------
      IF( ln_ironice ) THEN  
         !                                              
         CALL wrk_alloc( jpi, jpj, zironice )
         !                                              
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep    = rfact2 / e3t_n(ji,jj,1)
               zwflux  = fmmflx(ji,jj) / 1000._wp
               zfminus = MIN( 0._wp, -zwflux ) * trb(ji,jj,1,jpfer) * zdep
               zfplus  = MAX( 0._wp, -zwflux ) * icefeinput * zdep
               zironice(ji,jj) =  zfplus + zfminus
            END DO
         END DO
         !
         tra(:,:,1,jpfer) = tra(:,:,1,jpfer) + zironice(:,:) 
         ! 
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironice" ) )   &
            &   CALL iom_put( "Ironice", zironice(:,:) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! iron flux from ice
         !
         CALL wrk_dealloc( jpi, jpj, zironice )
         !                                              
      ENDIF

      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         CALL wrk_alloc( jpi, jpj,      zpdep, zsidep )
         CALL wrk_alloc( jpi, jpj, jpk, zirondep      )
         !                                              ! Iron and Si deposition at the surface
         IF( ln_solub ) THEN
            zirondep(:,:,1) = solub(:,:) * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ELSE
            zirondep(:,:,1) = dustsolub  * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ENDIF
         zsidep(:,:) = 8.8 * 0.075 * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 28.1 
         zpdep (:,:) = 0.1 * 0.021 * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 31. / po4r 
         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 270. * rday )
         DO jk = 2, jpkm1
            zirondep(:,:,jk) = dust(:,:) * mfrac * zwdust * rfact2 * EXP( -gdept_n(:,:,jk) / 540. )
         END DO
         !                                              ! Iron solubilization of particles in the water column
         tra(:,:,1,jppo4) = tra(:,:,1,jppo4) + zpdep   (:,:)
         tra(:,:,1,jpsil) = tra(:,:,1,jpsil) + zsidep  (:,:)
         tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + zirondep(:,:,:) 
         ! 
         IF( lk_iomput ) THEN
            IF( knt == nrdttrc ) THEN
                IF( iom_use( "Irondep" ) )   &
                &  CALL iom_put( "Irondep", zirondep(:,:,1) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! surface downward dust depo of iron
                IF( iom_use( "pdust" ) )   &
                &  CALL iom_put( "pdust"  , dust(:,:) / ( wdust * rday )  * tmask(:,:,1) ) ! dust concentration at surface
            ENDIF
         ELSE                                    
            IF( ln_diatrc )  &
              &  trc2d(:,:,jp_pcs0_2d + 11) = zirondep(:,:,1) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1)
         ENDIF
         CALL wrk_dealloc( jpi, jpj,      zpdep, zsidep )
         CALL wrk_dealloc( jpi, jpj, jpk, zirondep      )
         !                                              
      ENDIF
     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 1, nk_rnf(ji,jj)
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) +  rivdip(ji,jj) * rfact2
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) +  rivdin(ji,jj) * rfact2
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) +  rivdic(ji,jj) * 5.e-5 * rfact2
                  tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) +  rivdsi(ji,jj) * rfact2
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +  rivdic(ji,jj) * rfact2
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +  ( rivalk(ji,jj) - rno3 * rivdin(ji,jj) ) * rfact2
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         tra(:,:,1,jpno3) = tra(:,:,1,jpno3) + nitdep(:,:) * rfact2
         tra(:,:,1,jptal) = tra(:,:,1,jptal) - rno3 * nitdep(:,:) * rfact2
      ENDIF

      ! Add the external input of iron from sediment mobilization
      ! ------------------------------------------------------
      IF( ln_ironsed ) THEN
         tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
         !
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironsed" ) )   &
            &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! iron inputs from sediments
      ENDIF

      ! Add the external input of iron from hydrothermal vents
      ! ------------------------------------------------------
      IF( ln_hydrofe ) THEN
         tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + hydrofe(:,:,:) * rfact2
         !
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "HYDR" ) )   &
            &   CALL iom_put( "HYDR", hydrofe(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! hydrothermal iron input
      ENDIF

      ! OA: Warning, the following part is necessary, especially with Kriest
      ! to avoid CFL problems above the sediments
      ! --------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = e3t_n(ji,jj,ikt) / xstep
            zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) )
            zwscal (ji,jj) = MIN( 0.99 * zdep, wscal (ji,jj,ikt) )
            zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) )
         END DO
      END DO

      ! Computation of the sediment denitrification proportion: The metamodel from midlleburg (2006) is being used
      ! Computation of the fraction of organic matter that is permanently buried from Dunne's model
      ! -------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
           IF( tmask(ji,jj,1) == 1 ) THEN
              ikt = mbkt(ji,jj)
              zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) )  * 1E3 * 1E6 / 1E4
              zflx  = LOG10( MAX( 1E-3, zflx ) )
              zo2   = LOG10( MAX( 10. , trb(ji,jj,ikt,jpoxy) * 1E6 ) )
              zno3  = LOG10( MAX( 1.  , trb(ji,jj,ikt,jpno3) * 1E6 * rno3 ) )
              zdep  = LOG10( gdepw_n(ji,jj,ikt+1) )
              zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
              &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
              zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) )
              !
              zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) ) * 1E6
              zbureff(ji,jj) = 0.013 + 0.53 * zflx**2 / ( 7.0 + zflx )**2
           ENDIF
         END DO
      END DO 

      ! Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
      ! First, the total loss is computed.
      ! The factor for calcite comes from the alkalinity effect
      ! -------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tmask(ji,jj,1) == 1 ) THEN
               ikt = mbkt(ji,jj) 
               zwork1(ji,jj) = trb(ji,jj,ikt,jpgsi) * zwsbio4(ji,jj)
               zwork2(ji,jj) = trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj) + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) 
               ! For calcite, burial efficiency is made a function of saturation
               zfactcal      = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal      = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
               zwork3(ji,jj) = trb(ji,jj,ikt,jpcal) * zwscal(ji,jj) * 2.e0 * zfactcal
            ENDIF
         END DO
      END DO
      zsumsedsi  = glob_sum( zwork1(:,:) * e1e2t(:,:) ) * r1_rday
      zsumsedpo4 = glob_sum( zwork2(:,:) * e1e2t(:,:) ) * r1_rday
      zsumsedcal = glob_sum( zwork3(:,:) * e1e2t(:,:) ) * r1_rday

      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
      zrivsil =  1._wp - ( sumdepsi + rivdsiinput * r1_ryyss ) / ( zsumsedsi + rtrn )

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) 
            zws4 = zwsbio4(ji,jj) * zdep
            zwsc = zwscal (ji,jj) * zdep
            zsiloss = trb(ji,jj,ikt,jpgsi) * zwsc
            zcaloss = trb(ji,jj,ikt,jpcal) * zwsc
            !
            tra(ji,jj,ikt,jpgsi) = tra(ji,jj,ikt,jpgsi) - zsiloss
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
            tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
            zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
            zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
            zrivalk  =  1._wp - ( rivalkinput * r1_ryyss ) * zfactcal / ( zsumsedcal + rtrn )
            tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
            tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zrivalk
            zsedcal(ji,jj) = (1.0 - zrivalk) * zcaloss * e3t_n(ji,jj,ikt) 
            zsedsi (ji,jj) = (1.0 - zrivsil) * zsiloss * e3t_n(ji,jj,ikt) 
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) 
            zws4 = zwsbio4(ji,jj) * zdep
            zws3 = zwsbio3(ji,jj) * zdep
            zrivno3 = 1. - zbureff(ji,jj)
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trb(ji,jj,ikt,jpgoc) * zws4 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trb(ji,jj,ikt,jppoc) * zws3
            tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trb(ji,jj,ikt,jpbfe) * zws4
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trb(ji,jj,ikt,jpsfe) * zws3
            zwstpoc              = trb(ji,jj,ikt,jpgoc) * zws4 + trb(ji,jj,ikt,jppoc) * zws3

            ! The 0.5 factor in zpdenit is to avoid negative NO3 concentration after
            ! denitrification in the sediments. Not very clever, but simpliest option.
            zpdenit  = MIN( 0.5 * ( trb(ji,jj,ikt,jpno3) - rtrn ) / rdenit, &
                            zdenit2d(ji,jj) * zwstpoc * zrivno3 * (1. - nitrfac2(ji,jj,ikt) ) )
            z1pdenit = zwstpoc * zrivno3 - zpdenit
            zolimit = MIN( ( trb(ji,jj,ikt,jpoxy) - rtrn ) / o2ut, z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
            tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + z1pdenit - zolimit
            tra(ji,jj,ikt,jppo4) = tra(ji,jj,ikt,jppo4) + zpdenit + zolimit
            tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + zpdenit + zolimit
            tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) - rdenit * zpdenit
            tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - zolimit * o2ut
            tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + rno3 * (zolimit + (1.+rdenit) * zpdenit )
            tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zpdenit + zolimit
            sdenit(ji,jj) = rdenit * zpdenit * e3t_n(ji,jj,ikt) 
            zsedc(ji,jj)   = (1. - zrivno3) * zwstpoc * e3t_n(ji,jj,ikt) 
         END DO
      END DO

      ! Nitrogen fixation process
      ! Small source iron from particulate inorganic iron
      !-----------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !                      ! Potential nitrogen fixation dependant on temperature and iron
               zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
               IF( zlim <= 0.2 )   zlim = 0.01
               zfact = zlim * rfact2
               ztrfer = biron(ji,jj,jk)       / ( concfediaz + biron(ji,jj,jk)       )
               ztrpo4 = trb  (ji,jj,jk,jppo4) / ( concnnh4   + trb  (ji,jj,jk,jppo4) ) 
               zlight =  ( 1.- EXP( -etot_ndcy(ji,jj,jk) / diazolight ) ) 
               nitrpot(ji,jj,jk) =  MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) * r1_rday )   &
                 &         *  zfact * MIN( ztrfer, ztrpo4 ) * zlight
               zsoufer(ji,jj,jk) = zlight * 2E-11 / (2E-11 + biron(ji,jj,jk))
            END DO
         END DO
      END DO

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zfact = nitrpot(ji,jj,jk) * nitrfix
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) +             zfact
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3      * zfact
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2nit     * zfact 
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,jk,jppo4) ) &
               &                     * 0.002 * trb(ji,jj,jk,jpdoc) * xstep
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * xstep
           END DO
         END DO 
      END DO

      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r  !  conversion from molC/l/kt  to molC/m3/s
            IF( iom_use("Nfix"   ) ) CALL iom_put( "Nfix", nitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation 
            IF( iom_use("INTNFIX") ) THEN   ! nitrogen fixation rate in ocean ( vertically integrated )
               zwork1(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork1(:,:) = zwork1(:,:) + nitrpot(:,:,jk) * nitrfix * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX" , zwork1 ) 
            ENDIF
            IF( iom_use("SedCal" ) ) CALL iom_put( "SedCal", zsedcal(:,:) * zfact )
            IF( iom_use("SedSi" ) )  CALL iom_put( "SedSi",  zsedsi (:,:) * zfact )
            IF( iom_use("SedC" ) )   CALL iom_put( "SedC",   zsedc  (:,:) * zfact )
            IF( iom_use("Sdenit" ) ) CALL iom_put( "Sdenit", sdenit (:,:) * rno3 * zfact )
         ENDIF
      ELSE
         IF( ln_diatrc ) THEN 
            zfact = 1.e+3 * rfact2r  !  conversion from molC/l/kt  to molC/m3/s
            trc2d(:,:,jp_pcs0_2d + 12) = nitrpot(:,:,1) * nitrfix * rno3 * zfact * e3t_n(:,:,1) * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zdenit2d, zwork1, zwork2, zwork3, zbureff )
      CALL wrk_dealloc( jpi, jpj, zsedcal , zsedsi, zsedc )
      CALL wrk_dealloc( jpi, jpj, zwsbio3, zwsbio4, zwscal )
      CALL wrk_dealloc( jpi, jpj, jpk, zsoufer )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed')
      !
 9100  FORMAT(i8,3f10.5)
      !
   END SUBROUTINE p4z_sed


   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(jpi,jpj,jpk), sdenit(jpi,jpj), STAT=p4z_sed_alloc )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_warn('p4z_sed_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_sed_alloc



   !!======================================================================
END MODULE p4zsed
