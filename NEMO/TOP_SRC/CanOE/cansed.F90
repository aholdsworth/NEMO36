MODULE cansed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   CanOE Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
#if defined key_canoe
   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE cansink         !  vertical flux of particulate matter due to sinking
   USE canopt          !  optical model
   USE cansbc          !  External source of nutrients 
   USE canint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_sed  
   PUBLIC   can_sed_alloc
 

   !! * Module variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday
   REAL(wp), POINTER, DIMENSION(:,:) ::  zw2d 
   REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d   

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cansed.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_sed  ***
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
      REAL(wp) ::   zdenitot, znitrpottot, zlim, zfact, zfactcal
      REAL(wp) ::   zcaloss, zwsbio3, zwsbio4, zwscal, zdep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: znitrpot, zirondep, zafe, zbfe                       ! afe and bfe indicate aeolian and benthic Fe sources
      REAL(wp), POINTER, DIMENSION(:,:) :: zocdep, zicdep, zburial                                ! deposition and burial of POC and PIC
      REAL(wp) :: kni         = 0.1_wp     !: half-saturation for NO3 inhibition of diazotrophy
      INTEGER  :: rmtss                      !: number of seconds per month 
      !!---------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('can_sed')
      !
      IF( kt == nittrc000 .AND. knt == 1 )   r1_rday  = 1. / rday
      ! Number of seconds per month 
      rmtss =  nyear_len(1) * rday / raamo
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, znitrpot,  zafe, zbfe     )
      CALL wrk_alloc( jpi, jpj, zocdep, zicdep, zburial     )



      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         CALL wrk_alloc( jpi, jpj, jpk, zirondep      )
         !                                              ! Iron and Si deposition at the surface
! dust is in kgFe m^-2 month^-1; zirondep is in nmolFe m^-3 s^-1
         zirondep(:,:,:) = 0.e0          ! Initialisation of variables USEd to compute deposition
      !------------------------------------------------------------------
      ! Iron deposition at the surface
         zirondep(:,:,1) = dustsolub * dust(:,:) / ( 55.85 * rmtss * fse3t(:,:,1) ) * 1.E+12
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
      ! Iron solubilization of particles in the water column
      ! ----------------------------------------------------
          DO jk = 2, jpkm1
            zirondep(:,:,jk) = dust(:,:) / ( wdust * 55.85 * rmtss ) * 1.e-4 * EXP( -fsdept(:,:,jk) / 1000. ) * 1.E+12
          END DO
         !                                              ! Iron solubilization of particles in the water column
            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + zirondep(:,:,:)*rfact2 
         ! 
         IF( lk_iomput ) THEN
           IF( knt == nrdttrc ) THEN
              IF( iom_use("Irondep"   ) ) THEN  ! nitrogen fixation 
                  zafe(:,:,:)  =   zirondep(:,:,:) * 1.E-9                      * tmask(:,:,:)      ! zirondep and ironsed are in nmol m^-3 s^-1
                  CALL iom_put( "Irondep", zafe  )  ! surface downward net flux of iron
               ENDIF 
            ENDIF
         ELSE                                    
            IF( ln_diatrc )  &
              &  trc2d(:,:,jp_can0_2d + 11) = zirondep(:,:,1) * 1.e+3 *rfact2r * fse3t(:,:,1) * tmask(:,:,1)
         ENDIF
            CALL wrk_dealloc( jpi, jpj, jpk, zirondep      )
         !                                              
      ENDIF
 ! Add the external input of iron which is 3D distributed
      ! (dust, river and sediment mobilization)
      ! ------------------------------------------------------
      IF( ln_ironsed ) THEN
         tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
         !
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironsed" ) )   &
            &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e-9 * tmask(:,:,:) ) ! iron inputs from sediments
      ENDIF

    
      ! Add the external input of nutrients from river
!      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 1, nk_rnf(ji,jj)
!                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) +  rivdip(ji,jj) * rfact2
!                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) +  rivdin(ji,jj) * rfact2
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) +  rivdic(ji,jj) * 3.e-5 * rfact2
!                  tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) +  rivdsi(ji,jj) * rfact2
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +  rivdic(ji,jj) * 2.631*rfact2
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +  ( rivalk(ji,jj) - rno3 * rivdic(ji,jj) ) * rfact2
               ENDDO
            ENDDO
         ENDDO
      ENDIF




!      ! OA: Warning, the following part is necessary, especially with Kriest
!      ! to avoid CFL problems above the sediments
!      ! --------------------------------------------------------------------
!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            ikt  = mbkt(ji,jj)
!            zdep = fse3t(ji,jj,ikt) / xstep
!            zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) )
!            zwscal (ji,jj) = MIN( 0.99 * zdep, wscal (ji,jj,ikt) )
!            zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) )
!         END DO
!      END DO
!
      ! Loss of Caco3 in the sediments. 
      ! -------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwscal  = wscal (ji,jj,ikt) * zdep
            zcaloss = trb(ji,jj,ikt,jpcal) * zwscal
            
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
            zfactcal = FLOAT(FLOOR(MIN( 1.-excess(ji,jj,ikt), 1.5 )))       ! set burial fraction to 1 if Omega>1 and 0 otherwise
            tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zfactcal * 2.E-6
            tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zfactcal * 1.E-6
            zicdep(ji,jj) = tra(ji,jj,ikt,jpcal) * wscal(ji,jj,ikt)         ! deposition in mmol m^-2 s^-1
            zburial(ji,jj) = tra(ji,jj,ikt,jpcal) * wscal(ji,jj,ikt) * zfactcal

         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwsbio4 = wsbio4(ji,jj,ikt) * zdep
            zwsbio3 = wsbio3(ji,jj,ikt) * zdep
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trb(ji,jj,ikt,jpgoc) * zwsbio4
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trb(ji,jj,ikt,jppoc) * zwsbio3
            zocdep(ji,jj) = tra(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) + tra(ji,jj,ikt,jpgoc) * wsbio4(ji,jj,ikt)      ! deposition in mmol m^-2 s^-1
! all deposition of POC is returned to bottom layer as inorganic nutrients
            tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + (tra(ji,jj,ikt,jpgoc) * zwsbio4 + tra(ji,jj,ikt,jppoc) * zwsbio3) * 1.E-6
            tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - (tra(ji,jj,ikt,jpgoc) * zwsbio4 + tra(ji,jj,ikt,jppoc) * zwsbio3)
            tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + (tra(ji,jj,ikt,jpgoc) * zwsbio4 + tra(ji,jj,ikt,jppoc) * zwsbio3) * rr_n2c
            tra(ji,jj,ikt,jpfer) = tra(ji,jj,ikt,jpfer) + (tra(ji,jj,ikt,jpgoc) * zwsbio4 + tra(ji,jj,ikt,jppoc) * zwsbio3) * rr_fe2c
            tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + (tra(ji,jj,ikt,jpgoc) * zwsbio4 + tra(ji,jj,ikt,jppoc) * zwsbio3) * rr_n2c * 1.E-6
         END DO
      END DO

      ! Nitrogen fixation (simple parameterization). The total gain
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         tra(:,:,1,jpno3) = tra(:,:,1,jpno3) + nitdep(:,:) * rfact2
         tra(:,:,1,jptal) = tra(:,:,1,jptal) - rno3 * nitdep(:,:) * rfact2
      ENDIF
!


      ! Potential nitrogen fixation dependant on temperature and iron
      ! -------------------------------------------------------------

      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zlim = kni / (kni + trn(ji,jj,jk,jpno3) + trn(ji,jj,jk,jpnh4))        ! CMOC function for NO3-inhibition
               znitrpot(ji,jj,jk) =  MAX( 0.e0, ( tgfuncp(ji,jj,jk) - 0.773 ) ) * 1.962   &
                 &                 *  zlim * tra(ji,jj,jk,jpfer) / ( concfediaz + tra(ji,jj,jk,jpfer) ) &
                 &                 * ( 1.- EXP( -etot(ji,jj,jk) / diazolight ) )
            END DO
         END DO 
      END DO

      znitrpottot = glob_sum( znitrpot(:,:,:) * cvol(:,:,:) )

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zfact = znitrpot(ji,jj,jk) * nitrfix * xstep
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 1.e-6 * zfact
           END DO
         END DO 
      END DO
      !
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r
            CALL wrk_alloc( jpi, jpj, zw2d )  
            CALL wrk_alloc( jpi, jpj, jpk, zw3d )
            IF( iom_use("Nfix"   ) ) THEN  ! nitrogen fixation 
                zw3d(:,:,:)  =   znitrpot(:,:,:) * nitrfix * r1_rday * 0.001  * tmask(:,:,:)      ! znitrpot is n.d., nitrfix is in mmol m^-3 d^-1
               CALL iom_put( "Nfix"   , zw3d  )  ! nitrogen fixation 
            ENDIF
            IF( iom_use("OCdep"   ) ) THEN  ! nitrogen fixation 
                !zw2d(:,:)  =   zocdep(:,:) * r1_rday * 0.001* tmask(:,:,1)      ! zocdep, zicdep, and zburial are in mmol m^-2 d^-1
                zw2d(:,:)  =   tmask(:,:,1)      ! zocdep, zicdep, and zburial are in mmol m^-2 d^-1
                CALL iom_put( "OCdep"  , zw2d)  ! OC sedimentation to seafloor
            ENDIF 
            IF(lwp) THEN
                !zw2d(:,:)  =   zicdep(:,:) * r1_rday * 0.001* tmask(:,:,1)
                zw2d(:,:)  =   tmask(:,:,1)      ! zocdep, zicdep, and zburial are in mmol m^-2 d^-1
                CALL iom_put( "ICdep"  , zw2d)  ! CaCO3 sedimentation to seafloor
            ENDIF 
            IF( iom_use("Burial"   ) ) THEN  ! nitrogen fixation 
                zw2d(:,:) =   zburial(:,:) * r1_rday * 0.001* tmask(:,:,1)
                CALL iom_put( "Burial" , zw2d) ! CaCO3 burial
            ENDIF 
         CALL wrk_dealloc( jpi, jpj, zw2d )
         CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
         ENDIF
      ELSE
         IF( ln_diatrc ) THEN
            trc2d(:,:,jp_can0_2d + 12) = znitrpot(:,:,1) * nitrfix * zfact * fse3t(:,:,1) * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, znitrpot, zafe, zbfe )
      CALL wrk_dealloc( jpi, jpj, zocdep, zicdep, zburial     )
      !

      !
      IF( nn_timing == 1 )  CALL timing_stop('can_sed')
      !
 9100  FORMAT(i8,3f10.5)
      !
   END SUBROUTINE can_sed


   INTEGER FUNCTION can_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(jpi,jpj,jpk), sdenit(jpi,jpj), STAT=can_sed_alloc )
      !
      IF( can_sed_alloc /= 0 )   CALL ctl_warn('can_sed_alloc: failed to allocate arrays')
      !
   END FUNCTION can_sed_alloc


#else
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_sed                         ! Empty routine
   END SUBROUTINE can_sed
#endif 

   !!======================================================================
END MODULE  cansed
