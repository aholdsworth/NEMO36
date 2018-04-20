MODULE canlys
   !!======================================================================
   !!                         ***  MODULE canlys  ***
   !! TOP :   CanOE 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_lys        :   Compute the CaCO3 dissolution 
   !!   can_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  ! USE canche
   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_lys         ! called in trcsms_canoe.F90
   PUBLIC   can_lys_init    ! called in trcsms_canoe.F90

   !! * Shared module variables
   REAL(wp), PUBLIC :: kdca=7.4E-3_wp !: diss. rate constant calcite
   REAL(wp), PUBLIC :: nca =1.0_wp !: order of reaction for calcite dissolution

   !! * Module variables
   REAL(wp) :: calcon = 1.03E-2           !: mean calcite concentration [Ca2+] in sea water [mole/kg solution]
   REAL(wp) :: r1_rday = 1.0_wp/86400._wp !: 1 / rday
 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canlys.F90 3321 2012-03-05 17:10:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE can_lys( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zalk, zdic, zph, zah2
      REAL(wp) ::   zdispot, zfact, zcalcon, zalka, zaldi
      REAL(wp) ::   zomegaca, zexcess, zexcess0
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zco3, zcaldiss   
      !REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zco3, zcaldiss   
      !!---------------------------------------------------------------------
      !

      !
      CALL wrk_alloc( jpi, jpj, jpk, zco3, zcaldiss )
      
!
      zco3(:,:,:) = 0.
      zcaldiss(:,:,:) = 0.
      !     -------------------------------------------
      
      DO jn = 1, 5                               !  BEGIN OF ITERATION
         !
!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zfact = rhop(ji,jj,jk) / 1000. + rtrn
                  zph  = hi(ji,jj,jk) * tmask(ji,jj,jk) / zfact + ( 1.-tmask(ji,jj,jk) ) * 1.e-9 ! [H+]
                  zdic  = trb(ji,jj,jk,jpdic)/zfact
                  zalka = trb(ji,jj,jk,jptal) / zfact
                  ! CALCULATE [ALK]([CO3--], [HCO3-])
                  zalk  = zalka - ( akw3(ji,jj,jk) / zph - zph + borat(ji,jj,jk) / ( 1. + zph / akb3(ji,jj,jk) ) )
                  ! CALCULATE [H+] and [CO3--]
                  zaldi = zdic - zalk
                  zah2  = SQRT( zaldi * zaldi + 4.* ( zalk * ak23(ji,jj,jk) / ak13(ji,jj,jk) ) * ( zdic + zaldi ) )
                  zah2  = 0.5 * ak13(ji,jj,jk) / zalk * ( zaldi + zah2 )
                  !
                  zco3(ji,jj,jk) = zalk / ( 2. + zah2 / ak23(ji,jj,jk) ) * zfact
                  hi(ji,jj,jk)   = zah2 * zfact
               END DO
            END DO
         END DO
         !
      END DO 

      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! DEVIATION OF [CO3--] FROM SATURATION VALUE
               ! Salinity dependance in zomegaca and divide by rhop/1000 to have good units
               zcalcon  = calcon * ( tsn(ji,jj,jk,jp_sal) / 35._wp )
               zfact    = rhop(ji,jj,jk) / 1000._wp
               zomegaca = ( zcalcon * zco3(ji,jj,jk) * zfact ) / aksp(ji,jj,jk) 

               ! SET DEGREE OF UNDER-/SUPERSATURATION
               excess(ji,jj,jk) = 1._wp - zomegaca
            !#zexcess0 = MAX( 0., excess(ji,jj,jk) )
             !  zexcess  = zexcess0**nca

               ! AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
               !       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
               !       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
               !zdispot = kdca * zexcess * trb(ji,jj,jk,jpcal)
               zdispot = kdca * r1_rday * trb(ji,jj,jk,jpcal)
              !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
              !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
              zcaldiss(ji,jj,jk)  = zdispot  *rfact2                       ! calcite dissolution (first order, no saturation state dependence)
              zco3(ji,jj,jk)      = zco3(ji,jj,jk) + zcaldiss(ji,jj,jk) 
              !
              tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) -         zcaldiss(ji,jj,jk)
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2.E-6 * zcaldiss(ji,jj,jk)       
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + 1.E-6 * zcaldiss(ji,jj,jk) 
            END DO
        END DO
     END DO
      !

      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         IF( iom_use( "PH"     ) ) CALL iom_put( "PH"    , -1. * LOG10( hi(:,:,:) )          * tmask(:,:,:) )
         IF( iom_use( "CO3"    ) ) CALL iom_put( "CO3"   , zco3(:,:,:) * 1.e+3               * tmask(:,:,:) )
         IF( iom_use( "CO3sat" ) ) CALL iom_put( "CO3sat", aksp(:,:,:) * 1.e+3 / calcon      * tmask(:,:,:) )
         IF( iom_use( "DCAL"   ) ) CALL iom_put( "DCAL"  , zcaldiss(:,:,:) * 1.e-3 *rfact2r    * tmask(:,:,:) )
      ELSE
         trc3d(:,:,:,jp_can0_3d    ) = -1. * LOG10( hi(:,:,:) ) * tmask(:,:,:)
         trc3d(:,:,:,jp_can0_3d + 1) = zco3(:,:,:)              * tmask(:,:,:)
         trc3d(:,:,:,jp_can0_3d + 2) = aksp(:,:,:) / calcon     * tmask(:,:,:)
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zco3, zcaldiss )
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_lys')
      !
   END SUBROUTINE can_lys

   SUBROUTINE can_lys_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_lys_init  ***
      !!
      !! ** Purpose :   Initialization of CaCO3 dissolution parameters
      !!
      !! ** Method  :   Read the nampiscal namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiscal
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::  ios                 ! Local integer output status for namelist read

      NAMELIST/nampiscal/ kdca, nca
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampiscal in reference namelist : Pisces CaCO3 dissolution
      READ  ( numnatp_ref, nampiscal, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscal in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiscal in configuration namelist : Pisces CaCO3 dissolution
      READ  ( numnatp_cfg, nampiscal, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscal in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiscal )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for CaCO3 dissolution, nampiscal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    diss. rate constant calcite (s^-1)   kdca      =', kdca
         WRITE(numout,*) '    order of reaction for calcite dissolution nca       =', nca
      ENDIF

      !
   END SUBROUTINE can_lys_init

   !!======================================================================
END MODULE  canlys
