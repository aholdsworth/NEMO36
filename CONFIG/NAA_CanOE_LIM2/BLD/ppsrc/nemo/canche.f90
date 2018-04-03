MODULE canche
   !!======================================================================
   !!                         ***  MODULE canche  ***
   !! TOP :   CanOE Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_che      :  Sea water chemistry computed following OCMIP protocol
   !!----------------------------------------------------------------------
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE sms_canoe    !  CanOE Source Minus Sink variables
   USE lib_mpp       !  MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_che         !
   PUBLIC   can_che_alloc   !

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   chemc    ! Solubilities of O2 and CO2

   REAL(wp), PUBLIC ::   atcox  = 0.20946         ! units atm

   REAL(wp) ::   salchl = 1. / 1.80655    ! conversion factor for salinity --> chlorinity (Wooster et al. 1969)
   REAL(wp) ::   o2atm  = 1. / ( 1000. * 0.20946 )  

   REAL(wp) ::   akcc1  = -171.9065       ! coeff. for apparent solubility equilibrium
   REAL(wp) ::   akcc2  =   -0.077993     ! Millero et al. 1995 from Mucci 1983
   REAL(wp) ::   akcc3  = 2839.319        
   REAL(wp) ::   akcc4  =   71.595        
   REAL(wp) ::   akcc5  =   -0.77712      
   REAL(wp) ::   akcc6  =    0.00284263   
   REAL(wp) ::   akcc7  =  178.34        
   REAL(wp) ::   akcc8  =   -0.07711     
   REAL(wp) ::   akcc9  =    0.0041249   

   REAL(wp) ::   rgas   = 83.143         ! universal gas constants
   REAL(wp) ::   oxyco  = 1. / 22.4144   ! converts from liters of an ideal gas to moles

   REAL(wp) ::   bor1   = 0.00023        ! borat constants
   REAL(wp) ::   bor2   = 1. / 10.82

   REAL(wp) ::   ca0    = -162.8301      ! WEISS & PRICE 1980, units mol/(kg atm)
   REAL(wp) ::   ca1    =  218.2968
   REAL(wp) ::   ca2    =   90.9241
   REAL(wp) ::   ca3    =   -1.47696
   REAL(wp) ::   ca4    =    0.025695
   REAL(wp) ::   ca5    =   -0.025225
   REAL(wp) ::   ca6    =    0.0049867

   REAL(wp) ::   c10    = -3633.86        ! Coeff. for 1. dissoc. of carbonic acid (Dickson et al., 2007)
   REAL(wp) ::   c11    =    61.2172    
   REAL(wp) ::   c12    =    -9.67770  
   REAL(wp) ::   c13    =     0.011555     
   REAL(wp) ::   c14    =    -0.0001152

   REAL(wp) ::   c20    =  -471.78       ! coeff. for 2. dissoc. of carbonic acid (Dickson et al., 2007)   
   REAL(wp) ::   c21    =   -25.9290   
   REAL(wp) ::   c22    =     3.16967   
   REAL(wp) ::   c23    =     0.01781
   REAL(wp) ::   c24    =    -0.0001122


   REAL(wp) ::   cb0    = -8966.90      ! Coeff. for 1. dissoc. of boric acid 
   REAL(wp) ::   cb1    = -2890.53      ! (Dickson and Goyet, 1994)
   REAL(wp) ::   cb2    =   -77.942
   REAL(wp) ::   cb3    =     1.728
   REAL(wp) ::   cb4    =    -0.0996
   REAL(wp) ::   cb5    =   148.0248
   REAL(wp) ::   cb6    =   137.1942
   REAL(wp) ::   cb7    =     1.62142
   REAL(wp) ::   cb8    =   -24.4344
   REAL(wp) ::   cb9    =   -25.085
   REAL(wp) ::   cb10   =    -0.2474 
   REAL(wp) ::   cb11   =     0.053105

   REAL(wp) ::   cw0    = -13847.26     ! Coeff. for dissoc. of water (Dickson and Riley, 1979 )
   REAL(wp) ::   cw1    =    148.9652  
   REAL(wp) ::   cw2    =    -23.6521
   REAL(wp) ::   cw3    =    118.67 
   REAL(wp) ::   cw4    =     -5.977 
   REAL(wp) ::   cw5    =      1.0495  
   REAL(wp) ::   cw6    =     -0.01615

   !                                    ! volumetric solubility constants for o2 in ml/L  
   REAL(wp) ::   ox0    =  2.00856      ! from Table 1 for Eq 8 of Garcia and Gordon, 1992.
   REAL(wp) ::   ox1    =  3.22400      ! corrects for moisture and fugacity, but not total atmospheric pressure
   REAL(wp) ::   ox2    =  3.99063      !      Original CanOE code noted this was a solubility, but 
   REAL(wp) ::   ox3    =  4.80299      ! was in fact a bunsen coefficient with units L-O2/(Lsw atm-O2)
   REAL(wp) ::   ox4    =  9.78188e-1   ! Hence, need to divide EXP( zoxy ) by 1000, ml-O2 => L-O2
   REAL(wp) ::   ox5    =  1.71069      ! and atcox = 0.20946 to add the 1/atm dimension.
   REAL(wp) ::   ox6    = -6.24097e-3   
   REAL(wp) ::   ox7    = -6.93498e-3 
   REAL(wp) ::   ox8    = -6.90358e-3
   REAL(wp) ::   ox9    = -4.29155e-3 
   REAL(wp) ::   ox10   = -3.11680e-7 

   !                                    ! coeff. for seawater pressure correction : millero 95
   !                                    ! AGRIF doesn't like the DATA instruction
   REAL(wp) :: devk11  = -25.5
   REAL(wp) :: devk12  = -15.82
   REAL(wp) :: devk13  = -29.48
   REAL(wp) :: devk14  = -25.60
   REAL(wp) :: devk15  = -48.76
   !
   REAL(wp) :: devk21  = 0.1271
   REAL(wp) :: devk22  = -0.0219
   REAL(wp) :: devk23  = 0.1622
   REAL(wp) :: devk24  = 0.2324
   REAL(wp) :: devk25  = 0.5304
   !
   REAL(wp) :: devk31  = 0.
   REAL(wp) :: devk32  = 0.
   REAL(wp) :: devk33  = 2.608E-3
   REAL(wp) :: devk34  = -3.6246E-3
   REAL(wp) :: devk35  = 0.
   !
   REAL(wp) :: devk41  = -3.08E-3
   REAL(wp) :: devk42  = 1.13E-3
   REAL(wp) :: devk43  = -2.84E-3
   REAL(wp) :: devk44  = -5.13E-3
   REAL(wp) :: devk45  = -11.76E-3
   !
   REAL(wp) :: devk51  = 0.0877E-3
   REAL(wp) :: devk52  = -0.1475E-3     
   REAL(wp) :: devk53  = 0.
   REAL(wp) :: devk54  = 0.0794E-3      
   REAL(wp) :: devk55  = 0.3692E-3      

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
   !! $Id: canche.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_che
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_che  ***
      !!
      !! ** Purpose :   Sea water chemistry computed following OCMIP protocol
      !!
      !! ** Method  : - ...
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   ztkel, zt   , zt2   , zsal  , zsal2 , zbuf1 , zbuf2
      REAL(wp) ::   ztgg , ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zpres, ztc  , zcl   , zcpexp, zoxy  , zcpexp2
      REAL(wp) ::   zsqrt, ztr  , zlogt , zcek1
      REAL(wp) ::   zis  , zis2 , zsal15, zisqrt
      REAL(wp) ::   zckb , zck1 , zck2  , zckw  , zak1 , zak2  , zakb , zaksp0, zakw
      REAL(wp) ::    zaksp1
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_che')
      !
      ! CHEMICAL CONSTANTS - SURFACE LAYER
      ! ----------------------------------
!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !                             ! SET ABSOLUTE TEMPERATURE
            ztkel = tsn(ji,jj,1,jp_tem) + 273.15
            zt    = ztkel * 0.01
            zt2   = zt * zt
            zsal  = tsn(ji,jj,1,jp_sal) + ( 1.- tmask(ji,jj,1) ) * 35.
            zsal2 = zsal * zsal
            zlogt = LOG( zt )
            !                             ! LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1980)
            !                             !     AND FOR THE ATMOSPHERE FOR NON IDEAL GAS
            zcek1 = ca0 + ca1 / zt + ca2 * zlogt + ca3 * zt2 + zsal * ( ca4 + ca5 * zt + ca6 * zt2 )
            !                             ! LN(K0) OF SOLUBILITY OF O2 and N2 in ml/L (EQ. 8, GARCIA AND GORDON, 1992)
            ztgg  = LOG( ( 298.15 - tsn(ji,jj,1,jp_tem) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
            ztgg2 = ztgg  * ztgg
            ztgg3 = ztgg2 * ztgg
            ztgg4 = ztgg3 * ztgg
            ztgg5 = ztgg4 * ztgg
            zoxy  = ox0 + ox1 * ztgg + ox2 * ztgg2 + ox3 * ztgg3 + ox4 * ztgg4 + ox5 * ztgg5   &
                   + zsal * ( ox6 + ox7 * ztgg + ox8 * ztgg2 + ox9 * ztgg3 ) +  ox10 * zsal2

            !                             ! SET SOLUBILITIES OF O2 AND CO2 
            chemc(ji,jj,1) = EXP( zcek1 ) * 1.e-6 * rhop(ji,jj,1) / 1000.  ! mol/(L uatm)
            chemc(ji,jj,2) = ( EXP( zoxy  ) * o2atm ) * oxyco              ! mol/(L atm)
            !
         END DO
      END DO


      ! CHEMICAL CONSTANTS - DEEP OCEAN
      ! -------------------------------
!CDIR NOVERRCHK
      DO jk = 1, jpk
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! SET PRESSION
               zpres   = 1.025e-1 * gdept_n(ji,jj,jk)

               ! SET ABSOLUTE TEMPERATURE
               ztkel   = tsn(ji,jj,jk,jp_tem) + 273.15
               zsal    = tsn(ji,jj,jk,jp_sal) + ( 1.-tmask(ji,jj,jk) ) * 35.
               zsqrt  = SQRT( zsal )
               zsal15  = zsqrt * zsal
               zlogt  = LOG( ztkel )
               ztr    = 1. / ztkel
               zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
               zis2   = zis * zis
               zisqrt = SQRT( zis )
               ztc     = tsn(ji,jj,jk,jp_tem) + ( 1.- tmask(ji,jj,jk) ) * 20.

               ! CHLORINITY (WOOSTER ET AL., 1969)
               zcl     = zsal * salchl


               ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
               zckb    = ( cb0 + cb1 * zsqrt + cb2  * zsal + cb3 * zsal15 + cb4 * zsal * zsal ) * ztr   &
                  &    + ( cb5 + cb6 * zsqrt + cb7  * zsal )                                            &
                  &    + ( cb8 + cb9 * zsqrt + cb10 * zsal ) * zlogt + cb11 * zsqrt * ztkel             

               zck1    = c10 * ztr + c11 + c12 * zlogt + c13 * zsal + c14 * zsal * zsal
               zck2    = c20 * ztr + c21 + c22 * zlogt + c23 * zsal + c24 * zsal * zsal

               ! PKW (H2O) (DICKSON AND RILEY, 1979)
               zckw    = cw0 * ztr + cw1 + cw2 * zlogt + ( cw3 * ztr + cw4 + cw5 * zlogt ) * zsqrt + cw6 * zsal


               ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
               !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
               zaksp0  = akcc1 + akcc2 * ztkel + akcc3 * ztr + akcc4 * LOG10( ztkel )   &
                  &   + ( akcc5 + akcc6 * ztkel + akcc7 * ztr ) * zsqrt + akcc8 * zsal + akcc9 * zsal15

               ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
               zak1    = 10**(zck1)
               zak2    = 10**(zck2)
               zakb    = EXP( zckb  )
               zakw    = EXP( zckw )
               zaksp1  = 10**(zaksp0)

               ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
               !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
               !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
               !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
               !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
               !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
               !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
               !        & GIESKES (1970), P. 1285-1286 (THE SMALL
               !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
               !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))
               zcpexp  = zpres /(rgas*ztkel)
               zcpexp2 = zpres * zpres/(rgas*ztkel)

               ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
               !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
               !        (CF. BROECKER ET AL., 1982)

               zbuf1  = -     ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
               zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
               ak13(ji,jj,jk) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
               zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
               ak23(ji,jj,jk) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
               zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
               akb3(ji,jj,jk) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk14 + devk24 * ztc + devk34 * ztc * ztc )
               zbuf2  = 0.5 * ( devk44 + devk54 * ztc )
               akw3(ji,jj,jk) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )


               ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE 
               !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
               !        (P. 1285) AND BERNER (1976)
               zbuf1  =     - ( devk15 + devk25 * ztc + devk35 * ztc * ztc )
               zbuf2  = 0.5 * ( devk45 + devk55 * ztc )
               aksp(ji,jj,jk) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               ! TOTAL BORATE CONCENTR. [MOLES/L]
               borat(ji,jj,jk) = bor1 * zcl * bor2


            END DO
         END DO
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('can_che')
      !
   END SUBROUTINE can_che


   INTEGER FUNCTION can_che_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_che_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( chemc(jpi,jpj,2), STAT=can_che_alloc )
      !
      IF( can_che_alloc /= 0 )   CALL ctl_warn('can_che_alloc : failed to allocate arrays.')
      !
   END FUNCTION can_che_alloc


   !!======================================================================
END MODULE  canche
