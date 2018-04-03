MODULE canprod
   !!======================================================================
   !!                         ***  MODULE canprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   can_prod_init  :   Initialization of the parameters for growth
   !!   can_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE canopt          !  optical model
!   USE canlim          !  Co-limitations of differents nutrients
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_prod         ! called in canbio.F90
   PUBLIC   can_prod_init    ! called in trcsms_canoe.F90
   PUBLIC   can_prod_alloc

   !! * Shared module variables
   !! these are hardwired parameters
   REAL(wp), PUBLIC ::  mw_c       = 12._wp            !:
   REAL(wp), PUBLIC ::  mw_n       = 14._wp            !:
   REAL(wp), PUBLIC ::  mw_fe      = 55.845_wp         !:
   !! these are input parameters in namelist_pisces
   REAL(wp), PUBLIC ::  QNmax1     = 0.172_wp          !: Small phytoplankton max N quota
   REAL(wp), PUBLIC ::  QNmin1     = 0.04_wp           !: Small phytoplankton min N quota
   REAL(wp), PUBLIC ::  QNmax2     = 0.172_wp          !: Large phytoplankton max N quota
   REAL(wp), PUBLIC ::  QNmin2     = 0.04_wp           !: Large phytoplankton min N quota
   REAL(wp), PUBLIC ::  VCNref     = 0.6_wp            !: Reference rate of N uptake
   REAL(wp), PUBLIC ::  QFemax1    = 93.075_wp         !: Small phytoplankton max Fe quota
   REAL(wp), PUBLIC ::  QFemin1    = 4.65_wp           !: Small phytoplankton min Fe quota
   REAL(wp), PUBLIC ::  QFemax2    = 69.8063_wp        !: Large phytoplankton max Fe quota
   REAL(wp), PUBLIC ::  QFemin2    = 4.65_wp           !: Large phytoplankton min Fe quota
   REAL(wp), PUBLIC ::  VCFref     = 79._wp            !: Reference rate of Fe uptake
   REAL(wp), PUBLIC ::  PCref      = 3._wp             !: Reference rate of photosynthesis
   REAL(wp), PUBLIC ::  alphachl   = 1.08_wp           !: Initial slope of P-E curve
   REAL(wp), PUBLIC ::  kn1        = 0.1_wp            !: Small P half-saturation for NO3 uptake
   REAL(wp), PUBLIC ::  ka1        = 0.05_wp           !: Small P half-saturation for NH4 uptake
   REAL(wp), PUBLIC ::  kf1        = 100._wp           !: Small P half-saturation for Fe uptake
   REAL(wp), PUBLIC ::  kn2        = 0.5_wp            !: Large P half-saturation for NO3 uptake
   REAL(wp), PUBLIC ::  ka2        = 0.05_wp           !: Large P half-saturation for NH4 uptake
   REAL(wp), PUBLIC ::  kf2        = 200._wp           !: Large P half-saturation for Fe uptake
   REAL(wp), PUBLIC ::  thetamax   = 0.18_wp           !: Maximum chlorophyll/nitrogen ratio
   REAL(wp), PUBLIC ::  eta        = 2._wp             !: Metabolic cost of biosynthesis
   REAL(wp), PUBLIC ::  kexh       = 1.7_wp            !: exhudation of excess intracellular C


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prmax    !: optimal production = f(temperature)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   REAL(wp) :: tpp                    !: Total primary production
   
   REAL(wp) :: r1_rday                !: 1 / rday
   REAL(wp) :: texcret                !: 1 - excret 
   REAL(wp) :: texcret2               !: 1 - excret2        


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
   !! $Id: canprod.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zfact  
      REAL(wp) ::    ztn 
      REAL(wp) ::    zproreg, zproreg2
      REAL(wp) ::   zrfact2
! some local variables required by the revised model
      REAL(wp) :: phyc,phyn,phyfe,chl,Ni,Na,Fe,Tf,QN,qndep,qfedep,VCNmax,Alim,Nlim,VCN,QFe,VCFmax,VCF
      REAL(wp) :: PCmax,thetac,PCphot,rhochl,ei,xsphsyn, mwr_n2c, imw_n
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprdia, zprbio, zprdch, zprnch, zysopt,zw3d   
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprorcad, zprofed, zprofen, zpronew, zpronewd
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprocn, zprocd, zpronn, zprond
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_prod')
      !
      !  Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, zprdia, zprbio, zprdch, zprnch, zysopt            ) 
      CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprorcad, zprofed, zprofen, zpronew, zpronewd )
      CALL wrk_alloc( jpi, jpj, jpk, zprocn, zprocd, zpronn, zprond                                             ) 
      !
      zprorca (:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
      zprochld(:,:,:) = 0._wp
      zpronew (:,:,:) = 0._wp
      zpronewd(:,:,:) = 0._wp
      zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zprdch  (:,:,:) = 0._wp
      zprnch  (:,:,:) = 0._wp
      zysopt  (:,:,:) = 0._wp
      zprocn  (:,:,:) = 0._wp
      zprocd  (:,:,:) = 0._wp
      zpronn  (:,:,:) = 0._wp
      zprond  (:,:,:) = 0._wp
! precalculate some constant terms to minimize divisions
      mwr_n2c=mw_n/mw_c
      imw_n=1./mw_n
      !
      ! Computation of the various production terms 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot(ji,jj,jk) > 1.E-3 ) THEN

                      ztn = tsn(ji,jj,jk,jp_tem)
                      phyc = MAX(trb(ji,jj,jk,jpphy),0.)*mw_c
                      phyn = MAX(trb(ji,jj,jk,jpnn),0.)*mw_n
                      phyfe = MAX(trb(ji,jj,jk,jpnfe),0.)*mw_fe
                      chl = MAX(trb(ji,jj,jk,jpnch),0.)
                      Ni = MAX(trb(ji,jj,jk,jpno3),0.)
                      Na = MAX(trb(ji,jj,jk,jpnh4),0.)
                      Fe = MAX(trb(ji,jj,jk,jpfer),0.)                        ! Fe variables are in nmol m^-3, others in mmol m^-3
                      ei = etot(ji,jj,jk)*4.15                        ! convert to umol m^-2 s^-1

! this is modified from ~/mexfiles/vrm/fwd/bsource_vrm.f via chemo_2P2Z_gmk.F
! small phytoplankton

                      Tf = tgfuncp(ji,jj,jk)
                      QN = MIN(QNmax1,phyn/(phyc+rtrn))
                      QN = MAX(QNmin1,QN)
                      qndep = MAX((QNmax1-QN)/(QNmax1-QNmin1),0.)                  ! in principle this should be nonegative but if roundoff makes it even slightly negative the exponent could go NaN
                      VCNmax = VCNref*Tf*qndep**0.05
                      Alim  =  Na/(ka1+Na)
                      Nlim  =  Ni/(kn1+Ni)
                      VCN = VCNmax*((1.-Alim)*Nlim+Alim)

                      QFe = MIN(QFemax1,phyfe/(phyc+rtrn))
                      QFe = MAX(QFemin1,QFe)
                      qfedep = MAX((QFemax1-QFe)/(QFemax1-QFemin1),0.) 
                      VCFmax = VCFref*Tf*qfedep**0.05
                      VCF = VCFmax*Fe/(kf1+Fe)

                      !PCmax = PCref*Tf*MIN((QFe-QFemin1+rtrn*1.e6)/(QFemax1-QFemin1),(QN-QNmin1+rtrn)/(QNmax1-QNmin1))
                      PCmax = PCref*Tf*(QN-QNmin1+rtrn)/(QNmax1-QNmin1)

                      PCmax = MAX(PCmax,1.0e-10)
                      thetac = MAX(chl/(phyc+rtrn),0.001)
                      PCphot = PCmax*(1.-EXP(-alphachl*ei*thetac/PCmax))
                      rhochl = thetamax*(PCphot/(alphachl*thetac*MAX(ei,0.001)))

! calculate excess intracellular C for exhudation
                      xsphsyn=(phyc/(phyn+rtrn)*mwr_n2c-rr_c2n)*phyn*imw_n
                      xsphsyn=MAX(xsphsyn,0.)

                      zprocn(ji,jj,jk) = (PCphot-eta*VCN)*trb(ji,jj,jk,jpphy)*xstep-kexh*xsphsyn*xstep        ! C production rate (in molar units)
                      zpronn(ji,jj,jk) = VCN/QN*trb(ji,jj,jk,jpnn)*xstep                                      ! N uptake rate
                      zprofen(ji,jj,jk) = VCF/QFe*trb(ji,jj,jk,jpnfe)*xstep                                   ! Fe uptake rate
                      zprochln(ji,jj,jk) = rhochl*VCN/thetac*trb(ji,jj,jk,jpnch)*xstep                        ! Chl production rate
                      zpronew(ji,jj,jk) = zpronn(ji,jj,jk)*Nlim/(Alim+Nlim+rtrn)                              ! NO3 uptake
                      xlimnn(ji,jj,jk) = 1.-qndep 
                      xlimnfe(ji,jj,jk) = 1.-qfedep 

! large phytoplankton

                      phyc = MAX(trb(ji,jj,jk,jpdia),0.)*mw_c
                      phyn = MAX(trb(ji,jj,jk,jpdn),0.)*mw_n
                      phyfe = MAX(trb(ji,jj,jk,jpdfe),0.)*mw_fe
                      chl = MAX(trb(ji,jj,jk,jpdch),0.)

                      QN = MIN(QNmax2,phyn/(phyc+rtrn))
                      QN = MAX(QNmin2,QN)
                      qndep = MAX((QNmax2-QN)/(QNmax2-QNmin2),0.) 
                      VCNmax = VCNref*Tf*qndep**0.05
                      Alim  =  Na/(ka2+Na)
                      Nlim  =  Ni/(kn2+Ni)
                      VCN = VCNmax*((1.-Alim)*Nlim+Alim)

                      QFe = MIN(QFemax2,phyfe/(phyc+rtrn))
                      QFe = MAX(QFemin2,QFe)
                      qfedep = MAX((QFemax2-QFe)/(QFemax2-QFemin2),0.) 
                      VCFmax = VCFref*Tf*qfedep**0.05
                      VCF = VCFmax*Fe/(kf2+Fe)

                      !PCmax = PCref*Tf*MIN((QFe-QFemin2+rtrn)/(QFemax2-QFemin2),(QN-QNmin2+rtrn)/(QNmax2-QNmin2))
                      !PCmax = PCref*Tf*MIN((QFe-QFemin2+rtrn)/(QFemax2-QFemin2),(QN-QNmin2+rtrn)/(QNmax2-QNmin2))
                      PCmax = PCref*Tf*(QN-QNmin1+rtrn)/(QNmax1-QNmin1)
                      PCmax = MAX(PCmax,1.0e-10)
                      thetac = MAX(chl/(phyc+rtrn),0.001)
                      PCphot = PCmax*(1.-EXP(-alphachl*ei*thetac/PCmax))
                      rhochl = thetamax*(PCphot/(alphachl*thetac*MAX(ei,0.001)))

                      xsphsyn=(phyc/(phyn+rtrn)*mwr_n2c-rr_c2n)*phyn*imw_n
                      xsphsyn=MAX(xsphsyn,0.)

                      zprocd(ji,jj,jk) = (PCphot-eta*VCN)*trb(ji,jj,jk,jpdia)*xstep-kexh*xsphsyn*xstep       ! C production rate (in molar units)
                      zprond(ji,jj,jk) = VCN/QN*trb(ji,jj,jk,jpdn)*xstep                                     ! N uptake rate
                      zprofed(ji,jj,jk) = VCF/QFe*trb(ji,jj,jk,jpdfe)*xstep                                  ! Fe uptake rate
                      zprochld(ji,jj,jk) = rhochl*VCN/thetac*trb(ji,jj,jk,jpdch)*xstep                       ! Chl production rate
                      zpronewd(ji,jj,jk) = zprond(ji,jj,jk)*Nlim/(Alim+Nlim+rtrn)                            ! NO3 uptake
                      xlimdn(ji,jj,jk) = 1.-qndep 
                      xlimdfe(ji,jj,jk) = 1.-qfedep 

               ENDIF
            END DO
         END DO
      END DO
      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zproreg  = zpronn(ji,jj,jk) - zpronew(ji,jj,jk)
              zproreg2 = zprond(ji,jj,jk) - zpronewd(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronew(ji,jj,jk) - zpronewd(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprocn(ji,jj,jk)
              tra(ji,jj,jk,jpnn) = tra(ji,jj,jk,jpnn) + zpronn(ji,jj,jk)
              tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln(ji,jj,jk) 
              tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk)
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprocd(ji,jj,jk) 
              tra(ji,jj,jk,jpdn) = tra(ji,jj,jk,jpdn) + zprond(ji,jj,jk) 
              tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld(ji,jj,jk) 
              tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) 
! O2 production equals DIC reduction + an additional nitrate term based on Laws 1991; this term is set to conserve O2 globally at steady state, i.e. 0.301887 = 2/rr_c2n where 2 mol O2 / mol N is the O2 sink to nitrification
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + (zprocn(ji,jj,jk) + zprocd(ji,jj,jk)) &
                 &                    + 0.301887 * (zpronew(ji,jj,jk) + zpronewd(ji,jj,jk)) * rr_c2n   
              tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zprofen(ji,jj,jk) - zprofed(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - (zprocn(ji,jj,jk) + zprocd(ji,jj,jk))*1.E-6
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + (zpronew(ji,jj,jk)+zpronewd(ji,jj,jk))*1.E-6 &
                 &                                      - (zproreg+zproreg2)*1.E-6
          END DO
        END DO
     END DO


     ! Total primary production per year
     tpp = tpp + glob_sum( ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )
     IF( kt == nitend .AND. knt == nrdttrc ) THEN
        WRITE(numout,*) 'Total PP (Gtc) :'
        WRITE(numout,*) '-------------------- : ',tpp * 12. / 1.E12
        WRITE(numout,*) 
      ENDIF

    ! Total primary production per year
  !  IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
  !       & tpp = glob_sum( ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )

    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          CALL wrk_alloc( jpi, jpj, jpk, zw3d )
          zfact = 1.e-3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHY" ) .OR. iom_use( "PPPHY2" ) )  THEN
              zw3d(:,:,:) = zprocn (:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHY"  , zw3d )
              !
              zw3d(:,:,:) = zprocd(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHY2"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) )  THEN
              zw3d(:,:,:) = zpronew (:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatomes
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by  diatomes
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LNN" ) .OR. iom_use( "LDN" ) )  THEN
              zw3d(:,:,:) = xlimnn (:,:,:) * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "LNN"  , zw3d )
              !
              zw3d(:,:,:) = xlimdn(:,:,:) *  tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "LDN"  , zw3d )
          ENDIF
          IF( iom_use( "LNFe" ) .OR. iom_use( "LDFe" ) )  THEN
              zw3d(:,:,:) = xlimnfe(:,:,:) * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "LNFe"  , zw3d )
              !
              zw3d(:,:,:) = xlimdfe(:,:,:) *  tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "LDFe"  , zw3d )
          ENDIF
          CALL wrk_dealloc( jpi, jpj, jpk, zw3d )
       ENDIF
     ELSE
        IF( ln_diatrc ) THEN
           zfact = 1.e-3 * rfact2r
           trc3d(:,:,:,jp_can0_3d + 4)  = zprorca (:,:,:) * zfact * tmask(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 5)  = zprorcad(:,:,:) * zfact * tmask(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 6)  = zpronew (:,:,:) * zfact * tmask(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 7)  = zpronewd(:,:,:) * zfact * tmask(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 8)  = zprorcad(:,:,:) * zfact * tmask(:,:,:) * zysopt(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 9)  = zprofed (:,:,:) * zfact * tmask(:,:,:)
           trc3d(:,:,:,jp_can0_3d + 10) = zprofen (:,:,:) * zfact * tmask(:,:,:)
        ENDIF
     ENDIF

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
     !
      CALL wrk_dealloc( jpi, jpj, jpk, zprdia, zprbio, zprdch, zprnch, zysopt            ) 
      CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprorcad, zprofed, zprofen, zpronew, zpronewd )
      CALL wrk_dealloc( jpi, jpj, jpk, zprocn, zprocd, zpronn, zprond                                             ) 
     !
     IF( nn_timing == 1 )  CALL timing_stop('can_prod')
     !
   END SUBROUTINE can_prod


   SUBROUTINE can_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE can_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      NAMELIST/nampisprod/ QNmax1, QNmin1, QNmax2, QNmin2, VCNref, QFemax1,  &
         &                 QFemin1, QFemax2, QFemin2, VCFref, PCref, alphachl,     &
         &                 kn1, ka1, kf1, kn2, ka2, kf2, thetamax, eta, kexh
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, nampisprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, nampisprod, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Small phytoplankton max N quota          QNmax1       =', QNmax1
         WRITE(numout,*) '    Small phytoplankton min N quota          QNmin1       =', QNmin1
         WRITE(numout,*) '    Large phytoplankton max N quota          QNmax2       =', QNmax2
         WRITE(numout,*) '    Large phytoplankton min N quota          QNmin2       =', QNmin2
         WRITE(numout,*) '    Reference rate of N uptake               VCNref       =', VCNref
         WRITE(numout,*) '    Small phytoplankton max Fe quota         QFemax1      =', QFemax1
         WRITE(numout,*) '    Small phytoplankton min Fe quota         QFemin1      =', QFemin1
         WRITE(numout,*) '    Large phytoplankton max Fe quota         QFemax2      =', QFemax2
         WRITE(numout,*) '    Large phytoplankton min Fe quota         QFemin2      =', QFemin2
         WRITE(numout,*) '    Reference rate of Fe uptake              VCFref       =', VCFref
         WRITE(numout,*) '    Reference rate of photosynthesis         PCref        =', PCref
         WRITE(numout,*) '    Initial slope of P-E curve               alphachl     =', alphachl
         WRITE(numout,*) '    Small P half-saturation for NO3 uptake   kn1          =', kn1
         WRITE(numout,*) '    Small P half-saturation for NH4 uptake   ka1          =', ka1
         WRITE(numout,*) '    Small P half-saturation for Fe uptake    kf1          =', kf1
         WRITE(numout,*) '    Large P half-saturation for NO3 uptake   kn2          =', kn2
         WRITE(numout,*) '    Large P half-saturation for NH4 uptake   ka2          =', ka2
         WRITE(numout,*) '    Large P half-saturation for Fe uptake    kf2          =', kf2
         WRITE(numout,*) '    Maximum chlorophyll/carbon ratio         thetamax     =', thetamax
         WRITE(numout,*) '    Metabolic cost of biosynthesis           eta          =', eta
         WRITE(numout,*) '    Exudation rate of excess intracellular C kexh         =', kexh
      ENDIF
      !
      !
      !r1_rday   = 1._wp / rday 
      tpp       = 0._wp
      !
   END SUBROUTINE can_prod_init


   INTEGER FUNCTION can_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE can_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( prmax(jpi,jpj,jpk), quotan(jpi,jpj,jpk), quotad(jpi,jpj,jpk), STAT = can_prod_alloc )
      !
      IF( can_prod_alloc /= 0 ) CALL ctl_warn('can_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION can_prod_alloc


   !!======================================================================
END MODULE  canprod
