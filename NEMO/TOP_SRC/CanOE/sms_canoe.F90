MODULE sms_canoe   
   !!----------------------------------------------------------------------
   !!                     ***  sms_canoe.F90  ***  
   !! TOP :   CanOE Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
    !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------
#if defined key_canoe  
   !!----------------------------------------------------------------------
   !!   'key_canoe'                                         CanOE model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist canoe
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist canoe
   INTEGER ::   numonp      = -1           !! Logical unit for namelist canoe output

   !!*  Biological fluxes for light : variables shared by canoe & lobster
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   !

   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology
   REAL(wp) ::   ryyss             !: number of seconds per year 
   REAL(wp) ::   r1_ryyss          !: inverse number of seconds per year 


   !!*  Biological parameters 
   REAL(wp) ::   rr_c2n, rr_fe2n
   REAL(wp) ::   rr_n2c, rr_fe2c, rr_c2fe, rr_n2fe
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   wsbioc            !: ???

   REAL(wp) ::   rno3              !: ???

   !!* restoring
   LOGICAL  ::  ln_pisdmp          !: restoring or not of nutrients to a mean value
   INTEGER  ::  nn_pisdmp          !: frequency of relaxation or not of nutrients to a mean value
   LOGICAL  ::   ln_pisclo         !: Restoring or not of nutrients to initial value

   !!* Mass conservation
   LOGICAL  ::  ln_check_mass      !: Flag to check mass conservation

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdn    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnn    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: ???


   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing1   !: microzooplankton grazing
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing2   !: mesozooplankton grazing on phytoplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing3   !: mesozooplankton grazing on microzooplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr     !: denitrification
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zdnf       !: N2 fixation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitnh4   !: annamox
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nh4ox      !: nitrification
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpn    !: nanophytoplankton mortality 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpd    !: diatoms mortality
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprochln   !: nanophytoplankton chl production 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprochld   !: nanophytoplankton chl production    
   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???
   
   !!* Temperature dependence of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncp   !: Temp. dependence of phytoplankton growth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncz   !: Temp. dependence of microzooplankton resp
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncz2   !: Temp. dependence of mesozooplankton resp
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncr   !: Temp. dependence of remineralization


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sms_canoe.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_canoe_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_canoe_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(5)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for canoe & lobster
      ALLOCATE(  neln(jpi,jpj), heup(jpi,jpj),  STAT=ierr(1) )
      !
      !*  Biological fluxes for primary production
      ALLOCATE(  xlimnfe (jpi,jpj,jpk), xlimdfe (jpi,jpj,jpk),       &
         &      xlimnn(jpi,jpj,jpk), xlimdn(jpi,jpj,jpk),  STAT=ierr(2) ) 
         !
      !*  SMS for the organic matter
      ALLOCATE( nitrfac(jpi,jpj,jpk),  denitr(jpi,jpj,jpk),        &
         &      denitnh4(jpi,jpj,jpk), nh4ox(jpi,jpj,jpk),         &
         &      zdnf(jpi,jpj,jpk),     prodcal(jpi,jpj,jpk),       & 
         &      grazing1(jpi,jpj,jpk), grazing2(jpi,jpj,jpk),      &
         &      grazing3(jpi,jpj,jpk), xdiss  (jpi,jpj,jpk),   STAT=ierr(3) )  
         !

      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)    , ak13  (jpi,jpj,jpk) ,       &
         &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,       &
         &      akw3(jpi,jpj,jpk)    , borat (jpi,jpj,jpk) ,       &
         &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,     STAT=ierr(4) )
         !
      !* Temperature dependancy of SMS terms
      ALLOCATE( tgfuncp(jpi,jpj,jpk)  , tgfuncz(jpi,jpj,jpk) ,     &
         &      tgfuncz2(jpi,jpj,jpk) , tgfuncr(jpi,jpj,jpk) , STAT=ierr(5) )
         !
      ALLOCATE( zmortpn(jpi,jpj,jpk)  , zmortpd(jpi,jpj,jpk) ,     &
         &      zprochln(jpi,jpj,jpk) , zprochld(jpi,jpj,jpk) , STAT=ierr(7) )
      sms_canoe_alloc = MAXVAL( ierr )
      !
      IF( sms_canoe_alloc /= 0 )   CALL ctl_warn('sms_canoe_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_canoe_alloc

#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO CanOE model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================   
END MODULE sms_canoe    
