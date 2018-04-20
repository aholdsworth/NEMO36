MODULE sms_cmoc   
   !!----------------------------------------------------------------------
   !!                     ***  sms_cmoc.F90  ***  
   !! TOP :   CMOC Source Minus Sink variables
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
#if defined key_cmoc  
   !!----------------------------------------------------------------------
   !!   'key_cmoc'                                         CMOC model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp_ref = -1           !! Logical units for namelist cmoc
   INTEGER ::   numnatp_cfg = -1           !! Logical units for namelist cmoc
   INTEGER ::   numonp      = -1           !! Logical unit for namelist cmoc output

   !!*  Biological fluxes for light : variables shared by cmoc & lobster
!   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln  !: number of T-levels + 1 in the euphotic layer
!   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup  !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot  !: par (photosynthetic available radiation)
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  LOBSTER : zooplakton closure
   !                                                       !:  CMOC  : silicon dependant half saturation

   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology
   REAL(wp) ::   ryyss             !: number of seconds per year 
   REAL(wp) ::   r1_ryyss          !: inverse number of seconds per year 


   !!*  Biological parameters 
   REAL(wp) ::   rno3              !: ???
   REAL(wp) ::   o2ut              !: ???
   REAL(wp) ::   po4r              !: ???
   REAL(wp) ::   rdenit            !: ???
   REAL(wp) ::   rdenita           !: ???
   REAL(wp) ::   o2nit             !: ???
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   xkmort            !: ???
   ! <CMOC code OR 10/21/2015> CMOC block start
   !!*  Remineralization
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:  ) ::   oomask         ! Open ocean mask 
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::  redet          !: detritus remineralization
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)    ::  redettot       !: detritus remineralization integrated below the euphotic zone
   
   !!*  Biological fluxes for primary production
   !  iron limitation mask
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xlimnfecmoc!:
   ! <  rain ratio 
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xrcico     !:

   !!*  CMOC model parameters
   !  Phytoplankton Growth
   !REAL(wp)    :: apar_cmoc  
   !REAL(wp)    :: kw_cmoc    
   !REAL(wp)    :: kchl_cmoc  
   REAL(wp)    :: achl_cmoc  
   REAL(wp)    :: thm_cmoc   
   REAL(wp)    :: tau_cmoc   
   REAL(wp)    :: itau_cmoc
   REAL(wp)    :: ep_cmoc    
   REAL(wp)    :: tvm_cmoc   
   REAL(wp)    :: vm_cmoc    
   REAL(wp)    :: kn_cmoc    
   !  Mortality Growth
   REAL(wp)    :: mpd_cmoc   
   REAL(wp)    :: mpd2_cmoc 
   !  Zooplankton
   REAL(wp)    :: rm_cmoc    
   REAL(wp)    :: kp_cmoc    
   REAL(wp)    :: ga_cmoc    
   REAL(wp)    :: mzn_cmoc   
   REAL(wp)    :: mzd_cmoc   
   REAL(wp)    :: mz2_cmoc   
   !   Detritus and remineralization
   REAL(wp)    :: ed_cmoc    
   REAL(wp)    :: reref_cmoc 
   !   Calcite parametrization
   REAL(wp)    :: rmcico_cmoc 
   REAL(wp)    :: trcico_cmoc
   REAL(wp)    :: aci_cmoc   
   REAL(wp)    :: dci_cmoc   
   !   Dinitrogen fixation
   REAL(wp)    :: phinf_cmoc 
   REAL(wp)    :: phi0_cmoc  
   REAL(wp)    :: anf_cmoc   
   REAL(wp)    :: pnf_cmoc   
   REAL(wp)    :: inf_cmoc   
   REAL(wp)    :: tnfMa_cmoc 
   REAL(wp)    :: tnfmi_cmoc 
   !   Redfield ratio and euphotic zone
   REAL(wp)    :: cnrr_cmoc  
   REAL(wp)    :: ncrr_cmoc  
   INTEGER     :: jk_eud_cmoc   ! <CMOC code OR 01/23/2016> scale of the euphotic zone
   INTEGER     :: nk_bal_cmoc   ! <CMOC code OR 01/23/2016> open ocean criterion
   ! <CMOC code OR 10/21/2015> CMOC block end


   !!*  diagnostic parameters 
   REAL(wp) ::  tpp                !: total primary production
   REAL(wp) ::  t_oce_co2_exp      !: total carbon export
   REAL(wp) ::  t_oce_co2_flx      !: Total ocean carbon flux
   REAL(wp) ::  t_oce_co2_flx_cum  !: Cumulative Total ocean carbon flux
   REAL(wp) ::  t_atm_co2_flx      !: global mean of atmospheric pco2

   !!* restoring
   LOGICAL  ::  ln_pisdmp          !: restoring or not of nutrients to a mean value
   INTEGER  ::  nn_pisdmp          !: frequency of relaxation or not of nutrients to a mean value


   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???



   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sms_cmoc.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_cmoc_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_cmoc_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      ! <CMOC code OR 11/13/2015> removing user-defined DNF diagnostics, revert to CMOC diagnostics
      INTEGER ::   ierr(5)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for cmoc & lobster
      ALLOCATE( etot(jpi,jpj,jpk),  xksi(jpi,jpj), STAT=ierr(1) )
      !
      !*  Biological fluxes for primary production
      ALLOCATE( xlimnfecmoc(jpi,jpj),            STAT=ierr(2) ) !  iron limitation mask
      ALLOCATE( xrcico  (jpi,jpj),               STAT=ierr(3) ) !  rain ratio 
         !
      !
      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)    , ak13  (jpi,jpj,jpk) ,       &
         &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,       &
         &      akw3(jpi,jpj,jpk)    , borat (jpi,jpj,jpk) ,       &
         &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,     STAT=ierr(4) )
         !
      !*  Remineralization
      ALLOCATE(redet(jpi,jpj,jpk), redettot(jpi,jpj),             &
        &     oomask(jpi,jpj), STAT=ierr(5) ) 
         !
      !
      sms_cmoc_alloc = MAXVAL( ierr )
      !
      IF( sms_cmoc_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_cmoc_alloc

#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO CMOC model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================   
END MODULE sms_cmoc    
