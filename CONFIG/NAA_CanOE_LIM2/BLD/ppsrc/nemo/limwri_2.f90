MODULE limwri_2
   !!======================================================================
   !!                     ***  MODULE  limwri_2  ***
   !!         Ice diagnostics :  write ice output files
   !!======================================================================
   !! history :  2.0  ! 2003-08  (C. Ethe)      original code
   !!            2.0  ! 2004-10  (C. Ethe )     1D configuration
   !!             -   ! 2009-06  (B. Lemaire )  iom_put + lim_wri_state_2
   !!-------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_lim2'                                    LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   lim_wri_2       : write of the diagnostics variables in ouput file 
   !!   lim_wri_init_2  : initialization and namelist read
   !!   lim_wri_state_2 : write for initial state or/and abandon:
   !!                     > output.init.nc (if ninist = 1 in namelist)
   !!                     > output.abort.nc
   !!----------------------------------------------------------------------
   USE phycst
   USE dom_oce
   USE sbc_oce
   USE sbc_ice
   USE dom_ice_2
   USE ice_2

   USE dianam           ! build name of file (routine)
   USE lbclnk
   USE in_out_manager
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE iom
   USE ioipsl
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE




   PUBLIC   lim_wri_state_2   ! called by dia_wri_state 
   PUBLIC   lim_wri_alloc_2   ! called by nemogcm.F90

   INTEGER, PARAMETER                       ::   jpnoumax = 40   ! maximum number of variable for ice output
   INTEGER                                  ::   noumef          ! number of fields
   REAL(wp)           , DIMENSION(jpnoumax) ::   cmulti ,     &  ! multiplicative constant
      &                                          cadd            ! additive constant
   CHARACTER(len = 35), DIMENSION(jpnoumax) ::   titn            ! title of the field
   CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::   nam             ! name of the field
   CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::   uni             ! unit of the field
   INTEGER            , DIMENSION(jpnoumax) ::   nc              ! switch for saving field ( = 1 ) or not ( = 0 )

   INTEGER ::   nice, nhorid, ndim, niter, ndepid       ! ????
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: ndex51   ! ????

   REAL(wp) ::   epsi16 = 1.e-16_wp   ! constant values
   REAL(wp) ::   zzero  = 0._wp       !     -      -
   REAL(wp) ::   zone   = 1._wp       !     -      -

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zcmo      ! Workspace array for netcdf writer. 


   !! * Substitutions
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
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limwri_2.F90 4696 2014-06-26 13:10:44Z clem $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lim_wri_alloc_2()
      !!-------------------------------------------------------------------
      !!                  ***   ROUTINE lim_wri_alloc_2  ***
      !!-------------------------------------------------------------------
      ALLOCATE( ndex51(jpij), zcmo(jpi,jpj,jpnoumax), STAT=lim_wri_alloc_2)
      !
      IF( lk_mpp               )   CALL mpp_sum ( lim_wri_alloc_2 )
      IF( lim_wri_alloc_2 /= 0 )   CALL ctl_warn('lim_wri_alloc_2: failed to allocate array ndex51')
      !
   END FUNCTION lim_wri_alloc_2



   SUBROUTINE lim_wri_state_2( kt, kid, kh_i )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lim_wri_state_2  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :
      !!   2.0  !  2009-06  (B. Lemaire)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt               ! ocean time-step index)
      INTEGER, INTENT( in ) ::   kid , kh_i       
      !!----------------------------------------------------------------------

      CALL histdef( kid, "isnowthi", "Snow thickness"          , "m"      , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicethic", "Ice thickness"           , "m"      , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iiceprod", "Ice produced"            , "m/kt"   , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "ileadfra", "Ice concentration"       , "-"      , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicetemp", "Ice temperature"         , "K"      , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicevelu", "i-Ice speed (I-point)"   , "m/s"    , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicevelv", "j-Ice speed (I-point)"   , "m/s"    , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "isstempe", "Sea surface temperature" , "C"      , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "isssalin", "Sea surface salinity"    , "PSU"    , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicestru", "i-Wind stress over ice (I-pt)", "Pa", jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicestrv", "j-Wind stress over ice (I-pt)", "Pa", jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicesflx", "Solar flux over ice"     , "w/m2"   , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicenflx", "Non-solar flux over ice" , "w/m2"   , jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "isnowpre", "Snow precipitation"      , "kg/m2/s", jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 

      CALL histend( kid, snc4set )   ! end of the file definition

      CALL histwrite( kid, "isnowthi", kt, hsnif          , jpi*jpj, (/1/) )   
      CALL histwrite( kid, "iicethic", kt, hicif          , jpi*jpj, (/1/) )    
      CALL histwrite( kid, "iiceprod", kt, hicifp         , jpi*jpj, (/1/) )   
      CALL histwrite( kid, "ileadfra", kt, 1. - frld(:,:) , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicetemp", kt, sist(:,:) - rt0, jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicevelu", kt, u_ice          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicevelv", kt, v_ice          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "isstempe", kt, sst_m          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "isssalin", kt, sss_m          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicestru", kt, utau_ice       , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicestrv", kt, vtau_ice       , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicesflx", kt, qsr_ice(:,:,1) , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicenflx", kt, qns_ice(:,:,1) , jpi*jpj, (/1/) )
      CALL histwrite( kid, "isnowpre", kt, sprecip        , jpi*jpj, (/1/) )

    END SUBROUTINE lim_wri_state_2


   !!======================================================================
END MODULE limwri_2
