MODULE par_cmoc
   !!======================================================================
   !!                        ***  par_cmoc  ***
   !! TOP :   set the CMOC parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!           CMOC1  !  2017-08  (A. Holdsworth)      implemented cmoc in nemo
   !3.6
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_cmoc.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE


#if defined key_cmoc
   !!---------------------------------------------------------------------
   !!   'key_cmoc'   :                         standard CMOC bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_cmoc     = .TRUE.  !: CMOC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc     = 8      !: number of CMOC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_3d  = 11      !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_trd =  1      !: number of sms trends for CMOC

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  4    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  5    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  6    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch =  7    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 =  8    !: Nitrates Concentration

#else
   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_cmoc     = .FALSE.  !: CMOC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_3d  =  0       !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_cmoc_trd =  0       !: number of sms trends for CMOC
#endif

   ! Starting/ending CMOC do-loop indices (N.B. no CMOC : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_cm0     = 1                  !: First index of CMOC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cm1     = jp_cmoc          !: Last  index of CMOC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_cm0_2d  = 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_cm1_2d  = jp_cmoc_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_cm0_3d  = 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_cm1_3d  = jp_cmoc_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_cm0_trd = 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_cm1_trd = jp_cmoc_trd  !: Last  index of bio diag


   !!======================================================================
END MODULE par_cmoc
