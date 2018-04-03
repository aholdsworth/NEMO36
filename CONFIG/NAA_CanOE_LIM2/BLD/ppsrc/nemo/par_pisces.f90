MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_pisces.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   !!---------------------------------------------------------------------
   !!   'key_pisces'   :                         standard PISCES bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_p4z        = .TRUE.  !: p4z flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .FALSE. !: Kriest flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 24      !: number of PISCES passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 11      !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_trd =  1      !: number of sms trends for PISCES

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppo4 =  5    !: phosphate concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  6    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  7    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  8    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  9    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc = 10    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 11    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 12    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 13    !: Diatoms Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 14    !: Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpbfe = 15    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgoc = 16    !: big particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsfe = 17    !: Small iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 18    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgsi = 19    !: (big) Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnfe = 20    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch = 21    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 22    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 23    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 24    !: Ammonium Concentration


   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_2d  = 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_2d  = jp_pisces_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_3d  = 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_3d  = jp_pisces_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0_trd = 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1_trd = jp_pisces_trd  !: Last  index of bio diag


   !!======================================================================
END MODULE par_pisces
