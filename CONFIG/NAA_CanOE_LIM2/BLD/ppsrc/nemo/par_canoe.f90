MODULE par_canoe
   !!======================================================================
   !!                        ***  par_canoe  ***
   !! TOP :   set the CanOE parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
    !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_canoe.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   !!---------------------------------------------------------------------
   !!   'key_canoe'   :                         standard CanOE bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_canoe     = .TRUE.  !: CanOE flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_kriest     = .FALSE. !: Kriest flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe     = 19      !: number of CanOE passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_2d  = 13      !: additional 2d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_3d  = 11      !: additional 3d output 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_trd =  1      !: number of sms trends for CanOE

   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  5    !: small particulate organic phosphate concentration
   
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  6    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpnn  =  7    !: Nanophytoplankton N concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpnfe =  8    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch =  9   !: Nano Chlorophyll Concentration
   
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 10    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdn =  11    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 12    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 13    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpzoo = 14    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 15    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 16    !: dissolved Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgoc = 17    !: big organic carbon  Concentration
   
   
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 18    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 19    !: Ammonium Concentration


   ! Starting/ending CanOE do-loop indices (N.B. no PISCES : jpl_can < jpf_can the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_can0     = 1                  !: First index of CanOE tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_can1     = jp_canoe          !: Last  index of CanOE tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_can0_2d  = 1               !: First index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_can1_2d  = jp_canoe_2d    !: Last  index of 2D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_can0_3d  = 1               !: First index of 3D diag
   INTEGER, PUBLIC, PARAMETER ::   jp_can1_3d  = jp_canoe_3d    !: Last  index of 3d diag
   INTEGER, PUBLIC, PARAMETER ::   jp_can0_trd = 1              !: First index of bio diag
   INTEGER, PUBLIC, PARAMETER ::   jp_can1_trd = jp_canoe_trd  !: Last  index of bio diag


   !!======================================================================
END MODULE par_canoe
