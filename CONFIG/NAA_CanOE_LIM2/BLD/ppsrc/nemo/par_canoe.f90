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
   !!   Default                                   No CFC geochemical model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_canoe     = .FALSE.  !: CanOE flag 
   LOGICAL, PUBLIC, PARAMETER ::   lk_can        = .FALSE.  !: can flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe     =  0       !: No CFC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_2d  =  0       !: No CFC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_3d  =  0       !: No CFC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_canoe_trd =  0       !: number of sms trends for CanOE

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
