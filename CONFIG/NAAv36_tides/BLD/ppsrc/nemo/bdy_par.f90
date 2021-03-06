MODULE bdy_par
   !!======================================================================
   !!                      ***  MODULE bdy_par   ***
   !! Unstructured Open Boundary Cond. :   define related parameters
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D. Storkey and E. O'Dea) update for Shelf configurations
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC


   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy  = .TRUE.   !: Unstructured Ocean Boundary Condition flag



   INTEGER, PUBLIC, PARAMETER ::   jp_bdy  = 10       !: Maximum number of bdy sets
   INTEGER, PUBLIC, PARAMETER ::   jpbgrd  = 3	      !: Number of horizontal grid types used  (T, U, V)


   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdy_par.F90 4292 2013-11-20 16:28:04Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE bdy_par
