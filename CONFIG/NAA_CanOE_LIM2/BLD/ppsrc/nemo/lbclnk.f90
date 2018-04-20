MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! Ocean        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)     Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)     F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment  
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco) Add 'lbc_bdy_lnk' 
   !!                            and lbc_obc_lnk' routine to optimize  
   !!                            the BDY/OBC communications
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add a C1D case  
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lbc_lnk      : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_lnk_e    : generic interface for mpp_lnk_2d_e routine defined in lib_mpp
   !!   lbc_bdy_lnk  : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE lib_mpp          ! distributed memory computing library


   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE mpp_lnk_2d_9, mpp_lnk_2d_multiple
   END INTERFACE

   INTERFACE lbc_lnk
      MODULE PROCEDURE mpp_lnk_3d_gather, mpp_lnk_3d, mpp_lnk_2d
   END INTERFACE

   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE mpp_lnk_bdy_2d, mpp_lnk_bdy_3d
   END INTERFACE

   INTERFACE lbc_lnk_e
      MODULE PROCEDURE mpp_lnk_2d_e
   END INTERFACE

   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb
   END INTERFACE

   PUBLIC lbc_lnk       ! ocean lateral boundary conditions
   PUBLIC lbc_lnk_multi ! modified ocean lateral boundary conditions
   PUBLIC lbc_lnk_e
   PUBLIC lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC lbc_lnk_icb

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: lbclnk.F90 8113 2017-06-01 13:53:24Z lovato $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------


   !!======================================================================
END MODULE lbclnk

