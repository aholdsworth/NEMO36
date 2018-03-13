MODULE dom_ice
   !!======================================================================
   !!                   ***  MODULE  dom_ice  ***
   !! LIM-3 Sea Ice :   Domain  variables
   !!======================================================================
   !! History :  3.0  ! 2003-08  (M. Vancoppenolle)  LIM-3 original code
   !!            3.5  ! 2011-02  (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC dom_ice_alloc   ! Routine called by nemogcm.F90

   LOGICAL, PUBLIC ::   l_jeq = .TRUE.       !: Equator inside the domain flag

   INTEGER, PUBLIC ::   njeq , njeqm1        !: j-index of the equator if it is inside the domain

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   fcor   !: coriolis coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   wght   !: weight of the 4 neighbours to compute averages

   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: dom_ice.F90 5123 2015-03-04 16:06:03Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dom_ice_alloc()
      !!-------------------------------------------------------------------
      !!            *** Routine dom_ice_alloc ***
      !!-------------------------------------------------------------------
      INTEGER :: dom_ice_alloc
      !!-------------------------------------------------------------------
      !
      ALLOCATE( fcor(jpi,jpj), wght(jpi,jpj,2,2), STAT = dom_ice_alloc )
      !
      IF( dom_ice_alloc /= 0 )   CALL ctl_warn( 'dom_ice_alloc: failed to allocate arrays.' )
      !
   END FUNCTION dom_ice_alloc

   !!======================================================================
END MODULE dom_ice
