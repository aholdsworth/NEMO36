MODULE p4zsbc
   !!======================================================================
   !!                         ***  MODULE p4sbc  ***
   !! TOP :   PISCES surface boundary conditions of external inputs of nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sbc_init   :  Initialization of p4z_sbc
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE iom             !  I/O manager
   USE fldread         !  time interpolation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sbc
   PUBLIC   p4z_sbc_init   

   !! * Shared module variables
   LOGICAL , PUBLIC  :: ln_dust     !: boolean for dust input from the atmosphere
   LOGICAL , PUBLIC  :: ln_solub    !: boolean for variable solubility of atmospheric iron
   LOGICAL , PUBLIC  :: ln_river    !: boolean for river input of nutrients
   LOGICAL , PUBLIC  :: ln_ndepo    !: boolean for atmospheric deposition of N
   LOGICAL , PUBLIC  :: ln_ironsed  !: boolean for Fe input from sediments
   LOGICAL , PUBLIC  :: ln_hydrofe  !: boolean for Fe input from hydrothermal vents
   LOGICAL , PUBLIC  :: ln_ironice  !: boolean for Fe input from sea ice
   REAL(wp), PUBLIC  :: sedfeinput  !: Coastal release of Iron
   REAL(wp), PUBLIC  :: dustsolub   !: Solubility of the dust
   REAL(wp), PUBLIC  :: mfrac       !: Mineral Content of the dust
   REAL(wp), PUBLIC  :: icefeinput  !: Iron concentration in sea ice
   REAL(wp), PUBLIC  :: wdust       !: Sinking speed of the dust 
   REAL(wp), PUBLIC  :: nitrfix     !: Nitrogen fixation rate   
   REAL(wp), PUBLIC  :: diazolight  !: Nitrogen fixation sensitivty to light 
   REAL(wp), PUBLIC  :: concfediaz  !: Fe half-saturation Cste for diazotrophs 
   REAL(wp)          :: hratio      !: Fe:3He ratio assumed for vent iron supply

   LOGICAL , PUBLIC  :: ll_sbc

   !! * Module variables
   LOGICAL  ::  ll_solub

   INTEGER , PARAMETER  :: jpriv  = 7   !: Maximum number of river input fields
   INTEGER , PARAMETER  :: jr_dic = 1   !: index of dissolved inorganic carbon
   INTEGER , PARAMETER  :: jr_doc = 2   !: index of dissolved organic carbon
   INTEGER , PARAMETER  :: jr_din = 3   !: index of dissolved inorganic nitrogen
   INTEGER , PARAMETER  :: jr_don = 4   !: index of dissolved organic nitrogen
   INTEGER , PARAMETER  :: jr_dip = 5   !: index of dissolved inorganic phosporus
   INTEGER , PARAMETER  :: jr_dop = 6   !: index of dissolved organic phosphorus
   INTEGER , PARAMETER  :: jr_dsi = 7   !: index of dissolved silicate


   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_dust      ! structure of input dust
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_solub      ! structure of input dust
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_river  ! structure of input riverdic
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_ndepo     ! structure of input nitrogen deposition
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_ironsed   ! structure of input iron from sediment
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_hydrofe   ! structure of input iron from hydrothermal vents

   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file
   INTEGER  :: ntimes_dust, ntimes_riv, ntimes_ndep       ! number of time steps in a file
   INTEGER  :: ntimes_solub, ntimes_hydro                 ! number of time steps in a file

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:) :: dust, solub       !: dust fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivdic, rivalk    !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivdin, rivdip    !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivdsi    !: river input fields
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:) :: nitdep    !: atmospheric N deposition 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ironsed   !: Coastal supply of iron
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: hydrofe   !: Hydrothermal vent supply of iron

   REAL(wp), PUBLIC :: sumdepsi, rivalkinput, rivdicinput, nitdepinput
   REAL(wp), PUBLIC :: rivdininput, rivdipinput, rivdsiinput


   !!* Substitution
   !!----------------------------------------------------------------------
   !!                    ***  top_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose : Statement function file: to be include in all passive tracer modules
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec) new architecture
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! s* or z*-coordinate (3D + time dependency) + use of additional now arrays (..._n)







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n





   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   Default option :                         eiv: dummy variables
   !!----------------------------------------------------------------------
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 6312 2016-02-15 11:43:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c3d' :                 aht: 3D coefficient
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
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: top_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsbc.F90 7522 2017-01-02 10:06:49Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER  :: ji,jj 
      REAL(wp) :: zcoef, zyyss
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc')

      !
      ! Compute dust at nit000 or only if there is more than 1 time record in dust file
      IF( ln_dust ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_dust > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_dust )
            IF( nn_ice_tr == -1 .AND. .NOT. ln_ironice ) THEN
               dust(:,:) = sf_dust(1)%fnow(:,:,1)
            ELSE
               dust(:,:) = sf_dust(1)%fnow(:,:,1) * ( 1.0 - fr_i(:,:) )
            ENDIF
         ENDIF
      ENDIF

      IF( ll_solub ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_solub > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_solub )
            solub(:,:) = sf_solub(1)%fnow(:,:,1)
         ENDIF
      ENDIF

      ! N/P and Si releases due to coastal rivers
      ! Compute river at nit000 or only if there is more than 1 time record in river file
      ! -----------------------------------------
      IF( ln_river ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_riv > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_river )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zcoef = ryyss * e1e2t(ji,jj) * h_rnf(ji,jj) 
                  rivalk(ji,jj) =   sf_river(jr_dic)%fnow(ji,jj,1)                                    &
                     &              * 1.E3        / ( 12. * zcoef + rtrn )
                  rivdic(ji,jj) = ( sf_river(jr_dic)%fnow(ji,jj,1) + sf_river(jr_doc)%fnow(ji,jj,1) ) &
                     &              * 1.E3         / ( 12. * zcoef + rtrn )
                  rivdin(ji,jj) = ( sf_river(jr_din)%fnow(ji,jj,1) + sf_river(jr_don)%fnow(ji,jj,1) ) &
                     &              * 1.E3 / rno3 / ( 14. * zcoef + rtrn )
                  rivdip(ji,jj) = ( sf_river(jr_dip)%fnow(ji,jj,1) + sf_river(jr_dop)%fnow(ji,jj,1) ) &
                     &              * 1.E3 / po4r / ( 31. * zcoef + rtrn )
                  rivdsi(ji,jj) =   sf_river(jr_dsi)%fnow(ji,jj,1)                                    &
                     &              * 1.E3        / ( 28.1 * zcoef + rtrn )
               END DO
            END DO
         ENDIF
      ENDIF

      ! Compute N deposition at nit000 or only if there is more than 1 time record in N deposition file
      IF( ln_ndepo ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_ndep > 1 ) ) THEN
             zcoef = rno3 * 14E6 * ryyss
             CALL fld_read( kt, 1, sf_ndepo )
             nitdep(:,:) = sf_ndepo(1)%fnow(:,:,1) / zcoef / e3t_n(:,:,1) 
         ENDIF
         IF( lk_vvl ) THEN
           zcoef = rno3 * 14E6 * ryyss
           nitdep(:,:) = sf_ndepo(1)%fnow(:,:,1) / zcoef / e3t_n(:,:,1) 
         ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc')
      !
   END SUBROUTINE p4z_sbc

   SUBROUTINE p4z_sbc_init

      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ji, jj, jk, jm, ifpr
      INTEGER  :: ii0, ii1, ij0, ij1
      INTEGER  :: numdust, numsolub, numriv, numiron, numdepo, numhydro
      INTEGER  :: ierr, ierr1, ierr2, ierr3
      INTEGER  :: ios                 ! Local integer output status for namelist read
      INTEGER  :: ik50                !  last level where depth less than 50 m
      INTEGER  :: isrow             ! index for ORCA1 starting row
      REAL(wp) :: zexpide, zdenitide, zmaskt
      REAL(wp) :: ztimes_dust, ztimes_riv, ztimes_ndep 
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      REAL(wp), DIMENSION(:), ALLOCATABLE :: rivinput
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zdust, zndepo, zriver, zcmask
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(jpriv) ::  slf_river    ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_dust, sn_solub, sn_ndepo, sn_ironsed, sn_hydrofe   ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_riverdoc, sn_riverdic, sn_riverdsi   ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_riverdin, sn_riverdon, sn_riverdip, sn_riverdop
      !
      NAMELIST/nampissbc/cn_dir, sn_dust, sn_solub, sn_riverdic, sn_riverdoc, sn_riverdin, sn_riverdon,     &
        &                sn_riverdip, sn_riverdop, sn_riverdsi, sn_ndepo, sn_ironsed, sn_hydrofe, &
        &                ln_dust, ln_solub, ln_river, ln_ndepo, ln_ironsed, ln_ironice, ln_hydrofe,    &
        &                sedfeinput, dustsolub, icefeinput, wdust, mfrac, nitrfix, diazolight, concfediaz, hratio
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc_init')
      !
      !                            !* set file information
      REWIND( numnatp_ref )              ! Namelist nampissbc in reference namelist : Pisces external sources of nutrients
      READ  ( numnatp_ref, nampissbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampissbc in configuration namelist : Pisces external sources of nutrients
      READ  ( numnatp_cfg, nampissbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissbc )

      IF ( ( nn_ice_tr >= 0 ) .AND. ln_ironice ) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' ln_ironice incompatible with nn_ice_tr = ', nn_ice_tr
            WRITE(numout,*) ' Specify your sea ice iron concentration in nampisice instead '
            WRITE(numout,*) ' ln_ironice is forced to .FALSE. '
            ln_ironice = .FALSE.
         ENDIF
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampissbc '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    dust input from the atmosphere           ln_dust     = ', ln_dust
         WRITE(numout,*) '    Variable solubility of iron input        ln_solub    = ', ln_solub
         WRITE(numout,*) '    river input of nutrients                 ln_river    = ', ln_river
         WRITE(numout,*) '    atmospheric deposition of n              ln_ndepo    = ', ln_ndepo
         WRITE(numout,*) '    Fe input from sediments                  ln_ironsed  = ', ln_ironsed
         WRITE(numout,*) '    Fe input from seaice                     ln_ironice  = ', ln_ironice
         WRITE(numout,*) '    fe input from hydrothermal vents         ln_hydrofe  = ', ln_hydrofe
         WRITE(numout,*) '    coastal release of iron                  sedfeinput  = ', sedfeinput
         WRITE(numout,*) '    solubility of the dust                   dustsolub   = ', dustsolub
         WRITE(numout,*) '    Mineral Fe content of the dust           mfrac       = ', mfrac
         WRITE(numout,*) '    Iron concentration in sea ice            icefeinput  = ', icefeinput
         WRITE(numout,*) '    sinking speed of the dust                wdust       = ', wdust
         WRITE(numout,*) '    nitrogen fixation rate                   nitrfix     = ', nitrfix
         WRITE(numout,*) '    nitrogen fixation sensitivty to light    diazolight  = ', diazolight
         WRITE(numout,*) '    fe half-saturation cste for diazotrophs  concfediaz  = ', concfediaz
         WRITE(numout,*) '    Fe to 3He ratio assumed for vent iron supply hratio  = ', hratio
      END IF

      IF( ln_dust .OR. ln_river .OR. ln_ndepo ) THEN  ;  ll_sbc = .TRUE.
      ELSE                                            ;  ll_sbc = .FALSE.
      ENDIF

      IF( ln_dust .AND. ln_solub ) THEN               ;  ll_solub = .TRUE.
      ELSE                                            ;  ll_solub = .FALSE.
      ENDIF


      ! dust input from the atmosphere
      ! ------------------------------
      IF( ln_dust ) THEN 
         !
         IF(lwp) WRITE(numout,*) '    initialize dust input from atmosphere '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
         !
         ALLOCATE( dust(jpi,jpj) )    ! allocation
         !
         ALLOCATE( sf_dust(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_dust structure' )
         !
         CALL fld_fill( sf_dust, (/ sn_dust /), cn_dir, 'p4z_sed_init', 'Atmospheric dust deposition', 'nampissed' )
                                   ALLOCATE( sf_dust(1)%fnow(jpi,jpj,1)   )
         IF( sn_dust%ln_tint )     ALLOCATE( sf_dust(1)%fdta(jpi,jpj,1,2) )
         !
         IF( Agrif_Root() ) THEN   !  Only on the master grid
            ! Get total input dust ; need to compute total atmospheric supply of Si in a year
            CALL iom_open (  TRIM( sn_dust%clname ) , numdust )
            CALL iom_gettime( numdust, zsteps, kntime=ntimes_dust)  ! get number of record in file
            ALLOCATE( zdust(jpi,jpj,ntimes_dust) )
            DO jm = 1, ntimes_dust
               CALL iom_get( numdust, jpdom_data, TRIM( sn_dust%clvar ), zdust(:,:,jm), jm )
            END DO
            CALL iom_close( numdust )
            ztimes_dust = 1._wp / FLOAT( ntimes_dust ) 
            sumdepsi = 0.e0
            DO jm = 1, ntimes_dust
               sumdepsi = sumdepsi + glob_sum( zdust(:,:,jm) * e1e2t(:,:) * tmask(:,:,1) * ztimes_dust )
            ENDDO
            sumdepsi = sumdepsi / ( nyear_len(1) * rday ) * 12. * 8.8 * 0.075 * mfrac / 28.1 
            DEALLOCATE( zdust)
         ENDIF
      ELSE
         sumdepsi  = 0._wp
      END IF

      ! Solubility of dust deposition of iron
      ! Only if ln_dust and ln_solubility set to true (ll_solub = .true.)
      ! -----------------------------------------------------------------
      IF( ll_solub ) THEN
         !
         IF(lwp) WRITE(numout,*) '    initialize variable solubility of Fe '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
         !
         ALLOCATE( solub(jpi,jpj) )    ! allocation
         !
         ALLOCATE( sf_solub(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_solub structure' )
         !
         CALL fld_fill( sf_solub, (/ sn_solub /), cn_dir, 'p4z_sed_init', 'Solubility of atm. iron ', 'nampissed' )
                                   ALLOCATE( sf_solub(1)%fnow(jpi,jpj,1)   )
         IF( sn_solub%ln_tint )    ALLOCATE( sf_solub(1)%fdta(jpi,jpj,1,2) )
         ! get number of record in file
         CALL iom_open (  TRIM( sn_solub%clname ) , numsolub )
         CALL iom_gettime( numsolub, zsteps, kntime=ntimes_solub)  ! get number of record in file
         CALL iom_close( numsolub )
      ENDIF

      ! nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         !
         slf_river(jr_dic) = sn_riverdic  ;  slf_river(jr_doc) = sn_riverdoc  ;  slf_river(jr_din) = sn_riverdin 
         slf_river(jr_don) = sn_riverdon  ;  slf_river(jr_dip) = sn_riverdip  ;  slf_river(jr_dop) = sn_riverdop
         slf_river(jr_dsi) = sn_riverdsi  
         !
         ALLOCATE( rivdic(jpi,jpj), rivalk(jpi,jpj), rivdin(jpi,jpj), rivdip(jpi,jpj), rivdsi(jpi,jpj) ) 
         !
         ALLOCATE( sf_river(jpriv), rivinput(jpriv), STAT=ierr1 )           !* allocate and fill sf_river (forcing structure) with sn_river_
         rivinput(:) = 0.0

         IF( ierr1 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_irver structure' )
         !
         CALL fld_fill( sf_river, slf_river, cn_dir, 'p4z_sed_init', 'Input from river ', 'nampissed' )
         DO ifpr = 1, jpriv
                                          ALLOCATE( sf_river(ifpr)%fnow(jpi,jpj,1  ) )
            IF( slf_river(ifpr)%ln_tint ) ALLOCATE( sf_river(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         IF( Agrif_Root() ) THEN   !  Only on the master grid
            ! Get total input rivers ; need to compute total river supply in a year
            DO ifpr = 1, jpriv
               CALL iom_open ( TRIM( slf_river(ifpr)%clname ), numriv )
               CALL iom_gettime( numriv, zsteps, kntime=ntimes_riv)
               ALLOCATE( zriver(jpi,jpj,ntimes_riv) )
               DO jm = 1, ntimes_riv
                  CALL iom_get( numriv, jpdom_data, TRIM( slf_river(ifpr)%clvar ), zriver(:,:,jm), jm )
               END DO
               CALL iom_close( numriv )
               ztimes_riv = 1._wp / FLOAT(ntimes_riv) 
               DO jm = 1, ntimes_riv
                  rivinput(ifpr) = rivinput(ifpr) + glob_sum( zriver(:,:,jm) * tmask(:,:,1) * ztimes_riv ) 
               END DO
               DEALLOCATE( zriver)
            END DO
            ! N/P and Si releases due to coastal rivers
            ! -----------------------------------------
            rivdicinput = (rivinput(jr_dic) + rivinput(jr_doc) ) * 1E3 / 12._wp
            rivdininput = (rivinput(jr_din) + rivinput(jr_don) ) * 1E3 / rno3 / 14._wp
            rivdipinput = (rivinput(jr_dip) + rivinput(jr_dop) ) * 1E3 / po4r / 31._wp
            rivdsiinput = rivinput(jr_dsi) * 1E3 / 28.1_wp
            rivalkinput = rivinput(jr_dic) * 1E3 / 12._wp
            !
         ENDIF
      ELSE
         rivdicinput = 0._wp
         rivdininput = 0._wp
         rivdipinput = 0._wp
         rivdsiinput = 0._wp
         rivalkinput = 0._wp
      END IF 
      ! nutrient input from dust
      ! ------------------------
      IF( ln_ndepo ) THEN
         !
         IF(lwp) WRITE(numout,*) '    initialize the nutrient input by dust from ndeposition.orca.nc'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ALLOCATE( nitdep(jpi,jpj) )    ! allocation
         !
         ALLOCATE( sf_ndepo(1), STAT=ierr3 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr3 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_ndepo structure' )
         !
         CALL fld_fill( sf_ndepo, (/ sn_ndepo /), cn_dir, 'p4z_sed_init', 'Nutrient atmospheric depositon ', 'nampissed' )
                                   ALLOCATE( sf_ndepo(1)%fnow(jpi,jpj,1)   )
         IF( sn_ndepo%ln_tint )    ALLOCATE( sf_ndepo(1)%fdta(jpi,jpj,1,2) )
         !
         IF( Agrif_Root() ) THEN   !  Only on the master grid
            ! Get total input dust ; need to compute total atmospheric supply of N in a year
            CALL iom_open ( TRIM( sn_ndepo%clname ), numdepo )
            CALL iom_gettime( numdepo, zsteps, kntime=ntimes_ndep)
            ALLOCATE( zndepo(jpi,jpj,ntimes_ndep) )
            DO jm = 1, ntimes_ndep
               CALL iom_get( numdepo, jpdom_data, TRIM( sn_ndepo%clvar ), zndepo(:,:,jm), jm )
            END DO
            CALL iom_close( numdepo )
            ztimes_ndep = 1._wp / FLOAT( ntimes_ndep ) 
            nitdepinput = 0._wp
            DO jm = 1, ntimes_ndep
              nitdepinput = nitdepinput + glob_sum( zndepo(:,:,jm) * e1e2t(:,:) * tmask(:,:,1) * ztimes_ndep )
            ENDDO
            nitdepinput = nitdepinput / rno3 / 14E6 
            DEALLOCATE( zndepo)
         ENDIF
      ELSE
         nitdepinput = 0._wp
      ENDIF

      ! coastal and island masks
      ! ------------------------
      IF( ln_ironsed ) THEN     
         !
         IF(lwp) WRITE(numout,*) '    computation of an island mask to enhance coastal supply of iron'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ALLOCATE( ironsed(jpi,jpj,jpk) )    ! allocation
         !
         CALL iom_open ( TRIM( sn_ironsed%clname ), numiron )
         ALLOCATE( zcmask(jpi,jpj,jpk) )
         CALL iom_get  ( numiron, jpdom_data, TRIM( sn_ironsed%clvar ), zcmask(:,:,:), 1 )
         CALL iom_close( numiron )
         !
         ik50 = 5        !  last level where depth less than 50 m
         DO jk = jpkm1, 1, -1
            IF( gdept_1d(jk) > 50. )  ik50 = jk - 1
         END DO
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,*) ' Level corresponding to 50m depth ',  ik50,' ', gdept_1d(ik50+1)
         IF (lwp) WRITE(numout,*)
         DO jk = 1, ik50
            DO jj = 2, jpjm1
               DO ji = 1, jpi
                  IF( tmask(ji,jj,jk) /= 0. ) THEN
                     zmaskt = tmask(ji+1,jj,jk) * tmask(ji-1,jj,jk) * tmask(ji,jj+1,jk)    &
                        &                       * tmask(ji,jj-1,jk) * tmask(ji,jj,jk+1)
                     IF( zmaskt == 0. )   zcmask(ji,jj,jk ) = MAX( 0.1, zcmask(ji,jj,jk) ) 
                  END IF
               END DO
            END DO
         END DO
         !
         CALL lbc_lnk( zcmask , 'T', 1. )      ! lateral boundary conditions on cmask   (sign unchanged)
         !
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zexpide   = MIN( 8.,( gdept_0(ji,jj,jk) / 500. )**(-1.5) )
                  zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
                  zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
               END DO
            END DO
         END DO
         ! Coastal supply of iron
         ! -------------------------
         ironsed(:,:,jpk) = 0._wp
         DO jk = 1, jpkm1
            ironsed(:,:,jk) = sedfeinput * zcmask(:,:,jk) / ( e3t_0(:,:,jk) * rday )
         END DO
         DEALLOCATE( zcmask)
      ENDIF
      !
      ! Iron from Hydrothermal vents
      ! ------------------------
      IF( ln_hydrofe ) THEN
         !
         IF(lwp) WRITE(numout,*) '    Input of iron from hydrothermal vents '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ALLOCATE( hydrofe(jpi,jpj,jpk) )    ! allocation
         !
         CALL iom_open ( TRIM( sn_hydrofe%clname ), numhydro )
         CALL iom_get  ( numhydro, jpdom_data, TRIM( sn_hydrofe%clvar ), hydrofe(:,:,:), 1 )
         CALL iom_close( numhydro )
         !
         DO jk = 1, jpk
            hydrofe(:,:,jk) = ( hydrofe(:,:,jk) * hratio ) / ( e1e2t(:,:) * e3t_0(:,:,jk) * ryyss + rtrn ) / 1000._wp
         ENDDO
         !
      ENDIF
      ! 
      IF( ll_sbc ) CALL p4z_sbc( nit000 ) 
      !
      IF(lwp) THEN 
         WRITE(numout,*)
         WRITE(numout,*) '    Total input of elements from river supply'
         WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    N Supply   : ', rivdininput*rno3*1E3/1E12*14.,' TgN/yr'
         WRITE(numout,*) '    Si Supply  : ', rivdsiinput*1E3/1E12*28.1,' TgSi/yr'
         WRITE(numout,*) '    P Supply   : ', rivdipinput*1E3*po4r/1E12*31.,' TgP/yr'
         WRITE(numout,*) '    Alk Supply : ', rivalkinput*1E3/1E12,' Teq/yr'
         WRITE(numout,*) '    DIC Supply : ', rivdicinput*1E3*12./1E12,'TgC/yr'
         WRITE(numout,*) 
         WRITE(numout,*) '    Total input of elements from atmospheric supply'
         WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    N Supply   : ', nitdepinput*rno3*1E3/1E12*14.,' TgN/yr'
         WRITE(numout,*) 
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc_init')
      !
   END SUBROUTINE p4z_sbc_init


   !!======================================================================
END MODULE p4zsbc
