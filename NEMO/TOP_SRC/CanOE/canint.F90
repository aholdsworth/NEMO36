MODULE canint
   !!======================================================================
   !!                         ***  MODULE canint  ***
   !! TOP :   CanOE interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                     2017-10  (A. Holdsworth) CanOE v1.1 NEMO3.6  
   !!----------------------------------------------------------------------
#if defined key_canoe
   !!----------------------------------------------------------------------
   !!   'key_canoe'                                       CanOE bio-model
   !!----------------------------------------------------------------------
   !!   can_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe      !  CanOE Source Minus Sink variables
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   can_int  
   PUBLIC   can_int_init  
   REAL(wp), PUBLIC ::  T0C        = 273.15_wp         !:
   REAL(wp), PUBLIC ::  Tref       = 25._wp            !:
   REAL(wp), PUBLIC ::  AEP        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ2       = -4500._wp         !:
   REAL(wp), PUBLIC ::  AER        = -4500._wp         !:


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canint.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE can_int( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE can_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  :: ji, jj                 ! dummy loop indices
      REAL(wp) :: zvar                   ! local variable
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('can_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfuncp(:,:,:) = EXP(AEP*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncz(:,:,:) = EXP(AEZ*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncz2(:,:,:) = EXP(AEZ2*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncr(:,:,:) = EXP(AER*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))

      !
      IF( nn_timing == 1 )  CALL timing_stop('can_int')
      !
   END SUBROUTINE can_int

   SUBROUTINE can_int_init
      INTEGER :: ios                 ! Local integer output status for namelist read
      !
      NAMELIST/nampistf/ AEP, AEZ, AEZ2, AER
      !!----------------------------------------------------------------------
      REWIND( numnatp_ref )              ! Namelist nampisext in reference namelist : Pisces atm. conditions
      READ  ( numnatp_ref, nampistf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampistf in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisext in configuration namelist : Pisces atm. conditions
      READ  ( numnatp_cfg, nampistf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampistf in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampistf )
      !


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for temperature, nampistf'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Activation energy for phytoplankton      AEP          =', AEP
         WRITE(numout,*) '    Activation energy for microzooplankton   AEZ          =', AEZ
         WRITE(numout,*) '    Activation energy for mesozooplankton    AEZ2         =', AEZ2
         WRITE(numout,*) '    Activation energy for remineralization   AER          =', AER
      ENDIF
      !
   END SUBROUTINE can_int_init

#else
   !!======================================================================
   !!  Dummy module :                                   No CanOE bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE can_int                   ! Empty routine
      WRITE(*,*) 'can_int: You should not have seen this print! error?'
   END SUBROUTINE can_int
#endif 

   !!======================================================================
END MODULE  canint
