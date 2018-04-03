MODULE trcbbl
  !!======================================================================
   !!                       ***  MODULE  trcbbl  ***
   !! Ocean passive tracers physics :  advective and/or diffusive bottom boundary 
   !!                                  layer scheme
   !!======================================================================
   !!==============================================================================
   !! History :  OPA  !  1996-06  (L. Mortier)  Original code
   !!            8.0  !  1997-11  (G. Madec)    Optimization
   !!   NEMO     1.0  !  2002-08  (G. Madec)  free form + modules
   !!             -   !  2004-01  (A. de Miranda, G. Madec, J.M. Molines ) add advective bbl
   !!            3.3  !  2009-11  (G. Madec)  merge trabbl and trabbl_adv + style + optimization 
   !!             -   !  2010-04  (G. Madec)  Campin & Goosse advective bbl 
   !!             -   !  2010-06  (C. Ethe, G. Madec)  merge TRA-TRC
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_bbl( kt )              ! Empty routine
      WRITE(*,*) 'tra_bbl: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bbl

   !!======================================================================
END MODULE trcbbl
