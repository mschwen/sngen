***************************************************************************
* Filename: flavour_swap.f
***************************************************************************
*
*   Inputs:
*	o  neutrino type
*       o  neutrino energy
*   Outputs:
*	o  swapped neutrino type
*
***************************************************************************

      SUBROUTINE flavour_swap(type,E_nu,t,nu_type_swap)
   
      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'physics_constants.inc'

      INTEGER       type
      INTEGER       nu_type_swap
      REAL*8        E_nu
      REAL*8        t
      REAL*8        E_swap_nu
      REAL*8        E_swap_nu_bar      


      nu_type_swap = type

****************************************************************************
* This routine can either use time dependent values of the critical energies stored to file
* or the user can write the energies which will be used for all time below

      CALL E_swap(t,E_swap_nu,E_swap_nu_bar)

c       E_swap_nu = 0
c       E_swap_nu_bar = 0
****************************************************************************

c      WRITE(*,*) t, E_swap_nu, E_swap_nu_bar

      IF (type .EQ. id_nu_e) THEN
         IF (E_nu .GT. E_swap_nu) THEN
            nu_type_swap = id_nu_mu
         ENDIF
      ELSEIF (type .EQ. id_nu_mu) THEN
         IF (E_nu .GT. E_swap_nu) THEN
            nu_type_swap = id_nu_e
         ENDIF
      ELSEIF (type .EQ. id_nu_tau) THEN
         IF (E_nu .GT. E_swap_nu) THEN
            nu_type_swap = id_nu_e
         ENDIF
      ELSEIF (type .EQ. id_nu_ebar) THEN
         IF (E_nu .GT. E_swap_nu_bar) THEN
            nu_type_swap = id_nu_mubar
         ENDIF
      ELSEIF (type .EQ. id_nu_mubar) THEN
         IF (E_nu .GT. E_swap_nu_bar) THEN
            nu_type_swap = id_nu_ebar
         ENDIF
      ELSEIF (type .EQ. id_nu_taubar) THEN
         IF (E_nu .GT. E_swap_nu_bar) THEN
            nu_type_swap = id_nu_ebar
         ENDIF
      ENDIF

C      WRITE(*,*) nu_type_swap 

      RETURN
      END
