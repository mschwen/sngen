*****************************************************************************
* Filename: flux_intgral.f
*****************************************************************************
*
*   Inputs:
*       o interaction id (not array)
*
*   Outputs:
*       o data files with integrated supernova flux (including normal and
*          swapped)
*
*****************************************************************************

      SUBROUTINE flux_integral(interaction_id)

      IMPLICIT NONE

*****************************************************************************
     
      INTEGER interaction_id

      COMMON swap_flag
      INTEGER swap_flag

      INTEGER nu_type

* Energy loop constants
      
      REAL*8 E_initial, E_final
      REAL*8 E_step
      REAL*8 E
      REAL*8 sn_flux     ! function
    
      REAL*8 D1, D2, D3, D4
      REAL*8 E_initial2
      REAL*8 flux_t1, flux_t2

* Time loop constants
  
      REAL*8 t_initial, t_final
      REAL*8 t_step
      REAL*8 t

* Flux constants

      REAL*8 d_flux, Flux, Flux_dum

      REAL*8 d_flux1, d_flux2

* Averaging variables

      REAL*8 Flux_total, Flux_weight, E_average

      DOUBLEPRECISION Ran, ran_val

      REAL*8 P1, P2

* This program is simply used to create plots to verify that flavour swapping
* has been implemented correctly

      E_initial = 0
      E_final = 60
      E_step = 0.25
      t_initial = -9
      t_final = 10000
      t_step = 1
      Flux_dum = 0
      Flux = 0
      d_flux = 0

      Flux_total = 0
      Flux_weight = 0
      E_average = 0

* Flux swapping parameters
      d_flux1 = 0
      d_flux2 = 0
      D1 = 0
      D2 = 0
      D3 = 0
      D4 = 0

* Choose file name according to flux type and swap_flag

c      CALL get_nu_type(interaction_id,nu_type)

c      IF (swap_flag .EQ. 1) THEN

c         IF (nu_type .EQ. 30) THEN
c            OPEN (unit = 1, file = 'sn_flux/e_flux_Bu_swap.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 31) THEN
c            OPEN (unit = 1, file = 'sn_flux/e_bar_flux_Bu_swap.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 32 .OR. nu_type .EQ. 34) THEN
c            OPEN (unit = 1, file = 'sn_flux/x_flux_Bu_swap.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 33 .OR. nu_type .EQ. 35) THEN
c            OPEN (unit = 1, file = 'sn_flux/x_bar_flux_Bu_swap.out',
c     +         status = 'old')

c         ENDIF

c      ELSEIF (swap_flag .EQ. 0) THEN

c         IF (nu_type .EQ. 30) THEN
c            OPEN (unit = 1, file = 'sn_flux/e_flux_Bu.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 31) THEN
c            OPEN (unit = 1, file = 'sn_flux/e_bar_flux_Bu.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 32 .OR. nu_type .EQ. 34) THEN
c            OPEN (unit = 1, file = 'sn_flux/x_flux_Bu.out',
c     +         status = 'old')

c         ELSEIF (nu_type .EQ. 33 .OR. nu_type .EQ. 35) THEN
c            OPEN (unit = 1, file = 'sn_flux/x_bar_flux_Bu.out',
c     +         status = 'old')


c         ENDIF

c     ENDIF

* Loop through energy and time integrating flux over time

      
c      DO E = E_initial, E_final, E_step
         
c         Flux = 0
c         Flux_dum = 0
c         d_flux = 0

c         DO t = t_initial, t_final, t_step

c            d_flux =  sn_flux(interaction_id,E,t)

c            Flux_dum = Flux

c            Flux = Flux_dum + d_flux*t_step

c            d_flux = 0

c       END DO

c         WRITE (*,*) E, Flux

c          Flux_weight = Flux_weight + E*Flux
c          Flux_total = Flux_total + Flux

c          WRITE (1,45) Flux, E
c45        FORMAT (E14.8,F14.8)

c      END DO

c      CLOSE (unit = 1)      

c      E_average = Flux_weight/Flux_total

c      WRITE (*,*) E_average 
c      WRITE (*,*) 'Fluxes written to file'

*******************************
* Determine flavour swapping energy 
*******************************

      E_final = 170
      E_final = 100
      E_step = 0.05
      t = 1

      DO t = 500, 10000 , 100

50    Ran = ran_val(0.0d0, 40.0d0)
      E_initial = Ran

      D1 = 0
      D2 = 0
      D3 = 0
      D4 = 0

      flux_t1 = 0
      flux_t2 = 0

      Ran = ran_val(0.0d0, 7.d0)
      E_initial2 = Ran

      DO E = E_initial, E_final, E_step
        d_flux1 = sn_flux(57,E,t)*E_step 
        d_flux2 = sn_flux(66,E,t)*E_step
        D1 = D1 + d_flux1 - 2*d_flux2 
c        flux_t1 = flux_t1 + d_flux1
c        flux_t2 = flux_t2 + d_flux2    
c        WRITE(*,*) d_flux2
      END DO

c      WRITE (*,*) flux_t1, flux_t2

      DO E = E_initial2, E_final, E_step
        d_flux1 = sn_flux(61,E,t)*E_step
        d_flux2 = sn_flux(66,E,t)*E_step
        D2 = D2 + d_flux1 - 2*d_flux2
        flux_t1 = flux_t1 + d_flux1
        flux_t2 = flux_t2 + d_flux2
      END DO
        
c       WRITE(*,*) flux_t1, flux_t2

      DO E = E_initial, E_final, E_step
        d_flux1 = sn_flux(66,E,t)*E_step
        d_flux2 = sn_flux(57,E,t)*E_step
        D3 = D3 + d_flux1 - 0.5*d_flux2
      END DO

      DO E = E_initial2, E_final, E_step
        d_flux1 = sn_flux(66,E,t)*E_step
        d_flux2 = sn_flux(61,E,t)*E_step
        D4 = D4 + d_flux1 - 0.5*d_flux2
      END DO



      P1 = (D1-D2)/D1
      P1 = ABS(P1)
      P2 = (D3-D4)/D3
      P2 = ABS(P2)

c      WRITE(*,*) P1, P2
c      WRITE(*,*) E_initial, E_initial2

c       WRITE(*,*) D1, D2
c       WRITE(*,*) D3, D4

      IF (P1 .LT. 0.05 .AND. P2 .LT. 0.05) THEN
c       IF (P2 .LT. 0.05) THEN
c       WRITE(*,*) D1, D2
c       WRITE(*,*) D3, D4
       WRITE(*,*) E_initial, E_initial2, t
       ELSE
       GOTO 50
      ENDIF

      END DO

      RETURN
      END
