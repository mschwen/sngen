***************************************************************
* Filename: E_swap.f
***************************************************************
*
* Inputs: 
*    o time (ms)
*  
* Outputs:
*    o Critical energies for neutrino and antineutrino swapping
*
***************************************************************

      SUBROUTINE E_swap(t,E_swap_nu,E_swap_nu_bar)

      IMPLICIT NONE

      INCLUDE 'aspline_size.inc'
      INCLUDE 'energy_size.inc'
      INCLUDE 'time_size.inc'

      REAL*8   t ! Input
      REAL*8   last_time

      REAL*8   E_swap_nu, E_swap_nu_bar ! Output
      REAL*8   store(50000), store_bar(50000) ! Stores values

      COMMON           /sn_model/ 
     +                 sn_model_switch
      INTEGER          sn_model_switch

      INTEGER  E_data_file_size, index

      REAL*8 time(file_size_large), 
     +       E_c(file_size_large), E_c_bar(file_size_large)
      REAL*8 b_E_c(aspline_size), c_E_c(aspline_size),
     +       d_E_c(aspline_size)
      REAL*8 b_E_cbar(aspline_size), c_E_cbar(aspline_size),
     +       d_E_cbar(aspline_size)
      REAL*8 E_c_fit, E_c_bar_fit

      CHARACTER*40 E_data_file

      INTEGER  i, file_flag

      SAVE last_time

************************************************************************
*  Initializations
************************************************************************

      DATA file_flag              /0/
      DATA last_time              /0/

      DATA time                /file_size_large * 0./,
     +     E_c                 /file_size_large * 0./,
     +     E_c_bar             /file_size_large * 0./

      DATA b_E_c             /aspline_size * 0./,
     +     c_E_c             /aspline_size * 0./,
     +     d_E_c             /aspline_size * 0./,
     +     b_E_cbar          /aspline_size * 0./,
     +     c_E_cbar          /aspline_size * 0./,
     +     d_E_cbar          /aspline_size * 0./

****************************************************************
* Critical energy is dependent on what model is being used
****************************************************************

      index = NINT(t) + 1000 ! add 1000 to avoid negative indicies

* check to see if we need to preform calculation of can use stored values
* Makes sure that the time integration has gone through at least one cycle

c* Write critical energies to file
c      IF (t .GT. 480 .AND. file_flag .EQ. 0) THEN
c        OPEN (unit = 1, file = 'energy.out', status ='new')
c        DO i = -10, 50000
c        IF (store(i+1000) .NE. 0) THEN
c          WRITE (1,45) i,store(i+1000),store_bar(i+1000)
c45        FORMAT (I6.6,E14.8,F14.8)
c        ENDIF
c        END DO
c        file_flag = 1
c       ENDIF
c       CLOSE (unit = 1)

      IF (t .EQ. last_time) THEN

c      WRITE(*,*) 'Using stored critical energy values'

         E_swap_nu = store(index)
         E_swap_nu = store_bar(index)

      ELSE ! need to evaluate

      IF (sn_model_switch .EQ. 1) THEN  ! Burrows

          E_data_file = 'dat/E_swap/E_swap_Bur.dat'
          E_data_file_size = file_size_Bur 

      ELSEIF (sn_model_switch .EQ. 2) THEN ! Beacom & Vogel

          E_swap_nu = 8.1d0
          E_swap_nu_bar = 3.d0

          RETURN

      ELSEIF (sn_model_switch .EQ. 3) THEN ! Bruenn 15
 
          E_data_file = 'dat/E_swap/E_swap_Br15.dat'
          E_data_file_size = file_size_Br15

      ELSEIF (sn_model_switch .EQ. 4) THEN ! Bruenn 25

          E_data_file = 'dat/E_swap/E_swap_Br25.dat'
          E_data_file_size = file_size_Br25

      ELSEIF (sn_model_switch .EQ. 5) THEN ! Mezzacappa 13

          E_data_file = 'dat/E_swap/E_swap_Mez13.dat'
          E_data_file_size = file_size_Mez13

      ENDIF

      IF (sn_model_switch .NE. 2) THEN ! not Beacom and Vogel

      OPEN (unit = 50, file = E_data_file, status = 'old')

      DO 20 i = 1, E_data_file_size
         READ (50,*) time(i), E_c(i), E_c_bar(i)
         E_c(i) = dlog10(E_c(i))
         E_c_bar(i) = dlog10(E_c_bar(i))

20    CONTINUE

      CLOSE (unit = 50)


               CALL aspline(E_data_file_size,time,E_c,
     +                      b_E_c,c_E_c,d_E_c)

               CALL aspline(E_data_file_size,time,E_c_bar,
     +                      b_E_cbar,c_E_cbar,d_E_cbar)


            CALL aspline_eval(E_data_file_size,time,E_c,
     +                        b_E_c,c_E_c,d_E_c,t,E_c_fit)

            CALL aspline_eval(E_data_file_size,time,E_c_bar,
     +                        b_E_cbar,c_E_cbar,d_E_cbar,
     +                        t,E_c_bar_fit)


      E_swap_nu = 10.d0**(E_c_fit)
      E_swap_nu_bar = 10.d0**(E_c_bar_fit)

      store(index) = E_swap_nu
      store_bar(index) = E_swap_nu_bar


c      WRITE(*,*) t, store(index), store_bar(index)

      ENDIF ! calculation for models other than Beacom & Vogel

      ENDIF ! does calculation need to be preformed

* Write it to file

      IF (file_flag .EQ. 1 .AND. last_time .NE. t) THEN
      
       OPEN (unit = 1, file = 'energy.out', access = 'append',
     +       status = 'unknown')
       WRITE(1,45) t, store(index), store_bar(index)
       CLOSE (unit = 1)
45     FORMAT (F14.3,F14.3,F14.3)
  
      ENDIF

      last_time = t

      RETURN
      END

