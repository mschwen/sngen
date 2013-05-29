*FILE MEMBER=GET_ESDIF_RAD LIBRARY=SNOMAN LANGUAGE=FORTRAN77 DATE=26:Feb:2003
      double precision function GET_ESDIF_RAD(NUTYPE, Enu, Te)
*
*     Monte Carlo generation - return differential cross section for a given
*	neutrino type, neutrino energy, and electron total energy
*       with radiative correction.
*
*     Contact:  Y. Takeuchi (Queen's)
*     Verified:
*     Refereed:
*
*     Parameters:-
*     ==========
*
*     NUTYPE         in    neutrino type
*     ENU            in	   neutrino energy (MeV)
*     Te             in	   electron kinetic energy (MeV)
*     GET_ESDIF_RAD  out   diff. cross section IN UNITS OF 10^-42 CM^2/MEV
*
*
*     Common Block Access:-
*     ===================
*
*     Specification:-
*     =============
*
*   o Return differential ES cross section for neutrino type NUTYPE and
*     energy ENU and electron energy En with radiative correction.
*
*     Assumptions:-
*     ===========
*
*   o See Bahcall et al, PRD 51, 6146-6158 (1995).
*   o This routine only support for nu_e and nu_mu at moment.
*   o This routine use pre-calculated table
*   o Original table: $QSNO_ROOT/qphysics/parameter/pnue_dif.dat
*
*     Program Notes:-
*     =============
*
*     Revision History:-
*     ================
*
*     4.02   Y.Takeuchi  initial version ported from get_esxt_e.for
*     4.02   Y.Takeuchi  fix a bug in "fpz" on 26-JUN-2003
 
      IMPLICIT NONE
 
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'
      INCLUDE 'id_nu.inc'
 
*     Argument Declarations:-
*     =====================
 
      INTEGER          NUTYPE
      DOUBLE PRECISION Enu, Te
 
*     Local Variable Declarations:-
*     ===========================
 
c*** pre-calculate table
      integer ndata_table
      parameter (ndata_table = 1440760) ! number of data in file (< ndata**2)
      integer ndata
      parameter (ndata = 1700)       ! number of max. array
      double precision Enu_step
      parameter (Enu_step = 0.1)   ! Enu & Ee step in table
 
      integer ireaddif
      data    ireaddif /0/
      double precision xdif_nue(ndata,ndata), xdif_num(ndata,ndata)
      SAVE ireaddif, xdif_nue, xdif_num
 
c*** for calculation
      double precision E, fm, fpz, fpm, x, z, el, IT
      double precision gl, gr, kappa, pnc, dl
c      double precision fsinthetaW2, alpha, dpi, Gf, factor
      double precision fsinthetaW2, dpi, Gf, factor
c**************************************************************************
c constants consistent with sngen parameters
c**************************************************************************
      parameter (fsinthetaW2 = sin2theta)
      parameter (dpi = pi)
c      parameter (alpha = 0.007297352568) ! fine structure; updated to PDG2004
      parameter (Gf = G_fermi)
      parameter (factor = hbar_c2 / 1.d-42) ! in units of 10^-42cm**2
c**************************************************************************
 
c*** others
      integer i, j, k
      double precision TMX, e1, e2, t1, t2, xdif_nue1, xdif_num1
      double precision e1_t1, e1_t2, e2_t1, e2_t2, r1, r2
 
*ENDHEADER
 
 
c*** read table if it has not done yet
      if (ireaddif.eq.0) then
 
c*** initialize
         do i = 1, ndata
            do j = 1, ndata
               xdif_nue(i,j) = 0.0d0
               xdif_num(i,j) = 0.0d0
            enddo
         enddo

         OPEN (unit = 10, file = 'dat/xsect/xs_nu_e_es_dif_rad.dat',
     +         status='old')

c*** read table (table is in 10^-42 cm^2/MeV)
         do i = 1, ndata_table
            READ (10,*) e1, e2, xdif_nue1, xdif_num1

c            e1        = RCONS(ldesrd + 4 * i - 3)
c            e2        = RCONS(ldesrd + 4 * i - 2)
c            xdif_nue1 = RCONS(ldesrd + 4 * i - 1)
c            xdif_num1 = RCONS(ldesrd + 4 * i - 0)
 
            j = int((e1+Enu_step/4.0) / Enu_step)
            k = int((e2+Enu_step/4.0) / Enu_step)
            if (j.lt.1.or.j.gt.ndata .or. k.lt.1.or.k.gt.ndata) then
               write(*,*) "get_esdif_rad: ERROR1 j,k,e1,e2 = ",
     $              j, k, e1, e2
               stop
            endif
            if (xdif_nue(j,k).ne.0.0 .or. xdif_nue(j,k).ne.0.0) then
               write(*,*) "get_esdif_rad: ERROR2 j,k,e1,e2 = ",
     $              j, k, e1, e2
               stop
            endif
            xdif_nue(j,k) = xdif_nue1
            xdif_num(j,k) = xdif_num1
c        write(IQPRNT,*) 'dif:',j,k,e1,e2,xdif_nue(j,k),xdif_num(j,k)
c        write(*,*)   'dif:',j,k,e1,e2,xdif_nue(j,k),xdif_num(j,k)
         enddo
c         write(IQPRNT,*) 'get_esdif_rad: read data=', ndata_table
c         write(*,*)   'get_esdif_rad: read data=', ndata_table
         ireaddif = 1
      endif
 
c*** MAX ALLOWED ELECTRON K.E.
      TMX = 2.0d0*ENU**2/(2.0d0*ENU+m_electron)
 
c*** check Te range
      if (Te .lt. 0.0d0 .or. Te .gt. TMX) then
         get_esdif_rad = 0.0d0
         return
      endif
 
c*** check Te == 0.0 (is this correct?)
      if (Te .eq. 0.0d0) then
         Te = 0.000001
      endif
 
c*** set j, k, etc.
      j = int(Enu / Enu_step) +1
      k = int(Te  / Enu_step) +1
 
c*** read from table
      if (j.gt.1 .and. j.lt.ndata .and. k.gt.1 .and. k.lt.ndata) then
         e1_t1 = 0.0d0
         e1_t2 = 0.0d0
         e2_t1 = 0.0d0
         e2_t2 = 0.0d0
         IF (NUTYPE.EQ.ID_NU_E) THEN
            e1_t1 = xdif_nue(j-1,k-1)
            e1_t2 = xdif_nue(j-1,k)
            e2_t1 = xdif_nue(j,k-1)
            e2_t2 = xdif_nue(j,k)
         ELSEIF ((NUTYPE.EQ.ID_NU_MU).OR.(NUTYPE.EQ.ID_NU_TAU)) THEN
            e1_t1 = xdif_num(j-1,k-1)
            e1_t2 = xdif_num(j-1,k)
            e2_t1 = xdif_num(j,k-1)
            e2_t2 = xdif_num(j,k)
         ELSE
c*** others are not supported yet
            write(*,*) "get_esdif_rad: not supported!  nutype=",
     $           nutype
            stop
         ENDIF
 
c*** check table contents
         if  (e1_t1.gt.0.0d0 .and. e1_t2.gt.0.0d0 .and.
     $        e2_t1.gt.0.0d0 .and. e2_t2.gt.0.0d0) then
 
c*** use table
            e1 = Enu_step * dfloat(j-1)
            e2 = Enu_step * dfloat(j)
            t1 = Enu_step * dfloat(k-1)
            t2 = Enu_step * dfloat(k)
 
            r1 = (e1_t2 - e1_t1) / (t2 - t1) * (Te - t1) + e1_t1
            r2 = (e2_t2 - e2_t1) / (t2 - t1) * (Te - t1) + e2_t1
            get_esdif_rad = ((r2 - r1) / (e2 - e1) * (Enu - e1) + r1)
 
c            print *, "==1==", enu, te, j, k
c            print *, "==2==", t1,Te,t2,e1_t1, r1, e1_t2
c            print *, "==3==", t1,Te,t2,e2_t1, r2, e2_t2
c            print *, "==4==", e1,Enu,e2,r1,get_esdif_rad, r2
 
            return
         endif
      endif
 
c*** don't use table
 
      E = Te + m_electron
      x = dsqrt(1.0d0 + 2.0d0 * m_electron / Te)
      z = Te / Enu
      el = sqrt(E * E - m_electron * m_electron)
      IT = 1.0d0 / 6.0d0 * (1.0d0 / 3.0d0 + (3.0d0 - x * x)
     $     * (0.5d0 * x * dlog((x + 1.0d0) / (x - 1.0d0)) - 1.0d0))
 
c*** To avoid log(-1.0e15) in fm,fpz,fpm when Te = Temax
      if (1.0-z-m_electron/(E+el) .le. 0.0) then
         write(*,*)
     $        "get_esdif_rad: warning: 1.0-z-me/(E+el) = %e\n",
     $        1.0-z-m_electron/(E+el)
         get_esdif_rad = 0.0d0
         return
      endif
 
c*** just use central values for pnc and kappa (correct ?)
      pnc = 1.0126
 
      IF (NUTYPE.EQ.ID_NU_E) THEN
         kappa = 0.9791 + 0.0097 * IT
         gl =  pnc * (0.5 - kappa * fsinthetaW2) - 1.0
         gr = -pnc * kappa * fsinthetaW2
      ELSEIF ((NUTYPE.EQ.ID_NU_MU).OR.(NUTYPE.EQ.ID_NU_TAU)) THEN
         kappa = 0.9970 - 0.00037 * IT
         gl =  pnc * (0.5 - kappa * fsinthetaW2)
         gr = -pnc * kappa * fsinthetaW2
      ELSE
c*** others are not supported yet
         write(*,*) "get_esdif_rad: not supported!  nutype=",
     $        nutype
         stop
      ENDIF
 
      fm = (E/el * dlog((E+el)/m_electron)-1.0d0) *
     $     (2.0d0*dlog(1.0d0-z-m_electron/(E+el))
     $     - dlog(1.0d0-z) - 0.5d0*dlog(z) - 5.0d0/12.0d0) + 0.5d0 *
     $     (dl(z) - dl(el/E)) - 0.5d0 * dlog(1.0d0-z)**2 -
     $     (11.0d0/12.0d0 + 0.5d0 * z) * dlog(1.0d0-z)
     $     + z * (dlog(z) + 0.5d0 * dlog(2.0d0*Enu/m_electron))
     $     - (31.0d0/18.0d0 + 1.0d0/12.0d0 * dlog(z))*el/E
     $     - 11.0d0/12.0d0 * z + z**2/24.0d0
 
      fpz = (E/el * dlog((E+el)/m_electron)-1.0d0) *
     $     ((1.0d0-z)**2 * (2.0d0*dlog(1.0d0-z-m_electron/(E+el))
     $     - dlog(1.0d0-z) - 0.5d0*dlog(z) - 2.0d0/3.0d0) - 0.5d0 *
     $     (z**2 * dlog(z) + 1.0d0 - z))
     $     - 0.5d0 * (1.0d0-z)**2 * (dlog(1.0d0-z)**2 + el/E *
     $     (dl(1.0d0-z) - dlog(z) * dlog(1.0d0-z)))
     $     + dlog(1.0d0-z) * (z**2/2.0d0*dlog(z) + (1.0d0-z)/3.0d0*
     $     (2.0d0*z-1.0d0/2.0d0))
     $     - 0.5d0 * z**2 * dl(1.0d0-z) - z*(1.0d0-2.0d0*z)/3.0d0 *
     $     dlog(z) - z*(1.0d0-z)/6.0d0 - el/E/12.0d0*(dlog(z) +
     $     (1.0d0-z)*(115.0d0-109.0d0*z)/6.0d0)
 
      fpm = (E/el * dlog((E+el)/m_electron)-1.0d0) * 2.0d0 *
     $     dlog(1.0d0-z-m_electron/(E+el))
 
      get_esdif_rad = 2.0d0 * Gf**2 * m_electron / dpi * (
     $     gl**2 * (1.0d0 + alpha / dpi * fm)
     $     + gr**2 * ((1.0d0 - z)**2 + alpha / dpi * fpz)
     $     - gr * gl * m_electron / Enu * z 
     $     * (1.0d0 + alpha / dpi * fpm))
 
cc 26-JUN-2003 bug fixed by y.takeuchi
cc   factor (1.0-z)**2 is already included in the fpz
cc     $     + gr**2 * (1.0d0 - z)**2 * (1.0d0 + alpha / dpi * fpz)
 
      get_esdif_rad = get_esdif_rad * factor
 
      RETURN
      END
 
cccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function dl(x)
cccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision x, t, dstep
      integer istep, i
c
      istep = 1000
 
      dstep = x / dfloat(istep)
      dl = 0.0d0
      do i = 1, istep
         t = dstep * dfloat(i)
         dl = dl + dlog(dabs(1.0d0 - t))/t * dstep
      enddo
      return
      end
 
*ENDFILE MEMBER=GET_ESDIF_RAD
