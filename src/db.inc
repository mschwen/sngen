*FILE MEMBER=DB LIBRARY=SNOMAN LANGUAGE=INC DATE=16:Jul:1995
 
*     SNOMAN:  SDB System:  vector of key values for HEPDB interface
 
*     Contact:  R. Komar (UBC)
*		D. Cowen (UPenn)
 
*     Revision History:-
*     ================
 
*     2.09  R. Komar     First version
*     3_06_03  D. Cowen  Change kshl_pawc_size=2000000 -> 4000000
*              (Nick West's request).  13Mar98 DFC.
*     3_06_04  D. Cowen  Change kshl_pawc_size=4000000 -> 8000000
*              to deal with large ecal files.  20Apr98 DFC.
*     3_07     R. Komar  Change maxdate and maxtime from EST to GMT
 
 
*     Variables:-
*     =========
 
      integer kshl_pawc_size,nhbook
      parameter (kshl_pawc_size=8000000)
      integer nwp,ixpawc,ihdiv,ixhigz,ixku,fenc,lpawc,hcv
      common/pawc/nwp,ixpawc,ihdiv,ixhigz,ixku,fenc(5),
     *    lpawc,hcv(kshl_pawc_size-11)        !Zebra common store
      integer ihepdiv                   !HEPDB user division
      integer lhepdb, lsave, lprev, lnext  !structural links
      common/sdblnk/ihepdiv,lhepdb,lsave,lnext,lprev !link area for /PAWC/ store
      integer nnq
      parameter (nnq=kshl_pawc_size)              !dummy parameter size
      integer lqq(nnq),iqq(nnq)
      real*4 qq(nnq)
      equivalence (lqq(1),lpawc),(qq(1),iqq(1),lqq(9))
 
      integer nsyskeys
      parameter (nsyskeys=10)   !number of HEPDB system keys
      integer npairs
      parameter (npairs=1)	  !number of validity range pairs in keys
      integer ksdb_idx_first_usrkey  !index of first user key
      parameter (ksdb_idx_first_usrkey=nsyskeys+2*npairs+1)
      integer nusrkeys
      parameter(nusrkeys=6)      !number of user-defined keys
      integer usrkeys(nusrkeys)  !vector of user-defined keys
      character*8 chtag(nusrkeys)        !key labels
      character*10 chfor      !format of user keys
      integer*4 itag1           !first user key (data type)
      character*8 chtag2        !second user key (bank alias)
      integer itag3             !third user key (bank number)
      integer itag4             !fourth user key (task type)
      integer itag5             !fifth user key (format number)
      integer itag6             !sixth user key (source ID)
      integer ntotkeys
      parameter (ntotkeys=nsyskeys+2*npairs+nusrkeys)	!total number of keys
      integer keydbs(ntotkeys)     !keys vector
      integer imask(ntotkeys)      !mask for selecting keys for object search
      integer maxdate, maxtime     ! Maximum inst. of validity (don't touch)
      parameter (maxdate=20380517)
      parameter (maxtime=033319)
 
      character*4 chtop       !name of top directory
      character*80 chfile     !name of database file
      integer lch               !#characters
      integer lenocc            !CERN function to get #characters
 
*ENDFILE MEMBER=DB
