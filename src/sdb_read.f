*FILE MEMBER=SDB_READ LIBRARY=SNOMAN LANGUAGE=FORTRAN77 DATE=11:Jan:2000
      subroutine sdb_read(ipath,lbank,idiv,iretc)
 
*     SNOMAN: SDB System: Fetches a data bank from a database file.
 
*     Contact:  R. Komar (UBC)
 
*     Parameters:-
*     ==========
 
*     ipath          in      vector of key values for selecting the object
*     lbank          out     address of data bank retrieved from database
*     idiv           in      division number where bank is placed
*     iretc          out     error return code
 
*     To retrieve a bank based on its serial number, set the date in
*     the `ipath' vector to the negative value of the serial number.
*     All other keys (besides the bank name) are then ignored.
 
*     In addition to the imported arguments, the following options can be
*     specified externally by calling SDB_UTIL_PUT_OPTION(index,option):
*	index = ksdb_idx_which_bank
*	option = ksdb_opt_valid	(default)
*		 ksdb_opt_next
*		 ksdb_opt_prev
*	index = ksdb_idx_validity_inst
*	option = ksdb_opt_date (default)
*		 ksdb_opt_run
*	index = ksdb_idx_bank_type
*	option = ksdb_opt_data (default)
*		 ksdb_opt_format
*	index = ksdb_idx_verification
*	option = ksdb_opt_verify (default)
*		 ksdb_opt_noverify
*
*     Note that the default value for each option is used if a non-default
*     value isn't explicitely set.  Once an option is set, it remains that
*     way until SDB_UTIL_PUT_OPTION is called to change it again (ie. it
*     does not revert to the default value after a call to SDB_READ).
*
*     (see SDB_OPTION.INC for more details)
 
*     Common Block Access:-
*     ===================
*     (ignoring internal and environmental access)
 
*     /PAWC/,/Z/,/MT/	 out       retrieved bank is placed in one of these
*				   stores
 
 
*     Specification:-
*     =============
 
* o   Read database object and store in the right data structure.
 
*     Program Notes:-
*     =============
 
*     Link to keys bank is returned through the link area /sdblnk/
*
*     The `next' or `previous' bank modes should NOT be used when running
*     in supernova monitoring mode.
 
*     Revision History:-
*     ================
 
*     2_09  Rob Komar      First version.
*     2.09  G. Jonkmans    Change reference to the first usrkeys
*                          as ipath(1-ksdb_ipath_offset)
*           G. Jonkmans    Encode 4hKYCD, throught a call uctoh(),
*                          since VMS won't support bank name in 4H.
*           N. West        Bug fix TZFIND calls.
*                          Protect lbank (call array to force pass by reference)
*	    Rob Komar	   Drop titles and keys bank support. Drop insert
*			   and replace modes.  Add titles format support.
*     3_00  G. Jonkmans    Add effective validity range and fix 'previous'
*                          and 'next' options. Values are return in
*                          the ipath array.
*           N. West        Use revised error codes.  Don't report BANK_NOT_FOUND
*                          Declare local IPATH.
*	    R. Komar	   Add checksum verification of data bank.
*			   Remove calculation of offsets between stores
*			   (unnecessary after titles bank support was dropped).
*	    R. Komar	   Changed error return code from `iretc' to `irc' in
*			   the last calls to cdfree to decouple them from
*			   error codes returned when producing the output bank.
*	    R. Komar	   Verify checksum in this routine as SDB_UTIL_CHECKSUM
*			   was modified to return the checksum without
*			   verifying it.
*	    R. Komar	   General clean up of lsave,lprev,lnext usage,
*			   replaced calls to cdfree with calls to mzdrop,
*			   added translation of validity time from HHMMSSCC
*			   to HHMMSS in this routine, added support for
*			   retrieving a bank based on serial number.
*    3_02   R. Komar       Include `sdb_contents_type.inc' (which is used
*			   for handling the returned bank according to the
*			   type of its contents (defined in key 10 of the
*			   HEPDB system keys)).  Also, update the bank's
*			   database header with the effective validity
*			   range and database insertion time if required.
*           R. Komar       Modify character option to MZCOPY to support
*                          data structures.
*                          Fix bank name and number if wrong (temporary kludge)
*	    R. Komar	   Y2K fix for updating insertion date in bank's
*			   database header.  Add support for supernova
*			   monitoring switch.
*           N. West        Update formats to allow 10 digit bank numbers.
      implicit none
 
      include 'memory_size.inc'
      include 'zunit.inc'
      include 'mt.inc'
      include 'sdb.inc'
      include 'sdb_options.inc'
      include 'sdb_ipath.inc'
      include 'id_errors.inc'
      include 'sdb_contents_type.inc'
      include 'titles_dbhdr.inc'
 
*     Argument Declarations:-
*     =====================
 
      integer ipath(num_ipath_args)
      integer lbank(1),idiv
      integer iretc
 
 
*     Local Variable Declarations:-
*     ===========================
 
      integer idate,itime
      integer isel, isel_tmp
      integer isum(2)
      character*8 chopt
      integer lchopt
      character*4 chbank
      character*7 chalias
      character*100 chpath
      character*8 chtype
      integer lchtype
      integer lenocc,lch
      integer startval, endval
      integer start_gap_val, end_gap_val
      integer koffset
 
      integer usrkeys(nusrkeys)  !vector of user-defined keys
      integer iuser
 
c
c  Initialize the error return code.
c
      iretc=0
c
c  Initialize temporary links used internally.
c
      lhepdb=0
      lsave=0
      lnext=0
      lprev=0
	lbank(1)=0
c
c Initialize effective validity range counters
 
      startval = 0
      endval = 0
      start_gap_val = 0              ! To get the previous bank
      call sdb_util_pkts(maxdate,maxtime,end_gap_val)
c                                    ! To get the next bank
c
c  If no target division is given, set the division number to that
c  of the HEPDB user division.
c
      if(idiv .eq. 0) then
         idiv=ihepdiv
      endif
c
c  Encode the instant of validity.
c  (for now assume there is only one instant of validity, and that
c  it is either the packed date/time or the run-number).  If the
c  date is a negative number, then it is interpreted as the serial
c  number of the bank to be retrieved, and no instant of validity
c  is necessary.
c
      if(iopt(ksdb_idx_validity_inst) .eq. ksdb_opt_run) then
         isel=ipath(ksdb_idx_ipath_runnum)
      else
	 iopt(ksdb_idx_validity_inst)=ksdb_opt_date
         idate=ipath(ksdb_idx_ipath_date)
         itime=ipath(ksdb_idx_ipath_time)/100	! HHMMSSCC -> HHMMSS
	 if (idate .ge. 0) then
            call sdb_util_pkts(idate,itime,isel)
	 else
	    isel = idate			! select based on serial number
	    iopt(ksdb_idx_which_bank) = ksdb_opt_valid
	 endif
      endif
c
c  Set the user key values for the search.
c
      call ucopy(ipath(1-ksdb_ipath_offset),usrkeys,nusrkeys)
c
c  Set the options for the call to fetch either just the keys bank
c  or both the keys and data banks for each object.
c
      if(iopt(ksdb_idx_which_bank) .eq. ksdb_opt_next  .or.
     +   iopt(ksdb_idx_which_bank) .eq. ksdb_opt_prev) then
         chopt='FK'     !fetch only the keys from the database file
      else
         iopt(ksdb_idx_which_bank)=ksdb_opt_valid
         chopt='F'      !fetch keys and data from the database file
      endif
c
c  If the supernova monitoring switch is on, add the `N' option.
c
      if(iopt(ksdb_idx_SN_mon) .eq. ksdb_opt_SN_mon) then
	 lchopt=lenocc(chopt)
	 chopt=chopt(1:lchopt)//'N'
      endif
c
c  Get the pathname for the directory holding the bank
c
      call uhtoc(ipath(ksdb_idx_ipath_bankname),4,chbank,4)
      if(iopt(ksdb_idx_bank_type) .eq. ksdb_opt_format) then
	 lch=lenocc(chbank)
	 chalias=chbank(1:lch)//'fmt'
	 chtype=' format '
	 lchtype=8
      else
	 chalias=chbank
	 chtype=' '
	 lchtype=1
      endif
      call sdb_util_pathname(chalias,chpath,iretc)
      if(iretc.ne.0) then
	 write(iqlog,9002) chtype(1:lchtype),chbank
	 write(iqprnt,9002) chtype(1:lchtype),chbank
 9002    format('SDB_READ: Cannot fetch',a,'bank "',a4,'" due to ',
     +      'incorrect pathname.')
	 call ztell(iretc,0)
	 goto 100
      endif
c
c  Fetch the object from the database file.
c  The link to the keys bank is returned in 'lhepdb' via /sdblnk/
c
      call sdb_read_fetch(chpath,isel,usrkeys,chopt,iretc)
      if ( iretc .eq. IDZT_SDB_BANK_NOT_FOUND ) goto 100
      if( iretc.ne.0 ) then
         write(iqlog,9003) chtype(1:lchtype),chbank,
     +                ipath(ksdb_idx_ipath_banknum)
         write(iqprnt,9003) chtype(1:lchtype),chbank,
     +                ipath(ksdb_idx_ipath_banknum)
 9003    format('SDB_READ: No',a,'bank returned for "',a4,1x,i10,'".')
	 call ztell(iretc,0)
	 goto 100
      endif
 
c
c     Get the effective validity range of valid bank.  Links to the
c     keys banks for the next and previous banks are returned in `lnext'
c     and 'lprev', respectively, via /sdblnk/.
c     (this is not done for a bank retrieved using the serial
c     number because no instant of validity is provided around
c     which to base the calculation of the range).
c
 
      if (isel .ge. 0) then
	 isel_tmp = isel
c-->  If supernova monitoring is enabled, and the desired instant of
c-->  validity is outside the range of the returned bank, then temporarily
c-->  reset the instant of validity to be within the range.
         if(iopt(ksdb_idx_SN_mon) .eq. ksdb_opt_SN_mon) then
	    if (isel .lt. iqq(lhepdb+nsyskeys+1))
     +         isel_tmp = iqq(lhepdb+nsyskeys+1)
	    if (isel .gt. iqq(lhepdb+nsyskeys+2))
     +         isel_tmp = iqq(lhepdb+nsyskeys+2)
         endif
         call sdb_get_evalid(ipath,chpath,chbank,isel_tmp,usrkeys,
     +           startval,endval,start_gap_val,end_gap_val,iretc)
c-->  If supernova monitoring is enabled, and the desired instant of
c-->  validity is outside the effective range of the returned bank, then
c-->  reset the effective validity range to span the instant of validity.
         if(iopt(ksdb_idx_SN_mon) .eq. ksdb_opt_SN_mon) then
	    if (isel .lt. startval) startval = isel
	    if (isel .gt. endval)   endval   = isel + 60	!add 1 minute
         endif
      endif
 
c
c  If the next or previous data bank is required, get that one and
c  calculate its effective validity range.
c
      if(iopt(ksdb_idx_which_bank) .ne. ksdb_opt_valid) then
c-->  Drop the previously retrieved valid keys bank
	 call mzdrop(ihepdiv,lhepdb,'.')
	 lhepdb = 0
 
c-->  Set the instant of validity based on the next or previous keys bank
         if(iopt(ksdb_idx_which_bank) .eq. ksdb_opt_prev) then
	    if (lprev.eq.0) then
	       iretc = IDZT_SDB_BANK_NOT_FOUND
	       goto 100
	    endif
            isel=iqq(lprev+nsyskeys+2)
         else 				! set for next bank
	    if (lnext.eq.0) then
	       iretc = IDZT_SDB_BANK_NOT_FOUND
	       goto 100
	    endif
            isel=iqq(lnext+nsyskeys+1)
	 endif
 
c-->  Drop the next and previous keys banks for the valid bank if they exist
	 if (lnext.ne.0) then
	    call mzdrop(ihepdiv,lnext,'.')
	    lnext = 0
	 endif
         if (lprev.ne.0) then
	    call mzdrop(ihepdiv,lprev,'.')
	    lprev = 0
	 endif
 
c-->  Fetch the next or previous keys and data banks
         chopt = 'F'
         call sdb_read_fetch(chpath,isel,usrkeys,chopt,iretc)
         if ( iretc .eq. IDZT_SDB_BANK_NOT_FOUND ) goto 100
         if( iretc.ne.0 ) then
            write(iqlog,9003) chtype(1:lchtype),chbank,
     +                      ipath(ksdb_idx_ipath_banknum)
            write(iqprnt,9003) chtype(1:lchtype),chbank,
     +                      ipath(ksdb_idx_ipath_banknum)
            call ztell(iretc,0)
            goto 100
         endif
c
c-->  Get the effective validity range of the next or previous bank
 
         call sdb_get_evalid(ipath,chpath,chbank,isel,usrkeys,
     +              startval,endval,start_gap_val,end_gap_val,iretc)
 
      endif
c
c  Verify that the data is uncorrupted (ie. that the checksum for
c  the data part of the bank matches the one stored in the
c  database keys) if verification is enabled.
c
      if(iopt(ksdb_idx_verification).eq.ksdb_opt_verify .and.
     +	 iqq(lhepdb+8).ne.0) then
         call sdb_util_checksum(lqq(lhepdb-1),ihepdiv,isum,iretc)
	 if (iretc.ne.0) then
	    write(iqlog,9004) chtype(1:lchtype),chbank,
     +                  ipath(ksdb_idx_ipath_banknum)
	    write(iqprnt,9004) chtype(1:lchtype),chbank,
     +                   ipath(ksdb_idx_ipath_banknum)
 9004       format(' SDB_READ: Got an error calculating the ',
     + 'checksum for',a,'bank "',a4,1x,i10,'".',/,' Will use the bank ',
     + 'regardless.')
	    iretc=0
	 else
	    if (iqq(lhepdb+8).ne.isum(1) .or.
     +          iqq(lhepdb+9).ne.isum(2)) then
	       iretc=IDZT_SDB_DATA_CORRUPTED
	    endif
	 endif
      else
	 iretc=0
      endif
      if(iretc.eq.IDZT_SDB_DATA_CORRUPTED) then
         write(iqlog,9003) chtype(1:lchtype),chbank,
     +               ipath(ksdb_idx_ipath_banknum)
         write(iqprnt,9003) chtype(1:lchtype),chbank,
     +                ipath(ksdb_idx_ipath_banknum)
         call ztell(IDZT_SDB_DATA_CORRUPTED,0)
      else
c
c  Copy the data part of the retrieved object to the target
c  division as a stand-alone bank and drop the bank in the HEPDB
c  user division.  The data part of the object is pointed to by the
c  first down-link in the keys bank (ie. lqq(lhepdb-1)).
c
         call mzcopy(ihepdiv,lqq(lhepdb-1),idiv,lbank(1),2,'L')
      endif
 
c
c  Get the offset to the target division.
c
      call sdb_util_get_offset(idiv,koffset,iretc)
 
c
c  Store the contents of system key 10 (the bank contents type)
c  in the /key10/ common block.
c
      key10 = iqq(lhepdb+10)
 
c
c  Update the effective validity range and bank insertion date
c  elements of the IPATH array if required (converting from HEPDB
c  or SNODB time to SNOMAN time).  Do the same for the database
c  header (if it exists).
c
      if (isel .ge. 0) then
	 if (iopt(ksdb_idx_validity_inst) .eq. ksdb_opt_date) then
            call sdb_util_upts(ipath(ksdb_idx_ipath_eff_start_date),
     +                   ipath(ksdb_idx_ipath_eff_start_time),
     +                   startval)
            ipath(ksdb_idx_ipath_eff_start_time) = 100*
     +                   ipath(ksdb_idx_ipath_eff_start_time)
 
            call sdb_util_upts(ipath(ksdb_idx_ipath_eff_end_date),
     +                   ipath(ksdb_idx_ipath_eff_end_time),
     +                   endval)
            ipath(ksdb_idx_ipath_eff_end_time) = 100*
     +                   ipath(ksdb_idx_ipath_eff_end_time)
	 else
	    ipath(ksdb_idx_ipath_eff_start_date) = startval
	    ipath(ksdb_idx_ipath_eff_start_time) = 0
	    ipath(ksdb_idx_ipath_eff_end_date) = endval
	    ipath(ksdb_idx_ipath_eff_end_time) = 0
	 endif
	 if (key10 .le. ksdb_ctype_with_hdr_max .and.
     +      iopt(ksdb_idx_bank_type) .ne. ksdb_opt_format) then
	    iqq(koffset + lbank(1) + KDBHD_EFF_START_DATE) =
     +         ipath(ksdb_idx_ipath_eff_start_date)
	    iqq(koffset + lbank(1) + KDBHD_EFF_START_TIME) =
     +         ipath(ksdb_idx_ipath_eff_start_time)
	    iqq(koffset + lbank(1) + KDBHD_EFF_END_DATE) =
     +         ipath(ksdb_idx_ipath_eff_end_date)
	    iqq(koffset + lbank(1) + KDBHD_EFF_END_TIME) =
     +         ipath(ksdb_idx_ipath_eff_end_time)
	 endif
      endif
c
c  Here we have to use the hepdb unpacking routine as the entry
c  time is packed by hepdb keydbs(4)
 
      call cduptm(ipath(ksdb_idx_ipath_entry_date),
     +                      ipath(ksdb_idx_ipath_entry_time),
     +                      iqq(lhepdb+4),iretc)
      if (ipath(ksdb_idx_ipath_entry_date) .ge. 360000) then
         ipath(ksdb_idx_ipath_entry_date) = 19000000 +
     +                      ipath(ksdb_idx_ipath_entry_date)
      else
         ipath(ksdb_idx_ipath_entry_date) = 20000000 +
     +                      ipath(ksdb_idx_ipath_entry_date)
      endif
      ipath(ksdb_idx_ipath_entry_time) = 10000*
     +                      ipath(ksdb_idx_ipath_entry_time)
      if (key10 .le. ksdb_ctype_with_hdr_max .and.
     +   iopt(ksdb_idx_bank_type) .ne. ksdb_opt_format) then
	 iqq(koffset + lbank(1) + KDBHD_ENTRY_DATE) =
     +      ipath(ksdb_idx_ipath_entry_date)
	 iqq(koffset + lbank(1) + KDBHD_ENTRY_TIME) =
     +      ipath(ksdb_idx_ipath_entry_time)
      endif
 
c
c  Fix bank name and number if wrong (this can happen if the data
c  was compressed when stored).
c
      call uctoh('USER',iuser,4,4)
      if (iqq(koffset+lbank(1)-4) .eq. iuser) then
         iqq(koffset+lbank(1)-4) = ipath(ksdb_idx_ipath_bankname)
         iqq(koffset+lbank(1)-5) = ipath(ksdb_idx_ipath_banknum)
      endif
c
c  Clean up temporary banks.
c
 100  continue
      if (lhepdb .ne. 0) call mzdrop(ihepdiv,lhepdb,'.')
      if (lnext .ne. 0) call mzdrop(ihepdiv,lnext,'.')
      if (lprev .ne. 0) call mzdrop(ihepdiv,lprev,'.')
c
      return
      end
 
*ENDFILE MEMBER=SDB_READ
