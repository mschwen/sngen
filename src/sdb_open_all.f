*FILE MEMBER=SDB_OPEN_ALL LIBRARY=SNOMAN LANGUAGE=FORTRAN77 DATE=21:Apr:1997
      subroutine sdb_open_all(iretc)
c
c	This routine opens every database file.
c
      implicit none
c
      include 'db.inc'
      include 'db_io.inc'
      include 'db_misc.inc'
      include 'mzca.inc'
c
      integer iretc			!return code
      integer i,iprefix
      integer num_dbfiles,max_num_dbfiles
      parameter (max_num_dbfiles=30)
      character*2 chprefix(max_num_dbfiles),chpreftmp
      common /prefices/ num_dbfiles,chprefix
      common/slate/is(40)
      integer is,ldir,lfile,ltmp
      character*80 chdir
      character*99 chtmp
      character*99 list
      integer llist
      character*8 chopt
c
      nwp=kshl_pawc_size
      iretc=0
c
c  See if Zebra has been initialized.  If so, then assume that PAWC
c  store has also been initialized, and set the options for the call
c  to CDPAW accordingly.
c
      if (INKEYS(1).ne.MQKEYS(1) .or. INKEYS(2).ne.MQKEYS(2) .or.
     +    INKEYS(3).ne.MQKEYS(3)) then
	 chopt='ZPU'
      else
	 chopt='U'
      endif
c
c  Set up a user division 'ihepdiv' for HEPDB.
c
      call cdpaw(kshl_pawc_size,nhbook,ihepdiv,'HEP-DIV',500000,1500000,
     *          chopt,irc)
      lunfz=2
      luntmp=20
      luncd=50
      lrecl=1024      !use the default record size (1024 bytes)
c
c  Get all of the prefices from `hepdb.names'.  First, find the
c  names file.
c
      call getenvf('CDSERV',chdir)
      ldir = is(1)
*IF VMS
*$      chfile = chdir(1:ldir) // 'hepdb.names'
*$      lfile = ldir + 11
*ELSE
*IF NT
*$      chfile = chdir(1:ldir) // '\\hepdb.names'
*ELSE
      chfile = chdir(1:ldir) // '/hepdb.names'
*ENDIF
      lfile = ldir + 12
      call cutol(chfile(1:lfile))
*ENDIF
c
c  Loop through the names file, and open the file for each database
c  prefix found.  The names file has the following rough layout:
c
c :nick.config
c    :list.as av ch cm cp cr ...        ! list of database prefices
c  ...
c :nick.as                      ! info on database file with prefix `as'
c :file./snodb/dbs/assay.dbs
c ...
c
c  So, i) get the list of database prefices from the `:list' line
c  under the `:nick.config' line, and then ii) add each prefix
c  from the list to the `chprefix' array if there is a corresponding
c  `:nick.XX' entry further below in the names file.
c
      num_dbfiles=0
      open(unit=luntmp,file=chfile(1:lfile),status='old')
 1    read(luntmp,'(a99)',end=3) chtmp
      call cutol(chtmp)
      iprefix=index(chtmp(1:99),':nick.')
      if (iprefix.ne.0) then
         if (chtmp(iprefix+6:iprefix+12).eq.'config') then
            read(luntmp,'(a99)') chtmp
            ltmp=lenocc(chtmp)
            iprefix=index(chtmp(1:ltmp),':list.')
            if (iprefix.ne.0) then		!get list of prefices
               list=chtmp(iprefix+6:ltmp)
               llist=ltmp-iprefix-5
               call cltou(list)
            endif
         else
            chpreftmp=chtmp(iprefix+6:iprefix+7)
            call cltou(chpreftmp)
            iprefix=index(list(1:llist),chpreftmp)
            if (iprefix.eq.0) goto 1		!not in prefix list
            num_dbfiles=num_dbfiles+1
            if (num_dbfiles.gt.max_num_dbfiles) then
               write(*,2) max_num_dbfiles
 2     format(' More than ',i2,' prefices found in the',
     + ' hepdb.names file; ignoring the rest.',/,' Increase the',
     + ' maximum number allowed in SDB_OPEN_ALL.F')
               num_dbfiles=max_num_dbfiles
               goto 3
            endif
            chprefix(num_dbfiles)=chpreftmp
         endif
      endif
      goto 1
 3    continue
      close(unit=luntmp)
c
c  For each prefix in the `chprefix' array, obtain the top directory
c  name and full file name for each database file.  If successful,
c  open the database file (and, thus, load its alias dictionary into
c  the HEPDB system division).
c
      do 10 i=1,num_dbfiles
         call cdpref(luntmp,chprefix(i)(1:2),chtop,chfile,irc)
         if(irc.le.4) then
            call cdopen(luncd,lunfz,chtop,chfile,lrecl,ihepdiv,' ',irc)
         else
            write(*,*)' sdb_open_all: error from cdpref for ',
     +                chprefix(i)(1:2)
	    write(*,*)' cdpref returned a code of:',irc
	    iretc=irc
         endif
         luncd=luncd+1
 10   continue
c
      return
      end
*ENDFILE MEMBER=SDB_OPEN_ALL
