	subroutine read_titles

	include 'z.inc'
	REAL              Z_ARRAY(100)
	EQUIVALENCE       ( IZSTOR , Z_ARRAY )

	include 'mt.inc'
	integer              mt_array(100)
	equivalence       ( mt_fence , mt_array )

	include 'memory_size.inc'
	include 'sdb.inc'
	include 'sdb_ipath.inc'
	include 'sdb_options.inc'
	include 'bank_mctk.inc'
	include 'id_particles.inc'
	include 'id_en_types.inc'
        include 'titles_ccxs.inc'
        include 'titles_ccxb.inc'
        include 'titles_ccxt.inc'
	include 'titles_ncxs.inc'
	include 'titles_ncxb.inc'
	include 'titles_ncxt.inc'


	integer irc
	integer ipath(num_ipath_args)
	character*4 chbank
	integer name(5)

	call define_memory_size

**    Initialise ZEBRA.

      CALL MZEBRA(-1)

      IQLOG  = 7
      IQPRNT = 6

      CALL MZVERS

      CALL MZSTOR( IZSTOR ,
     +             '/Z/' ,
     +             ' ' ,
     +             Z_FENCE ,
     +             Z_STR_LINKS ,
     +             Z_REF_LINKS ,
     +             Z_WORK ,
     +             Z_ARRAY(ISHL_Z_SIZE_LIM2) ,
     +             Z_ARRAY(ISHL_Z_SIZE) )

**    Create phony TCKP stand-alone bank with no links and one data word.

	chbank='TCKP'
	call uctoh(chbank,name(1),4,4)	! name TCKP
	name(2)=0			! no links
	name(3)=0			! no structural links
	name(4)=1			! 1 data word
	name(5)=2			! integer data

	call mzlift(izstor, lmtckp, 0, 2, name, 0)

**    Initialize MT

        call mzstor( mtstor ,
     +               '/MT/' ,
     +               '.' ,
     +               mt_fence ,
     +               mt_str_links ,
     +               mt_ref_links ,
     +               mt_work ,
     +               mt_array(ishl_mt_size_lim2) ,
     +               mt_array(ishl_mt_size) )
 
        mtdiv = mtstor + 2

**    Initialize SDB

	call hlimit(-ishl_pawc_size)
	call sdb_open_all(irc)
	call mzlink(ixpawc,'/sdblnk/',lhepdb,lprev,lhepdb)

**    Get titles banks.

	call sdb_util_put_option(ksdb_idx_validity_inst,ksdb_opt_date,
     + irc)
	call sdb_util_put_option(ksdb_idx_which_bank,ksdb_opt_valid,
     + irc)
	call sdb_util_put_option(ksdb_idx_bank_type,ksdb_opt_data,
     +  irc)
	ipath(ksdb_idx_ipath_date)=19990601
	ipath(ksdb_idx_ipath_time)=12000000
	ipath(ksdb_idx_ipath_formatnum)=-1
	ipath(ksdb_idx_ipath_datatype)=0
	ipath(ksdb_idx_ipath_tasktype)=-1
	ipath(ksdb_idx_ipath_source_id)=-1
	ipath(ksdb_idx_ipath_banknum)=1

*****************************************************************
* cross section titles
* (see SNOMAN Companion titles banks for mneumonics)
*****************************************************************

* CC total cross sections (not used to calculate total SN events)
c        chbank='CCNA'   ! NSAPMGK
        chbank='CCXT'   ! KN
	call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
	call sdb_read(ipath, lccxt, mtdiv, irc)
	ldccxt = lccxt + 30

* CC differential cross sections (nu + d)
c        chbank='CCNP'   ! NSAPMGK (extend to Enu=170 MeV)
        chbank='CCXS'   ! KN
        call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
        call sdb_read(ipath, lccxs, mtdiv, irc)
        ldccxs = lccxs + 30

* CC differential cross sections (nu_bar + d)
c        chbank='CCNB'  ! NSAPMGK (extend to Enu=170 MeV)
        chbank='CCXB'  ! KN
	call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
	call sdb_read(ipath, lccxb, mtdiv, irc)
	ldccxb = lccxb + 30

* NC total cross sections (not used to calculate total SN events)
c        chbank='NCNA'  ! NSAPMGK
	chbank='NCXT'  ! KN
	call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
	call sdb_read(ipath, lncxt, mtdiv, irc)
	ldncxt = lncxt + 30

* NC differential cross sections (nu + d)
c        chbank='NCNP'  ! NSAPMGK (do not extend past Enu=50 MeV)
	chbank='NCXB'  ! KN
	call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
	call sdb_read(ipath, lncxb, mtdiv, irc)
	ldncxb = lncxb + 30

* NC differential cross sections (nu_bar + d)
c        chbank='NCNB'  ! NSAPMGK (do not extend past Enu=50 MeV)
	chbank='NCXS'  ! KN
	call uctoh(chbank,ipath(ksdb_idx_ipath_bankname),4,4)
	call sdb_read(ipath, lncxs, mtdiv, irc)
	ldncxs = lncxs + 30

**    Close database files.

	call cdend(' ','AC',irc)

	return
	end

