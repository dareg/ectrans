program transform_test
	use parkind1,only: jpim,jprb,jprd
	use oml_mod,only: oml_max_threads
	use mpl_mpif
	use mpl_module
#ifdef GRIB_API
	use grib_api
#endif
	use yomgstats,only: jpmaxstat,ylstats=>lstats

	implicit none

	integer(kind=jpim),parameter :: nfldx=412,nrgri=1280
	integer(kind=jpim) :: nerr,nlin,insf,nsmax,ndgl,nq,ii,nout,noutdump,nspec2,ngptot,&
		ngptotg,iflds,jroc,jb,itag,nspec2g,i,ibl,iret
	integer :: jstep,rgri(nrgri)
	integer(kind=jpim) :: nstats_mem,ntrace_stats,nprnt_stats,nprintnorms,niter,&
		nmax_resol,nprintlev,npromatr,ncombflen,nproc,nthread,nprgpns,nprgpew,nprtrv,nprtrw,&
		nspecresmin,mysetv,mysetw,myseta,mysetb,mp_type,mbx_size,ivsetsc(1),npsp,nflevg,&
		nflevl,isqr,nproma,ngpblks,iprtrv,iprtrw,iprsp,ilevpp,irest,ilev,jl,ilastlev,&
		ndimgmv,ndimgmvs,nuv,myproc,jj
	integer(kind=jpim),allocatable :: nloen(:),ito(:),nprcids(:),numll(:),ivset(:),&
		npsurf(:)
	logical :: lstack,ldone,luserpnm,lkeeprpnm,luseflt,ltrace_stats,lstats_omp,&
		lstats_comms,lstats_mpl,lstats,lbarrier_stats,lbarrier_stats2,ldetailed_stats,&
		lstats_alloc,lsyncstats,lstatscpu,lstats_mem,lxml_stats,lfftw,lmpoff,lsync_trans,&
		leq_regions,llinfo,lvorgp,ldivgp,luvder,lscders,lresetgp
	character(len=1) :: ctypeg
	character(len=16) :: cgrid
	character(len=127) :: cinsf,crt,crtable,cfname
	real(kind=jprb) :: znormsp(1),znormsp1(1),zmaxerr(4)
	real(kind=jprd) :: ztinit,ztloop,timef,ztstepmax,ztstepmin,ztstepavg,&
		ztstepmax1,ztstepmin1,ztstepavg1,ztstepmax2,ztstepmin2,ztstepavg2,&
		zaveave(0:jpmaxstat)
	real(kind=jprb),parameter :: zra=6371229
	real(kind=jprd),allocatable :: ztstep(:),ztstep1(:),ztstep2(:)
	real(kind=jprb),allocatable :: znormdiv(:),znormdiv1(:),&
		znormvor(:),znormvor1(:),znormt(:),znormt1(:)
	real(kind=jprb),allocatable,target :: zgmv(:,:,:,:),zgmvs(:,:,:),sp3d(:,:,:),zsp(:,:)
	real(kind=jprb),pointer :: zvor(:,:),zdiv(:,:),zt(:,:,:),zuv(:,:,:,:),zgpt(:,:,:,:),&
		zwind(:,:,:,:)

	namelist/namrgri/ndgl,nlin,nq,nsmax,rgri
	!namelist/namtrans/nprintnorms,niter,nproma,npromatr,nspecresmin,mbx_size
	namelist/namstats/lstats,lstack,lbarrier_stats,lbarrier_stats2,ltrace_stats,&
		ldetailed_stats,lstats_alloc,lstats_comms,lstats_mpl,lstats_omp,lstats_mem,&
		lxml_stats,lsyncstats,lstatscpu,nprnt_stats
	namelist/namtrans/lfftw,luserpnm,lkeeprpnm,luseflt,lsync_trans,leq_regions,&
		lvorgp,ldivgp,luvder,lscders,nprintlev,nprintnorms

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "dist_spec.h"
#include "gath_grid.h"
#include "trans_inq.h"
#include "gpnorm_trans.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"

	ztinit = timef()

	call mpl_init()

	itag = 123456

	nproc = mpl_nproc()
	myproc = mpl_myrank()
	nthread = oml_max_threads()

	if (nproc < 1) call abor1("Error: nproc=0")
	if (nproc == 1) lmpoff = .true.

	nout = 6
	if (myproc /= 1) open(nout,file="/dev/null")

	ztinit = (timef()-ztinit)/1000
	write(nout,"('Startup time:',F9.3,'s')") ztinit
	ztinit = timef()

	write(nout,*) "MPI/OMP setting:",nproc,"x",nthread

	nerr = 0
	noutdump = 7
	ncombflen = 1800000
	lmpoff = .false.
	mp_type = 2
	mbx_size = 150000000

	nstats_mem = 0
	ntrace_stats = 0

	write(nout,*) "Read namstat"
	lstats = .false.
	lstack = .false.
	lbarrier_stats = .false.
	lbarrier_stats2 = .false.
	ldetailed_stats = .false.
	lstats_alloc = .false.
	lstats_comms = .false.
	lstats_mpl = .false.
	lstats_omp = .false.
	lstats_mem = .false.
	lxml_stats = .false.
	lsyncstats = .false.
	lstatscpu = .false.
	ltrace_stats = .false.
	nprnt_stats = 1
	open(unit=4,file="stats.nam",form="formatted")
	read(4,namstats)
	close(4)
	write(nout,*) "lstats/ldetailed_stats:",lstats,ldetailed_stats
	write(nout,*) "ltrace_stats/lstatscpu/nprnt_stats:",ltrace_stats,lstatscpu,nprnt_stats

	if (ldetailed_stats) then
		write(nout,*) "Info: detailed stats asked for, forces lstats_omp, lstats_comms, &
			&lstats_mpl, lstatscpu on and prints switched on on all tasks"
		lstats_omp = .true.
		lstats_comms = .true.
		lstats_mpl = .true.
		lstatscpu = .true.
		nprnt_stats = nproc
	end if

	write(nout,*) "Read namtrans"
	luserpnm = .false.
	lkeeprpnm = .false.
	luseflt = .false.
	lfftw = .true.
	leq_regions = .true.
	lsync_trans = .true.
	lvorgp = .false.
	ldivgp = .true.
	luvder = .false.
	lscders = .true.
	nprintlev = 1
	nprintnorms = 1
	npromatr = 0
	nspecresmin = 0
	open(unit=4,file="spectraltrans.nam",form="formatted")
	read(4,namtrans)
	close(4)
	write(nout,*) "lfftw/l(use/keep)rpnm/luseflt:",lfftw,luserpnm,lkeeprpnm,luseflt
	write(nout,*) "leq_regions/npromatr/nspecresmin:",leq_regions,npromatr,nspecresmin
	write(nout,*) "lvorgp/ldivgp/luvder/lscders:",lvorgp,ldivgp,luvder,lscders

	! defaults: cubic octahedral grid,...
	cgrid = ""
	call parseargs(niter,nproma,nflevg,cgrid)

	if (niter <= 0) call abor1("Error: niter <= 0")
	if (nflevg <= 0) call abor1("Error: nflevg <= 0")

	write(nout,*) "Set grid and levels"
	inquire(file="rgri.nam",exist=llinfo)
	if (llinfo) then
		write(nout,*) "Read namrgri"
		nsmax = 0
		ndgl = 0
		nlin = 1
		nq = 0
		rgri(:) = 0

		open(unit=4,file="rgri.nam",form="formatted")
		read(4,namrgri)
		close(4)
		write(nout,*) "nsmax/ndgl/nlin/nq:",nsmax,ndgl,nlin,nq

		if (ndgl <= 0.or.ndgl > nrgri) call abor1("Error: ndgl <= 0 or ndgl > nrgri")

		if (len_trim(cgrid) == 0) then
			if (ndgl <= 0) call abor1("Error: ndgl <= 0 and no grid tag")
			if (any(rgri(1:ndgl) <= 0)) call abor1("Error: rgri <= 0")

			allocate(nloen(ndgl))
			nloen(:) = rgri(1:ndgl)
		else
			if (ndgl > 0) write(nout,*) "Info: resetting ndgl as to grid tag ",trim(cgrid)
			call parse_grid(cgrid,ndgl,nloen)
		end if
#ifdef GRIB_API
	else
		if (myproc == 1) then
			cinsf = ""
			if (cinsf == "") cinsf = "fort.11"

			write(nout,*) "Read GRIB file info:",cinsf
			inquire(file=cinsf,exist=llinfo)
			if (.not.llinfo) call abor1("Error: rgri <= 0")

			call grib_open_file(insf,cinsf,"R",iret)
			if (iret /= grib_success) call abor1("ERROR OPENING FILE INPUT SPECTRAL FILE")

			call readgrib_resol(insf,nsmax)
			call grib_close_file(insf)

			if (ndgl == 0) call setndgl(nlin,nq,nsmax,ndgl)
			call check_ndgl(nlin,ndgl)
		end if

		if (nproc > 1) call mpl_broadcast(ndgl,itag,1,cdstring='transform_test:')

		allocate(nloen(ndgl))

		ctypeg = "r"
		write(nout,*) "--> grid type",ctypeg
		if (.not.(ctypeg == "r".or.ctypeg == "f")) call abor1("WRONG GRID TYPE")

		if (myproc == 1) then
			if (ctypeg == "f") then
				nloen(:) = 2*ndgl
			else
				if (nlin == 1) then
					crt = "/rtablel_2"
				else if (nq == 1) then
					crt = "/rtable_3"
				else if (nq == 2) then
					crt = "/rtable_4"
				else
					crt = "/rtable_2"
				end if

				crtable = "."
				ii = index(crtable," ")
				if (nsmax < 1000) then
					write(crtable(ii:ii+len_trim(crt)+2),"(A,I3.3)") trim(crt),nsmax
				else
					write(crtable(ii:ii+len_trim(crt)+3),"(A,I4.4)") trim(crt),nsmax
				end if

				write(nerr,*) "Read namrgri in nloen:",crtable
				open(15,file=crtable,form="FORMATTED",action="READ")
				read(15,namrgri)
				close(15)
			end if
		end if

		if (nproc > 1) then
			call mpl_broadcast(nsmax,itag,1,cdstring='transform_test:')
			call mpl_broadcast(nloen,itag,1,cdstring='transform_test:')
		end if
#endif
	end if

	write(nout,*) "Set up distribution (tasks and sets)"
	allocate(nprcids(nproc))
	do jj=1,nproc
		nprcids(jj) = jj
	end do

	write(nout,*) "Set up distribution: A/B sets"
	if (nproc == 1) then
		nprgpns = 1
		nprgpew = 1
	else
		call setab(nproc,nprgpns,nprgpew)
		if (nprgpns*nprgpew /= nproc) call abor1("Error: no partition found for A/B sets")
	end if

	if (nspecresmin == 0) nspecresmin = nproc

	write(nout,*) "Set trans distribution: W/V sets"
	if (nproc == 1) then
		nprtrv = 1
		nprtrw = 1
	else
		do i=2,nproc-1
			nprtrv = i
			nprtrw = nproc/i
			if (nprtrv*nprtrw /= nproc) cycle

			if (nprtrv > nprtrw) exit
			if (nprtrw > nspecresmin) cycle
			if (nprtrw <= nspecresmin/(2*oml_max_threads())) exit
		end do

		if (nprtrv*nprtrw /= nproc.or.nprtrw > nspecresmin.or.nprtrv > nprtrw) then
			call setab(nproc,nprtrw,nprtrv)
		end if

		if (nprtrw*nprtrv /= nproc) call abor1("NPRTRW*NPRTRV /= NPROC")
		if (nprtrw > nspecresmin) call abor1("NPRTRW > NSPECRESMIN")
	end if

	write(nout,*) "Distribution values (NS EW W V):",nprgpns,nprgpew,nprtrw,nprtrv

	if (lmpoff) then
		mysetw = (myproc-1)/nprtrv+1
		mysetv = mod(myproc-1,nprtrv)+1
	else
		call mpl_groups_create(nprtrw,nprtrv)
		call mpl_cart_coords(myproc,mysetw,mysetv)
		iprtrv = mod(myproc-1,nprtrv)+1
		iprtrw = (myproc-1)/nprtrv+1
		if (iprtrv /= mysetv.or.iprtrw /= mysetw)&
			call abor1("inconsistency when computing MYSETW and MYSETV")

		llinfo = myproc == 1
		call mpl_buffer_method(kmp_type=mp_type,kmbx_size=mbx_size,kprocids=nprcids,&
			ldinfo=llinfo)
	end if

	allocate(numll(nprtrv+1),npsurf(nprtrv))

	ilevpp = nflevg/nprtrv
	irest = nflevg-ilevpp*nprtrv
	do jroc=1,nprtrv
		if (jroc <= irest) then
			numll(jroc) = ilevpp+1
		else
			numll(jroc) = ilevpp
		end if
	end do

	iprsp = min(nflevg+1,nprtrv)
	numll(iprsp+1:nprtrv+1) = 0

	nflevl = numll(mysetv)

	npsurf(1:iprsp) = 0
	npsurf(iprsp) = 1
	npsp = npsurf(mysetv)
	ivsetsc(1) = iprsp

	write(nout,*) "Setup transforms",nproc
	nmax_resol = 37
	call setup_trans0(kout=nout,kerr=nerr,kprintlev=nprintlev,kmax_resol=nmax_resol,&
		kpromatr=npromatr,kprgpns=nprgpns,kprgpew=nprgpew,kprtrw=nprtrw,kcombflen=ncombflen,&
		ldmpoff=lmpoff,ldsync_trans=lsync_trans,ldeq_regions=leq_regions,prad=zra,&
		ldalloperm=.true.)

	call setup_trans(nsmax,ndgl,kloen=nloen,ldsplit=.true.,kflev=nflevl,&
		ldusefftw=lfftw,lduserpnm=luserpnm,ldkeeprpnm=lkeeprpnm,lduseflt=luseflt)

	call trans_inq(kspec2=nspec2,kspec2g=nspec2g,kgptot=ngptot,kgptotg=ngptotg)

	if (nproma == 0) nproma = ngptot

	ngpblks = (ngptot-1)/nproma+1

	allocate(sp3d(nflevl,nspec2,3),zsp(1,nspec2))

	write(nout,*) "Initialize spectral data"
	call initspec(nsmax,zsp,sp3d)

	zvor => sp3d(:,:,1)
	zdiv => sp3d(:,:,2)
	zt => sp3d(:,:,3:3)

	write(nout,*) "START OF RUNTIME PARAMETERS"
	write(nout,'("NLIN/NQ/NSMAX/NDGL:",4(x,I10))') nlin,nq,nsmax,ndgl
	write(nout,'("NPROC/NTHREAD/NPRGPNS/NPRGPEW:",4(x,I10))') nproc,nthread,nprgpns,nprgpew
	write(nout,'("NPRTRW/NPRTRV:",2(x,I10))') nprtrw,nprtrv
	write(nout,'("NPROMA/NGPBLKS:",2(x,I10))') nproma,ngpblks
	write(nout,'("NGPTOT/NGPTOTG/NFLEVG:",3(x,I10))') ngptot,ngptotg,nflevg
	write(nout,'("LUSEFLT=",L10)') luseflt
	write(nout,'("NSPEC2/NSPEC2G:",3(x,I10))') nspec2,nspec2g

	allocate(ivset(nflevg))

	ilev = 0
	do jb=1,nprtrv
		ivset(ilev+1:ilev+numll(jb)) = jb
		ilev = ilev+numll(jb)
	end do

#ifdef GRIB_API
	inquire(file=cinsf,exist=llinfo)
	if (llinfo) then
		call grib_open_file(insf,cinsf,"R",iret)
		if (iret /= grib_success) call abor1("ERROR OPENING FILE INPUT SPECTRAL FILE")

		call readgrib(insf,nflevg,nspec2g,ivset,ivsetsc,zsp,zvor,zdiv,zt)
	end if
#endif

	write(nout,*) "allocate local GMV data"
	! reset all GP fields, even those not active on SP fields (ie unused in dir_trans)
	lresetgp = .false.
	ndimgmvs = 1
	if (lscders)  ndimgmvs = 3
	allocate(zgmvs(nproma,ndimgmvs,ngpblks))
	nuv = 4
	if (luvder) nuv = 8
	if (.true.) then
		ndimgmv = 5
		allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
		allocate(zwind(nproma,nflevg,nuv,ngpblks))
		if (lscders) then
			zgpt => zgmv(:,:,1:3,:)
		else
			zgpt => zgmv(:,:,1:1,:)
		end if
	else
		ndimgmv = 5+nuv
		allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
		zwind => zgmv(:,:,1:nuv,:)
		if (lscders) then
			zgpt => zgmv(:,:,5:7,:)
		else
			zgpt => zgmv(:,:,5:5,:)
		end if
	end if

	zgmv(:,:,:,:) = -2
	zwind(:,:,:,:) = -1
	zgmvs(:,:,:) = -3

	if (lvorgp.and.ldivgp) then
		zuv => zwind(:,:,3:4,:)
	else if (lvorgp.or.ldivgp) then
		zuv => zwind(:,:,2:3,:)
	else
		zuv => zwind(:,:,1:2,:)
	end if

	allocate(znormvor(nflevg),znormvor1(nflevg))
	allocate(znormdiv(nflevg),znormdiv1(nflevg))
	allocate(znormt(nflevg),znormt1(nflevg))

	if (nprintnorms > 0) then
		write(nout,*) "Spectral norms at init:"
		call spnorm(1,zsp,ivsetsc,znormsp,"SP")
		call spnorm(nflevg,zvor,ivset,znormvor,"VOR")
		call spnorm(nflevg,zdiv,ivset,znormdiv,"DIV")
		call spnorm(nflevg,zt(:,:,1),ivset,znormt,"T")
	end if

	ztinit = (timef()-ztinit)/1000
	write(nout,"('Initialisation time:',F9.3,'s')") ztinit

	allocate(ztstep(niter),ztstep1(niter),ztstep2(niter))

	write(nout,*) "START OF SPEC TRANSFORMS"

	if (lstats) then
		call gstats(0,0)
		call gstats_setup(nproc,myproc,nprcids,lstats,lstatscpu,lsyncstats,ldetailed_stats,&
			lbarrier_stats,lbarrier_stats2,lstats_omp,lstats_comms,lstats_mem,nstats_mem,&
			lstats_alloc,ltrace_stats,ntrace_stats,nprnt_stats,lxml_stats)
		call gstats_psut
		!call gstats_label_ifs
	end if

	ztloop = timef()

	ylstats = .false.

	do jstep=1,niter
		write(nout,*) ". iter",jstep
		if (jstep > 1) ylstats = .true.

		ztstep(jstep) = timef()
		ztstep1(jstep) = timef()

		call inv_trans(pspvor=zvor,pspdiv=zdiv,pspsc2=zsp,pspsc3a=zt,ldscders=lscders,&
			ldvorgp=lvorgp,lddivgp=ldivgp,lduvder=luvder,kresol=1,kproma=nproma,&
			kvsetuv=ivset,kvsetsc2=ivsetsc,kvsetsc3a=ivset,pgpuv=zwind,pgp2=zgmvs,pgp3a=zgpt)
		ztstep1(jstep) = (timef()-ztstep1(jstep))/1000

		if (nprintnorms > 0) then
			call gpnorms2d(zgmvs)
			call gpnorms3d(zgpt,zwind,zuv)
		end if

		if (jstep == 1) then
			write(nout,*) "--> step 1, reset GMV to constant values"
         zgmvs(:,:,:) = 0
			do ibl=1,ngpblks
				zgmvs(:,1,ibl) = 11+mod(ibl+myproc,2)/10.
			end do
			if (lresetgp.and.lscders) then
				zgmvs(:,2,:) = .5+mod(myproc,2)/100.
				zgmvs(:,3,:) = .25+mod(myproc,2)/200.
			end if
         zgmv(:,:,:,:) = 0
         zwind(:,:,:,:) = 0
			do jl=1,nflevg
				do ibl=1,ngpblks
					zuv(:,jl,1,ibl) = 1+mod(ibl+myproc,2)+jl
					zuv(:,jl,2,ibl) = mod(ibl+myproc,2)+sign(jl,(nflevg+1)/2-jl)
					zgpt(:,jl,1,ibl) = 280+2*mod(ibl+myproc,2)+10*(jl-(nflevg+1)/2)
					if (lresetgp) zgpt(:,jl,2,ibl) = ibl+50+mod(myproc,2)+2*(jl-(nflevg+1)/2)
					if (lresetgp) zgpt(:,jl,3,ibl) = ibl+mod(myproc,2)+(jl-(nflevg+1)/2)
				end do
			end do

			if (lresetgp) then
				if (lvorgp) zwind(:,:,1,:) = .02*mod(myproc,2)
				if (ldivgp) then
					if (lvorgp) then
						zwind(:,:,2,:) = .01*mod(myproc,2)
					else
						zwind(:,:,1,:) = .01*mod(myproc,2)
					end if
				end if
			end if
			call gpnorms2d(zgmvs)
			call gpnorms3d(zgpt,zwind,zuv)
		end if

		if (.false.) then
		call dump_gp_field(noutdump,jstep,myproc,nproma,ngpblks,"P",zgmvs(:,1,:))
		call dump_gp_field(noutdump,jstep,myproc,nproma,ngpblks,"U",zuv(:,nflevg,1,:))
		call dump_gp_field(noutdump,jstep,myproc,nproma,ngpblks,"V",zuv(:,nflevg,2,:))
		call dump_gp_field(noutdump,jstep,myproc,nproma,ngpblks,"T",zgmv(:,nflevg,1,:))
		end if

		ztstep2(jstep) = timef()
		call dir_trans(pspvor=zvor,pspdiv=zdiv,pspsc2=zsp,pspsc3a=zt,kresol=1,&
			kproma=nproma,kvsetuv=ivset,kvsetsc2=ivsetsc,kvsetsc3a=ivset,&
			pgpuv=zuv,pgp2=zgmvs(:,1:1,:),pgp3a=zgpt(:,:,1:1,:))
		ztstep2(jstep) = (timef()-ztstep2(jstep))/1000

		ztstep(jstep) = (timef()-ztstep(jstep))/1000

		write(nout,'(". time at step ",I6,":", F8.4)') jstep,ztstep(jstep)

		if (nprintnorms > 0) then
			call spnorm(1,zsp,ivsetsc,znormsp,"SP")
			call spnorm(nflevg,zvor,ivset,znormvor,"VOR")
			call spnorm(nflevg,zdiv,ivset,znormdiv,"DIV")
			call spnorm(nflevg,zt(:,:,1),ivset,znormt,"T")

			if (jstep == 1) then
				znormsp1 = znormsp
				znormdiv1(:) = znormdiv(:)
				znormvor1(:) = znormvor(:)
				znormt1(:) = znormt(:)
			end if

			if (myproc == 1) then
				zmaxerr(1) = sperror(1,znormsp,znormsp1)
				zmaxerr(2) = sperror(nflevg,znormdiv,znormdiv1)
				zmaxerr(3) = sperror(nflevg,znormvor,znormvor1)
				zmaxerr(4) = sperror(nflevg,znormt,znormt1)

				write(nout,"('. max error:',4(x,e10.3))") zmaxerr(:)
			end if
		end if
	end do

	ztloop = (timef()-ztloop)/1000

	if (nprintnorms == 0) then
		write(nout,*) "Spectral norms:"
		call spnorm(1,zsp,ivsetsc,znormsp,"SP")
		call spnorm(nflevg,zvor,ivset,znormvor,"VOR")
		call spnorm(nflevg,zdiv,ivset,znormdiv,"DIV")
		call spnorm(nflevg,zt(:,:,1),ivset,znormt,"T")
	end if

	if (myproc == 1) then
		zmaxerr(1) = sperror(1,znormsp,znormsp1)
		zmaxerr(2) = sperror(nflevg,znormdiv,znormdiv1)
		zmaxerr(3) = sperror(nflevg,znormvor,znormvor1)
		zmaxerr(4) = sperror(nflevg,znormt,znormt1)

		write(nout,'("SURFACE PRESSURE MAX ERROR=",E10.3)') zmaxerr(1)
		write(nout,'("DIVERGENCE MAX ERROR=",E10.3)') zmaxerr(2)
		write(nout,'("VORTICITY MAX ERROR=",E10.3)') zmaxerr(3)
		write(nout,'("TEMPERATURE MAX ERROR=",E10.3)') zmaxerr(4)
		write(nout,'("GLOBAL MAX ERROR=",E10.3)') maxval(zmaxerr)
	end if

	write(nout,*) "Statistics for iterations"
	call mpl_allreduce(ztloop,"SUM",ldreprod=.false.)

	call stepstat(ztstep,ztstepavg,ztstepmin,ztstepmax)
	call stepstat(ztstep1,ztstepavg1,ztstepmin1,ztstepmax1)
	call stepstat(ztstep2,ztstepavg2,ztstepmin2,ztstepmax2)

	if (myproc == 1) then
		ztloop = ztloop/nproc
		ztstepavg = ztstepavg/nproc/niter
		ztstepavg1 = ztstepavg1/nproc/niter
		ztstepavg2 = ztstepavg2/nproc/niter
		ztstep(:) = ztstep(:)/nproc
		ztstep1(:) = ztstep1(:)/nproc
		ztstep2(:) = ztstep2(:)/nproc

		call sort(ztstep,niter)
		call sort(ztstep1,niter)
		call sort(ztstep2,niter)

		write(nout,"(/,'Time step statistics (min/avg/max)')")
		write(nout,"(' Inverse transforms:',3(x,f11.4))") ztstepmin1,ztstepavg1,ztstepmax1
		write(nout,"(' Direct transforms: ',3(x,f11.4))") ztstepmin2,ztstepavg2,ztstepmax2
		write(nout,"(' Total transforms:  ',3(x,f11.4))") ztstepmin,ztstepavg,ztstepmax
		write(nout,"(' Iteration loop :',f11.4)") ztloop
	end if

	if (lstack) call printstack(nproc,nprcids)

	if (lstats) then
		call gstats(0,1)
		call gstats_print(nout,zaveave,jpmaxstat)
	end if

	if (myproc /= 1) close(nout)

	deallocate(zgmv,zgmvs)
	if (ndimgmv == 5) deallocate(zwind)

	call mpl_barrier()
	call mpl_end()
contains
	subroutine parseargs(niter,nproma,nflevg,cgrid)
		integer,intent(out) :: niter,nproma,nflevg
		character(len=16),intent(inout) :: cgrid

		integer :: iarg
		character(len=128) :: carg

		niter = 0
		nproma = 0
		nflevg = 0
		iarg = 1

		do while (iarg <= command_argument_count())
			call get_command_argument(iarg,carg)

			select case(carg)
			case("-n")
				iarg = iarg+1
				niter = get_int_value(iarg)
			case("--nproma")
				iarg = iarg+1
				nproma = get_int_value(iarg)
			case("-l","--nlev")
				iarg = iarg+1
				nflevg = get_int_value(iarg)
			case("-g","--grid")
				iarg = iarg+1
				call get_command_argument(iarg,cgrid)
			case default
				stop("Unrecognised argument: "//trim(carg))
			end select

			iarg = iarg+1
		end do
	end subroutine

	subroutine parse_grid(cgrid,ndgl,nloen)
		character(len=*) :: cgrid
		integer,intent(out) :: ndgl
		integer,intent(inout),allocatable :: nloen(:)

		integer :: ios,gaussian_number,ndlon,ndgnh
		real,parameter :: rpi=2*asin(1.)

		read(cgrid(2:len_trim(cgrid)),*,iostat=ios) gaussian_number
		if (ios /= 0) stop("ERROR: Unsupported grid: "//trim(cgrid))

		if (cgrid(1:1) == "F") then
			ndgl = 2*gaussian_number
			allocate(nloen(ndgl))

			nloen(:) = 4*gaussian_number
		else if (cgrid(1:1) == "O") then
			ndgl = 2*gaussian_number
			allocate(nloen(ndgl))

			do i=1,ndgl/2
				nloen(i) = 20+4*(i-1)
				nloen(ndgl-i+1) = nloen(i)
			end do
		else if (cgrid(1:1) == "N") then
			ndgl = (gaussian_number+1)/2*2
			ndlon = 2*ndgl
			ndgnh = ndgl/2
			allocate(nloen(ndgl))

			do i=1,ndgnh
				nloen(i) = 10+int(.5+(ndlon-10)*cos(rpi/2*(1-i/real(ndgnh+1))))/2*2
				nloen(ndgl-i+1) = nloen(i)
			end do
		else
			stop("ERROR: Unsupported grid: "//trim(cgrid))
		end if
	end subroutine

	integer function get_int_value(iarg)
		integer,intent(in) :: iarg

		character(len=128) :: carg

		call get_command_argument(iarg,carg)
		read(carg,"(i)") get_int_value
	end function

	subroutine setab(n,na,nb)
		integer(kind=jpim),intent(in) :: n
		integer(kind=jpim),intent(out) :: na,nb

		integer(kind=jpim) :: ja,ib,isqr

		isqr = sqrt(real(n))

		do ja=isqr,n
			ib = n/ja
			if (ja*ib == n) exit
		end do

		if (ja > n) then
			na = 0
			nb = 0
			return
		end if

		na = max(ja,ib)
		nb = min(ja,ib)
	end subroutine

	subroutine setndgl(nlin,nq,nsmax,ndgl)
		integer(kind=jpim),intent(in) :: nlin,nq,nsmax
		integer(kind=jpim),intent(out) :: ndgl

		integer(kind=jpim) :: i
		integer(kind=jpim),parameter :: nsmaxc(*)=(/79,95,127,159,199,255,319,399,511,639,&
			799,1023,1279,1599,1999,3999,7999/)
		integer(kind=jpim),parameter :: nsmaxq(*)=(/21,42,63,106,213,341,426,533,682,853,&
			1364,1706/)
		integer(kind=jpim),parameter :: nsmaxl(*)=(/63,95,127,159,191,199,255,319,399,511,&
			639,799,1023,1279,2047,3999,7999/)

		if (nlin == 1) then
			i = findloc(nsmaxl,nsmax,1)
			if (i == 0) call abor1("UNSUPPORTED SPECTRAL RESOLUTION-LIN. GRID")

			ndgl = nsmax+1
		else if (nq == 1.or.nq == 2) then
			i = findloc(nsmaxc,nsmax,1)
			if (i == 0) call abor1("UNSUPPORTED SPECTRAL RESOLUTION-CUBIC GRID ")

			ndgl = 2*(nsmax+1)
		else
			i = findloc(nsmaxq,nsmax,1)
			if (i == 0) call abor1("UNSUPPORTED SPECTRAL RESOLUTION-QUAD. GRID ")

			if (nsmax == 63) then
				ndgl = 96 ! 64+32
			else if (nsmax == 1364) then
				ndgl = 2048
			else
				! previous cases would give 95 and 2047
				ndgl = 1.5*nsmax+1.4
			end if
		end if
	end subroutine

	subroutine check_ndgl(nlin,ndgl)
		integer(kind=jpim),intent(in) :: nlin,ndgl

		integer(kind=jpim),parameter :: igl0(*)=(/32,64,96,160,320,512,640,800,1024,1280/),&
			igl1(*)=(/32,64,96,128,160,256,320,400,512,640,800,1024/)
		integer(kind=jpim) :: i

		if (nlin == 0) then
			i = findloc(igl0,ndgl,1)
			if (i == 0) call abor1("WRONG SPECTRAL RESOLUTION, QUAD. GRID")
		else if (nlin == 1) then
			i = findloc(igl1,ndgl,1)
			if (i == 0) call abor1("WRONG SPECTRAL RESOLUTION, LIN. GRID")
		else
			call abor1("WRONG NLIN")
		end if
	end subroutine

#ifdef GRIB_API
	subroutine readgrib_resol(insf,nsmax)
		integer(kind=jpim),intent(in) :: insf
		integer(kind=jpim),intent(out) :: nsmax

		integer(kind=jpim) :: iret,igrib(1)
		character(len=127) :: cgridtype

		call grib_new_from_file(insf,igrib(1),iret)
		if (iret /= grib_success) call abor1("ERROR GRIB_NEW_FROM_FILE")

		call grib_get(igrib(1),"gridType",cgridtype,iret)
		if (cgridtype /= "sh") call abor1("INPUT DATA NOT IN SPECTRAL FORM")

		call grib_get(igrib(1),"pentagonalResolutionParameterJ",nsmax)

		call grib_release(igrib(1))
	end subroutine

	subroutine readgrib(insf,nflevg,nspec2g,ivset,ivsetsc,zsp,zvor,zdiv,zt)
		integer(kind=jpim),intent(in) :: insf,nflevg,nspec2g,ivset(:),ivsetsc(:)
		real(jprb),intent(out) :: zsp(:,:),zvor(:,:),zdiv(:,:),zt(:,:,:)

		integer(kind=jpim) :: iflds,icode,jroc,itag,iret,igribout,ioutsf
		integer(kind=jpim) :: iparam(1),igrib(1),iedition(1),icurlev(1)
		real(kind=jprb),allocatable :: zvorg(:,:),zdivg(:,:),zspg(:,:),ztg(:,:,:),zfpdat(:)

		itag = 54321
		igrib(1) = 0
		icode = 0
		ilastlev = 0

		if (myproc == 1) then
			write(nout,*) "allocate full-grid spectral data",nflevg,nspec2g
			allocate(zfpdat(nspec2g))
			allocate(zvorg(nflevg,nspec2g))
			allocate(zdivg(nflevg,nspec2g))
			allocate(ztg(nflevg,nspec2g,1))
			allocate(zspg(1,nspec2g))
		end if

		do
			if (myproc == 1) then
				write(nout,*) "read spectral data from GRIB file"
				call grib_new_from_file(insf,igrib(1),iret)
				if (iret == grib_end_of_file) exit

				if (iflds < nfldx) then
					if (iret /= grib_success) call abor1("ERROR GRIB_NEW_FROM_FILE")

					call grib_get(igrib(1),"edition",iedition(1),iret)
					call grib_get(igrib(1),"paramId",iparam(1),iret)

					if (iparam(1) /= icode.and.icode /= 0) then
						print*,"FIELD ",iparam(1)," NOT TRANSFORMED"
						igribout = 0

						call grib_clone(igrib(1),igribout,iret)
						if (iret /= grib_success) call abor1("ERROR GRIB_CLONE")

						call grib_write(igribout,ioutsf,iret)
						if (iret /= grib_success) call abor1("ERROR GRIB_WRITE")

						call grib_release(igribout)

						cycle
					end if

					call grib_get(igrib(1),"level",icurlev,iret)
					call grib_get(igrib(1),"shortName",cfname,iret)
					call grib_get(igrib(1),"values",zfpdat,iret)

					call grib_release(igrib(1))

					ilev = icurlev(1)
					if (ilev > nflevg) call abor1("NFLEVG < ILASTLEV")
					ilastlev = max(ilev,ilastlev)

					if (cfname == "lnsp") then
						zspg(1,:) = zfpdat(:)
						iflds = iflds+1
					else if (cfname == "vo") then
						zvorg(ilev,:) = zfpdat(:)
						iflds = iflds+1
					else if (cfname == "d") then
						zdivg(ilev,:) = zfpdat(:)
						iflds = iflds+1
					else if (cfname == "t") then
						ztg(ilev,:,1) = zfpdat(:)
						iflds = iflds+1
					end if
				end if

				if (nproc > 1) then
					do jroc=2,nproc
						call mpl_send(ilev,nprcids(jroc),itag)
					end do
				end if
			else
				call mpl_recv(ilev,nprcids(1),itag)
			end if

			if (nproc > 1) then
				if (ilev > 0) call mpl_broadcast(iparam(1),itag,1,cdstring='transform_test:')
				call mpl_barrier()
			end if

			if (ilev == 0) exit
		end do

		if (myproc == 1) call grib_close_file(insf)

		if (nproc > 1) then
			call mpl_broadcast(iflds,itag,1,cdstring='transform_test:')
			call mpl_broadcast(ilastlev,itag,1,cdstring='transform_test:')
		end if

		if (ilastlev < 1) call abor1("ILASTLEV < 1")

		if (ilastlev < nflevg) then
			if (myproc == 1) then
				write(nout,*) "fill sp data for non read levels, from",ilastlev

				do ilev=ilastlev+1,nflevg
					zvorg(ilev,:) = zvorg(mod(ilev-1,ilastlev)+1,:)
					zdivg(ilev,:) = zdivg(mod(ilev-1,ilastlev)+1,:)
					ztg(ilev,:,1) = ztg(mod(ilev-1,ilastlev)+1,:,1)
				end do
			end if
		end if

		write(nout,'("SPECTRAL FIELDS HAVE BEEN SUCCESSFULY READ, IFLDS=",I3)') iflds

		allocate(ito(iflds))
		ito(:) = 1

		call dist_spec(pspecg=zspg,kfdistg=1,kfrom=ito,pspec=zsp,kvset=ivsetsc)
		call dist_spec(zvorg,kfdistg=nflevg,kfrom=ito,pspec=zvor,kvset=ivset)
		call dist_spec(zdivg,kfdistg=nflevg,kfrom=ito,pspec=zdiv,kvset=ivset)
		call dist_spec(ztg(:,:,1),kfdistg=nflevg,kfrom=ito,pspec=zt(:,:,1),kvset=ivset)
		call dist_spec(zspg,kfdistg=1,kfrom=ito,pspec=zsp,kvset=ivsetsc)

		if (myproc == 1) deallocate(ito,zfpdat,zvorg,zdivg,ztg,zspg)
	end subroutine
#endif

	subroutine initspec(nsmax,zsp,sp3d)
		integer,intent(in) :: nsmax
		real(kind=jprb),intent(inout) :: zsp(:,:),sp3d(:,:,:)

		integer(kind=jpim) :: nfield,jf,jk,ioff,im,nump,m_num=4
		integer,allocatable :: myms(:),nasm0(:)
		logical :: linit

		if (m_num > nsmax) call abor1("Error: mnum > nsmax")

		allocate(nasm0(0:nsmax))
		call trans_inq(kasm0=nasm0)

		call trans_inq(knump=nump)

		allocate(myms(nump))
		call trans_inq(kmyms=myms)

		!linit = true.any(myms(:) == m_num)
		linit = .true.
		nfield = size(sp3d,3)

		write(nout,"(/,'SP dims (nsmax/nump):',3(x,i0),l2)") nsmax,nump,nfield,linit
		write(nout,*) "zsp:",shape(zsp)
		write(nout,"('myms:',10(x,i5))") myms(:)
		write(nout,"('nasm0:',10(x,i5))") nasm0(myms(:))

		zsp(1,:) = 0
		if (linit) then
			do jk=1,nump
				im = myms(jk)
				if (im /= 0) cycle

				ioff = nasm0(im)
				if (ioff+2*(nsmax-im) > size(zsp,2)) then
					write(0,*) "--> proc/jk/im/ioff/nspec2:",myproc,jk,im,ioff,size(zsp,2)
					cycle
				end if

				!zsp(1,ioff:ioff+2*(nsmax-im)) = 1/8.
				zsp(1,ioff:ioff+2*(nsmax-im)/5) = 1/16.
			end do

			do jf=1,nfield
				sp3d(:,:,jf) = 0

				do jk=1,nump
					im = myms(jk)
					if (im /= 0) cycle

					ioff = nasm0(im)
					if (ioff+2*(nsmax-im) > size(sp3d,2)) cycle

					! +2*(jf-1): some offset for fields
					sp3d(:,ioff+2*(jf-1):ioff+2*(nsmax-im),jf) = jf/64.
				end do
			end do
		end if

		deallocate(nasm0,myms)
	end subroutine

	subroutine sort(a,n)
		integer(kind=jpim),intent(in) :: n
		real(kind=jprd),intent(inout) :: a(n)

		integer(kind=jpim) :: i,j
		real(kind=jprd) :: x

		do i=2,n
			x = a(i)

			do j=i-1,1,-1
				if (a(j) <= x) exit

				a(j+1) = a(j)
			end do

			a(j+1) = x
		end do
	end subroutine

	subroutine dump_gp_field(iunit,jstep,myproc,nproma,ngpblks,fldchar,fld)
		integer(kind=jpim),intent(in) :: iunit,jstep,myproc,nproma,ngpblks
		character,intent(in) :: fldchar
		real(kind=jprb),intent(in) :: fld(nproma*ngpblks)

		character(len=14) :: filename

		write(filename,"(A1,I3.3,I4.4,'.dat')") fldchar,jstep,myproc

		open(iunit,file=filename,form="unformatted")
		write(iunit) fld(:)
		close(iunit)
	end subroutine

	subroutine spnorm(nl,zsp,ivset,pm,cname)
		integer(jpim),intent(in) :: nl,ivset(:)
		character(len=*),intent(in) :: cname
		real(jprb),intent(in) :: zsp(:,:)
		real(jprb),intent(out) :: pm(nl)

		integer(jpim) :: i

		call specnorm(zsp,ivset,pnorm=pm)

		do i=1,nl
			write(nout,"('SPEC',a3,':',x,i3,x,g22.15)") trim(cname),i,pm(i)
		end do
	end subroutine

	real function sperror(n,znorm,znorm1)
		integer(kind=jpim),intent(in) :: n
		real(kind=jprb),intent(in) :: znorm(n),znorm1(n)

		sperror = maxval(abs(znorm1/znorm-1),1)
	end function

	function mnx(z) result(zmnx)
		real(kind=jprb),intent(in) :: z(:,:)

		real(kind=jprb) :: zmnx(3)

		zmnx(:) = (/minval(z),sum(z)/size(z),maxval(z)/)
	end function

	subroutine stepstat(zts,ztsavg,ztsmin,ztsmax)
		real(kind=jprd),intent(inout) :: zts(:)
		real(kind=jprd),intent(out) :: ztsavg,ztsmin,ztsmax

		ztsavg = sum(zts)
		ztsmin = minval(zts,1)
		ztsmax = maxval(zts,1)

		call mpl_allreduce(zts,"SUM",ldreprod=.false.)
		call mpl_allreduce(ztsavg,"SUM",ldreprod=.false.)
		call mpl_allreduce(ztsmin,"MIN",ldreprod=.false.)
		call mpl_allreduce(ztsmax,"MAX",ldreprod=.false.)
	end subroutine

	subroutine printstack(nproc,nprcids)
		integer(kind=jpim),intent(in) :: nproc,nprcids(:)

		integer(kind=jpim) :: i,istack,getstackusage

		istack = getstackusage()

		if (myproc == 1) then
			print *,"Stack Utilisation Information"
			print*,"task size(bytes)"
			print "(11x,i10)",istack

			do i=2,nproc
				call mpl_recv(istack,nprcids(i),i,cdstring='transform_test:')
				print "('task',x,i3,7x,i10) ",i,istack
			end do
		else
			call mpl_send(istack,nprcids(1),myproc,cdstring='transform_test:')
		end if
	end subroutine

	subroutine gpnorm(cname,zgp,pm,pn,px,psd)
		character(len=*),intent(in) :: cname
		real(jprb),intent(in) :: zgp(:,:,:)
		real(jprb),intent(out) :: pm(nflevg),pn(nflevg),px(nflevg),psd(nflevg)

		integer(jpim) :: i
		real(jprb) :: zdev(nproma,nflevg,ngpblks),zn(nflevg),zx(nflevg)

		call gpnorm_trans(zgp,nflevg,nproma,pm,pn,px,.false.,kresol=1)

		!do i=1,nflevg
		!	zdev(:,i,:) = (zgp(:,i,:)-pm(i))**2
		!end do

		!call gpnorm_trans(zdev,nflevg,nproma,psd,zn,zx,.false.,kresol=1)
		!psd(:) = sqrt(psd)

      psd(:) = 0
		do i=1,nflevg
			write(nout,"('GPUV ',a,':',x,i3,4(x,g22.15))") trim(cname),i,pn(i),pm(i),px(i),&
				psd(i)
		end do
	end subroutine

	subroutine gpnorms2d(zgmvs)
		real(jprb),intent(in) :: zgmvs(:,:,:)

		integer(jpim) :: nf
		real(jprb) :: z3(nproma,3,ngpblks),zm(3),zn(3),zx(3),zsdm(3),zsdn(3),zsdx(3)

		nf = 1
		if (lscders) nf = 3
		call gpnorm_trans(zgmvs,nf,nproma,zm,zn,zx,.false.,kresol=1)
		!do i=1,nf
		!	z3(:,i,:) = (zgmvs(:,i,:)-zm(i))**2
		!end do
		!call gpnorm_trans(z3,nf,nproma,zsdm,zsdn,zsdx,.false.,kresol=1)

      zsdm(:) = 0
		write(nout,"('GMVS P *: 1',4(x,g22.15))") zn(1),zm(1),zx(1),zsdm(1)
		if (lscders) then
			write(nout,"('GMVS Pl: 1',4(x,g22.15))") zn(2),zm(2),zx(2),zsdm(2)
			write(nout,"('GMVS Pm: 1',4(x,g22.15))") zn(3),zm(3),zx(3),zsdm(3)
		end if
	end subroutine

	subroutine gpnorms3d(zgpt,zvd,zuv)
		real(jprb),intent(in) :: zgpt(:,:,:,:),zvd(:,:,:,:),zuv(:,:,:,:)

		integer(jpim) :: i
		real(jprb) :: zm(nflevg),zn(nflevg),zx(nflevg),zsd(nflevg)

		if (lvorgp) call gpnorm("VOR *",zvd(:,:,1,:),zm,zn,zx,zsd)

		if (ldivgp) then
			i = 1
			if (lvorgp) i = 2
			call gpnorm("DIV *",zvd(:,:,i,:),zm,zn,zx,zsd)
		end if

		call gpnorm("U *",zuv(:,:,1,:),zm,zn,zx,zsd)
		call gpnorm("V *",zuv(:,:,2,:),zm,zn,zx,zsd)

		call gpnorm("T *",zgpt(:,:,1,:),zm,zn,zx,zsd)
		if (lscders) then
			call gpnorm("Tl?",zgpt(:,:,2,:),zm,zn,zx,zsd)
			call gpnorm("Tm?",zgpt(:,:,3,:),zm,zn,zx,zsd)
		end if
	end subroutine
end program
