       Program vulcan_paraview
C***********************************************************************
C                                                                      *
C      INTERFACE PROGRAM BETWEEN VULCAN & PARAVIEW                     *
C                                                                      *
C***********************************************************************
      implicit none
c------------------------------------------------------------- variables
      integer    maxlon,itwo
      parameter (maxlon=81835000,itwo=1)
      character  name*6,
     .           subti*64,title*64,inputfile*300,yesno*1,
     .           outfil*300,outfilres*300,comment*75,typeprob*2
      integer    i,
     .           numstp,iiter,istep,itime,
     .           kprob,ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,ngrup,
     .           nmats,npoin,nprel,nprop,nstr1,ndata,nstat,ksgau,
     .           idata(11),istat(4),ipoin,npoin1,
     .           pointpos,
     .           flush(10),idimension,itype,
     .           FindLastPoint,nnuin,npoic,nnupc,nnuxx,nporo,npoxx,
     .           large,largex,
     .           kdyna,naccer,
     .           icheckop,nskip
      real*8       x,y,z
      real       ttime
      real,allocatable :: coord(:,:), gpcod(:,:,:), proel(:,:),
     .                    props(:,:), disto(:,:), elvar(:,:),
     .                    stren(:,:), fisno(:,:), strep(:,:),
     .                    fpcha(:,:), porot(:,:), displ(:,:),
     .                    react(:,:), distr(:,:), dista(:,:)
      integer,allocatable :: lnods(:,:), matno(:)
      logical    finish, yes_forever
c--------------------------------------------------type of problem------
c
c     sg and sun
c
c     type '(''Enter type of problem(2m,2t,3m,3t): '',$)'
c     accept '(A2)',typeprob
c
c     linux
c
      i = IARGC()
      if(i.ne.2) then
      write(*,*)'Error: Bad number of arguments'
      write(*,*)'Usage interpara problem_type pos_filename'
      stop
      end if
      call getarg(1, typeprob)
      call getarg(2, Inputfile)
c      write(6,*) 'Enter type of problem(2m,2t,3m,3t):' !!
c      read(5,*) typeprob  !!
c
      idimension=0
      if(typeprob.eq.'2m') then
       idimension=2
       itype=1
      end if
      if(typeprob.eq.'2t') then
       idimension=2
       itype=2
      end if
      if(typeprob.eq.'3m') then
       idimension=3
       itype=1
      end if
      if(typeprob.eq.'3t') then
       idimension=3
       itype=2
      end if
      if(idimension.eq.0) call error('incorrect Type of problem.')
c------------------------------------------------checking operations----
      if(itype.eq.2) then
c
c     sg
c
c      type '(''Do you want checking operations ? {y)es, n)o}: '',$)'
c      accept '(A1)',yesno
c
c     linux
c
c      write(6,*)'(''Do you want checking operations? {y)es, n)o}:'',$)'
c      read(5,*) yesno
c
       yesno='n'                                        ! avoid question
c
       icheckop=0
       if(yesno.eq.'Y'.or.yesno.eq.'y') icheckop=1
      endif
c-------------------------------------------------------- Open data file
c
c     sg
c
c   1 type '(''Enter post-process file name( with extension ): '',$)'
c     accept '(A300)',inputfile
c
c     linux
c
    1  continue
c    1 write(6,*)'Enter post-process file name (with extension):' !!
c      read(5,*) Inputfile !!
c
c     pointpos = index(Inputfile,'.') ! Find '.' into Inputfile
      pointpos = FindLastPoint(Inputfile)
c
c
      outfil = Inputfile(1:pointpos-1)//'-pv.vtk'  ! Append '.neu'
      open(unit=11,file=outfil,access='sequential',form='formatted')
      
c      outfilres = Inputfile(1:pointpos-1)//'-pv_000.vtt'  ! Append '.neu'
c      open(unit=12,file=outfilres,access='sequential',form='formatted')
c     sg
c
c     open (unit=7, file=inputfile, status='old', form='unformatted', 
c    .     err=1, readonly)
c
c     linux
c     unit 7 file . pos
       open(unit=7,file=Inputfile,status='old',form='unformatted',err=1)
c
c------------------------------------------------ Read general variables
      read(7,err=2,end=2)  ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,
     .                     ngrup,nmats,npoin,nprel,nprop,nstr1,title,
     .                     ndata,
     .                     nstat,ksgau,nnuin,nnuxx,npoxx,large,kdyna
c
c     nnuxx=npoic for mechanical problem
c     nnuxx=nnupc for thermal problem
c
c     npoxx=0     for mechanical problem
c     npoxx=nporo for thermal problem
c
c     large not used for mechanical problem
c     large is a "deformed shape indicator" for thermal problem
c
      if(itype.eq.1) then   ! mechanical
       npoic=nnuxx
       nnupc=0
       nporo=0
       largex=0
       naccer=1                             ! velocities & accelerations
      endif
      if(itype.eq.2) then   ! thermal
       npoic=0
       nnupc=nnuxx
       nporo=npoxx
       largex=0
       if(large.gt.0) largex=1
       naccer=0                             ! temperature rates
      endif
c
c     sg
c
c     type '(//,''Problem title: '',a64,/)', title
c
c     linux
c
c     write(6,*)'(//,''Problem title: '',/)', title      ! avoid writing
c
      write (11,'(A26)')'# vtk DataFile Version 2.0'

c     write (11,14) title  ! Ricardo: Nuevo formato no acepta titulo ni nombre de archivo
c     write (11,14) inputfile
c
c     sg
c
c     type '(''Enter first comment ): '',$)'
c     accept '(A75)',comment
c
c     linux
c
c     write(6,*) '(''Enter first comment ): '',$)'
c     read(5,*) comment
      comment='a'                                       ! avoid question
c
      write (11,'(A1)') comment
c
c     sg
c
c     type '(''Enter second comment ): '',$)'
c     accept '(A75)',comment
c
c     linux
c
c     write(6,*) '(''Enter second comment ): '',$)'
c     read(5,*) comment
c
      write (11,'(A5)')"ASCII"
c
c     sg
c
c     type '(''Enter third comment ): '',$)'
c     accept '(A75)',comment
c
c     linux
c
c     write(6,*) '(''Enter third comment ): '',$)'
c     read(5,*) comment
c
c     comment='a'                                       ! avoid question
c
c     write (11,14) comment
c
c     sg
c
c     write (11,14) 'c'
   14 format (a75)
c
c     linux
c
c     write (11,*) 'c'
c------------------------------------------------------ Memory dimension
      allocate(lnods(nnode,nelem))
      allocate(matno(nelem))
      allocate(proel(nprel,ngrup))
      allocate(props(nprop,nmats))
      allocate(coord(ndime,npoin))
      allocate(gpcod(ndime,ngaus,nelem))
      allocate(disto(ndofc,npoin))
      allocate(elvar(nstat,nelem))
      allocate(stren(nstr1,npoin))
      allocate(fisno(nstr1,npoin))
      allocate(strep(nnuin,npoin))
      allocate(fpcha(nnupc,npoin))
      allocate(porot(nporo,npoin))
      if (largex.ne.0) then
        allocate(displ(ndime,npoin))
      else
        allocate(displ(1,1))
      end if
      allocate(react(ndofc,npoin))
      if (kdyna.ne.0) then
        allocate(distr(ndofc,npoin))
        if (naccer.ne.0) then
          allocate(dista(ndofc,npoin))
        else
          allocate(dista(1,1))
        end if
      else
        allocate(distr(1,1))
      end if
c
      open(unit=8,file='Messages',access='sequential',form='formatted')
c--------------------------------------------------------- Read Geometry
      call readgeom(ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,npoic,
     .              nprel,nprop,idata,istat,
     .              coord,gpcod,lnods,
     .              matno,proel,props)
c-------------------------------------------------------- Write Geometry
      call writegeom(lnods,matno,proel,
     .               props,coord,
     .               ndime,nnode,nelem,ngrup,nmats,npoin,npoic,
     .               nprel,nprop,name,npoin1)
c------------------------------------------Write Gauss Points-----------
c      call writegauss(gpcod,ndime,ngaus,nelem)
c------------------------------------- Prepare to read and write results
      numstp=0
      finish= .false.
      yes_forever=.false.
      nskip=0

      do while(.not.finish) ! ------------------- Read and Write results
        

        call readres(disto,elvar,fisno,iiter,istep,
     .              itime,ndime,ndofc,nelem,npoin,nstat,nstr1,subti,
     .              title,ttime,stren,ksgau,finish,nnuin,
     .              strep,fpcha,porot,
     .              displ,react,distr,dista,
     .              npoic,nnupc,nporo,large,itype,istat,kdyna)
       if(.not.finish) then                               ! Deal with it
c
c     sg
c
c       type '(/,a64)', subti
c       type '(''Interval='',i5,'', Step='',i5,'', Iteration='',i5,
c    .         '', Time='',f12.5)', itime, istep, iiter, ttime
c
c     linux
c
c       write(6,*) subti
c       write(6,*) '(''Interval='','', Step='','', Iteration='',
c    . '', Time='')', itime, istep, iiter, ttime
        write(*,*) ttime       
        write(6,345) itime, istep, iiter, ttime
345   format('Interval:', I5, 'Step:',I5, 'Iteration:',I5,'Time:',e15.6)
c 
        if (yes_forever) then
         yesno='y'
c        i=flush(6)
        else
c
c     sg
c
c        type '(''Must it be loaded? {y)es,n)o,s)top,a)ll,sk)ip }:'',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*)'Must it be loaded? {y)es,n)o,s)top,a)ll,sk)ip }:'
         read(5,*) yesno
c
        end if
        if(yesno.eq.'S'.or.yesno.eq.'s') then
         finish=.true.
c
c     sg
c
c        type '(//,''Please wait...'')'
c
c     linux
c
         write(6,*)'Please wait...'
c
c        i=flush(6)
        end if
        if(yesno.eq.'A'.or.yesno.eq.'a') then
         yes_forever=.true.
         yesno='y'
        end if
        if(yesno.eq.'K'.or.yesno.eq.'k') then
c
c     sg
c
c        type '(''Skipping steps:'',$)'
c        accept nskip
c
c     linux
c
         write(6,*)'Skipping steps:'
         read(5,*) nskip
c
         yes_forever=.true.
        endif
        if(nskip.ne.0) then
         if(numstp.eq.1) then
          yesno='y'
         else
          yesno='n'
          if(mod((numstp-1),nskip).eq.0) yesno='y'
         end if
        end if
        if(yesno.eq.'Y'.or.yesno.eq.'y') then ! ---------- Write results

       write(outfilres,'(A,I3.3,A)')Inputfile(1:pointpos-1)//'-pv_',
     .       numstp,'.vtt'
       if(.not.finish) then
       open(unit=12,file=outfilres,access='sequential',form='formatted')
       endif
 600  format('POINT_DATA ', i7)
      write(12, 600) npoin1
         call wrtres(disto,elvar,fisno,gpcod,
     .               lnods,matno,proel,props,
     .               stren,coord,strep,fpcha,
     .               porot,displ,react,distr,
     .               dista,
     .               kprob,ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,
     .               nmats,nnode,npoin,nprel,nprop,nstat,nstr1,iiter,
     .               istat,istep,itime,ttime,ksgau,numstp,name,
     .               yesno,npoin1,itype,nnuin,npoic,nnupc,nporo,
     .               large,kdyna,
     .               icheckop)
        end if
       end if                                         ! not finish
       numstp=numstp+1
       close(unit=12)
c       outfilres  = Inputfile(1:pointpos-1)//'-pv_'//numstp//'.vtt'
      end do                                          ! while not finish
    3 close(unit=8,status='delete')
      call exit
    2 call error('read general variab.')
      end 
c
c------------------- S U B R O U T I N E S -----------------------------
c
c--------------------------------------------------------- Read Geometry
      subroutine readgeom(ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,
     .                    npoic,
     .                    nprel,nprop,idata,istat,coord,gpcod,lnods,
     .                    matno,proel,props)
      implicit none
      integer   ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,npoic,
     .          nprel,nprop,ielem,idime,igaus,
     .          lnods(nnode,nelem),      matno(nelem),
     .          idata(11),               istat(4)
      real      coord(ndime,npoin-npoic), gpcod(ndime,ngaus,nelem),
     .          proel(nprel,ngrup),      props(nprop,nmats)
c
c     do ielem=1,nelem                     ! Coordinates at Gauss points
c      read(7,err=1,end=1) ((gpcod(idime,igaus,ielem),
c    .                       idime=1,ndime),igaus=1,ngaus)
c     end do
c
      read(7,err=1,end=1) coord,lnods,matno,proel,props,idata,istat
      return
    1 call error('read geometry  read ')
      end                                                     ! readgeom
c-------------------------------------------------------- Write Geometry
      subroutine writegeom(lnods,matno,proel,props,coord,ndime,
     .                     nnode,nelem,ngrup,nmats,npoin,npoic,
     .                     nprel,nprop,name,npoin1)
      implicit none
      integer   ndime,nnode,nelem,ngrup,nmats,npoin,nprel,nprop,
     .          numstp,ipoin,isystm,ielem,
     .          lgrup,nnodp,imats,ltype,ntype,i,npoin1,
     .          lnods(nnode,nelem),matno(nelem),npoic,tmppara(nnode),
     .          j,nnodl,nnodx,ielim,ntoble,nelem3,ielem3,ioptmesh
      real      proel(nprel,ngrup),props(nprop,nmats),
     .          coord(ndime,npoin),x,y,z
      character name*6
c-------------------------------------- Write model elements  (Topology)
c
      nnode=int(proel(2,1))  !number of nodes of 1st material's elements
      npoin1=npoin-npoic
      ielim=0
      
c**** Header Description
      write (11,'(A25)')'DATASET UNSTRUCTURED_GRID'
      write (11, 500),npoin1
 500  format("POINTS ",i8," float")

c**** End Header Description
c


      do ielem=1,nelem
       lgrup=matno(ielem)
       nnodl=int(proel(2,lgrup))
c
       if(nnodl.lt.nnode) then
       do j=nnodl+1,nnode
        lnods(j,ielem)=lnods(nnodl,ielem)    ! ?? R: WTF?!
       end do
       endif
c
       ltype=int(proel(5,lgrup))

c**** delirio para cambiar la conectividad de los elementos de contacto
c
c      if(ltype.eq.104) then
c       i=lnods(3,ielem)
c       lnods(3,ielem)=lnods(4,ielem)
c       lnods(4,ielem)=i
c      endif
c
c**** delirio para que los elem 101 pasen a tener 4 nodos (2D) u 8 (3D)
c
c      if(ltype.eq.101) then                               ! old version
c       ipoin=npoin1+1
c       coord(1,ipoin)=coord(1,(lnods(1,ielem)))
c       coord(2,ipoin)=coord(2,(lnods(1,ielem)))
c       if (ndime.eq.3) coord(3,ipoin)=coord(3,(lnods(1,ielem)))
c       lnods(3,ielem)=ipoin
c       ipoin=ipoin+1
c       coord(1,ipoin)=coord(1,(lnods(2,ielem)))
c       coord(2,ipoin)=coord(2,(lnods(2,ielem)))
c       if (ndime.eq.3) coord(3,ipoin)=coord(3,(lnods(2,ielem)))
c       lnods(4,ielem)=ipoin
c       npoin1=ipoin
c      endif
c
       if(ltype.eq.101.or.ltype.eq.32) then                ! new version
        do j=1,nnode/2
         lnods(j+nnode/2,ielem)=lnods(j,ielem)
        enddo
       endif
c
c**** mas delirio: Elemento 3 o sea el gap mecanico
c
c      if(ltype.eq.4) then
c       do j=3,nnode
c        lnods(j,ielem)=lnods(2,ielem)
c       end do
c      endif
c
c**** algoritmo para eliminar los elementos de contacto para lograr que
c     los contornos salgan bien en FLAVIA (ioptmesh=0)
c
c     estos elementos de contacto no deben eliminarse en mallas con
c     numeracion elemental optimizada para ser usada con el "frontal
c     solver" (ioptmesh=1)
c
c      ioptmesh=0                                    ! better as input !
       ioptmesh=1                                    ! better as input !
       if(ioptmesh.eq.0) then
        ltype=int(proel(5,lgrup))
        if(ltype.eq.101.or.ltype.eq.104.or.
     .     ltype.eq.4.or.ltype.eq.32) then
         ielim=ielim+1
        endif
       endif
      end do                ! ielem=1,nelem
c
      nelem=nelem-ielim     ! only solid elements
c
c----------------------------------------------------------- Write model
c
c  nnodx=12 => hexahedra
c        25 => hexahedra 20-noded
c        10 => tetrahedra
c        24 => tetrahedra 10-noded
c        13 => toblerone  >>  transformed to 3 tetrahedra
c        5  => triangle
c        22 => triangle 6-noded
c        9  => quadrilateral
c        23 => quadrilateral 8-noded
c        3  => line
c
      nnodx=0
      if(nnode.eq. 8.and.ndime.eq.3) nnodx=12     ! hexahedra
      if(nnode.eq.20.and.ndime.eq.3) nnodx=25     ! hexahedra 20-noded
      if(nnode.eq. 4.and.ndime.eq.3) nnodx=10     ! tetrahedra
      if(nnode.eq.10.and.ndime.eq.3) nnodx=24     ! tetrahedra 10-noded
      if(nnode.eq. 6.and.ndime.eq.3) nnodx=13     ! toblerone
      if(nnode.eq. 3.and.ndime.eq.2) nnodx=5      ! triangle
      if(nnode.eq. 6.and.ndime.eq.2) nnodx=22     ! triangle 6-noded
      if(nnode.eq. 4.and.ndime.eq.2) nnodx=9      ! quadrilateral
      if(nnode.eq. 8.and.ndime.eq.2) nnodx=23     ! quadrilateral 8-noded
      if(nnode.eq. 2) nnodx=3                     ! line
c
c**** deals with toblerone
c
      ntoble=0
      if(nnode.eq.6.and.ndime.eq.3) ntoble=1
c
      if(ntoble.eq.0) then
c      write(11,*) nelem, npoin1, nnodx
      else
       nelem3=nelem*3
c      write(11,*) nelem3, npoin1, nnodx
      endif
c
c     write(11,*) 'c'
      do ipoin=1,npoin1
         x=coord(1,ipoin)
         y=coord(2,ipoin)
         if (ndime.eq.3) then
            z=coord(3,ipoin)
            write(11,901) x,y,z      ! formatted
         else
            write(11,901) x,y,0.0        ! formatted
         endif
      end do
  901 format(3e15.6)
c     write(11,*) 'c'
c**** Cells Header
      write(11, 501) nelem, (nnode+1)*nelem
  501 format('CELLS ',i12,2x,i12)
c**** End Cells Header
c
      if(ntoble.eq.0) then
       do ielem=1,nelem
c       write(11,*) ielem,(lnods(i,ielem),i=1,nnode),matno(ielem)
c       write(11,900) ielem,(lnods(i,ielem),i=1,nnode),matno(ielem) !FALLA CON TOBLERONE!!
        do i=1, nnode
          tmppara(i) = lnods(i,ielem) - 1
        enddo
        write(11,900) nnode, (tmppara(j),j=1,nnode)
       enddo
  900  format(i1,27i8)
      else
       do ielem=1,nelem
        do j=1,3
         ielem3=(ielem-1)*3+j
         if(j.eq.1) write(11,900) ielem3,lnods(1,ielem),lnods(2,ielem),
     .                                   lnods(3,ielem),lnods(4,ielem),
     .                                   matno(ielem)
         if(j.eq.2) write(11,900) ielem3,lnods(3,ielem),lnods(4,ielem),
     .                                   lnods(6,ielem),lnods(2,ielem),
     .                                   matno(ielem)
         if(j.eq.3) write(11,900) ielem3,lnods(4,ielem),lnods(2,ielem),
     .                                   lnods(5,ielem),lnods(6,ielem),
     .                                   matno(ielem)
        enddo
       enddo
      endif           ! ntoble.eq.0
      write(11,'(A11,i8)')'CELL_TYPES ', nelem
      do ielem=1,nelem
        write(11,*) nnodx
      enddo
      end ! Write geometry
c-----------------------------------------------------Write Gauss Points
      subroutine writegauss(gpcod,ndime,ngaus,nelem)
      implicit none
      integer ndime,ngaus,nelem,idime,igaus,ielem,icounter
      real gpcod(ndime,ngaus,nelem)
      if(ndime.eq.2) then
c        write(12,*) 'Gauss Points   ',0,ngaus,0,0,0
         icounter=0
         do ielem=1,nelem
          do igaus=1,ngaus
           icounter=icounter+1
c          write(12,*) icounter,(gpcod(idime,igaus,ielem),idime=1,ndime)
          end do
         end do
      end if
      end
c---------------------------------------------------------- Read results
      subroutine readres(disto,elvar,fisno,iiter,istep,itime,
     .                   ndime,ndofc,nelem,npoin,nstat,nstr1,
     .                   subti,title,ttime,stren,ksgau,finish,nnuin,
     .                   strep,fpcha,porot,
     .                   displ,react,distr,dista,
     .                   npoic,nnupc,nporo,large,itype,istat,kdyna)
      implicit none
      integer   iiter,istep,itime,ndime,ndofc,nelem,npoin,
     .          nstat,nstr1,ksgau,ielem,istat1,ipoin,idime,indofc,nnuin,
     .          npoic,
     .          nnupc,nporo,large,itype,istat(4),kdyna
      real      disto(ndofc,npoin-npoic), elvar(nstat,nelem),
     .          stren(nstr1,npoin-npoic), fisno(nstr1,npoin-npoic),
     .          strep(nnuin,npoin-npoic), fpcha(nnupc,npoin),
     .          porot(nporo,npoin),       displ(ndime,npoin),
     .          react(ndofc,npoin),       distr(ndofc,npoin-npoic),
     .          dista(ndofc,npoin-npoic)
      real ttime
c
c**** stren: nodal stresses & nodal fluxes
c     fisno: nodal strains
c     fpcha: nodal (macroscopic) phase-change function (thermal prob.)
c     strep: nodal internal variables
c     porot: nodal porosity criteria (thermal problem)
c     displ: nodal displacements (thermal problem)
c     react: nodal reactions (mechanical problem)
c     distr: nodal velocities or temperature rates
c                                        (mechanical or thermal problem)
c     dista: nodal accelerations (mechanical problem)
c
      character subti*64,title*64
      logical finish
c
c**** displacements or temperatures (mechanical or thermal)
c
    1 read(7,err=1,end=2) title,subti,itime,istep,iiter,ttime,disto
c
c
c     ttime=ttime+2820   ! useful when restart is used (better as input)
c
c
c**** velocities & accelerations (only mechanical for dynamic problems)
c
      if(itype.eq.1) then
       if(kdyna.gt.0) then
        read(7) distr
        read(7) dista
       endif
      endif
c
c**** temperature rates (only thermal for transient problems; itype=2)
c
      if(itype.eq.2) then
       if(kdyna.gt.0) then
        read(7) distr
       endif
      endif
c
c**** displacements (only thermal for large strain problems; itype=2)
c
      if(itype.eq.2) then
       if(large.gt.0) then
        read(7) displ
       endif
      endif
c
c**** reactions (only mechanical; itype=1)
c
      if(itype.eq.1) then
       read(7) react
      endif
c
c**** phase-change function (only thermal; itype=2)
c
      if(nnupc.gt.0) then                   ! nnupc > 0 only for itype=2
       read(7) fpcha
      endif
c
c**** internal variables (mechanical or thermal)
c
      if(itype.eq.1) then
       if(ksgau.ne.0.and.nnuin.ne.0) then
        read(7) stren,fisno,strep
       end if
      endif
c
      if(itype.eq.2) then
       if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
        read(7) stren
       end if
       if(ksgau.ne.0.and.nnuin.ne.0) then
        read(7) strep
       end if
      endif
c
c**** porosity criteria (only thermal; itype=2)
c
      if(ksgau.ne.0.and.nporo.ne.0) then    ! nporo > 0 only for itype=2
       read(7) porot
      end if
c
      return
 2    finish= .true.
c
c     sg
c
c     type '(//,''End of data. Please wait...'')'
c
c     linux
c
      write(6,*) '(//,''End of data. Please wait...'')'
c
      end ! read results
c--------------------------------------------------------- Write results
      subroutine wrtres(disto,elvar,fisno,gpcod,lnods,matno,proel,props,
     .                  stren,Coord,strep,fpcha,porot,displ,react,distr,
     .                  dista,
     .                  kprob,ndime,ndofc,ndofn,nelem,
     .                  ngaus,ngrup,nhist,nmats,nnode,npoin,nprel,
     .                  nprop,nstat,nstr1,iiter,
     .                  istat,istep,itime,ttime,ksgau,numstp,name,
     .                  yesno,npoin1,itype,nnuin,npoic,nnupc,
     .                  nporo,large,kdyna,
     .                  icheckop)
      implicit none
      integer nelem,nnode
      integer numstp,lnods(nnode,nelem),matno(nelem),
     .        kprob,ndime,ndofc,ndofn,ngaus,ngrup,nhist,nmats,
     .        npoin,nprel,nprop,nstat,nstr1,iiter,istat(4),
     .        istep,itime,ksgau,ncrit,npoin1,itype,
     .        ipoin,ielem,indofc,inode,inhist,igaus,instr1,icounter,
     .        idime,nnupc,nporo,inupc,iporo,
     .        vechist(20),ivec,nvec,nnuin,inuin,npoic,icheckop,large,
     .        kdyna
      real disto(ndofc,npoin1),elvar(nstat,nelem),fisno(nstr1,npoin),
     .     gpcod(*),proel(nprel,ngrup),props(nprop,nmats),
     .     stren(nstr1,npoin),
     .     Coord(ndime,npoin),
     .     strep(nnuin,npoin), fpcha(nnupc,npoin), porot(nporo,npoin),
     .     pmean, pdevi(6), pvonm,tempmin,tempmax, displ(ndime,npoin),
     .     react(ndofc,npoin), distr(ndofc,npoin1),dista(ndofc,npoin1),
     .     ai1j2
      real ttime
      character  name*6, yesno*1
c
      save tempmin,tempmax
c
      data nvec /10/
      data vechist /1,2,3,4,5,6,7,8,9,10,
     .              0,0,0,0,0,0,0,0,0,0/
c
c**** displacements or temperatures
c**** Displacements: Vectors
c**** Temperatures:  Scalars
c      write(12, 600) npoin1
c  600 format('POINT_DATA ', i7)
  601 format(3e15.6)
      if(itype.eq.1) then
        write(12, '(A33)')'VECTORS Nodal_Displacements float'
      else
        write(12, '(A34)')'SCALARS Nodal_Temperatures float 1'
        write(12, '(A20)')'LOOKUP_TABLE default'
c
       if(icheckop.eq.1) then
        if(ttime.eq.0.0) then
         tempmin= 10.0E+10
         tempmax=-10.0E+10
         do indofc=1,ndofc
          do ipoin=1,npoin1
           if(disto(indofc,ipoin).lt.tempmin)
     .      tempmin=disto(indofc,ipoin)
           if(disto(indofc,ipoin).gt.tempmax)
     .      tempmax=disto(indofc,ipoin)
          enddo
         enddo
        else
         do indofc=1,ndofc
          do ipoin=1,npoin1
           if(disto(indofc,ipoin).lt.tempmin)
     .      disto(indofc,ipoin)=tempmin
           if(disto(indofc,ipoin).gt.tempmax)
     .      disto(indofc,ipoin)=tempmax
          enddo
         enddo
        endif
       endif
c
      end if
      if(ndime.lt.3.and.itype.eq.1) then
        do ipoin=1,npoin1
          write(12,601) (disto(indofc,ipoin),indofc=1,ndofc),0.0
        end do
      else 
        do ipoin=1,npoin1
          write(12,601) (disto(indofc,ipoin),indofc=1,ndofc)
        end do
      end if

c
c**** velocities & accelerations
c**** Both: Vectors
      if(itype.eq.1) then
       if(kdyna.gt.0) then
         write(12, '(A28)')'VECTORS Nodal_Velocity float'
        if(ndime.lt.3) then
          do ipoin=1,npoin1
            write(12,601) (disto(indofc,ipoin),indofc=1,ndofc),0.0
          end do
        else
          do ipoin=1,npoin1
            write(12,601) (disto(indofc,ipoin),indofc=1,ndofc)
          end do
        end if
        write(12, '(A32)')'VECTORS Nodal_Acceleration float'
        if(ndime.lt.3) then
          do ipoin=1,npoin1
            write(12,601) (disto(indofc,ipoin),indofc=1,ndofc),0.0
          end do
        else
          do ipoin=1,npoin1
            write(12,601) (disto(indofc,ipoin),indofc=1,ndofc)
          end do
        end if
       end if
      end if
c
c**** temperature rates
c**** Scalars
      if(itype.eq.2) then
       if(kdyna.gt.0) then
        write(12, '(A33)')'SCALARS Temperature_Rates float 1'
        write(12, '(A20)')'LOOKUP_TABLE default'
        do ipoin=1,npoin1
         write(12,601) (distr(indofc,ipoin),indofc=1,ndofc)
        end do
       end if
      end if
c
c**** phase-change function UNIMPLEMENTED
c**** Scalars
      if(itype.eq.2) then
       if(nnupc.gt.0) then
        do inupc=1,nnupc
         write(12,1902) inupc
         do ipoin=1,npoin1
           write(12,103) fpcha(inupc,ipoin)
         end do
        enddo
       endif
      endif
 1902 format('SCALARS Phase-Change',I2.2,' float 1',/,
     . 'LOOKUP_TABLE default')
c
c**** porosity criteria 
c**** Scalars
      if(itype.eq.2) then
       if(nporo.gt.0.and.ksgau.ne.0) then
        do iporo=1,nporo
         if(iporo.eq.1) write(12,1903)
         if(iporo.eq.2) write(12,1904)
         if(iporo.eq.3) write(12,1905)
         if(iporo.eq.4) write(12,1906)
         if(iporo.eq.5) write(12,1907)
         if(iporo.eq.6) write(12,1908)
         do ipoin=1,npoin1
           write(12,103) porot(iporo,ipoin)
         end do
        enddo
       endif
      endif
 1903 format('SCALARS Solidific_Time float 1',/,'LOOKUP_TABLE default')
 1904 format('SCALARS T.Grad_(f=0) float 1',/,'LOOKUP_TABLE default')
 1905 format('SCALARS T.Grad*Sol.time float 1',/,'LOOKUP_TABLE default')
 1906 format('SCALARS T.Grad/T.rate float 1',/,'LOOKUP_TABLE default')
 1907 format('SCALARS Solidus_veloc float 1',/,'LOOKUP_TABLE default')
 1908 format('SCALARS Gs_ts2/3_/_Vs float 1',/,'LOOKUP_TABLE default')
c
c**** stresses
c**** Tensor
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,'(A26)')'TENSORS Nodal_Stress float'
       endif
      endif
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        do ipoin=1,npoin1
c        write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)  ! old
         if(ndime.eq.2) then                                       ! new
          write(12,103) stren(1,ipoin),stren(3,ipoin),
     .                  0.0,stren(3,ipoin),
     .                  stren(2,ipoin),0.0,0.0,0.0,0.0
         endif
         if(ndime.eq.3) then                  ! ordered according to GiD
          write(12,103) stren(1,ipoin),stren(3,ipoin),
     .                  stren(5,ipoin),stren(3,ipoin),
     .                  stren(2,ipoin),stren(6,ipoin),
     .                  stren(5,ipoin),stren(6,ipoin),
     .                  stren(4,ipoin)
         endif
        enddo
       endif
      endif
 1909 format(3e15.6,2x,3e15.6,2x,3e15.6)
c
c**** Hoop stress (FLAVIA 2D does not plot it, but GiD does!)
c
c     if(itype.eq.1) then
c      if(ksgau.ne.0) then
c       if(ndime.eq.2) then
c        write(12,912) ttime
c       endif
c      endif
c     endif
c 912 format('Nodal Hoop St.',' 1 ',e15.6,' 1 1 0')
c
c     if(itype.eq.1) then
c      if(ksgau.ne.0) then
c       if(ndime.eq.2) then
c        do ipoin=1,npoin1
c         write(12,103) ipoin,stren(4,ipoin)
c        enddo
c       endif
c      endif
c     endif
c
c**** Von Mises stress
c**** Scalar
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,'(A32)')'SCALARS Von_Mises_Stress float 1'
        write(12,'(A20)')'LOOKUP_TABLE default'
       endif
      endif
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        if(ndime.eq.2) then
         do ipoin=1,npoin1
          pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
          pdevi(1)=stren(1,ipoin)-pmean
          pdevi(2)=stren(2,ipoin)-pmean
          pdevi(3)=stren(3,ipoin)
          pdevi(4)=stren(4,ipoin)-pmean
          pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .      2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)
          pvonm=sqrt(3.0/2.0*pvonm)
           write(12,601) pvonm
         enddo
        endif
c
        if(ndime.eq.3) then
         do ipoin=1,npoin1
          pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
          pdevi(1)=stren(1,ipoin)-pmean
          pdevi(2)=stren(2,ipoin)-pmean
          pdevi(3)=stren(3,ipoin)
          pdevi(4)=stren(4,ipoin)-pmean
          pdevi(5)=stren(5,ipoin)
          pdevi(6)=stren(6,ipoin)
          pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .      2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)+
     .      2.0*pdevi(5)*pdevi(5)+2.0*pdevi(6)*pdevi(6)
          pvonm=sqrt(3.0/2.0*pvonm)
          write(12,601) pvonm
         enddo
        endif
       endif
      endif
c
c**** I1/J2 (Mean stress/Von Mises stress)
c**** Scalar
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,'(A33)')'SCALARS I1/J2_Mean_Stress float 1'
        write(12,'(A20)')'LOOKUP_TABLE default'
       end if
      end if
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        if(ndime.eq.2) then
         do ipoin=1,npoin1
          pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
          pdevi(1)=stren(1,ipoin)-pmean
          pdevi(2)=stren(2,ipoin)-pmean
          pdevi(3)=stren(3,ipoin)
          pdevi(4)=stren(4,ipoin)-pmean
          pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .      2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)
          pvonm=sqrt(3.0/2.0*pvonm)
c         ai1j2=0.0
c         if(pvonm.gt.1.0e-08) ai1j2=pmean/pvonm
          ai1j2=pmean
          write(12,601) ai1j2
         enddo
        endif
c
        if(ndime.eq.3) then
         do ipoin=1,npoin1
          pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
          pdevi(1)=stren(1,ipoin)-pmean
          pdevi(2)=stren(2,ipoin)-pmean
          pdevi(3)=stren(3,ipoin)
          pdevi(4)=stren(4,ipoin)-pmean
          pdevi(5)=stren(5,ipoin)
          pdevi(6)=stren(6,ipoin)
          pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .      2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)+
     .      2.0*pdevi(5)*pdevi(5)+2.0*pdevi(6)*pdevi(6)
          pvonm=sqrt(3.0/2.0*pvonm)
c         ai1j2=0.0
c         if(pvonm.gt.1.0e-08) ai1j2=pmean/pvonm
          ai1j2=pmean
          write(12,601) ai1j2
         enddo
        endif
       endif
      endif
c
c**** heat fluxes
c**** Vector
      if(itype.eq.2) then
       if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
        write(12,'(A29)')'VECTORS Nodal_Heat_Flux float'
       end if
      end if
c
      if(itype.eq.2) then
       if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
         if(ndime.lt.3) then 
           do ipoin=1,npoin1
c            write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)
c                                                            ! standard
             write(12,601) (stren(instr1,ipoin),instr1=1,nstr1),0.0
c          
           end do                                          ! see below
         else
           do ipoin=1,npoin1
c            write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)
c                                                            ! standard
             write(12,601) (stren(instr1,ipoin),instr1=1,nstr1)
c                                                            ! see below
           end do
         end if
       endif
      endif
c
c**** strains
c**** Tensor
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,'(A26)')'TENSORS Nodal_Strain float'
       end if
      end if
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        do ipoin=1,npoin1
c        write(12,103) ipoin,(fisno(instr1,ipoin),instr1=1,nstr1)  ! old
         if(ndime.eq.2) then                                       ! new
           write(12,103) fisno(1,ipoin),fisno(3,ipoin),
     .                  0.0,fisno(3,ipoin),
     .                  fisno(2,ipoin),0.0,0.0,0.0,0.0
         endif
         if(ndime.eq.3) then                  ! ordered according to GiD
            write(12,103) fisno(1,ipoin),fisno(3,ipoin),
     .                  fisno(5,ipoin),fisno(3,ipoin),
     .                  fisno(2,ipoin),fisno(6,ipoin),
     .                  fisno(5,ipoin),fisno(6,ipoin),
     .                  fisno(4,ipoin)
         
         endif
        enddo
       endif
      endif
c
c**** internal variables (for mechanical or thermal)
c**** Scalar
      if(ksgau.ne.0.and.nnuin.ne.0) then
       do inuin=1,nnuin
         write(12,'(A18,i3.3,A8)')'SCALARS Nodal_Var_',inuin,' float 1'
         write(12,'(A21)')'LOOKUP_TABLE default'
        do ipoin=1,npoin1
c        write(12,103) ipoin,strep(inuin,ipoin)              ! standard
         write(12,104) strep(inuin,ipoin)              ! see below commented by Ricardo
        enddo
       end do
      endif
  104 format(15e15.6)                 ! allows to write large numbers
c
c**** displacements (thermal problem)
c**** Vector
      if(itype.eq.2) then
       if(large.gt.0) then
         write(12,'(A20)')'VECTORS Nodal_Displacements float'
         if(ndime.lt.3) then
           do ipoin=1,npoin1
             write(12,103) (displ(indofc,ipoin),indofc=1,ndime)
           end do
         else
           do ipoin=1,npoin1
             write(12,103) (displ(indofc,ipoin),indofc=1,ndime)
           end do
         end if
       endif
      endif

  103 format(15f15.5)
 1103 format(i7,15f25.5)
      end 
c--------------------------------------------------------- Out of memory
      subroutine outmem(lastmem)
      implicit none
      integer lastmem
 1    format('Insuficient memory.',/,'Dimension array mem to more than',
     . i10)
      open(unit=8,file='Messages',access='sequential',form='formatted')
      write(8,1) lastmem
      call error('Memory dimension    ')
      end ! out of memory
c----------------------------------------------------------------- Error
      subroutine error(message)
      implicit none
      character message*20
 1    format(' Error in: ',a20,', Code:',i5)
      write(8,1) message
c
c     sg
c
c     type '(//,''Interfla: ** Error was found. LOOK MESSAGE FILE **'')'
c
c     linux
c
      write(6,*)'(//,''Interfla:**Error was found.LOOK MESSAGE FILE*'')'
c
      call exit
      end ! error
c------------------------------------------------------------- uppercase
c
c     sg or linux (Red Hat 4.0)
c
c     subroutine uppercase(name)
c     implicit none
c     integer ipos,iocv
c     character name*6
c
c     do ipos=1,6
c      iocv=ichar(name(ipos:ipos))
c      if('141'O.le.iocv.and.'172'O.ge.iocv) then
c       iocv=iocv-'40'O
c       name(ipos:ipos)=char(iocv)
c      end if
c     end do
c     end                                                    ! uppercase
c--------------------------------------------------------- FindLastPoint
      function FindLastPoint(name)
      integer FindLastPoint
      character name(300)*1
      do i=300,1,-1
       if(name(i).eq.'.') then
        FindLastPoint=i
        return
       end if
      end do
      end
c========================= E N D =======================================
