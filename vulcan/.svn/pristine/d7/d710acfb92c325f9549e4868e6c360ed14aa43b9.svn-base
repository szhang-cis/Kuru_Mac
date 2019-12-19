       Program vulcan_flavia
C***********************************************************************
C                                                                      *
C      INTERFACE PROGRAM BETWEEN VULCAN & FLAVIA                       *
C                                                                      *
C***********************************************************************
      implicit none
c------------------------------------------------------------- variables
      character  name*6,
     .           subti*64,title*64,inputfile*300,yesno*1,
     .           outfil1*300,outfil2*300,comment*75,typeprob*2,
     .           argument*75
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
     .           icheckop,nskip,
     .           nargument,iargument,IOS,leng
      real       x,y,z,ttime
      real,allocatable :: coord(:,:), gpcod(:,:,:), proel(:,:),
     .                    props(:,:), disto(:,:), elvar(:,:),
     .                    stren(:,:), fisno(:,:), strep(:,:),
     .                    fpcha(:,:), porot(:,:), displ(:,:),
     .                    react(:,:), distr(:,:), dista(:,:)
      integer,allocatable :: lnods(:,:), matno(:)
      logical    finish, yes_forever,newprocedure,existt,argall
c--------------------------------------------------type of problem------
c
c     sg and sun
c
c     type '(''Enter type of problem(2m,2t,3m,3t): '',$)'
c     accept '(A2)',typeprob
c
c     linux
c
      argall=.false.

      nargument=command_argument_count()
      if(nargument.EQ.0) then
       newprocedure=.false.
      else
       newprocedure=.true.
       iargument=1
      end if

      if (newprocedure) then
       call get_command_argument(iargument,argument)
       if ((argument.eq.'--help').or.(argument.eq.'-?')
     .     .or.(nargument.eq.1).or.(nargument.gt.3)) then
        write(6,*) "[Zero argument '-a' or '--all' for all time]",
     .                                                " optional"
        write(6,*) 'First argument enter type of problem(2m,2t,3m,3t)'
        write(6,*) 'Second argument enter post-process file name ',
     .                                       '( with extension )'
        STOP
       end if
      end if


      if (newprocedure) then
       if ((argument.eq.'-a').or.(argument.eq.'--all')) then
        argall=.true.
        iargument=iargument+1
       end if
       call get_command_argument(iargument,argument)
       typeprob=argument(1:2)
       iargument=iargument+1
      else
       write(6,*) '(''Enter type of problem(2m,2t,3m,3t): '',$)'
       read(5,*) typeprob
      endif
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
c     if(idimension.eq.0) call error('incorrect Type of problem.',0,0)
      if(idimension.eq.0) then
       if(newprocedure) then
        write(*,*) 'incorrect Type of problem [',typeprob,']'
        write(*,*) 'types of problems (2m,2t,3m,3t)'
        STOP
       else
        call error('incorrect Type of problem.')
       endif
      endif
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
      
      if (newprocedure) then
       call get_command_argument(iargument,Inputfile,leng)
       inquire(FILE=Inputfile,EXIST=existt)
       if (.not.existt) then
        write(*,*) "The file ",Inputfile(1:leng)," don't exist"
        STOP
       endif
        open (unit=7, file=Inputfile, status='old', form='unformatted',
     .      IOSTAT=IOS)
        if(IOS.ne.0) then
         write(*,*) "Ocurrio un error al abrir el archivo"
         STOP
        endif
      else
    1  write(6,*)'(''Enter post-process file name(with extension):'',$)'
       read(5,*) inputfile
c
c     sg
c
c     open (unit=7, file=inputfile, status='old', form='unformatted', 
c    .     err=1, readonly)
c
c     linux
c
       open (unit=7, file=inputfile, status='old', form='unformatted',
     .      err=1)
      endif
      write(*,*) "1"
c
c     pointpos = index(Inputfile,'.') ! Find '.' into Inputfile
      pointpos = FindLastPoint(Inputfile)
c
      if(idimension.eq.2) then
       outfiL1  = Inputfile(1:pointpos-1)//'-flavia.dat' ! Append '.neu'
      else
       outfiL1  = Inputfile(1:pointpos-1)//'-flavia.msh' ! Append '.neu'
      end if
      outfiL2  = Inputfile(1:pointpos-1)//'-flavia.res'  ! Append '.neu'
      open(unit=11,file=outfil1,access='sequential',form='formatted')
      open(unit=12,file=outfil2,access='sequential',form='formatted')

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
      write (11,14) title
      write (11,14) inputfile
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
      write (11,14) comment
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
      comment='a'                                       ! avoid question
c
      write (11,14) comment
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
      comment='a'                                       ! avoid question
c
      write (11,14) comment
c
c     sg
c
c     write (11,14) 'c'
   14 format (a75)
c
c     linux
c
      write (11,*) 'c'
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
        if(naccer.ne.0) then
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
     .                  props,coord,
     .                  ndime,nnode,nelem,ngrup,nmats,npoin,npoic,
     .                  nprel,nprop,name,npoin1)
c------------------------------------------Write Gauss Points-----------
c      call writegauss(gpcod),ndime,ngaus,nelem)
c------------------------------------- Prepare to read and write results
      numstp=1
      finish= .false.
      if(newprocedure) then
       yes_forever=argall
      else
       yes_forever=.false.
      endif
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
        write(6,*) '(''Interval='',''Step='',''Iteration='',''Time='')',
     .                itime, istep, iiter, ttime            ! new option
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
         write(6,*)
     .        '(''Must it be loaded? {y)es,n)o,s)top,a)ll,sk)ip }:'',$)'
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
         write(6,*)  '(//,''Please wait...'')'
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
         write(6,*)'(''Skipping steps:'',$)'
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
      end do                                          ! while not finish
    3 close(unit=8,status='delete')
      call exit
c   2 call error('read general variab.',0,0)                  ! old form
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
     .          lnods(nnode,nelem),matno(nelem),npoic,
     .          j,nnodl,nnodx,ielim,ntoble,nelem3,ielem3,ioptmesh
      real      proel(nprel,ngrup),props(nprop,nmats),
     .          coord(ndime,npoin),x,y,z
      character name*6
c-------------------------------------- Write model elements  (Topology)
c
      nnode=int(proel(2,1))  !number of nodes of 1st material's elements
      npoin1=npoin-npoic
      ielim=0
      do ielem=1,nelem
       lgrup=matno(ielem)
       nnodl=int(proel(2,lgrup))
c
       if(nnodl.lt.nnode) then
       do j=nnodl+1,nnode
        lnods(j,ielem)=lnods(nnodl,ielem)    ! ??
       end do
       endif
c
       ltype=int(proel(5,lgrup))
c
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
c  nnodx=1 => hexahedra
c        1 => hexahedra 20-noded
c        3 => tetrahedra
c        3 => tetrahedra 10-noded
c        3 => toblerone  >>  transformed to 3 tetrahedra
c        3 => triangle
c        3 => triangle 6-noded
c        4 => quadrilateral
c        4 => quadrilateral 8-noded
c       11 => line
c
      nnodx=0
      if(nnode.eq. 8.and.ndime.eq.3) nnodx=1     ! hexahedra
      if(nnode.eq.20.and.ndime.eq.3) nnodx=1     ! hexahedra 20-noded
      if(nnode.eq. 4.and.ndime.eq.3) nnodx=3     ! tetrahedra
      if(nnode.eq.10.and.ndime.eq.3) nnodx=3     ! tetrahedra 10-noded
      if(nnode.eq. 6.and.ndime.eq.3) nnodx=3     ! toblerone
      if(nnode.eq. 3.and.ndime.eq.2) nnodx=3     ! triangle
      if(nnode.eq. 6.and.ndime.eq.2) nnodx=3     ! triangle 6-noded
      if(nnode.eq. 4.and.ndime.eq.2) nnodx=4     ! quadrilateral
      if(nnode.eq. 8.and.ndime.eq.2) nnodx=4     ! quadrilateral 8-noded
      if(nnode.eq. 2) nnodx=11                   ! line
c
c**** deals with toblerone
c
      ntoble=0
      if(nnode.eq.6.and.ndime.eq.3) ntoble=1
c
      if(ntoble.eq.0) then
       write(11,*) nelem, npoin1, nnodx
      else
       nelem3=nelem*3
       write(11,*) nelem3, npoin1, nnodx
      endif
c
      write(11,*) 'c'
      do ipoin=1,npoin1
         x=coord(1,ipoin)
         y=coord(2,ipoin)
         if (ndime.eq.3) then
            z=coord(3,ipoin)
            write(11,901) ipoin,x,y,z      ! formatted
         else
            write(11,901) ipoin,x,y        ! formatted
         endif
      end do
  901 format(i8,3e15.6)
      write(11,*) 'c'
c
      if(ntoble.eq.0) then
       do ielem=1,nelem
c       write(11,*) ielem,(lnods(i,ielem),i=1,nnode),matno(ielem)
        write(11,900) ielem,(lnods(i,ielem),i=1,nnode),matno(ielem)
       enddo
  900  format(i8,27i8,i5)
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
      end ! Write geometry
c-----------------------------------------------------Write Gauss Points
      subroutine writegauss(gpcod,ndime,ngaus,nelem)
      implicit none
      integer ndime,ngaus,nelem,idime,igaus,ielem,icounter
      real gpcod(ndime,ngaus,nelem)
      if(ndime.eq.2) then
         write(12,*) 'Gauss Points   ',0,ngaus,0,0,0
         icounter=0
         do ielem=1,nelem
          do igaus=1,ngaus
           icounter=icounter+1
           write(12,*) icounter,(gpcod(idime,igaus,ielem),idime=1,ndime)
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
     .          dista(ndofc,npoin-npoic),
     .          ttime
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
     .     ttime,
     .     gpcod(*),proel(nprel,ngrup),props(nprop,nmats),
     .     stren(nstr1,npoin),
     .     Coord(ndime,npoin),
     .     strep(nnuin,npoin), fpcha(nnupc,npoin), porot(nporo,npoin),
     .     pmean, pdevi(6), pvonm,tempmin,tempmax, displ(ndime,npoin),
     .     react(ndofc,npoin), distr(ndofc,npoin1),dista(ndofc,npoin1),
     .     ai1j2
      character  name*6, yesno*1
c
      save tempmin,tempmax
c
      data nvec /10/
      data vechist /1,2,3,4,5,6,7,8,9,10,
     .              0,0,0,0,0,0,0,0,0,0/
c
c**** displacements or temperatures
c
      if(itype.eq.1) then
       write(12,901) ttime
  901  format('Nodal Displac  ',' 1 ',e15.6,' 2 1 0')
      else
       write(12,900) ttime
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
  900 format('Nodal Temperat ',' 1 ',e15.6,' 1 1 0')
      end if
      do ipoin=1,npoin1
       write(12,103) ipoin,(disto(indofc,ipoin),indofc=1,ndofc)
      end do
c
c**** velocities & accelerations
c
      if(itype.eq.1) then
       if(kdyna.gt.0) then
        write(12,1911) ttime
        do ipoin=1,npoin1
         write(12,103) ipoin,(distr(indofc,ipoin),indofc=1,ndofc)
        end do
        write(12,1912) ttime         ! accelerations can be very large !
        do ipoin=1,npoin1
         write(12,1103) ipoin,(dista(indofc,ipoin),indofc=1,ndofc)
        end do
       end if
      end if
 1911 format('Nodal Velocity ',' 1 ',e15.6,' 2 1 0')
 1912 format('Nodal Accelerat',' 1 ',e15.6,' 2 1 0')
c
c**** temperature rates
c
      if(itype.eq.2) then
       if(kdyna.gt.0) then
        write(12,911) ttime
        do ipoin=1,npoin1
         write(12,103) ipoin,(distr(indofc,ipoin),indofc=1,ndofc)
        end do
       end if
      end if
  911 format('Nodal TempRate ',' 1 ',e15.6,' 1 1 0')
c
c**** phase-change function
c
      if(itype.eq.2) then
       if(nnupc.gt.0) then
        do inupc=1,nnupc
         write(12,1902) inupc,ttime
         do ipoin=1,npoin1
          write(12,103) ipoin,fpcha(inupc,ipoin)
         end do
        enddo
       endif
      endif
 1902 format('Phase-Change ',i2,' 1 ',e15.6,' 1 1 0')
c
c**** porosity criteria
c
      if(itype.eq.2) then
       if(nporo.gt.0.and.ksgau.ne.0) then
        do iporo=1,nporo
         if(iporo.eq.1) write(12,1903) ttime
         if(iporo.eq.2) write(12,1904) ttime
         if(iporo.eq.3) write(12,1905) ttime
         if(iporo.eq.4) write(12,1906) ttime
         if(iporo.eq.5) write(12,1907) ttime
         if(iporo.eq.6) write(12,1908) ttime
         do ipoin=1,npoin1
          write(12,1103) ipoin,porot(iporo,ipoin)
         end do
        enddo
       endif
      endif
 1903 format('Solidific. time',' 1 ',e15.6,' 1 1 0')
 1904 format('T.Grad (f=0)   ',' 1 ',e15.6,' 1 1 0')
 1905 format('T.Grad*Sol.time',' 1 ',e15.6,' 1 1 0')
 1906 format('T.Grad/T.rate  ',' 1 ',e15.6,' 1 1 0')
 1907 format('Solidus veloc. ',' 1 ',e15.6,' 1 1 0')
 1908 format('Gs ts2/3 / Vs  ',' 1 ',e15.6,' 1 1 0')
c
c**** stresses
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,902) ttime
  902   format('Nodal Stresses',' 1 ',e15.6,' 3 1 0')
       endif
      endif
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        do ipoin=1,npoin1
c        write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)  ! old
         if(ndime.eq.2) then                                       ! new
          write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)
         endif
         if(ndime.eq.3) then                  ! ordered according to GiD
          write(12,103) ipoin,stren(1,ipoin),stren(2,ipoin),
     .                        stren(4,ipoin),stren(3,ipoin),
     .                        stren(6,ipoin),stren(5,ipoin)
         endif
        enddo
       endif
      endif
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
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,903) ttime
  903   format('Von M Stresses',' 1 ',e15.6,' 1 1 0')
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
          write(12,103) ipoin,pvonm
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
          write(12,103) ipoin,pvonm
         enddo
        endif
       endif
      endif
c
c**** I1/J2 (Mean stress/Von Mises stress)
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,906) ttime
c 906   format('I1/J2         ',' 1 ',e15.6,' 1 1 0')
  906   format('I1            ',' 1 ',e15.6,' 1 1 0')
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
          write(12,103) ipoin,ai1j2
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
          write(12,103) ipoin,ai1j2
         enddo
        endif
       endif
      endif
c
c**** heat fluxes
c
      if(itype.eq.2) then
       if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
        write(12,2902) ttime
 2902   format('Nodal Heat Flu',' 1 ',e15.6,' 2 1 0')
       end if
      end if
c
      if(itype.eq.2) then
       if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
        do ipoin=1,npoin1
c        write(12,103) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)
c                                                            ! standard
         write(12,104) ipoin,(stren(instr1,ipoin),instr1=1,nstr1)
c                                                            ! see below
        end do
       endif
      endif
c
c**** strains
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        write(12,904) ttime
  904   format('Nodal strains ',' 1 ',e15.6,' 3 1 0')
       end if
      end if
c
      if(itype.eq.1) then
       if(ksgau.ne.0) then
        do ipoin=1,npoin1
c        write(12,103) ipoin,(fisno(instr1,ipoin),instr1=1,nstr1)  ! old
         if(ndime.eq.2) then                                       ! new
          write(12,103) ipoin,(fisno(instr1,ipoin),instr1=1,nstr1)
         endif
         if(ndime.eq.3) then                  ! ordered according to GiD
          write(12,103) ipoin,fisno(1,ipoin),fisno(2,ipoin),
     .                        fisno(4,ipoin),fisno(3,ipoin),
     .                        fisno(6,ipoin),fisno(5,ipoin)
         endif
        enddo
       endif
      endif
c
c**** internal variables (for mechanical or thermal)
c
      if(ksgau.ne.0.and.nnuin.ne.0) then
       do inuin=1,nnuin
        write(12,905) inuin,ttime
        do ipoin=1,npoin1
c        write(12,103) ipoin,strep(inuin,ipoin)              ! standard
         write(12,104) ipoin,strep(inuin,ipoin)              ! see below
        enddo
       end do
      endif
  905 format('Nodal IntVar',i3,' 1 ',e15.6,' 1 1 0')
  104 format(i7,15e15.6)                 ! allows to write large numbers
c
c**** displacements (thermal problem)
c
      if(itype.eq.2) then
       if(large.gt.0) then
        write(12,901) ttime
        do ipoin=1,npoin1
         write(12,103) ipoin,(displ(indofc,ipoin),indofc=1,ndime)
        end do
       endif
      endif

  103 format(i7,15f15.5)
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
