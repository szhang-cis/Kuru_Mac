      Program vulcants_gnuplot
C***********************************************************************
C                                                                      *
C**** INTERFACE PROGRAM BETWEEN VULCAN-TS & GNUPLOT                    *
C                                                                      *
C     Note: change maxlon to increase the dimension                    *
C                                                                      *
C***********************************************************************
      implicit none
c------------------------------------------------------------- variables
c      integer    maxlon,itwo
c      parameter (maxlon=367500000,itwo=1)
c      integer    mem(maxlon),                      ! Memory to read data
c     .           plnods,pmatno,pproel,pprops,pgpcod,pcoord,   ! Pointers
c     .           pdisto,pelvar,pstren,pfisno,pstrep,pfpcha,pporot,
c     .           pdispl,pdistr,
c     .           lastmem
      character  name*6,
     .           subti*64,title*64,inputfile*30,inputfilem*30,
     .           yesno*1,yesnot*1,yesnom*1,yesnop*1,yesnor*1,yesnomic*1,
     .           outfil1*30,outfil2*30,comment*75,
     .           curva*5,curvam*5,curvap*5,curvamp*5,curvar*5,curvamic*5
      integer    i,
     .           ico,icom,icop,icomp,icor,icomic,
     .           nod,nodm,nodm1,nodp,nodp1,nodmp,nodm1p,nodr,nodmic,
     .           numstp,iiter,istep,itime,
     .           ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,ngrup,
     .           nmats,npoin,nprel,nprop,nstr1,ndata,nstat,ksgau,
     .           idata(10),istat(10),ipoin,
     .           pointpos,
     .           flush,cont,
     .           nodos(10),nodosm(10),nodosm1(10),nodosp(10),
     .           nodosp1(10),nodosmp(10),nodosm1p(10),nodosr(10),
     .           nodosmic(10),
     .           nnuin,nnupc,nporo,large,kdyna,
     .           icheckop,ioptionm,kplas,nnum4,nfami,ipointn,ipointr
      real       x,y,z,anfami(20),arfami(20),ttime
      logical    finish, yes_forever
      integer,allocatable:: lnods(:,:),matno(:)
      real,allocatable:: proel(:,:),props(:,:),gpcod(:,:,:),
     .                   coord(:,:),disto(:,:),elvar(:,:)  ,
     .                   stren(:,:),fisno(:,:),strep(:,:)  ,
     .                   fpcha(:,:),porot(:,:),displ(:,:)  ,
     .                   distr(:,:)
     
c-------------------------------------------------------- Open data file
c
c     sg
c
c   1 type '(''Enter post-process file name( with extension ): '',$)'
c     accept '(A30)',inputfile
c
c     linux
c
    1 write(6,*) '(''Enter post-process file name(with extension):'',$)'
      read(5,*) inputfile
c
      pointpos = index(Inputfile,'.')          ! Find '.' into Inputfile
c
c     sg
c
c     open (unit=7, name=inputfile, status='old', form='unformatted', 
c    .      err=1, readonly)
c
c     linux
c
      open (unit=7, file=inputfile, status='old', form='unformatted',
     .      err=1)
c------------------------------------------------checking operations----
c
c     sg
c
c     type '(''Do you want checking operations ? {y)es, n)o}: '',$)'
c     accept '(A1)',yesno
c
c     linux
c
c     write(6,*) '(''Do you want checking operations? {y)es, n)o}:'',$)'
c     read(5,*) yesno
c
      yesno='N'                ! avoid question
c
      icheckop=0
      if(yesno.eq.'Y'.or.yesno.eq.'y') icheckop=1
c------------------------------------------------ Read general variables
      read(7,err=2,end=2)  ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,
     .                     ngrup,nmats,npoin,nprel,nprop,nstr1,title,
     .                     ndata,
     .                     nstat,ksgau,nnuin,nnupc,nporo,large,kdyna
c------------------------------------------------------ Memory dimension
c      plnods=1                               !lnods
      allocate(lnods(nnode,nelem))
c      pmatno=plnods+nelem*nnode              !matno
      allocate(matno(nelem))
c      pproel=pmatno+nelem                    !proel
      allocate(proel(nprel,ngrup))
c      pprops=pproel+nprel*ngrup*itwo         !props
      allocate(props(nprop,nmats))
c      pgpcod=pprops+nprop*nmats*itwo         !gpcod
      allocate(gpcod(ndime,ngaus,nelem))
c      pcoord=pgpcod+ndime*ngaus*nelem*itwo   !coord
      allocate(coord(ndime,npoin))
 
c      lastmem=pcoord+npoin*ndime*itwo
c      if(lastmem.gt.maxlon) call outmem(lastmem)

c      pdisto =pcoord                         !disto
      allocate(disto(ndofc,npoin))
c      pelvar =pdisto +ndofc*npoin*itwo       !elvar
      allocate(elvar(nstat,nelem))
c      pstren =pelvar +nstat*nelem*itwo       !stren
      allocate(stren(nstr1,npoin))
c      pfisno =pstren +nstr1*npoin*itwo       !fisno
      allocate(fisno(nstr1,npoin))
c      pstrep =pfisno +npoin*itwo*nstr1       !strep
      allocate(strep(nnuin,npoin))
c      pfpcha =pstrep +npoin*itwo*nnuin       !fpcha
      allocate(fpcha(nnupc,npoin))
c      pporot =pfpcha +npoin*itwo*nnupc       !porot
      allocate(porot(nporo,npoin))
c      pdispl =pporot +npoin*itwo*nporo       !displ
      if (large.ne.0) then
        allocate(displ(ndime,npoin))
      else
        allocate(displ(1,1))
      end if
     
c      pdistr =pdispl +ndime*npoin*itwo*large !distr
      if (kdyna.ne.0) then
        allocate(distr(ndofc,npoin))
      else
        allocate(distr(1,1))
      end if
c      lastmem=pdistr +ndofc*npoin*itwo*kdyna
c      if(lastmem.gt.maxlon) call outmem(lastmem)
c
      open(unit=8,file='Messages',access='sequential',form='formatted')
c--------------------------------------------------------- Read Geometry
      call readgeom(ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,
     .              nprel,nprop,idata,istat,
     .              coord,gpcod,lnods,
     .              matno,proel,props)
c-------------------------------------recibe numeros de nodos requeridos
      yesno='Y'
      yesnot='Y'
      yesnom='Y'
      yesnop='Y'
      yesnor='Y'
      yesnomic='Y'
      ico=0
      icom=0
      icop=0
      icomp=0
      icor=0
      icomic=0
      curva(1:3)='cur'
      curvam(1:3)='mic'
      curvap(1:3)='pha'
      curvamp(1:3)='por'
      curvar(1:3)='rat'
      curvamic(1:3)='fam'
c
c     sg
c
c     type '(''temperature evolution ? { y)es, n)o }: '',$)'
c     accept '(a1)',yesnot
c
c     linux
c
      write(6,*) '(''temperature evolution ? { y)es, n)o }: '',$)'
      read(5,*) yesnot
c
      if(((yesnot.eq.'Y').OR.(yesnot.eq.'y')).AND.(ico.le.10)) then
       do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(ico.le.10))
        ico=ico+1
c
c     sg
c
c       type '(''enter node number: '',$)'
c       accept '(i3)' ,nod
c
c     linux
c
        write(6,*) '(''enter node number: '',$)'
        read(5,*) nod
c
        nodos(ico)=nod
        if (ico.le.10)then
         write(curva(4:5),'(I1)')ico
        else
         write(curva(4:5),'(I2)')ico
        endif
c
c     sg
c
c       open(unit=(100+ico),name=curva,
c    .       access='sequential',form='formatted')
c       write((100+ico),10) inputfile
c       write((100+ico),11) nodos(ico)
c
c     linux
c
        open(unit=(10+ico),file=curva,
     .       access='sequential',form='formatted')
        write((10+ico),10) inputfile
        write((10+ico),11) nodos(ico)
c
   10   format('#  example:',a30)
   11   format('#  node number:' ,i8)
c
c     sg
c
c       type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c       accept '(a1)',yesno
c
c     linux
c
        write(6,*) '(''do you wish another curve ? { y)es, n)o }:'',$)'
        read(5,*) yesno
c
       enddo
      endif
C
      if(kdyna.gt.0) then
       yesno='Y'
c
c     sg
c
c      type '(''temperature rate evolution ? { y)es, n)o }: '',$)'
c      accept '(a1)',yesnor
c
c     linux
c
       write(6,*) '(''temperature rate evolution ? { y)es, n)o }: '',$)'
       read(5,*) yesnor
c
       if(((yesnor.eq.'Y').OR.(yesnor.eq.'y')).AND.(icor.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(icor.le.10))
         icor=icor+1
c
c     sg
c
c        type '(''enter node number: '',$)'
c        accept '(i3)',nodr
c
c     linux
c
         write(6,*) '(''enter node number: '',$)'
         read(5,*) nodr
c
         nodosr(icor)=nodr
         if(icor.le.10)then
          write(curvar(4:5),'(I1)')icor
         else
          write(curvar(4:5),'(I2)')icor
         endif
c
c     sg
c
c        open(unit=(150+icor),name=curvar,
c    .        access='sequential',form='formatted')
c        write((150+icor),10) inputfile
c        write((150+icor),11) nodosr(icor)
c
c     linux
c
         open(unit=(20+icor),file=curvar,
     .        access='sequential',form='formatted')
         write((20+icor),10) inputfile
         write((20+icor),11) nodosr(icor)
c
c     sg
c
c        type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*) '(''do you wish another curve ? { y)es, n)o }:'',$)'
         read(5,*) yesno
c
        enddo
       endif
      endif                        ! kdyna.gt.0
C
      if(ksgau.ne.0.and.nnuin.ne.0) then
       yesno='Y'
c
c     sg
c
c      type'(''microscopical variables evolution ? {y)es, n)o}:'',$)'
c      accept '(a1)',yesnom
c
c     linux
c
       write(6,*)
     .          '(''microscopical variables evolution?{y)es, n)o}:'',$)'
       read(5,*) yesnom
c
       if(((yesnom.eq.'Y').OR.(yesnom.eq.'y')).AND.(ico.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(ico.le.10))
         icom=icom+1
c
c     sg
c
c        type '(''enter node number: '',$)'
c        accept '(i3)' ,nodm
c
c     linux
c
         write(6,*) '(''enter node number: '',$)'
         read(5,*) nodm
c
         nodosm(icom)=nodm
c
c     sg
c
c        type '(''enter microscopical variable number: '',$)'
c        accept '(i3)' ,nodm1
c
c     linux
c
         write(6,*) '(''enter microscopical variable number: '',$)'
         read(5,*) nodm1
c
         nodosm1(icom)=nodm1
         if (icom.le.10)then
          write(curvam(4:5),'(I1)')icom
         else
          write(curvam(4:5),'(I2)')icom
         endif
c
c     sg
c
c        open(unit=(200+icom),name=curvam,
c    .        access='sequential',form='formatted')
c        write((200+icom),20) inputfile
c        write((200+icom),21) nodosm(icom)
c        write((200+icom),22) nodosm1(icom)
c
c     linux
c
         open(unit=(30+icom),file=curvam,
     .        access='sequential',form='formatted')
         write((30+icom),20) inputfile
         write((30+icom),21) nodosm(icom)
         write((30+icom),22) nodosm1(icom)
c
   20    format('#  example:',a30)
   21    format('#  node number:' ,i4)
   22    format('#  microscopical variable Nro:' ,i4)
c
c     sg
c
c        type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*) '(''do you wish another curve? { y)es, n)o }: '',$)'
         read(5,*) yesno
c
        enddo
       endif
C
C     SGCI: Spheroidal Graphite Cast Iron
c
       yesno='Y'
c
c     sg
c
c      type'(''nuclei computation (only for SGCI)? {y)es, n)o}:'',$)'
c      accept '(a1)',yesnomic
c
c     linux
c
       write(6,*)
     .       '(''nuclei computation (only for SGCI)? {y)es, n)o}:'',$)'
       read(5,*) yesnomic
c
       if(((yesnomic.eq.'Y').OR.(yesnomic.eq.'y')).AND.
     .                                              (icomic.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(icomic.le.10))
         icomic=icomic+1
c
c     sg
c
c        type '(''enter node number: '',$)'
c        accept '(i3)' ,nodmic
c
c     linux
c
         write(6,*) '(''enter node number: '',$)'
         read(5,*) nodmic
c
         nodosmic(icomic)=nodmic
         if (icomic.le.10)then
          write(curvamic(4:5),'(I1)')icomic
         else
          write(curvamic(4:5),'(I2)')icomic
         endif
c
         nodosmic(icomic)=nodmic
c
c     sg
c
c   4    type '(''Enter file with SGCI model parameters: '',$)'
c        accept '(A30)',inputfilem
c
c     linux
c
    4    write(6,*) '(''Enter file with SGCI model parameters:'',$)'
         read(5,*) inputfilem
c
c     sg
c
c        open (unit=9, name=inputfilem, status='old', form='formatted',
c    .         err=4, readonly)
c
c     linux
c
         open (unit=9, file=inputfilem, status='old', form='formatted',
     .         err=4)
c
         read(9,*) kplas,nfami
         if(nfami.gt.20) then
c
c     sg
c
c         type '(''Number of families gt 20 - Stop'',$)'
c
c     linux
c
          write(6,*) '(''Number of families gt 20 - Stop'',$)'
c
         endif
c
         do i=1,nfami+1
          read(9,*) arfami(i)
         enddo
c
         if(kplas.eq.4) then                ! Boeri multinodular model
          nnum4=(nnuin-4)/2
          ipointn=4
          ipointr=4+nnum4
         endif
         if(kplas.eq.5) then                ! Su uninodular model
          nnum4=(nnuin-4)/3
          ipointn=4
          ipointr=4+nnum4
         endif
         if(kplas.eq.9) then                ! Dardati multinodular model
          nnum4=(nnuin-10)/5
          ipointn=10                        ! to be revised !!!!
          ipointr=10+nnum4                  ! to be revised !!!!
         endif
c
c                   Esta parte estoy seguro que funciona bien para los
c                   modelos de Boeri y Su.
c                   No estoy seguro que funcione bien para el modelo
c                   Dardati ya que no me acuerdo en este modelo Dardati 
c                   cual es la diferencia entre RNU y RNUZ1. Por lo tanto, 
c                   no me queda claro cual de estos radios hay que tener 
c                   en cuenta para calcular el numero de nucleos (cada 
c                   familia queda determinada por el intervalo de radios
c                   (minimo y maximo por cada familia)
c                   definidos en el archivo adicional que hay que entregarle
c                   cuando se corre este programa intergplos).
c                   En resumen, no tengo claro cual es el valor de los
c                   punteros ipointn e ipointr en este caso.
c                   Tampoco tengo claro cual de las variables del modelo
c                   Carazo hay que calcular.
c
         if(kplas.eq.10) then               ! Carazo eutectoid model
          nnum4=(nnuin-2)/3
          ipointn=2
          ipointr=2+nnum4
         endif

         if(kplas.eq.910) then              ! Dardati + Carazo models
          nnum4=(nnuin-12)/8
          ipointn=2                         ! a definir !!!!
          ipointr=2+nnum4                   ! a definir !!!!
         endif


c
c     sg
c
c        open(unit=(350+icomic),name=curvamic,
c    .        access='sequential',form='formatted')
c        write((350+icomic),10) inputfilem
c        write((350+icomic),21) nodosmic(icomic)
c        write((350+icomic),22) nfami,arfami(1),arfami(nfami+1)
c
c     linux
c
         open(unit=(60+icomic),file=curvamic,
     .        access='sequential',form='formatted')
         write((60+icomic),10) inputfilem
         write((60+icomic),21) nodosmic(icomic)
         write((60+icomic),12) nfami,arfami(1),arfami(nfami+1)
c
   12    format('#  number of nuclei families:',i4,' - min/max radii:',
     .                                                      e10.3,e10.3)
c
c     sg
c
c        type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*) '(''do you wish another curve ? { y)es, n)o }:'',$)'
         read(5,*) yesno
c
        enddo
       endif
      endif       ! ksgau.ne.0.and.nnuin.ne.0
C
      if(nnupc.gt.0) then
       yesno='Y'
c
c     sg
c
c      type '(''macroscopical f_pc evolution ? { y)es, n)o }: '',$)'
c      accept '(a1)',yesnop
c
c     linux
c
       write(6,*)'(''macroscopical f_pc evolution? { y)es, n)o }: '',$)'
       read(5,*) yesnop
c
       if(((yesnop.eq.'Y').OR.(yesnop.eq.'y')).AND.(ico.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(ico.le.10))
         icop=icop+1
c
c     sg
c
c        type '(''enter node number: '',$)'
c        accept '(i3)' ,nodp
c
c     linux
c
         write(6,*) '(''enter node number: '',$)'
         read(5,*) nodp
c
         nodosp(icop)=nodp
c
c     sg
c
c        type'(''enter macroscopical phase-change function number:'',$)'
c        accept '(i3)' ,nodp1
c
c     linux
c
         write(6,*)
     .          '(''enter macroscopical phase-change funct.number:'',$)'
         read(5,*) nodp1
c
         nodosp1(icop)=nodp1
         if (icop.le.10)then
          write(curvap(4:5),'(I1)')icop
         else
          write(curvap(4:5),'(I2)')icop
         endif
c
c     sg
c
c        open(unit=(250+icop),name=curvap,
c    .        access='sequential',form='formatted')
c        write((250+icop),30) inputfile
c        write((250+icop),31) nodosp(icop)
c        write((250+icop),32) nodosp1(icop)
c
c     linux
c
         open(unit=(40+icop),file=curvap,
     .        access='sequential',form='formatted')
         write((40+icop),30) inputfile
         write((40+icop),31) nodosp(icop)
         write((40+icop),32) nodosp1(icop)
c
   30    format('#  example:',a30)
   31    format('#  node number:' ,i4)
   32    format('#  macroscopical phase-change function Nro:' ,i4)
c
c     sg
c
c        type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*) '(''do you wish another curve? { y)es, n)o }: '',$)'
         read(5,*) yesno
c
        enddo
       endif
      endif       ! nnupc.gt.0
C
      if(ksgau.ne.0.and.nporo.ne.0) then
       yesno='Y'
c
c     sg
c
c      type'(''porosity criteria evolution ? {y)es, n)o}:'',$)'
c      accept '(a1)',yesno
c
c     linux
c
       write(6,*)'(''porosity criteria evolution ? {y)es, n)o}:'',$)'
       read(5,*) yesno
c
       if(((yesnom.eq.'Y').OR.(yesnom.eq.'y')).AND.(ico.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(ico.le.10))
         icomp=icomp+1
c
c     sg
c
c        type '(''enter node number: '',$)'
c        accept '(i3)' ,nodmp
c
c     linux
c
         write(6,*) '(''enter node number: '',$)'
         read(5,*) nodmp
c
         nodosmp(icomp)=nodmp
c
c     sg
c
c        type '(''enter porosity criterion number: '',$)'
c        accept '(i3)' ,nodm1p
c
c     linux
c
         write(6,*) '(''enter porosity criterion number: '',$)'
         read(5,*) nodm1p
c
         nodosm1p(icomp)=nodm1p
         if (icomp.le.10)then
          write(curvamp(4:5),'(I1)')icomp
         else
          write(curvamp(4:5),'(I2)')icomp
         endif
c
c     sg
c
c        open(unit=(300+icomp),name=curvamp,
c    .        access='sequential',form='formatted')
c        write((300+icomp),40) inputfile
c        write((300+icomp),41) nodosmp(icomp)
c        write((300+icomp),42) nodosm1p(icomp)
c
c     linux
c
         open(unit=(50+icomp),file=curvamp,
     .        access='sequential',form='formatted')
         write((50+icomp),40) inputfile
         write((50+icomp),41) nodosmp(icomp)
         write((50+icomp),42) nodosm1p(icomp)
c
   40    format('#  example:',a30)
   41    format('#  node number:' ,i4)
   42    format('#  porosity criterion Nro:' ,i4)
c
c     sg
c
c        type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c        accept '(a1)',yesno
c
c     linux
c
         write(6,*)'(''do you wish another curve ? { y)es, n)o }: '',$)'
         read(5,*) yesno
c
        enddo
       endif
      endif      ! ksgau.ne.0.and.nporo.ne.0
c------------------------------------- Prepare to read and write results
      numstp=1
      finish= .false.
      do while(.not.finish) ! ------------------- Read and Write results
       call readres(disto,elvar,fisno,iiter,istep,
     .              itime,ndime,ndofc,nelem,npoin,nstat,nstr1,subti,
     .              title,ttime,stren,ksgau,finish,nnuin,
     .              strep,fpcha,porot,displ,
     .              distr,
     .              nporo,nnupc,icheckop,large,istat,kdyna)
       if(.not.finish) then                               ! Deal with it
        call wrtres(disto,elvar,fisno,gpcod,
     .              lnods,matno,proel,props,
     .              stren,coord,strep,fpcha,
     .              porot,distr,
     .              ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,nmats,
     .              nnode,npoin,nprel,nprop,nstat,nstr1,iiter,
     .              istat,istep,itime,ttime,ksgau,numstp,name,
     .              yesno,
     .              ico,icom,icop,icomp,icor,icomic,
     .              nodos,nodosm,nodosm1,nodosp,nodosp1,nodosmp,
     .              nodosm1p,nodosr,nodosmic,
     .              nnupc,nnuin,nporo,
     .              nnum4,nfami,ipointn,ipointr,anfami,arfami)
        numstp=numstp+1
       end if                                         ! not finish
      end do                                          ! while not finish
    3 close(unit=8,status='delete')
      call exit
c   2 call error('read general variab.',0,0)
    2 call error('read general variab.')
      end
c
c----------------------- S U B R O U T I N E S -------------------------
c
c--------------------------------------------------------- Read Geometry
      subroutine readgeom(ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,
     .                    nprel,nprop,idata,istat,coord,gpcod,lnods,
     .                    matno,proel,props)
      implicit none
      integer   ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,
     .          nprel,nprop,ielem,idime,igaus,
     .          lnods(nnode,nelem),      matno(nelem),
     .          idata(8),                istat(4)
      real      coord(ndime,npoin),      gpcod(ndime,ngaus,nelem),
     .          proel(nprel,ngrup),      props(nprop,nmats)
c
c     do ielem=1,nelem                                        ! old form
c      read(7,err=1,end=1) ((gpcod(idime,igaus,ielem),
c    .                       idime=1,ndime),igaus=1,ngaus)
c     end do
c
      read(7,err=1,end=1) coord,lnods,matno,proel,props,idata,istat
      igaus=1
      return
c   1 call error('read geometry  read ',0,0)
    1 call error('read geometry  read ')
      end ! readgeom
c---------------------------------------------------------- Read results
      subroutine readres(disto,elvar,fisno,iiter,istep,itime,
     .                   ndime,ndofc,nelem,npoin,nstat,nstr1,
     .                   subti,title,ttime,stren,ksgau,finish,nnuin,
     .                   strep,fpcha,porot,displ,
     .                   distr,
     .                   nporo,nnupc,icheckop,large,istat,kdyna)
      implicit none
      integer   iiter,istep,itime,ndime,ndofc,nelem,npoin,
     .          nstat,nstr1,ksgau,ielem,jstat,nnuin,nnupc,icheckop,
     .          indofc,ipoin,nporo,large,istat(4),kdyna
      real      disto(ndofc,npoin),       elvar(nstat,nelem),
     .          stren(nstr1,npoin),       fisno(ndime,npoin),
     .          strep(nnuin,npoin),       fpcha(nnupc,npoin),
     .          porot(nporo,npoin),       displ(ndime,npoin),
     .          distr(ndofc,npoin),
     .          ttime, tempmin, tempmax
c
c**** stren: nodal heat fluxes (read but not write!!!)
c     fisno: nodal temp. gradients (not implemented yet !)
c     fpcha: nodal (macroscopic) phase-change function
c     strep: nodal (microscopic) internal variables
c     porot: nodal porosity criteria
c     displ: nodal displacements (thermal problem)
c     distr: nodal temperature rate
c
      character subti*64,title*64
      logical finish

      save tempmin,tempmax

    1 read(7,err=1,end=2) title,subti,itime,istep,iiter,ttime,disto
c
c     do ielem=1,nelem                                        ! old form
c      read(7) (elvar(jstat,ielem),jstat=1,nstat)
c     end do
c
      if(icheckop.eq.1) then
       if(ttime.eq.0.0) then
        tempmin= 10.0E+10
        tempmax=-10.0E+10
        do indofc=1,ndofc
         do ipoin=1,npoin
          if(disto(indofc,ipoin).lt.tempmin)
     .     tempmin=disto(indofc,ipoin)
          if(disto(indofc,ipoin).gt.tempmax)
     .     tempmax=disto(indofc,ipoin)
         enddo
        enddo
       else
        do indofc=1,ndofc
         do ipoin=1,npoin
          if(disto(indofc,ipoin).lt.tempmin)
     .     disto(indofc,ipoin)=tempmin
          if(disto(indofc,ipoin).gt.tempmax)
     .     disto(indofc,ipoin)=tempmax
         enddo
        enddo
       endif
      endif
c
c**** temperature rates (only transient problems)
c
      if(kdyna.gt.0) then
       read(7) distr
      endif
c
c**** displacements (only coupled problems)
c
      if(large.gt.0) then
       read(7) displ
      endif
c
c**** phase-change function
c
      if(nnupc.gt.0) then
       read(7) fpcha
      endif
c
c**** internal variables
c
      if(ksgau.ne.0.and.(istat(3).ne.istat(4))) then
       read(7) stren
      end if
c
      if(ksgau.ne.0.and.nnuin.ne.0) then
       read(7) strep
      end if
c
c**** porosity criteria
c
      if(ksgau.ne.0.and.nporo.ne.0) then
       read(7) porot
      end if
c
      return
    2 finish= .true.
c
c     sg
c
c     type '(//,''End of data. Please wait...'')'
c
c     linux
c
      write(6,*) '(//,''End of data. Please wait...'')'
c
      end                                                 ! read results
c--------------------------------------------------------- Write results
      subroutine wrtres(disto,elvar,fisno,gpcod,
     .                  lnods,matno,proel,props,
     .                  stren,Coord,strep,fpcha,
     .                  porot,distr,
     .                  ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,nmats,
     .                  nnode,npoin,nprel,nprop,nstat,nstr1,iiter,
     .                  istat,istep,itime,ttime,ksgau,numstp,name,
     .                  yesno,
     .                  ico,icom,icop,icomp,icor,icomic,
     .                  nodos,nodosm,nodosm1,nodosp,nodosp1,nodosmp,
     .                  nodosm1p,nodosr,nodosmic,
     .                  nnupc,nnuin,nporo,
     .                  nnum4,nfami,ipointn,ipointr,anfami,arfami)
      implicit none
      integer numstp,lnods(*),matno(*),
     .        i,ii,j,jj,k,
     .        ico,icom,icop,icomp,icor,icomic,
     .        nodos(*),nodosm(*),nodosm1(*),nodosp(*),nodosp1(*),
     .        nodosmp(*),nodosm1p(*),nodosr(*),nodosmic(*),
     .        ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,nmats,
     .        nnode,npoin,nprel,nprop,nstat,nstr1,iiter,istat(4),
     .        istep,itime,ksgau,ncrit,nnupc,nnuin,nporo,
     .        nnum4,nfami,ipointn,ipointr
      real    disto(*),elvar(*),fisno(*),ttime,
     .        gpcod(*),proel(*),props(*),stren(*),Coord(ndime,*),
     .        strep(nnuin,npoin),fpcha(nnupc,npoin),porot(nporo,npoin),
     .        distr(*),anfami(*),arfami(*)
      character name*6, yesno*1
c-----------------------------------------------------------------------
      if(ico.gt.0) then
       do i=1,ico
c
c     sg
c
c       write((100+i),*) ttime, disto(nodos(i))
c
c     linux
c
        write((10+i),*) ttime, disto(nodos(i))
c
       end do
      endif
c-----------------------------------------------------------------------
      if(icor.gt.0) then
       do i=1,icor
c
c     sg
c
c       write((150+i),*) ttime, distr(nodosr(i))
c
c     linux
c
        write((20+i),*) ttime, distr(nodosr(i))
c
       end do
      endif
c-----------------------------------------------------------------------
      if(icom.gt.0) then
       do i=1,icom
c
c     sg
c
c       write((200+i),*) ttime, strep(nodosm1(i),nodosm(i))
c
c     linux
c
        write((30+i),*) ttime, strep(nodosm1(i),nodosm(i))
c
       end do
      endif
c-----------------------------------------------------------------------
      if(icop.gt.0) then
       do i=1,icop
c
c     sg
c
c       write((250+i),*) ttime, fpcha(nodosp1(i),nodosp(i))
c
c     linux
c
        write((40+i),*) ttime, fpcha(nodosp1(i),nodosp(i))
c
       end do
      endif
c-----------------------------------------------------------------------
      if(icomp.gt.0) then
       do i=1,icomp
c
c     sg
c
c       write((300+i),*) ttime, porot(nodosm1p(i),nodosmp(i))
c
c     linux
c
        write((50+i),*) ttime, porot(nodosm1p(i),nodosmp(i))
c
       end do
      endif
c-----------------------------------------------------------------------
      if(icomic.gt.0) then
       do i=1,icomic
c
        ii=nodosmic(i)
        do k=1,nfami
         anfami(k)=0.0d0
        enddo
c
        do j=ipointr+1,ipointr+nnum4
         jj=ipointn+j-ipointr
         do k=1,nfami
          if((strep(j,ii).gt.arfami(k)).and.
     .       (strep(j,ii).le.arfami(k+1)))
     .                                  anfami(k)=anfami(k)+strep(jj,ii)
         enddo
        enddo
c
c     sg
c
c       write((350+i),*) ttime,(anfami(k),k=1,nfami)
c
c     linux
c
c       write((60+i),*) ttime,(anfami(k),k=1,nfami)
        write((60+i),10) ttime,(anfami(k),k=1,nfami)
c
       end do
      endif
c
   10 format(f10.3,1x,20e10.3)
c
      end                                                       ! wrtres
c--------------------------------------------------------- Out of memory
      subroutine outmem(lastmem)
      implicit none
      integer lastmem
    1 format('Insuficient memory.',/,'Dimension array mem to more than',
     . i10)
      open(unit=8,file='Messages',access='sequential',form='formatted')
      write(8,1) lastmem
c     call error('Memory dimension    ',0,0)
      call error('Memory dimension    ')
      end                                                ! out of memory
c----------------------------------------------------------------- Error
      subroutine error(message)
      implicit none
      character message*20
 1    format(' Error in: ',a20,', Code:',i5)
      write(8,1) message
c
c     sg
c
c     type '(//,''Intgplos: ** Error was found. LOOK MESSAGE FILE **'')'
c
c     linux
c
      write(6,*)'(//,''Intgplos:**Error was found.LOOK MESSAGE FILE*'')'
c
      call exit
      end                                                        ! error
c------------------------------------------------------------- uppercase
c
c     SG or Linux Red Hat 4.0
c
c     subroutine uppercase(name)
c     implicit none
c     integer ipos,iocv
c     character name*6
c     do ipos=1,6
c      iocv=ichar(name(ipos:ipos))
c      if('141'O.le.iocv.and.'172'O.ge.iocv) then
c       iocv=iocv-'40'O
c       name(ipos:ipos)=char(iocv)
c      end if
c     end do
c     end                                                    ! uppercase
c============================== E N D ==================================
