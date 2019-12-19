      Program vulcanm_gnuplot
C***********************************************************************
C                                                                      *
C**** INTERFACE PROGRAM BETWEEN VULCAN-M & GNUPLOT                     *
C                                                                      *
C     Note: change maxlon to increase the dimension                    *
C                                                                      *
C***********************************************************************
      implicit none
c------------------------------------------------------------- variables
c      integer    maxlon,itwo
c      parameter (maxlon=36000000,itwo=1)
c      integer    mem(maxlon),                      ! Memory to read data
c     .           plnods,pmatno,pproel,pprops,pgpcod,pcoord,   ! Pointers
c     .           pdisto,pelvar,pstren,pfisno,pstrep,preact,
c     .           preac1,preac2,pdistr,pdista,lastmem
      character  opcode*1,ctype*1
      character  name*6,text*20,
     .           subti*64,title*64,inputfile*30,inputfiler*30,
     .           yesno*1,yesnot*1,yesnop*1,yesnopo*1,yesnor*1,
     .           outfil1*30,gstra*1,gstre*1,comment*75,
     .           curva*5,curvap*5,curvamp*5,curvar*5,curvav*5,
     .           curvare*5,curvaret*5
      integer    i,key,atype,
     .           ico,icop,icomp,icor,
     .           nod,nodp,nodp1,nodmp,nodm1p,
     .           numstp,iiter,istep,itime,
     .           kprob,ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,ngrup,
     .           nmats,npoin,nprel,nprop,nstr1,ndata,nstat,ksgau,
     .           idata(10),istat(10),ipoin,npoin1,
     .           pointpos,
     .           flush(10),cont,ictype,nnuin,npoic,nporo,large,
     .           nodos(40),nodosp(40),nodosp1(40),nodosmp(40),
     .           nodosm1p(40),nodosd(40),
     .           nodos1(40),nodos2(40),nodosr(40),nodosr1(40),
     .           nodosrt(40),
     .           icheckop,indexg,ntens,ninte,ndisp,ioptionr,numberr,
     .           kdyna,nskip
      real       x,y,z,stpval,xc,yc,ttime
      logical    finish, yes_forever
      integer, allocatable :: lnods(:,:),matno(:),reac1(:),reac2(:)
      real, allocatable :: proel(:,:),props(:,:),gpcod(:,:,:),
     .                     coord(:,:),disto(:,:),elvar(:,:)  ,
     .                     stren(:,:),fisno(:,:),strep(:,:)  ,
     .                     react(:,:),distr(:,:),dista(:,:)
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
c------------------------------------------------ Read general variables
      read(7,err=2,end=2)  ndime,ndofc,ndofn,ngaus,nhist,nnode,nelem,
     .                     ngrup,nmats,npoin,nprel,nprop,nstr1,title,
     .                     ndata,
     .                     nstat,ksgau,nnuin,npoic,nporo,large,kdyna
c------------------------------------------------------ Memory dimension
c
c     Note: nporo not used (only to compatibilize thermal & mechanical
c           problems for interfla.f)
c
c      plnods=1                                   !lnods
      allocate(lnods(nnode,nelem))
c      pmatno=plnods+nelem*nnode                  !matno
      allocate(matno(nelem))
c      pproel=pmatno+nelem                        !proel
      allocate(proel(nprel,ngrup))
c      pprops=pproel+nprel*ngrup*itwo             !props
      allocate(props(nprop,nmats))
c      pgpcod=pprops+nprop*nmats*itwo             !gpcod
      allocate(gpcod(ndime,ngaus,nelem))
c      pcoord=pgpcod+ndime*ngaus*nelem*itwo       !coord 
      allocate(coord(ndime,npoin))

c      lastmem=pcoord+npoin*ndime*itwo
c      if(lastmem.gt.maxlon) call outmem(lastmem)
      
c      pdisto =pcoord +npoin*ndime*itwo           !disto
      allocate(disto(ndofc,npoin))
c      pelvar =pdisto +ndofc*npoin*itwo           !elvar
      allocate(elvar(nstat,nelem))
c      pstren =pelvar +nstat*nelem*itwo           !stren
      allocate(stren(nstr1,npoin))
c      pfisno =pstren +nstr1*npoin*itwo           !fisno
      allocate(fisno(nstr1,npoin))
c      pstrep =pfisno +npoin*itwo*nstr1           !strep
      allocate(strep(nstr1,npoin))
c      preact =pstrep +npoin*itwo*nstr1           !react
      allocate(react(ndofc,npoin))
c      preac1 =preact +ndofc*npoin*itwo           !nodosrt1
      allocate(reac1(npoin))
c      preac2 =preac1 +npoin                      !nodosrt2
      allocate(reac2(npoin))
c      pdistr =preac2 +npoin                      !distr
      if (kdyna.ne.0) then
        allocate(distr(ndofc,npoin))
        allocate(dista(ndofc,npoin))
      else
        allocate(distr(1,1))
        allocate(dista(1,1))
      endif
c      pdista =pdistr +ndofc*npoin*itwo*kdyna     !dista
c      lastmem=pdista +ndofc*npoin*itwo*kdyna
c      if(lastmem.gt.maxlon) call outmem(lastmem)
c
      open(unit=8,file='Messages',access='sequential',form='formatted')
c--------------------------------------------------------- Read Geometry
      call readgeom(ndime,nnode,nelem,ngaus,ngrup,nmats,npoin,npoic,
     .              nprel,nprop,idata,istat,
     .              coord,gpcod,lnods,
     .              matno,proel,props)
      npoin1=npoin-npoic
c-------------------------------------recibe numeros de nodos requeridos
      yesno='Y'
      yesnot='Y'
      yesnop='Y'
      yesnopo='Y'
      yesnor='Y'
      ico=0
      icop=0
      icomp=0
      icor=0
      curva(1:3)='des'
      curvap(1:3)='ten'
      curvamp(1:3)='int'
      curvar(1:3)='dto'
      curvav(1:3)='von'
      curvare(1:3)='rea'
      curvaret(1:3)='ret'
c
c     sg
c
c     type '(''Total reaction calculation ? (y)es,(n)o: '',$)'
c     accept '(A30)',yesno
c
c     linux
c
      write(6,*) '(''Total reaction calculation ? (y)es,(n)o:'',$)'
      read(5,*) yesno
c
      ioptionr=1
      if(yesno.eq.'Y'.or.yesno.eq.'y') ioptionr=2
c
      if(ioptionr.eq.2) then
c
c     sg
c
c      type '(''Enter file with a node list( with extension ): '',$)'
c      accept '(A30)',inputfiler
c
c     linux
c
       write(6,*) '(''Enter file with a node list(with extension):'',$)'
       read(5,*) inputfiler
c
c     sg
c
c      open (unit=9, name=inputfiler, status='old', form='formatted',
c    .       err=1, readonly)
c
c     linux
c
       open (unit=9, file=inputfiler, status='old', form='formatted',
     .       err=1)
c
       read(9,*) numberr
       do i=1,numberr
        read(9,*) reac1(i),reac2(i)
       enddo
c
c      particular cases
c
       if(reac2(1).eq.4) then ! forces & torque for 2D rolling prob.
        read(9,*) xc,yc          ! coord. of roll center
       endif
      endif                      ! ioptionr.eq.2
c
      indexg=0
 9999 continue
      indexg=indexg+1            ! total number of curves per node
c
c     sg
c
c     type '(''Evolution of variables - Enter node number: '',$)'
c     accept '(i8)' ,nod
c
c     linux
c
      write(6,*) '(''Evolution of variables - Enter node number: '',$)'
      read(5,*) nod
c
c**** displacements (components)
c
      yesno='Y'
      ndisp=ndofn
      if(((yesnot.eq.'Y').OR.(yesnot.eq.'y')).AND.(indexg.le.10)) then
       do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(indexg.le.10))
        ico=ico+1
        nodos(ico)=nod
        nodosd(ico)=ico-(indexg-1)*ndisp
        if (ico.lt.10)then
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
c       write((100+ico),12) nodosd(ico)
c
c     linux
c
        open(unit=(10+ico),file=curva,
     .       access='sequential',form='formatted')
        write((10+ico),10) inputfile
        write((10+ico),11) nodos(ico)
        write((10+ico),12) nodosd(ico)
c
   10   format('#  example:',a30)
   11   format('#  node number:' ,i8)
   12   format('#  Displacement component (1=X, 2=Y, 3=Z):' ,i4)
c
        yesno='y'
        if((ico/indexg).eq.ndisp) then
         yesno='n'
         yesnot='n'
        endif
c
       end do
      endif
c
c**** total displacement
c
      nodos1(indexg)=nod
      if(indexg.lt.10)then
       write(curvar(4:5),'(I1)')indexg
      else
       write(curvar(4:5),'(I2)')indexg
      endif
c
c     sg
c
c     open(unit=(50+indexg),name=curvar,
c    .     access='sequential',form='formatted')
c     write((50+indexg),20) inputfile
c     write((50+indexg),21) nodos1(indexg)
c     write((50+indexg),22)
c
c     linux
c
      open(unit=(20+indexg),file=curvar,
     .     access='sequential',form='formatted')
      write((20+indexg),20) inputfile
      write((20+indexg),21) nodos1(indexg)
      write((20+indexg),22)
c
   20 format('#  example:',a30)
   21 format('#  node number:' ,i8)
   22 format('#  Total displacement:')
c
c**** reactions
c
      yesno='Y'
      if(((yesnor.eq.'Y').OR.(yesnor.eq.'y')).AND.(indexg.le.10)) then
       do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(indexg.le.10))
        icor=icor+1
        nodosr(icor)=nod
        nodosr1(icor)=icor-(indexg-1)*ndisp
        if(icor.lt.10)then
         write(curvare(4:5),'(I1)')icor
        else
         write(curvare(4:5),'(I2)')icor
        endif
c
c     sg
c
c       open(unit=(80+icor),name=curvare,
c    .       access='sequential',form='formatted')
c       write((80+icor),80) inputfile
c       write((80+icor),81) nodosr(icor)
c       write((80+icor),82) nodosr1(icor)
c
c     linux
c
        open(unit=(80+icor),file=curvare,
     .       access='sequential',form='formatted')
        write((80+icor),80) inputfile
        write((80+icor),81) nodosr(icor)
        write((80+icor),82) nodosr1(icor)
c
   80   format('#  example:',a30)
   81   format('#  node number:',i8)
   82   format('#  Reaction component (1=X, 2=Y, 3=Z):' ,i4)
c
        yesno='y'
        if((icor/indexg).eq.ndisp) then
         yesno='n'
         yesnor='n'
        endif
c
       end do
      endif
c
c**** total reactions
c
c     two options: 1) total reaction at a particular node (uses nodosrt)
c                  2) total reaction of a list of nodes (uses nodosrt1 &
c                     nodosrt2 with exact dimension a priori unknown but
c                     with a maximum dimension of npoin)
c
      if(ioptionr.eq.1) then
       nodosrt(indexg)=nod
      endif                                 ! ioptionr.eq.1
c
      if(indexg.lt.10)then
       write(curvaret(4:5),'(I1)')indexg
      else
       write(curvaret(4:5),'(I2)')indexg
      endif
c
c     sg
c
c     open(unit=(90+indexg),name=curvaret,
c    .     access='sequential',form='formatted')
c     write((90+indexg),20) inputfile
c     if(ioptionr.eq.1) then
c      write((90+indexg),21) nodosrt(indexg)
c     else
c      write((90+indexg),22)
c     endif
c     write((90+indexg),93)
c
c     linux
c
      open(unit=(90+indexg),file=curvaret,
     .     access='sequential',form='formatted')
      write((90+indexg),90) inputfile
      if(ioptionr.eq.1) then
       write((90+indexg),91) nodosrt(indexg)
      else
       write((90+indexg),92) inputfiler
      endif
      if(reac2(1).le.3) write((90+indexg),93)
      if(reac2(1).eq.4) write((90+indexg),94)
c
   90 format('#  example:',a30)
   91 format('#  node number:',i8)
   92 format('#  list of nodes:',a30)
   93 format('#  Total reaction')
   94 format('#  Total reaction y, Torque, x_coordinate')
c
c**** stresses & strains (components)
c
      if(ndime.eq.2) ntens=4                         ! stress components
      if(ndime.eq.3) ntens=6                         ! stress components
      if(ntens.gt.0) then
       if(ksgau.ne.0) then
        yesno='Y'
c
        if(((yesnop.eq.'Y').OR.(yesnop.eq.'y')).AND.(indexg.le.10)) then
         do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(indexg.le.10))
          icop=icop+1
          nodosp(icop)=nod                        ! former
          nodosp1(icop)=icop-(indexg-1)*ntens     ! always=1
c
          if (icop.lt.10)then
           write(curvap(4:5),'(I1)')icop
          else
           write(curvap(4:5),'(I2)')icop
          endif
c
c     sg
c
c         open(unit=(300+icop),name=curvap,
c    .         access='sequential',form='formatted')
c         write((300+icop),30) inputfile
c         write((300+icop),31) nodosp(icop)
c         write((300+icop),32) nodosp1(icop)
c
c     linux
c
          open(unit=(30+icop),file=curvap,
     .         access='sequential',form='formatted')
          write((30+icop),30) inputfile
          write((30+icop),31) nodosp(icop)
          write((30+icop),32) nodosp1(icop)
c
   30     format('#  example:',a30)
   31     format('#  node number:' ,i4)
   32     format('#  Stress & strain component:' ,i4)
c
          yesno='y'
          if((icop/indexg).eq.ntens) then
           yesno='n'
           yesnop='n'
          endif
c
         enddo
        endif
       endif       ! ksgau.ne.0
      endif        ! ntens.gt.0
c
c**** Von Mises stress (vms=sqrt(3J2)) & I1/vms parameter
c
      if(ksgau.ne.0) then
       nodos2(indexg)=nod
       if(indexg.lt.10)then
        write(curvav(4:5),'(I1)')indexg
       else
        write(curvav(4:5),'(I2)')indexg
       endif
c
c     sg
c
c      open(unit=(70+indexg),name=curvav,
c    .      access='sequential',form='formatted')
c      write((70+indexg),50) inputfile
c      write((70+indexg),51) nodos2(indexg)
c      write((70+indexg),52)
c
c     linux
c
       open(unit=(40+indexg),file=curvav,
     .      access='sequential',form='formatted')
       write((40+indexg),50) inputfile
       write((40+indexg),51) nodos2(indexg)
       write((40+indexg),52)
c
      endif                               ! ksgau.ne.0
c
   50 format('#  example:',a30)
   51 format('#  node number:' ,i8)
c  52 format('#  Von Mises stresses (vms) & I1/vms:')
   52 format('#  Von Mises stresses (vms) & I1:')
c
c**** Internal variables
c
      if(ksgau.ne.0.and.nnuin.ne.0) then
       ninte=nnuin
       yesno='Y'
c
       if(((yesnopo.eq.'Y').OR.(yesnopo.eq.'y')).AND.
     .                                              (indexg.le.10)) then
        do while(((yesno.eq.'Y').OR.(yesno.eq.'y')).AND.(indexg.le.10))
         icomp=icomp+1
         nodosmp(icomp)=nod                                     ! former
         nodosm1p(icomp)=icomp-(indexg-1)*ninte
c
         if(icomp.lt.10)then
          write(curvamp(4:5),'(I1)')icomp
         else
          write(curvamp(4:5),'(I2)')icomp
         endif
c
c    sg
c
c        open(unit=(250+icomp),name=curvamp,
c    .        access='sequential',form='formatted')
c        write((250+icomp),40) inputfile
c        write((250+icomp),41) nodosmp(icomp)
c        write((250+icomp),42) nodosm1p(icomp)
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
c  42    format('#  Effective plastic strain:' ,i4)
   42    format('#  Internal variable:' ,i4)
c
         yesno='y'
         if((icomp/indexg).eq.ninte) then
          yesno='n'
          yesnopo='n'
         endif
c
        enddo
       endif
      endif                         ! ksgau.ne.0.and.nnuin.ne.0
c
c     avoid questions
c
c
c     sg
c
c     type '(''do you wish another curve ? { y)es, n)o }: '',$)'
c     accept '(a1)',yesno
c
c     linux
c
      write(6,*) '(''do you wish another curve ? { y)es, n)o }:'',$)'
      read(5,*) yesno
c
      if(yesno.eq.'y'.or.yesno.eq.'Y') then
       yesnot='y'
       yesnor='y'
       yesnop='y'
       yesnopo='y'
       go to 9999
      endif
c------------------------------------- Prepare to read and write results
      numstp=1
      finish= .false.
      yes_forever=.false.
      nskip=0
      do while(.not.finish) ! ------------------- Read and Write results
       call readres(disto,elvar,fisno,iiter,istep,
     .              itime,ndime,ndofc,nelem,npoin,nstat,nstr1,subti,
     .              title,ttime,stren,ksgau,finish,nnuin,
     .              strep,react,distr,dista,
     .              npoic,kdyna)
       if(.not.finish) then ! Deal with it
        call wrtres(disto,elvar,fisno,gpcod,
     .              lnods,matno,proel,props,
     .              stren,coord,strep,react,
     .              reac1,reac2,distr,dista,
     .              kprob,ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,
     .              nmats,nnode,npoin,nprel,nprop,nstat,nstr1,iiter,
     .              istat,istep,itime,ttime,ksgau,numstp,name,
     .              yesno,npoin1,atype,text,gstre,gstra,nnuin,
     .              npoic,
     .              ico,icop,icomp,icor,
     .              nodos,nodosp,nodosp1,nodosmp,nodosm1p,nodosd,
     .              nodos1,nodos2,nodosr,nodosr1,nodosrt,
     .              ioptionr,numberr,indexg,kdyna,xc,yc)
        numstp=numstp+1
       end if                                               ! not finish
      end do                                          ! while not finish
    3 close(unit=8,status='delete')
      call exit
c   2 call error('read general variab.',0,0)                  ! old form
    2 call error('read general variab.')
      end                               ! Main program Vulcan -> gnuplot
c
c----------------------- S U B R O U T I N E S -------------------------
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
     .          idata(9),                istat(4)
      real      coord(ndime,npoin-npoic), gpcod(ndime,ngaus,nelem),
     .          proel(nprel,ngrup),       props(nprop,nmats)
c
c     do ielem=1,nelem                     ! coordinates at Gauss points
c      read(7,err=1,end=1) ((gpcod(idime,igaus,ielem),
c    .                        idime=1,ndime),igaus=1,ngaus)
c     end do
c
      read(7,err=1,end=1) coord,lnods,matno,proel,props,idata,istat
      return
c   1 call error('read geometry  read ',0,0)                  ! old form
    1 call error('read geometry  read ')
      end                                                     ! readgeom
c---------------------------------------------------------- Read results
      subroutine readres(disto,elvar,fisno,iiter,istep,itime,
     .                   ndime,ndofc,nelem,npoin,nstat,nstr1,
     .                   subti,title,ttime,stren,ksgau,finish,nnuin,
     .                   strep,react,distr,dista,
     .                   npoic,kdyna)
      implicit none
      integer   iiter,istep,itime,ndime,ndofc,nelem,npoin,
     .          nstat,nstr1,ksgau,ielem,i,nnuin,npoic,kdyna
      real      disto(ndofc,npoin-npoic), elvar(nstat,nelem),
     .          stren(nstr1,npoin-npoic), fisno(nstr1,npoin-npoic),
     .          strep(nnuin,npoin-npoic), react(ndofc,npoin-npoic),
     .          distr(ndofc,npoin-npoic), dista(ndofc,npoin-npoic),
     .          ttime
c
      character subti*64,title*64
      logical finish
c
c**** stren: nodal stresses
c     fisno: nodal strains
c     strep: nodal internal variables (nnuin > 0 always)
c     react: nodal reactions
c     distr: nodal velocities
c     dista: nodal accelerations
c
c**** displacements
c
    1 read(7,err=1,end=2) title,subti,itime,istep,iiter,ttime,disto
c
c**** velocities & accelerations
c
      if(kdyna.eq.1) then
       read(7,err=1,end=2) distr
       read(7,err=1,end=2) dista
      endif
c
c     do ielem=1,nelem              ! Internal variables at Gauss points
c      read(7) (elvar(i,ielem),i=1,nstat)
c     end do
c
c**** reactions
c
      read(7) react
c
c**** nodal stresses, strains and internal variables
c
      if(ksgau.ne.0) then
       read(7) stren,fisno,strep
      end if
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
      subroutine wrtres(disto,elvar,fisno,gpcod,lnods,matno,proel,props,
     .                  stren,
     .                  coord,strep,react,nodosrt1,nodosrt2,distr,dista,
     .                  kprob,ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,
     .                  nmats,
     .                  nnode,npoin,nprel,nprop,nstat,nstr1,iiter,
     .                  istat,istep,itime,ttime,ksgau,numstp,name,
     .                  yesno,npoin1,atype,text,gstre,gstra,nnuin,
     .                  npoic,
     .                  ico,icop,icomp,icor,
     .                  nodos,nodosp,nodosp1,nodosmp,nodosm1p,nodosd,
     .                  nodos1,nodos2,nodosr,nodosr1,nodosrt,
     .                  ioptionr,numberr,indexg,kdyna,xc,yc)
      implicit none
      integer numstp,lnods(*),matno(*),istat(4),
     .        kprob,ndime,ndofc,ndofn,nelem,ngaus,ngrup,nhist,nmats,
     .        nnode,npoin,nprel,nprop,nstat,nstr1,iiter,i,idofc,
     .        istep,itime,ksgau,ncrit,npoin1,atype,key,
     .        ncomps,irtype,menu,ictype,icind1,icind2,iexist,
     .        ipoin,imats,nsrf,istr1,numb,nnuin,npoic,ioptionr,numberr,
     .        kdyna
c
      integer ico,icop,icomp,icor,nodos(*),
     .        nodosp(*),nodosp1(*),nodosmp(*),nodosm1p(*),nodosd(*),
     .        indexg,nodos1(*),nodos2(*),nodosr(*),nodosr1(*),
     .        nodosrt(*),nodosrt1(*),nodosrt2(*)
c
      real    disto(ndofc,npoin1),elvar(*),fisno(nstr1,npoin1),ttime,
     .        gpcod(*),proel(*),props(*),stren(nstr1,npoin1),
     .        coord(ndime,npoin1),
     .        strep(nnuin,npoin1),react(ndofc,npoin1),
     .        distr(ndofc,npoin1),dista(ndofc,npoin1)
c
      real    pmean,pdevi(6),pvonm,totdi,ai1j2,
     .        totdix,totdiy,tottox,tottoy,tottor,xc,yc,xa,ya,
     .        totnor,tottan,amu,rroll
c
      character name*6, yesno*1, opcode*1, text*20, namet*8,
     .          icname*8, gstra*1, gstre*1
c--------------------------------------------------------- Displacements
      do i=1,ico
c
c     sg
c
c      write((100+i),*) ttime, disto(nodosd(i),nodos(i))
c
c     linux
c
       if(kdyna.eq.0) then
        write((10+i),*) ttime, disto(nodosd(i),nodos(i))
       else
        write((10+i),*) ttime, disto(nodosd(i),nodos(i)),
     .                         distr(nodosd(i),nodos(i)),
     .                         dista(nodosd(i),nodos(i))
       endif
c
      end do
c-----------------------------------------------------------------------
c--------------------------------------------------- Total displacements
      do i=1,indexg
       if(ndime.eq.2) then
        ipoin=nodos1(i)
        totdi=disto(1,ipoin)*disto(1,ipoin)+
     .        disto(2,ipoin)*disto(2,ipoin)
        totdi=sqrt(totdi)
c
c     sg
c
c       write((50+i),*) ttime, totdi
c
c     linux
c
        write((20+i),*) ttime, totdi
c
       endif
c
       if(ndime.eq.3) then
        ipoin=nodos1(i)
        totdi=disto(1,ipoin)*disto(1,ipoin)+
     .        disto(2,ipoin)*disto(2,ipoin)+
     .        disto(3,ipoin)*disto(3,ipoin)
        totdi=sqrt(totdi)
c
c     sg
c
c       write((50+i),*) ttime, totdi
c
c     linux
c
        write((20+i),*) ttime, totdi
c
       endif
      end do
c-----------------------------------------------------------------------
c------------------------------------------------------------- Reactions
      do i=1,icor
c
c     sg
c
c      write((80+i),*) ttime, react(nodosr1(i),nodosr(i))
c
c     linux
c
       write((80+i),*) ttime, react(nodosr1(i),nodosr(i))
c
      end do
c-----------------------------------------------------------------------
c------------------------------------------------------- Total reactions
      if(ioptionr.eq.1) then
       do i=1,indexg
        if(ndime.eq.2) then
         ipoin=nodosrt(i)
         totdi=react(1,ipoin)*react(1,ipoin)+
     .         react(2,ipoin)*react(2,ipoin)
         totdi=sqrt(totdi)
c
c     sg
c
c        write((90+i),*) ttime, totdi
c
c     linux
c
         write((90+i),*) ttime, totdi
c
        endif
c
        if(ndime.eq.3) then
         ipoin=nodosrt(i)
         totdi=react(1,ipoin)*react(1,ipoin)+
     .         react(2,ipoin)*react(2,ipoin)+
     .         react(3,ipoin)*react(3,ipoin)
         totdi=sqrt(totdi)
c
c     sg
c
c        write((90+i),*) ttime, totdi
c
c     linux
c
         write((90+i),*) ttime, totdi
c
        endif
       end do
      endif                                  ! ioptionr.eq.1
C
      if(ioptionr.eq.2) then
       if(nodosrt2(1).le.3) then             ! standard reaction
        totdi=0.0
        do i=1,numberr
         ipoin=nodosrt1(i)
         idofc=nodosrt2(i)
         totdi=totdi+react(idofc,ipoin)
        enddo
c
c     sg
c
c       write((90+1),*) ttime, totdi         ! approx: only one reaction
c
c     linux
c
        write((90+1),*) ttime, totdi         ! approx: only one reaction
       endif
c
       if(nodosrt2(1).eq.4) then             ! part. case: 2D rolling
        totdix=0.0
        totdiy=0.0
        tottox=0.0
        tottoy=0.0
        xa=0.0                               ! just to avoid NAN
        ya=0.0
        do i=1,numberr
         ipoin=nodosrt1(i)
         totdix=totdix+react(1,ipoin)             ! force x
         totdiy=totdiy+react(2,ipoin)             ! force y
         tottox=tottox+react(1,ipoin)*(coord(2,ipoin)+disto(2,ipoin)-yc)
         tottoy=tottoy+react(2,ipoin)*(coord(1,ipoin)+disto(1,ipoin)-xc)
        enddo
        tottor=-tottox+tottoy
c       if(totdiy.ne.0.0) xa=tottoy/totdiy+xc     ! absolute coordinates
c       if(totdix.ne.0.0) ya=tottox/totdix+yc
        if(totdiy.ne.0.0) xa=tottoy/totdiy        ! relative coordinates
        if(totdix.ne.0.0) ya=tottox/totdix
c       rroll=sqrt((xa-xc)*(xa-xc)+               ! roll radius (abs.c.)
c    .             (ya-yc)*(ya-yc))
        rroll=sqrt(xa*xa+ya*ya)                   ! roll radius (rel.c.)
        tottan=0.0
        if(rroll.ne.0.0) tottan=tottor/rroll      ! tangential force
        totnor=sqrt(totdix*totdix+totdiy*totdiy-  ! normal force
     .              tottan*tottan)
        amu=0.0
        if(totnor.ne.0.0) amu=tottan/totnor       ! global friction
c
c     sg
c
c       write((90+1),*) ttime, totdix, totdiy, tottor, xa, ya
c
c     linux
c
c       write((90+1),*) ttime, totdix, totdiy, tottor, xa, ya
        write((90+1),*) ttime, totdiy, tottor, xa
c       write((90+1),*) ttime, tottor
c       write((90+1),*) ttime, xa
c       write((90+1),*) ttime, ya
       endif
      endif                                  ! ioptionr.eq.2
c-----------------------------------------------------------------------
c-------------------------------------------------------------- Stresses
      do i=1,icop
       ipoin=nodosp(i)
       istr1=nodosp1(i)
c
c     sg
c
c      write((300+i),*) ttime, stren(istr1,ipoin), fisno(istr1,ipoin)
c
c     linux
c
       write((30+i),*) ttime, stren(istr1,ipoin), fisno(istr1,ipoin)
c
      end do
c-----------------------------------------------------------------------
c--------------------------------------------------------------Von Mises
c
c**** Von Mises stress (vms=sqrt(3J2)) & I1/vms parameter
c
      if(icomp.gt.0) then
       do i=1,indexg
        if(ndime.eq.2) then
         ipoin=nodos2(i)
         pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
         pdevi(1)=stren(1,ipoin)-pmean
         pdevi(2)=stren(2,ipoin)-pmean
         pdevi(3)=stren(3,ipoin)
         pdevi(4)=stren(4,ipoin)-pmean
         pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .     2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)
         pvonm=sqrt(3.0/2.0*pvonm)
c
c        ai1j2=0.0
c        if(pvonm.gt.1.0e-08) ai1j2=pmean/pvonm
         ai1j2=pmean
c
c     sg
c
c        write((70+i),*) ttime, pvonm,ai1j2
c
c     linux
c
         write((40+i),*) ttime, pvonm,ai1j2
c
        endif
c
        if(ndime.eq.3) then
         ipoin=nodos2(i)
         pmean=(stren(1,ipoin)+stren(2,ipoin)+stren(4,ipoin))/3.0
         pdevi(1)=stren(1,ipoin)-pmean
         pdevi(2)=stren(2,ipoin)-pmean
         pdevi(3)=stren(3,ipoin)
         pdevi(4)=stren(4,ipoin)-pmean
         pdevi(5)=stren(5,ipoin)
         pdevi(6)=stren(6,ipoin)
         pvonm=pdevi(1)*pdevi(1)+pdevi(2)*pdevi(2)+
     .     2.0*pdevi(3)*pdevi(3)+pdevi(4)*pdevi(4)+
     .     2.0*pdevi(5)*pdevi(5)+2.0*pdevi(6)*pdevi(6)
         pvonm=sqrt(3.0/2.0*pvonm)
c
c        ai1j2=0.0
c        if(pvonm.gt.1.0e-08) ai1j2=pmean/pvonm
         ai1j2=pmean
c
c     sg
c
c        write((70+i),*) ttime, pvonm,ai1j2
c
c     linux
c
         write((40+i),*) ttime, pvonm,ai1j2
c
        endif
       enddo
      endif
c-----------------------------------------------------------------------
c---------------------------------------------------- Internal variables
      if(icomp.gt.0) then
       do i=1,icomp
c
c     sg
c
c       write((250+i),*) ttime, strep(nodosm1p(i),nodosmp(i))
c
c     linux
c
        write((50+i),*) ttime, strep(nodosm1p(i),nodosmp(i))
c
       end do
      endif
c
      end                                                       ! wrtres
c-----------------------------------------------------------------------
c--------------------------------------------------------- Out of memory
      subroutine outmem(lastmem)
      implicit none
      integer lastmem
c
    1 format('Insuficient memory.',/,'Dimension array mem to more than',
     . i10)
      open(unit=8,file='Messages',access='sequential',form='formatted')
      write(8,1) lastmem
c     call error('Memory dimension    ',0,0)                  ! old form
      call error('Memory dimension    ')
      end                                                ! out of memory
c----------------------------------------------------------------- Error
      subroutine error(message)
      implicit none
      character message*20
c
c   1 format(' Error in: ',a20,', Code:',i5)                  ! old form
    1 format(' Error in: ',a20,'')
      write(8,1) message
c
c     sg
c
c     type '(//,''Intfgvm: ** Error was found. LOOK MESSAGE FILE **'')'
c
c
c
      write(6,*)'(//,''Intfgvm:* Error was found.LOOK MESSAGE FILE *'')'
c
      call exit
      end                                                        ! error
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
c============================== E N D ==================================
