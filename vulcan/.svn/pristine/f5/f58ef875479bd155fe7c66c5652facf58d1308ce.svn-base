      program fantout
c
c***********************************************************************
c                                                                      *
c                          F A N T O U T                               *
c                                                                      *
c      This program reads FANTOM's neutral file and writes the files   *
c      to be used by FEMVIEW, FLAVIA (GiD) or PATRAN.                  *
c      In the case of FEMVIEW or FLAVIA (GiD) some elements have to be *
c      subdivided. The translation of topology is made according to the*
c      following correspondence:                                       *
c                                                                      *
c      FEMVIEW:  TRI   NNODE = 4, 7      ---> 3, 4                     *
c                QUA   NNODE = 5, 9      ---> 3, 4                     *
c                TET   NNODE = 5, 11, 15 ---> 4, 4, 8                  *
c                HEX   NNODE = 9, 27     ---> 4, 8                     *
c                                                                      *
c      FLAVIA:   TRI   NNODE = 4, 7      ---> 3, 4                     *
c      (GiD)     QUA   NNODE = 5         ---> 3                        *
c                TET   not implemented yet                             *
c                HEX   not implemented yet                             *
c                                                                      *
c***********************************************************************
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      parameter    (mxnod=27, mxelm=60000, mxpoi=60000, mxdat=11,
     .              mxtot=mxpoi*mxdat)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      common/precon/kprec 
      character*74  title
      character*20  filen
      character*6   s_tit
      character*5   wopos(mxdat)
      character*4   ext
      dimension     lnods(mxnod,mxelm), lnodw(12,12*mxelm),
     .              coord(3,mxpoi), bridge(mxtot)
c
      call opfile(nin,nog,nor,i_pos,filen,ext,length,s_tit)
c
c...  read and write control data and geometry
c
      read(nin) ndime,nnode,nelty,nelem,npoin,ktemp,kprec,title
c
      if(nelem.gt.mxelm) then                                 ! controls
       write(6,*) 'WARNING: increase mxelm !'
       stop
      endif
      if(npoin.gt.mxpoi) then
       write(6,*) 'WARNING: increase mxpoi !'
       stop
      endif
c
      call geomet(lnods,lnodw,coord,s_tit)
c
c...  read and write nodal results
c
      call result(bridge,wopos,filen,ext,length)
c
      end
      subroutine f_comp(fname,filen,fifix,length,ext,len,ifext)
c
      character*20 fname,filen,fifix
      character*4  ext
      integer      length,len,ifext
c
      if(fname(1:3).eq.'   ') then
        filen = fifix(1:len)//ext
      elseif(ifext.eq.0) then
        filen = fname(1:len)//ext
      else
        filen = fname(1:length)
      endif
c
      end
      subroutine f_leng(fname,length,len,lfix,ifext)
c
      character*20 fname
      integer      length,len,lfix,ifext
c
      do i = 1,20
        if (fname(i:i).ne.' ') length = length+1
        if (fname(i:i).eq.'.') len = i-1
      enddo
c
      if(lfix.eq.0.and.len.eq.0) then
        len = length
        length = length+4
        ifext = 0
      elseif(length.eq.0) then
        len = lfix
        length = lfix+4
        ifext = 0
      elseif(len.eq.0) then
        len = length
        length = length+4
        ifext = 0
      endif
c
      end
      subroutine f_read(filen,len,ext,fname)
c
      character*20 filen,fname
      character*4  ext
      integer      len
c
      if(ext.eq.'.pos') then
        write(6,10) ext
        read (5,'(a)') fname
      elseif(ext.eq.'.fem') then
        write(6,20) filen(1:len),ext
        read (5,'(a)') fname
      elseif(ext.eq.'.dat') then
        write(6,50)
        write(6,30) filen(1:len),ext
        read (5,'(a)') fname
      elseif(ext.eq.'.res') then
        write(6,40) filen(1:len),ext
        read (5,'(a)') fname
      endif
c
   10 format('neutral file name [default = *',a,']: ',$)
   20 format('post-process file name [default = ',a,a,']: ',$)
   30 format('post-process geometry file name [default = ',a,a,']: ',$)
   40 format('post-process results file name  [default = ',a,a,']: ',$)
   50 format('* * WARNING: for 3D problems, use extension "msh" * *')
c
      end
      subroutine geofem(lnods,lnodw,coord,s_tit)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*6   s_tit
      dimension     lnods(nnode,nelem), lnodw(12,*), coord(ndime,npoin)
c
      write(nog,10) 1,'C',s_tit
      write(nog,10) 2,'C','Coord.'
      do ipoin  = 1,npoin
        write(nog,20) -1,ipoin,(coord(i,ipoin),i=1,ndime)
      enddo
      write(nog,30) -3
      write(nog,10) 3,'C','Conne.'
c
      nnodw = 0
      if(nnode.eq.4.and.nelty.eq.2.and.ndime.eq.2) then
        call top004(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 7
      endif
      if(nnode.eq.5.and.ndime.eq.2) then
        call top05a(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 7
      endif
      if(nnode.eq.5.and.ndime.eq.3) then
        call top05b(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 3
      endif
      if(nnode.eq.7) then
        call top007(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 9
      endif
      if(nnode.eq.9.and.ndime.eq.2) then
        call top09a(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 9
      endif
      if(nnode.eq.9.and.ndime.eq.3) then
        call top09b(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 3
      endif
      if(nnode.eq.11) then
        call top011(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 3
      endif
      if(nnode.eq.15) then
        call top015(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 1
      endif
      if(nnode.eq.27) then
        call top027(lnods,nelem,nnode,lnodw,melem,nnodw)
        iflag = 1
      endif
c
      if(nnodw.eq.0) then
        if(nnode.eq.3)                iflag = 7
        if(nnode.eq.4.and.ndime.eq.2) iflag = 9
        if(nnode.eq.4.and.ndime.eq.3) iflag = 3
        if(nnode.eq.6)                iflag = 8
        if(nnode.eq.8.and.ndime.eq.2) iflag = 10
        if(nnode.eq.8.and.ndime.eq.3) iflag = 1
        if(nnode.eq.10)               iflag = 6
        if(nnode.eq.20)               iflag = 17
        do ielem = 1,nelem
          write(nog,40) -1,ielem,iflag,1,1
          write(nog,40) -2,(lnods(inode,ielem),inode=1,nnode)
        end do
      else 
        do ielem = 1,melem
          write(nog,40) -1,ielem,iflag,1,1
          write(nog,40) -2,(lnodw(inode,ielem),inode=1,nnodw)
        end do
      end if
      write(nog,30) -3
c
   10 format(1x,i4,a1,a6)
   20 format(1x,i2,i5,3e12.5)
   30 format(1x,i2)
   40 format(1x,i2,20i5)
c
      end 
      subroutine geofla(lnods,lnodw,coord,s_tit)
c
c     ********************* only for 2-D **************************
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*6   s_tit
      dimension     lnods(nnode,nelem), lnodw(12,*), coord(ndime,npoin)
c
      nnodw = 0
      if(nnode.eq.4.and.nelty.eq.2)
     .               call top004(lnods,nelem,nnode,lnodw,melem,nnodw)
      if(nnode.eq.5) call top05a(lnods,nelem,nnode,lnodw,melem,nnodw)
c     if(nnode.eq.6) call top006(lnods,nelem,nnode,lnodw,melem,nnodw)
      if(nnode.eq.7) call top007(lnods,nelem,nnode,lnodw,melem,nnodw)
c     if(nnode.eq.8) call top008(lnods,nelem,nnode,lnodw,melem,nnodw)
c     if(nnode.eq.9) call top09a(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      if(nnode.eq.3.or.nnodw.eq.3) then
        write(nog,10) s_tit,nelem,npoin
      else
        if(ndime.eq.2) write(nog,20) s_tit,nelem,npoin,nnode
        if(ndime.eq.3) then                    ! 3D case: see interfla.f
          nnodx=0
          if(nnode.eq.4) nnodx=3
          if(nnode.eq.8) nnodx=1
          if(nnodx.eq.0) then
            write(6,*) 'ERROR: this element is not implemented !'
            stop
          endif 
          write(nog,20) s_tit,nelem,npoin,nnodx
        endif
      endif
c
c...  write coordinates 
c
      do ipoin=1,npoin
        write(nog,30) ipoin,(coord(idime,ipoin),idime=1,ndime)
      enddo
c
      write(nog,40)
c
c...  write nodal set for elements
c
      if(nnodw.eq.0) then
        do ielem=1,nelem
          if(ndime.eq.2) then
            if(nnode.eq.6.or.nnode.eq.8.or.nnode.eq.9) 
     .        call nodren(nnode,lnods(1,ielem))
          endif
          write(nog,50) ielem,(lnods(inode,ielem),inode=1,nnode)
        end do
      else
        do ielem=1,melem
          write(nog,50) ielem,(lnodw(inode,ielem),inode=1,nnodw)
        end do
      endif
c
      close (nog)
c
   10 format(a6/'version 1.0',4(/),'nelem  npoin'/2i5/'coordenadas')
   20 format(a6/'version 1.0',4(/),'nelem  npoin  nnode'/3i7/
     .       'coordenadas')
c  30 format(7x,i7,3e15.7)
   30 format(i8,3e15.6)                ! idem interfla
   40 format('conectividades')
c  50 format(i7,5x,8i7)
   50 format(i8,8i8)                   ! idem interfla
c
      end 
      subroutine geomet(lnods,lnodw,coord,s_tit)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*6   s_tit
      dimension     lnods(nnode,nelem), lnodw(*), coord(ndime,npoin)
c
      read(nin) ((coord(idime,ipoin),idime=1,ndime),ipoin=1,npoin),
     .            lnods
      if(i_pos.eq.1) then
        call geofem(lnods,lnodw,coord,s_tit)
      elseif(i_pos.eq.2.or.i_pos.eq.4) then
        call geofla(lnods,lnodw,coord,s_tit)
      elseif(i_pos.eq.3) then
        call geopat(lnods,coord,s_tit)
      endif
c
      end
      subroutine geopat(lnods,coord,s_tit)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*6   s_tit
      dimension     lnods(nnode,nelem), coord(ndime,npoin)
c
      write(nog,10) s_tit
c
      if(nelty.eq.1) then
        ishap = 4
        if(ndime.eq.3) ishap = 8
      elseif(nelty.eq.2) then
        ishap = 3
        if(ndime.eq.3) ishap = 5
      endif
c
      kc = 1+int((nnode+9)/10)
c
      iline = int(nnode/10)
      if(mod(nnode,10).gt.0.or.iline.eq.0) iline = iline+1
c
c...  write total number of nodes and elements
c
      write(nog,20) npoin,nelem,0,0,0
      write(nog,*)
c
c...  write coordinates
c
      veloz = 0.0d0
      do ipoin=1,npoin
        write(nog,30) ipoin
        if(ndime.eq.2) then
          write(nog,40) (coord(idime,ipoin),idime=1,ndime),veloz
        else
          write(nog,40) (coord(idime,ipoin),idime=1,ndime)
        endif
        write(nog,50)
      enddo
c
c...  write nodal connectivities
c
      do ielem = 1,nelem
c
c...    renumber nodal connectivities for PATRAN
c
        if(nnode.eq.15) then 
          call nodren(nnode,lnods(11,ielem))
        elseif(nnode.eq.27) then 
          call nodren(nnode,lnods(21,ielem))
        endif
c
c...    write nodal connectivities
c
        write(nog,60) ielem,ishap,kc,0,0,0,0,0
        write(nog,70) nnode,0,0,0,0,0,0
        k = 0
        do j = 1,iline
          if(j.lt.iline) then
            write(nog,80) (lnods(k+inode,ielem),inode=1,10)
          else
            write(nog,80) (lnods(inode,ielem),inode=k+1,nnode)
          endif
          k = k+10
        enddo
      end do
c
      write(nog,90) 0,0,1,0,0,0,0,0
c
   10 format('25',2(7x,'0'),7x,'1',/,a6)
   20 format('26',2(7x,'0'),7x,'1',5i8)
   30 format(' 1',i8,7x,'0',7x,'2')
   40 format(3e16.9)
   50 format('1G',7x,'6',2(7x,'0'),'  000000')
   60 format(' 2',8i8)
   70 format(4i8,3e16.9)
   80 format(10i8)
   90 format('99',8i8)
c
      end 
      subroutine nodren(nnode,nrnode)
c
      dimension nrnode(7)
c
      if(nnode.eq.6) then
        nrsave = nrnode(2)
        nrnode(2) = nrnode(4)
        nrnode(4) = nrnode(5)
        nrnode(5) = nrnode(3)
        nrnode(3) = nrsave
      elseif(nnode.eq.8.or.nnode.eq.9) then
        nrsave = nrnode(2)
        nrnode(2) = nrnode(5)
        nrnode(5) = nrnode(3)
        nrnode(3) = nrsave
        nrsave = nrnode(4)
        nrnode(4) = nrnode(6)
        nrnode(6) = nrnode(7)
        nrnode(7) = nrsave
      elseif(nnode.eq.15) then
        nrsave = nrnode(1)
        nrnode(1) = nrnode(5)
        nrnode(5) = nrnode(2)
        nrnode(2) = nrsave
      elseif(nnode.eq.27) then
        nrsave = nrnode(1)
        nrnode(1) = nrnode(7)
        nrnode(7) = nrnode(4)
        nrnode(4) = nrnode(5)
        nrnode(5) = nrnode(3)
        nrnode(3) = nrnode(6)
        nrnode(6) = nrnode(2)
        nrnode(2) = nrsave
      endif
c
      end
      subroutine opfile(nin,nog,nor,i_pos,filer,exr,lengtr,s_tit)
c      
      character*4  ext(4), exr
      character*6  s_tit
      character*20 fname(4) , filen(4), filer
      dimension    length(4), len(4)  ,ifext(4)
c      
      nin = 20
      nog = 21
      nor = 22
c
c...  conventions for file extention (if not specified)
c
      ext(1) = '.pos'
      ext(2) = '.fem'
      ext(3) = '.dat'
      ext(4) = '.res'
c
c...  initialize some variables
c
      do i = 1,4
        length(i) = 0
        len(i)    = 0
        ifext(i)  = 1
      enddo
c
c...  read name of post-processor
c
      write(6,'(a40,a23,$)') 'post-processor: (1) FEMVIEW, (2) FLAVIA,',
     .                       ' (3) PATRAN, (4) GiD ? '
      read (5,'(i1)') i_pos
c
c...  read name of input-file
c
      call f_read(filen(1),len(1),ext(1),fname(1))
c
c...  calculate length with and without extention
c
      call f_leng(fname(1),length(1),len(1),len(1),ifext(1))
c
c...  compose filename if necessary
c
      call f_comp(fname(1),filen(1),fname(1),length(1),ext(1),len(1),
     .           ifext(1))
c
c...  read other filenames, calculate lengths & compose final filenames
c
      if(i_pos.eq.1) then
        call f_read(filen(1),len(1),ext(2),fname(2))
        call f_leng(fname(2),length(2),len(2),len(1),ifext(2))
        call f_comp(fname(2),filen(2),fname(1),length(2),ext(2),len(2),
     .             ifext(2))
        write(6,10) filen(2)(1:6)
        read (5,'(a)') s_tit
        if(s_tit(1:6).eq.'      ') s_tit = filen(2)(1:6)
      else
        call f_read(filen(1),len(1),ext(3),fname(3))
        call f_read(filen(1),len(1),ext(4),fname(4))
        call f_leng(fname(3),length(3),len(3),len(1),ifext(3))
        call f_leng(fname(4),length(4),len(4),len(1),ifext(4))
        call f_comp(fname(3),filen(3),fname(1),length(3),ext(3),len(3),
     .             ifext(3))
        call f_comp(fname(4),filen(4),fname(1),length(4),ext(4),len(4),
     .             ifext(4))
        write(6,10) filen(3)(1:6)
        read (5,'(a)') s_tit
        if(s_tit(1:6).eq.'      ') s_tit = filen(3)(1:6)
      endif
c
c...  open files
c
      open(nin,file=filen(1),status='old',form='unformatted')
      if (i_pos.eq.1) then
        open(nog,file=filen(2),status='unknown')
      else
        open(nog,file=filen(3),status='unknown')
      endif
c
      filer  = filen(4)
      exr    = ext(4)
      lengtr = length(4) 
c
   10 format('title (a6) [default = ',a,']: ',$)
c
      end
      subroutine opresu(fname,ext,len,istep,nor)
c
      integer istep,len,l4,nor
      character*20 fname,filen
      character*4  ext
      character*4  cstep4
      character*3  cstep3
      character*2  cstep2
      character*1  cstep1
c      
      l4 = len-4
c     if (istep.lt.10) then                 ! asi estaba y no funco (??)
c       encode(4,'(i1)',cstep1) istep
c       filen = fname(1:l4)//'_'//cstep1//ext
c
      if (istep.eq.1) then ! en .res irian las cond.inic.(a implementar)
        filen = fname(1:l4)//'_'//'1'//ext
c
c     sg
c
c     elseif (istep.lt.10) then
c       encode(4,'(i1)',cstep1) istep
c       filen = fname(1:l4)//'_'//cstep1//ext
c
c     elseif (istep.lt.100) then
c       encode(4,'(i2)',cstep2) istep
c       filen = fname(1:l4)//'_'//cstep2//ext
c     elseif (istep.lt.1000) then
c       encode(4,'(i3)',cstep3) istep
c       filen = fname(1:l4)//'_'//cstep3//ext
c     else
c       encode(4,'(i4)',cstep4) istep
c       filen = fname(1:l4)//'_'//cstep4//ext
c     endif
c
c     linux
c
      elseif (istep.eq.2) then
        filen = fname(1:l4)//'_'//'2'//ext
      elseif (istep.eq.3) then
        filen = fname(1:l4)//'_'//'3'//ext
      elseif (istep.eq.4) then
        filen = fname(1:l4)//'_'//'4'//ext
      elseif (istep.eq.5) then
        filen = fname(1:l4)//'_'//'5'//ext
      elseif (istep.eq.6) then
        filen = fname(1:l4)//'_'//'6'//ext
      elseif (istep.eq.7) then
        filen = fname(1:l4)//'_'//'7'//ext
      elseif (istep.eq.8) then
        filen = fname(1:l4)//'_'//'8'//ext
      elseif (istep.eq.9) then
        filen = fname(1:l4)//'_'//'9'//ext
      elseif (istep.eq.10) then
        filen = fname(1:l4)//'_'//'10'//ext
      elseif (istep.eq.11) then
        filen = fname(1:l4)//'_'//'11'//ext
      elseif (istep.eq.12) then
        filen = fname(1:l4)//'_'//'12'//ext
      elseif (istep.eq.13) then
        filen = fname(1:l4)//'_'//'13'//ext
      elseif (istep.eq.14) then
        filen = fname(1:l4)//'_'//'14'//ext
      elseif (istep.eq.15) then
        filen = fname(1:l4)//'_'//'15'//ext
      elseif (istep.eq.16) then
        filen = fname(1:l4)//'_'//'16'//ext
      elseif (istep.eq.17) then
        filen = fname(1:l4)//'_'//'17'//ext
      elseif (istep.eq.18) then
        filen = fname(1:l4)//'_'//'18'//ext
      elseif (istep.eq.19) then
        filen = fname(1:l4)//'_'//'19'//ext
      elseif (istep.eq.20) then
        filen = fname(1:l4)//'_'//'20'//ext
      elseif (istep.eq.21) then
        filen = fname(1:l4)//'_'//'21'//ext
      elseif (istep.eq.22) then
        filen = fname(1:l4)//'_'//'22'//ext
      elseif (istep.eq.23) then
        filen = fname(1:l4)//'_'//'23'//ext
      elseif (istep.eq.24) then
        filen = fname(1:l4)//'_'//'24'//ext
      elseif (istep.eq.25) then
        filen = fname(1:l4)//'_'//'25'//ext
      elseif (istep.eq.26) then
        filen = fname(1:l4)//'_'//'26'//ext
      elseif (istep.eq.27) then
        filen = fname(1:l4)//'_'//'27'//ext
      elseif (istep.eq.28) then
        filen = fname(1:l4)//'_'//'28'//ext
      elseif (istep.eq.29) then
        filen = fname(1:l4)//'_'//'29'//ext
      elseif (istep.eq.30) then
        filen = fname(1:l4)//'_'//'30'//ext
      elseif (istep.eq.31) then
        filen = fname(1:l4)//'_'//'31'//ext
      elseif (istep.eq.32) then
        filen = fname(1:l4)//'_'//'32'//ext
      elseif (istep.eq.33) then
        filen = fname(1:l4)//'_'//'33'//ext
      elseif (istep.eq.34) then
        filen = fname(1:l4)//'_'//'34'//ext
      elseif (istep.eq.35) then
        filen = fname(1:l4)//'_'//'35'//ext
      elseif (istep.eq.36) then
        filen = fname(1:l4)//'_'//'36'//ext
      elseif (istep.eq.37) then
        filen = fname(1:l4)//'_'//'37'//ext
      elseif (istep.eq.38) then
        filen = fname(1:l4)//'_'//'38'//ext
      elseif (istep.eq.39) then
        filen = fname(1:l4)//'_'//'39'//ext
      elseif (istep.eq.40) then
        filen = fname(1:l4)//'_'//'40'//ext
      elseif (istep.eq.41) then
        filen = fname(1:l4)//'_'//'41'//ext
      elseif (istep.eq.42) then
        filen = fname(1:l4)//'_'//'42'//ext
      elseif (istep.eq.43) then
        filen = fname(1:l4)//'_'//'43'//ext
      elseif (istep.eq.44) then
        filen = fname(1:l4)//'_'//'44'//ext
      elseif (istep.eq.45) then
        filen = fname(1:l4)//'_'//'45'//ext
      elseif (istep.eq.46) then
        filen = fname(1:l4)//'_'//'46'//ext
      elseif (istep.eq.47) then
        filen = fname(1:l4)//'_'//'47'//ext
      elseif (istep.eq.48) then
        filen = fname(1:l4)//'_'//'48'//ext
      elseif (istep.eq.49) then
        filen = fname(1:l4)//'_'//'49'//ext
      elseif (istep.eq.50) then
        filen = fname(1:l4)//'_'//'50'//ext
      elseif (istep.eq.51) then
        filen = fname(1:l4)//'_'//'51'//ext
      elseif (istep.eq.52) then
        filen = fname(1:l4)//'_'//'52'//ext
      elseif (istep.eq.53) then
        filen = fname(1:l4)//'_'//'53'//ext
      elseif (istep.eq.54) then
        filen = fname(1:l4)//'_'//'54'//ext
      elseif (istep.eq.55) then
        filen = fname(1:l4)//'_'//'55'//ext
      elseif (istep.eq.56) then
        filen = fname(1:l4)//'_'//'56'//ext
      elseif (istep.eq.57) then
        filen = fname(1:l4)//'_'//'57'//ext
      elseif (istep.eq.58) then
        filen = fname(1:l4)//'_'//'58'//ext
      elseif (istep.eq.59) then
        filen = fname(1:l4)//'_'//'59'//ext
      elseif (istep.eq.60) then
        filen = fname(1:l4)//'_'//'60'//ext
      elseif (istep.eq.61) then
        filen = fname(1:l4)//'_'//'61'//ext
      elseif (istep.eq.62) then
        filen = fname(1:l4)//'_'//'62'//ext
      elseif (istep.eq.63) then
        filen = fname(1:l4)//'_'//'63'//ext
      elseif (istep.eq.64) then
        filen = fname(1:l4)//'_'//'64'//ext
      elseif (istep.eq.65) then
        filen = fname(1:l4)//'_'//'65'//ext
      elseif (istep.eq.66) then
        filen = fname(1:l4)//'_'//'66'//ext
      elseif (istep.eq.67) then
        filen = fname(1:l4)//'_'//'67'//ext
      elseif (istep.eq.68) then
        filen = fname(1:l4)//'_'//'68'//ext
      elseif (istep.eq.69) then
        filen = fname(1:l4)//'_'//'69'//ext
      elseif (istep.eq.70) then
        filen = fname(1:l4)//'_'//'70'//ext
      elseif (istep.eq.71) then
        filen = fname(1:l4)//'_'//'71'//ext
      elseif (istep.eq.72) then
        filen = fname(1:l4)//'_'//'72'//ext
      elseif (istep.eq.73) then
        filen = fname(1:l4)//'_'//'73'//ext
      elseif (istep.eq.74) then
        filen = fname(1:l4)//'_'//'74'//ext
      elseif (istep.eq.75) then
        filen = fname(1:l4)//'_'//'75'//ext
      elseif (istep.eq.76) then
        filen = fname(1:l4)//'_'//'76'//ext
      elseif (istep.eq.77) then
        filen = fname(1:l4)//'_'//'77'//ext
      elseif (istep.eq.78) then
        filen = fname(1:l4)//'_'//'78'//ext
      elseif (istep.eq.79) then
        filen = fname(1:l4)//'_'//'79'//ext
      elseif (istep.eq.80) then
        filen = fname(1:l4)//'_'//'80'//ext
      elseif (istep.eq.81) then
        filen = fname(1:l4)//'_'//'81'//ext
      elseif (istep.eq.82) then
        filen = fname(1:l4)//'_'//'82'//ext
      elseif (istep.eq.83) then
        filen = fname(1:l4)//'_'//'83'//ext
      elseif (istep.eq.84) then
        filen = fname(1:l4)//'_'//'84'//ext
      elseif (istep.eq.85) then
        filen = fname(1:l4)//'_'//'85'//ext
      elseif (istep.eq.86) then
        filen = fname(1:l4)//'_'//'86'//ext
      elseif (istep.eq.87) then
        filen = fname(1:l4)//'_'//'87'//ext
      elseif (istep.eq.88) then
        filen = fname(1:l4)//'_'//'88'//ext
      elseif (istep.eq.89) then
        filen = fname(1:l4)//'_'//'89'//ext
      elseif (istep.eq.90) then
        filen = fname(1:l4)//'_'//'90'//ext
      elseif (istep.eq.91) then
        filen = fname(1:l4)//'_'//'91'//ext
      elseif (istep.eq.92) then
        filen = fname(1:l4)//'_'//'92'//ext
      elseif (istep.eq.93) then
        filen = fname(1:l4)//'_'//'93'//ext
      elseif (istep.eq.94) then
        filen = fname(1:l4)//'_'//'94'//ext
      elseif (istep.eq.95) then
        filen = fname(1:l4)//'_'//'95'//ext
      elseif (istep.eq.96) then
        filen = fname(1:l4)//'_'//'96'//ext
      elseif (istep.eq.97) then
        filen = fname(1:l4)//'_'//'97'//ext
      elseif (istep.eq.98) then
        filen = fname(1:l4)//'_'//'98'//ext
      elseif (istep.eq.99) then
        filen = fname(1:l4)//'_'//'99'//ext
      elseif (istep.eq.100) then
        filen = fname(1:l4)//'_'//'100'//ext
      elseif (istep.eq.101) then
        filen = fname(1:l4)//'_'//'101'//ext
      elseif (istep.eq.102) then
        filen = fname(1:l4)//'_'//'102'//ext
      elseif (istep.eq.103) then
        filen = fname(1:l4)//'_'//'103'//ext
      elseif (istep.eq.104) then
        filen = fname(1:l4)//'_'//'104'//ext
      elseif (istep.eq.105) then
        filen = fname(1:l4)//'_'//'105'//ext
      elseif (istep.eq.106) then
        filen = fname(1:l4)//'_'//'106'//ext
      elseif (istep.eq.107) then
        filen = fname(1:l4)//'_'//'107'//ext
      elseif (istep.eq.108) then
        filen = fname(1:l4)//'_'//'108'//ext
      elseif (istep.eq.109) then
        filen = fname(1:l4)//'_'//'109'//ext
      elseif (istep.eq.110) then
        filen = fname(1:l4)//'_'//'110'//ext
      elseif (istep.eq.111) then
        filen = fname(1:l4)//'_'//'111'//ext
      elseif (istep.eq.112) then
        filen = fname(1:l4)//'_'//'112'//ext
      elseif (istep.eq.113) then
        filen = fname(1:l4)//'_'//'113'//ext
      elseif (istep.eq.114) then
        filen = fname(1:l4)//'_'//'114'//ext
      elseif (istep.eq.115) then
        filen = fname(1:l4)//'_'//'115'//ext
      elseif (istep.eq.116) then
        filen = fname(1:l4)//'_'//'116'//ext
      elseif (istep.eq.117) then
        filen = fname(1:l4)//'_'//'117'//ext
      elseif (istep.eq.118) then
        filen = fname(1:l4)//'_'//'118'//ext
      elseif (istep.eq.119) then
        filen = fname(1:l4)//'_'//'119'//ext
      elseif (istep.eq.120) then
        filen = fname(1:l4)//'_'//'120'//ext
      elseif (istep.eq.121) then
        filen = fname(1:l4)//'_'//'121'//ext
      elseif (istep.eq.122) then
        filen = fname(1:l4)//'_'//'122'//ext
      elseif (istep.eq.123) then
        filen = fname(1:l4)//'_'//'123'//ext
      elseif (istep.eq.124) then
        filen = fname(1:l4)//'_'//'124'//ext
      elseif (istep.eq.125) then
        filen = fname(1:l4)//'_'//'125'//ext
      elseif (istep.eq.126) then
        filen = fname(1:l4)//'_'//'126'//ext
      elseif (istep.eq.127) then
        filen = fname(1:l4)//'_'//'127'//ext
      elseif (istep.eq.128) then
        filen = fname(1:l4)//'_'//'128'//ext
      elseif (istep.eq.129) then
        filen = fname(1:l4)//'_'//'129'//ext
      elseif (istep.eq.130) then
        filen = fname(1:l4)//'_'//'130'//ext
      elseif (istep.eq.131) then
        filen = fname(1:l4)//'_'//'131'//ext
      elseif (istep.eq.132) then
        filen = fname(1:l4)//'_'//'132'//ext
      elseif (istep.eq.133) then
        filen = fname(1:l4)//'_'//'133'//ext
      elseif (istep.eq.134) then
        filen = fname(1:l4)//'_'//'134'//ext
      elseif (istep.eq.135) then
        filen = fname(1:l4)//'_'//'135'//ext
      elseif (istep.eq.136) then
        filen = fname(1:l4)//'_'//'136'//ext
      elseif (istep.eq.137) then
        filen = fname(1:l4)//'_'//'137'//ext
      elseif (istep.eq.138) then
        filen = fname(1:l4)//'_'//'138'//ext
      elseif (istep.eq.139) then
        filen = fname(1:l4)//'_'//'139'//ext
      elseif (istep.eq.140) then
        filen = fname(1:l4)//'_'//'140'//ext
      elseif (istep.eq.141) then
        filen = fname(1:l4)//'_'//'141'//ext
      elseif (istep.eq.142) then
        filen = fname(1:l4)//'_'//'142'//ext
      elseif (istep.eq.143) then
        filen = fname(1:l4)//'_'//'143'//ext
      elseif (istep.eq.144) then
        filen = fname(1:l4)//'_'//'144'//ext
      elseif (istep.eq.145) then
        filen = fname(1:l4)//'_'//'145'//ext
      elseif (istep.eq.146) then
        filen = fname(1:l4)//'_'//'146'//ext
      elseif (istep.eq.147) then
        filen = fname(1:l4)//'_'//'147'//ext
      elseif (istep.eq.148) then
        filen = fname(1:l4)//'_'//'148'//ext
      elseif (istep.eq.149) then
        filen = fname(1:l4)//'_'//'149'//ext
      elseif (istep.eq.150) then
        filen = fname(1:l4)//'_'//'150'//ext
      elseif (istep.eq.151) then
        filen = fname(1:l4)//'_'//'151'//ext
      elseif (istep.eq.152) then
        filen = fname(1:l4)//'_'//'152'//ext
      elseif (istep.eq.153) then
        filen = fname(1:l4)//'_'//'153'//ext
      elseif (istep.eq.154) then
        filen = fname(1:l4)//'_'//'154'//ext
      elseif (istep.eq.155) then
        filen = fname(1:l4)//'_'//'155'//ext
      elseif (istep.eq.156) then
        filen = fname(1:l4)//'_'//'156'//ext
      elseif (istep.eq.157) then
        filen = fname(1:l4)//'_'//'157'//ext
      elseif (istep.eq.158) then
        filen = fname(1:l4)//'_'//'158'//ext
      elseif (istep.eq.159) then
        filen = fname(1:l4)//'_'//'159'//ext
      elseif (istep.eq.160) then
        filen = fname(1:l4)//'_'//'160'//ext
      elseif (istep.eq.161) then
        filen = fname(1:l4)//'_'//'161'//ext
      elseif (istep.eq.162) then
        filen = fname(1:l4)//'_'//'162'//ext
      elseif (istep.eq.163) then
        filen = fname(1:l4)//'_'//'163'//ext
      elseif (istep.eq.164) then
        filen = fname(1:l4)//'_'//'164'//ext
      elseif (istep.eq.165) then
        filen = fname(1:l4)//'_'//'165'//ext
      elseif (istep.eq.166) then
        filen = fname(1:l4)//'_'//'166'//ext
      elseif (istep.eq.167) then
        filen = fname(1:l4)//'_'//'167'//ext
      elseif (istep.eq.168) then
        filen = fname(1:l4)//'_'//'168'//ext
      elseif (istep.eq.169) then
        filen = fname(1:l4)//'_'//'169'//ext
      elseif (istep.eq.170) then
        filen = fname(1:l4)//'_'//'170'//ext
      elseif (istep.eq.171) then
        filen = fname(1:l4)//'_'//'171'//ext
      elseif (istep.eq.172) then
        filen = fname(1:l4)//'_'//'172'//ext
      elseif (istep.eq.173) then
        filen = fname(1:l4)//'_'//'173'//ext
      elseif (istep.eq.174) then
        filen = fname(1:l4)//'_'//'174'//ext
      elseif (istep.eq.175) then
        filen = fname(1:l4)//'_'//'175'//ext
      elseif (istep.eq.176) then
        filen = fname(1:l4)//'_'//'176'//ext
      elseif (istep.eq.177) then
        filen = fname(1:l4)//'_'//'177'//ext
      elseif (istep.eq.178) then
        filen = fname(1:l4)//'_'//'178'//ext
      elseif (istep.eq.179) then
        filen = fname(1:l4)//'_'//'179'//ext
      elseif (istep.eq.180) then
        filen = fname(1:l4)//'_'//'180'//ext
      elseif (istep.eq.181) then
        filen = fname(1:l4)//'_'//'181'//ext
      elseif (istep.eq.182) then
        filen = fname(1:l4)//'_'//'182'//ext
      elseif (istep.eq.183) then
        filen = fname(1:l4)//'_'//'183'//ext
      elseif (istep.eq.184) then
        filen = fname(1:l4)//'_'//'184'//ext
      elseif (istep.eq.185) then
        filen = fname(1:l4)//'_'//'185'//ext
      elseif (istep.eq.186) then
        filen = fname(1:l4)//'_'//'186'//ext
      elseif (istep.eq.187) then
        filen = fname(1:l4)//'_'//'187'//ext
      elseif (istep.eq.188) then
        filen = fname(1:l4)//'_'//'188'//ext
      elseif (istep.eq.189) then
        filen = fname(1:l4)//'_'//'189'//ext
      elseif (istep.eq.190) then
        filen = fname(1:l4)//'_'//'190'//ext
      elseif (istep.eq.191) then
        filen = fname(1:l4)//'_'//'191'//ext
      elseif (istep.eq.192) then
        filen = fname(1:l4)//'_'//'192'//ext
      elseif (istep.eq.193) then
        filen = fname(1:l4)//'_'//'193'//ext
      elseif (istep.eq.194) then
        filen = fname(1:l4)//'_'//'194'//ext
      elseif (istep.eq.195) then
        filen = fname(1:l4)//'_'//'195'//ext
      elseif (istep.eq.196) then
        filen = fname(1:l4)//'_'//'196'//ext
      elseif (istep.eq.197) then
        filen = fname(1:l4)//'_'//'197'//ext
      elseif (istep.eq.198) then
        filen = fname(1:l4)//'_'//'198'//ext
      elseif (istep.eq.199) then
        filen = fname(1:l4)//'_'//'199'//ext
      elseif (istep.eq.200) then
        filen = fname(1:l4)//'_'//'200'//ext
      else
       write(6,*) 'istep gt 200 - stop *****  istep=',istep
       stop
      endif
c
      open(nor,file=filen,status='unknown')
c
      end
      subroutine opresi(fname,ext,len,nor)
c
      integer istep,len,l4,nor
      character*20 fname,filen
      character*4  ext
      character*4  cstep4
      character*3  cstep3
      character*2  cstep2
      character*1  cstep1
c
      l4 = len-4
      filen = fname(1:l4)//ext
c
      open(nor,file=filen,status='unknown')
c
      end
      subroutine othfem(bridge,wopos)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      dimension     bridge(*)
      character*5   wopos
c
      read(nin) istep,ttime,(bridge(itotv),itotv=1,npoin)
c
      write(nog,10) 100,'C','Step  ',ttime,1,istep
      write(nog,20) -4,wopos,1,1,0
      write(nog,30) -5,wopos,1,1
      do ipoin=1,npoin
        write(nog,40) -1,ipoin,bridge(ipoin)
      enddo
      write(nog,50) -3
c
   10 format(1x,i4,a1,a6,e12.5,32x,i2,i5)
   20 format(1x,i2,2x,a5,3x,3i5)
   30 format(1x,i2,2x,a5,3x,2i5)
   40 format(1x,i2,i5,e12.5)
   50 format(1x,i2)
c
      end
      subroutine othfla(bridge,wopos)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      dimension     bridge(*)
      character*5   wopos
c
      write(nor,10) wopos
c
      do ipoin = 1,npoin
        write(nor,20) ipoin,bridge(ipoin)
      enddo
c
   10 format(a5)
   20 format(1x,i5,2x,e12.5)
c
      end
      subroutine othgid(bridge,wopos,ttime)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      dimension     bridge(*)
      character*5   wopos
c
c...  write pressure
c
      if(wopos.eq.'PRESS') then
        write(nor,10) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c...  write vorticity
c
      if(wopos.eq.'VORTI') then
        write(nor,30) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c...  write streamlines
c
      if(wopos.eq.'STREA') then
        write(nor,50) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c...  write viscosities
c
      if(wopos.eq.'VISCO') then
        write(nor,70) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c...  write densities
c
      if(wopos.eq.'DENSI') then
        write(nor,80) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c...  write material front (pseudo-concentration)
c
      if(wopos.eq.'MATER') then
        write(nor,90) ttime
        do ipoin = 1,npoin
          write(nor,40) ipoin,bridge(ipoin)
        enddo
      endif
c
c  10 format('Pressure        1      ',e10.3,' 1 1 0')
c  30 format('Vorticities     1      ',e10.3,' 1 1 0')
c  40 format(1x,i5,2x,e12.5)
c  50 format('Streamlines     1      ',e10.3,' 1 1 0')
c  70 format('Viscosities     1      ',e10.3,' 1 1 0')
c  80 format('Densities       1      ',e10.3,' 1 1 0')
c  90 format('Material front  1      ',e10.3,' 1 1 0')
   10 format('Pressure        1      ',e10.5,' 1 1 0')
   30 format('Vorticities     1      ',e10.5,' 1 1 0')
   40 format(1x,i5,2x,e12.5)
   50 format('Streamlines     1      ',e10.5,' 1 1 0')
   70 format('Viscosities     1      ',e10.5,' 1 1 0')
   80 format('Densities       1      ',e10.5,' 1 1 0')
   90 format('Material front  1      ',e10.5,' 1 1 0')
c
      end
      subroutine resfem(bridge)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*5   wopos
      dimension     bridge(*)
c
   10 read(nin,end=20) wopos
      if(wopos.eq.'UNKNO') then
        call unkfem(bridge)
      else
        call othfem(bridge,wopos)
      end if
      go to 10
   20 write(nog,'(i5)') 9999
c
      end
      subroutine resfla(bridge,filen,ext,length)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      common/precon/kprec 
      character*20  filen
      character*5   wopos
      character*4   ext
      dimension     bridge(*)
c
      ntotv  = npoin*ndime+npoin*ktemp
      ndite  = ndime+ktemp
      istold =-1
c
c...  read nodal results
c
   10 read(nin,end=20) wopos
c
      if(wopos.eq.'UNKNO') then
	if(kprec.eq.1) then
        read(nin) istep,ttime
	do ipoin=1,npoin
	do idime=1,ndime
	itotv=ndime*(ipoin-1)+idime 
        read(nin) bridge(itotv)
	end do 
	end do
	if(ktemp.eq.1) then
        ktotv=ndime*npoin+1
        read(nin) (bridge(itotv),itotv=ktotv,ntotv)
	end if 
	else 
        read(nin) istep,ttime, (bridge(itotv),itotv=1,ntotv)
	end if 
      else
        read(nin) istep,ttime, (bridge(itotv),itotv=1,npoin)
      endif
c
c...  open results file and write velocities and mach numbers
c
      if(istep.ne.istold) then
        call opresu(filen,ext,length,istep,nor)
        if(wopos.ne.'UNKNO') then
          write(nor,'(a10)') 'VELOCITIES'
          do i = 1,npoin
            if(ndime.eq.2) then
              write(nor,'(1x,i5,2(2x,e12.5))') i,0.,0.
            else
              write(nor,'(1x,i5,3(2x,e12.5))') i,0.,0.,0.
            endif
          enddo
          write(nor,'(a11)') 'MACH NUMBER'
          do i = 1,npoin
            write(nor,'(1x,i5,2x,e12.5)') i,0.
          enddo
        endif
        istold = istep
      endif
c
c...  write nodal results
c
      if(wopos.eq.'UNKNO') then
        call unkfla(bridge)
      else
        call othfla(bridge,wopos)
      endif
      go to 10
   20 continue
c
      end
      subroutine respat(bridge,wopos,filen,ext,length)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*20  filen
      character*5   wopos(*)
      character*4   ext
      dimension     bridge(*)
c
      ntotv  = npoin*ndime+npoin*ktemp
      ndite  = ndime+ktemp
      istold =-1
      iend   = 0
  110 continue
      icoun = 0
  120 continue
      ipos = icoun*npoin
      icoun = icoun+1
      read(nin,end=130) wopos(icoun)
      read(nin,end=130) istep
      backspace(nin)
      if(icoun.eq.1) istold = istep
c
c...  read nodal results
c
      if(istep.eq.istold) then
        if(wopos(icoun).eq.'UNKNO') then
          read(nin) istep,ttime,(bridge(ipos+itotv),itotv=1,ntotv)
          icoun = icoun+ndite-1
        else
          read(nin) istep,ttime,(bridge(ipos+itotv),itotv=1,npoin)
        endif
        goto 120
      endif
      backspace(nin)
      goto 140
  130 continue
      iend = 1
  140 continue
      icoun = icoun-1
c
c...  open results file
c
      call opresu(filen,ext,length,istold,nor)
c
c...  write heading information of results file
c
      write(nor,10) 'PATRAN ','RESULTS'
      write(nor,20) npoin,0,0,0,icoun
      if(wopos(1).eq.'UNKNO') then
        if(ndime.eq.2.and.ktemp.eq.0) then
          write(nor,10) ' VELO_X',' VELO_Y'
        elseif(ndime.eq.3.and.ktemp.eq.0) then
          write(nor,10) ' VELO_X',' VELO_Y',' VELO_Z'
        elseif(ndime.eq.2.and.ktemp.eq.1) then
          write(nor,10) ' VELO_X',' VELO_Y',' TEMPER'
        elseif(ndime.eq.3.and.ktemp.eq.1) then
          write(nor,10) ' VELO_X',' VELO_Y',' VELO_Z',' TEMPER'
        endif
        if(icoun.gt.ndite) then
          write(nor,10) (wopos(i),i=ndite+1,icoun)
        else
          write(nor,10) ' '
        endif
      else
          write(nor,10) ' '
          write(nor,10) (wopos(i),i=1,icoun)
      endif
c
c...  write nodal results
c
      if(wopos(1).eq.'UNKNO') then
        do ipoin = 1,npoin
          iposi = ndime*(ipoin-1)+1
          write(nor,30) ipoin,(bridge(iposi+idime),idime=1,ndime),
     .                  (bridge(icolu*npoin+ipoin),icolu=ndime,icoun-1)
        enddo
      else
        do ipoin = 1,npoin
          write(nor,30) ipoin,(bridge(icolu*npoin+ipoin),
     .                                icolu=0,icoun-1)
        enddo
      endif
c
      if(iend.eq.0) goto 110
   10 format(7a7)
   20 format(2i9,e15.6,2i9)
c  30 format(i8,<icoun>(e13.7))
   30 format(i8,20(e13.7))
c
      end
      subroutine resgid(bridge,filen,ext,length)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      common/precon/kprec
      character*20  filen
      character*5   wopos
      character*4   ext
      dimension     bridge(*)
c
      ntotv  = npoin*ndime+npoin*ktemp
      ndite  = ndime+ktemp
      istold =-1
c
c...  open results file
c
      call opresi(filen,ext,length,nor)
c
c...  read nodal results
c
   10 read(nin,end=20) wopos
c
      if(wopos.eq.'UNKNO') then
        if(kprec.eq.1) then
        read(nin) istep,ttime
        do ipoin=1,npoin
        do idime=1,ndime
        itotv=ndime*(ipoin-1)+idime
        read(nin) bridge(itotv)
        end do
        end do
        if(ktemp.eq.1) then
        ktotv=ndime*npoin+1
        read(nin) (bridge(itotv),itotv=ktotv,ntotv)
        end if
        else
        read(nin) istep,ttime, (bridge(itotv),itotv=1,ntotv)
        end if
      else
        read(nin) istep,ttime, (bridge(itotv),itotv=1,npoin)
      endif
c
c...  write nodal results
c
      if(wopos.eq.'UNKNO') then
        call unkgid(bridge,ttime)
      else
        call othgid(bridge,wopos,ttime)
      endif
      go to 10
   20 continue
c
      end
      subroutine result(bridge,wopos,filen,ext,length)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      character*20  filen
      character*5   wopos(*)
      character*4   ext
      dimension     bridge(*)
c
      if(i_pos.eq.1) then
        call resfem(bridge)
      elseif(i_pos.eq.2) then
        call resfla(bridge,filen,ext,length)
      elseif(i_pos.eq.3) then
        call respat(bridge,wopos,filen,ext,length)
      elseif(i_pos.eq.4) then
        call resgid(bridge,filen,ext,length)
      endif
c
      end
      subroutine top004(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 3
      melem = 3*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(3,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(4,ielem)
        lnodw(3,kelem) = lnods(3,ielem)
      end do
c
      end
      subroutine top007(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 4
      melem = 3*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(4,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(6,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(4,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(7,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(6,ielem)
        lnodw(2,kelem) = lnods(7,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(3,ielem)
      end do
c
      end
      subroutine top011(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 4
      melem = 12*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(8,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(6,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(3,ielem)
        lnodw(2,kelem) = lnods(7,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(10,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(4,ielem)
        lnodw(2,kelem) = lnods(10,ielem)
        lnodw(3,kelem) = lnods(9,ielem)
        lnodw(4,kelem) = lnods(8,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(8,ielem)
        lnodw(2,kelem) = lnods(9,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(10,ielem)
        lnodw(3,kelem) = lnods(8,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(6,ielem)
        lnodw(2,kelem) = lnods(9,ielem)
        lnodw(3,kelem) = lnods(10,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(8,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(5,ielem)
        lnodw(2,kelem) = lnods(9,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(6,ielem)
        lnodw(2,kelem) = lnods(10,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(8,ielem)
        lnodw(2,kelem) = lnods(10,ielem)
        lnodw(3,kelem) = lnods(9,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
      end do
c
      end
      subroutine top015(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 8
      melem = 4*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(11,ielem)
        lnodw(4,kelem) = lnods(7,ielem)
        lnodw(5,kelem) = lnods(8,ielem)
        lnodw(6,kelem) = lnods(12,ielem)
        lnodw(7,kelem) = lnods(15,ielem)
        lnodw(8,kelem) = lnods(14,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(5,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        lnodw(5,kelem) = lnods(12,ielem)
        lnodw(6,kelem) = lnods(9,ielem)
        lnodw(7,kelem) = lnods(13,ielem)
        lnodw(8,kelem) = lnods(15,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(11,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(3,ielem)
        lnodw(5,kelem) = lnods(14,ielem)
        lnodw(6,kelem) = lnods(15,ielem)
        lnodw(7,kelem) = lnods(13,ielem)
        lnodw(8,kelem) = lnods(10,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(8,ielem)
        lnodw(2,kelem) = lnods(12,ielem)
        lnodw(3,kelem) = lnods(15,ielem)
        lnodw(4,kelem) = lnods(14,ielem)
        lnodw(5,kelem) = lnods(4,ielem)
        lnodw(6,kelem) = lnods(9,ielem)
        lnodw(7,kelem) = lnods(13,ielem)
        lnodw(8,kelem) = lnods(10,ielem)
      end do
c
      end
      subroutine top027(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 8
      melem = 8*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(9,ielem)
        lnodw(3,kelem) = lnods(21,ielem)
        lnodw(4,kelem) = lnods(12,ielem)
        lnodw(5,kelem) = lnods(13,ielem)
        lnodw(6,kelem) = lnods(22,ielem)
        lnodw(7,kelem) = lnods(27,ielem)
        lnodw(8,kelem) = lnods(25,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(9,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(10,ielem)
        lnodw(4,kelem) = lnods(21,ielem)
        lnodw(5,kelem) = lnods(22,ielem)
        lnodw(6,kelem) = lnods(14,ielem)
        lnodw(7,kelem) = lnods(23,ielem)
        lnodw(8,kelem) = lnods(27,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(21,ielem)
        lnodw(2,kelem) = lnods(10,ielem)
        lnodw(3,kelem) = lnods(3,ielem)
        lnodw(4,kelem) = lnods(11,ielem)
        lnodw(5,kelem) = lnods(27,ielem)
        lnodw(6,kelem) = lnods(23,ielem)
        lnodw(7,kelem) = lnods(15,ielem)
        lnodw(8,kelem) = lnods(24,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(12,ielem)
        lnodw(2,kelem) = lnods(21,ielem)
        lnodw(3,kelem) = lnods(11,ielem)
        lnodw(4,kelem) = lnods(4,ielem)
        lnodw(5,kelem) = lnods(25,ielem)
        lnodw(6,kelem) = lnods(27,ielem)
        lnodw(7,kelem) = lnods(24,ielem)
        lnodw(8,kelem) = lnods(16,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(13,ielem)
        lnodw(2,kelem) = lnods(22,ielem)
        lnodw(3,kelem) = lnods(27,ielem)
        lnodw(4,kelem) = lnods(25,ielem)
        lnodw(5,kelem) = lnods(5,ielem)
        lnodw(6,kelem) = lnods(17,ielem)
        lnodw(7,kelem) = lnods(26,ielem)
        lnodw(8,kelem) = lnods(20,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(22,ielem)
        lnodw(2,kelem) = lnods(14,ielem)
        lnodw(3,kelem) = lnods(23,ielem)
        lnodw(4,kelem) = lnods(27,ielem)
        lnodw(5,kelem) = lnods(17,ielem)
        lnodw(6,kelem) = lnods(6,ielem)
        lnodw(7,kelem) = lnods(18,ielem)
        lnodw(8,kelem) = lnods(26,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(27,ielem)
        lnodw(2,kelem) = lnods(23,ielem)
        lnodw(3,kelem) = lnods(15,ielem)
        lnodw(4,kelem) = lnods(24,ielem)
        lnodw(5,kelem) = lnods(26,ielem)
        lnodw(6,kelem) = lnods(18,ielem)
        lnodw(7,kelem) = lnods(7,ielem)
        lnodw(8,kelem) = lnods(19,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(25,ielem)
        lnodw(2,kelem) = lnods(27,ielem)
        lnodw(3,kelem) = lnods(24,ielem)
        lnodw(4,kelem) = lnods(16,ielem)
        lnodw(5,kelem) = lnods(20,ielem)
        lnodw(6,kelem) = lnods(26,ielem)
        lnodw(7,kelem) = lnods(19,ielem)
        lnodw(8,kelem) = lnods(8,ielem)
      end do
c
      end
      subroutine top05a(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 3
      melem = 4*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(3,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(3,ielem)
        lnodw(2,kelem) = lnods(4,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
      end do
c
      end
      subroutine top05b(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 4
      melem = 4*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(4,ielem)
        lnodw(3,kelem) = lnods(2,ielem)
        lnodw(4,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(4,ielem)
        lnodw(3,kelem) = lnods(3,ielem)
        lnodw(4,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(3,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        lnodw(4,kelem) = lnods(5,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(3,ielem)
        lnodw(4,kelem) = lnods(5,ielem)
      end do
c
      end
      subroutine top09a(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 4
      melem = 4*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(9,ielem)
        lnodw(4,kelem) = lnods(8,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(5,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(6,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(9,ielem)
        lnodw(2,kelem) = lnods(6,ielem)
        lnodw(3,kelem) = lnods(3,ielem)
        lnodw(4,kelem) = lnods(7,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(8,ielem)
        lnodw(2,kelem) = lnods(9,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(4,ielem)
      end do
c
      end
      subroutine top09b(lnods,nelem,nnode,lnodw,melem,nnodw)
c
      dimension lnods(nnode,nelem), lnodw(12,*)
c
      nnodw = 4
      melem = 12*nelem
      kelem = 0
      do ielem = 1,nelem
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(2,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(5,ielem)
        lnodw(2,kelem) = lnods(6,ielem)
        lnodw(3,kelem) = lnods(2,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(6,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(3,ielem)
        lnodw(3,kelem) = lnods(2,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(3,ielem)
        lnodw(2,kelem) = lnods(7,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(8,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(4,ielem)
        lnodw(2,kelem) = lnods(8,ielem)
        lnodw(3,kelem) = lnods(5,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(4,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(1,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(6,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(7,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(7,ielem)
        lnodw(2,kelem) = lnods(5,ielem)
        lnodw(3,kelem) = lnods(8,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(1,ielem)
        lnodw(2,kelem) = lnods(2,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
        kelem = kelem+1
        lnodw(1,kelem) = lnods(2,ielem)
        lnodw(2,kelem) = lnods(3,ielem)
        lnodw(3,kelem) = lnods(4,ielem)
        lnodw(4,kelem) = lnods(9,ielem)
      end do
c
      end
      subroutine unkfem(bridge)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      common/precon/kprec 
      dimension     bridge(*)
c
      ntotv = npoin*ndime+npoin*ktemp
	if(kprec.eq.1) then
        read(nin) istep,ttime
	do ipoin=1,npoin
	do idime=1,ndime
	itotv=ndime*(ipoin-1)+idime 
        read(nin) bridge(itotv)
	end do 
	end do
	if(ktemp.eq.1) then
        ktotv=ndime*npoin+1
        read(nin) (bridge(itotv),itotv=ktotv,ntotv)
	end if 
	else 
        read(nin) istep,ttime, (bridge(itotv),itotv=1,ntotv)
	end if 
c
      write(nog,10) 100,'C','Step  ',ttime,1,istep
      write(nog,20) -4,'VELOCITY',4,1,0
      write(nog,30) -5,'X-COMP  ',1,2,1,0,0
      write(nog,30) -5,'Y-COMP  ',1,2,2,0,0
      write(nog,30) -5,'Z-COMP  ',1,2,3,0,0
      write(nog,80) -5,'ALL     ',1,2,0,0,1,'all     '
      veloz = 0.0d0
      do ipoin = 1,npoin
        iposi = ndime*(ipoin-1)
        if(ndime.eq.2) then
          write(nog,50) -1,ipoin,(bridge(iposi+idime),idime=1,ndime),
     .                            veloz
        else
          write(nog,50) -1,ipoin,(bridge(iposi+idime),idime=1,ndime)
        endif
      enddo
      write(nog,70) -3
      if(ktemp.eq.1) then
        write(nog,10) 100,'C','Step  ',ttime,1,istep
        write(nog,20) -4,'TEMPE   ',1,1,0
        write(nog,40) -5,'TEMPE   ',1,1
        do ipoin = 1,npoin
          itotv = ndime*npoin+ipoin
          write(nog,60) -1,ipoin,bridge(itotv)
        enddo
        write(nog,70) -3
      end if
c
   10 format(1x,i4,a1,a6,e12.5,32x,i2,i5)
   20 format(1x,i2,2x,a8,3i5)
   30 format(1x,i2,2x,a8,5i5)
   40 format(1x,i2,2x,a8,2i5)
   50 format(1x,i2,i5,3e12.5)
   60 format(1x,i2,i5,e12.5)
   70 format(1x,i2)
   80 format(1x,i2,2x,a8,5i5,a8)
c
      end
      subroutine unkfla(bridge)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      dimension     bridge(*)
c
c...  write velocities
c
      write(nor,10)
      do ipoin = 1,npoin
        iposi = ndime*(ipoin-1)
        write(nor,20) ipoin,(bridge(iposi+idime),idime=1,ndime)
      enddo
c
c...  write mach numbers (=0)
c
      write(nor,30)
      do ipoin = 1,npoin
        write(nor,40) ipoin,0.
      enddo
c
c...  write temperatures if present
c
      if(ktemp.eq.1) then
        write(nor,50)
        do ipoin = 1,npoin
          itotv = ndime*npoin+ipoin
          write(nor,40) ipoin,bridge(itotv)
        enddo
      end if
c
   10 format('VELOCITIES')
   20 format(1x,i5,3(2x,e12.5))
   30 format('MACH NUMBER')
   40 format(1x,i5,2x,e12.5)
   50 format('TEMPERATURE')
c
      end
      subroutine unkgid(bridge,ttime)
c
      implicit      real*4 (a-h,o-z)
c     implicit      real*8 (a-h,o-z)
      common/contro/ndime,nnode,nelty,nelem,npoin,ktemp,nin,nog,nor,
     .              i_pos
      dimension     bridge(*)
c
c...  write velocities
c
      write(nor,10) ttime
      do ipoin = 1,npoin
        iposi = ndime*(ipoin-1)
        write(nor,20) ipoin,(bridge(iposi+idime),idime=1,ndime)
      enddo
c
c...  write mach numbers (=0)
c
c     write(nor,30) ttime
c     do ipoin = 1,npoin
c       write(nor,40) ipoin,0.
c     enddo
c
c...  write temperatures if present
c
      if(ktemp.eq.1) then
        write(nor,50) ttime
        do ipoin = 1,npoin
          itotv = ndime*npoin+ipoin
          write(nor,40) ipoin,bridge(itotv)
        enddo
      end if
c
c  10 format('Velocities      1      ',e10.3,' 2 1 0')
c  20 format(1x,i5,3(2x,e12.5))
c  30 format('Mach number     1      ',e10.3,' 1 1 0')
c  40 format(1x,i5,2x,e12.5)
c  50 format('Temperature     1      ',e10.3,' 1 1 0')
   10 format('Velocities      1      ',e10.5,' 2 1 0')
   20 format(1x,i5,3(2x,e12.5))
   30 format('Mach number     1      ',e10.5,' 1 1 0')
   40 format(1x,i5,2x,e12.5)
   50 format('Temperature     1      ',e10.5,' 1 1 0')
c
      end
