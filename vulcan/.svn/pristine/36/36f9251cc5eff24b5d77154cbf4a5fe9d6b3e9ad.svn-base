      program malla3d
c***********************************************************************
c
c     programa para hacer una malla 3D de hexaedros a partir de una 2D
c     de cuadrilateros.
c
c***********************************************************************
      implicit double precision(a-h,o-z)

      parameter(mpoin=100000, mdime=3, melem=100000, mnode=27,
     .          melemz=40)

      character*40 filein1,filein2,fileout1,fileout2,fileout3

c      dimension lnods(melem,mnode),imat(melem),
c     .          ln(melem,mnode,melemz+1),ipres(mpoin,mdime+2),
c     .          ifunc(mpoin)
      integer, allocatable :: lnods(:,:),imat(:), 
     .                        ln(:,:,:),ipres(:,:), ifunc(:)
c      dimension coord(mpoin,mdime),pres(mpoin,mdime+2)
      real(kind=8),allocatable :: coord(:,:),pres(:,:)

      twopi=6.283185307179586d0

      write(6,*)'problema mecanico o termico (=1) o fluido (=2)='
      read(5,*) itfm

      if(itfm.eq.2) then
       write(6,*)'hay frente material (1:no; 2:si)'
       read(5,*) ifront
      endif

      write(6,*)'Reading from VULCAN (1) or GiD (2) geometry file ?'
      read(5,*) igid

      if(igid.eq.1) then
       write(6,*)'nombre del archivo de la malla 2d (*.geo - VULCAN)='
       read(5,*) filein1
      endif
      if(igid.eq.2) then
       write(6,*)'nombre del archivo de la malla 2d (*.dat - GiD)='
       read(5,*) filein1
      endif

      write(6,*)'numero total de elementos (2D)='
      read(5,*) nelem

      if(nelem.gt.melem) then
       write(6,*) 'nelem gt melem',nelem,melem
       stop
      endif

      write(6,*)'numero de nodos por elemento (2D)='
      read(5,*) nnode

      inoco=1
      if(itfm.eq.1) then              ! to be extended for flow problems
       write(6,*)'hay elem. contacto con nodos NO coinc. (1:no; 2:si)='
       read(5,*) inoco
      endif
      if(inoco.eq.2) then
       write(6,*)'numero de material de contacto (menor numeracion)'
       write(6,*)'(warning: non-coincident elements must be at the end'
       write(6,*)'          of the element list in the geo file)'
       read(5,*) imatnc
      endif

      write(6,*)'numero total de nodos (2D)='
      read(5,*) npoin

      if(npoin.gt.mpoin) then
       write(6,*) 'npoin gt mpoin',npoin,mpoin
       stop
      endif

      ndime=3

      write(6,*)'Writing for VULCAN (1), GiD (2) or both (3) ?'
      read(5,*) igidw

      if(igidw.eq.1.or.igidw.eq.3) then
       write(6,*)'nombre del archivo de la malla 3D (*.geo - VULCAN)='
       read(5,*) fileout1
       write(6,*)'nombre del archivo de la malla 3D (*.fix - VULCAN)='
       read(5,*) fileout3
      endif
      if(igidw.eq.2.or.igidw.eq.3) then
       write(6,*)'NAME FILEOUT2 (.msh - GiD)='
       read(5,*) fileout2
      endif

      if(igid.eq.1.and.(igidw.eq.1.or.igidw.eq.3)) then
       write(6,*)'considerar cond. de contorno 2D en 3D (1:no; 2:si)='
       read(5,*) ifixi
       if(ifixi.eq.2) then
        write(6,*)'nombre del archivo de la malla 2D (*.fix)='
        read(5,*) filein2
       endif
      endif

      write(6,*)'(Warning: the 2D mesh is assumed in the z=0 plane)'
      write(6,*)' '
      write(6,*)'3D Cartesian (1) or cylindrical (2) mesh ?'
      read(5,*) icacy

      if(icacy.eq.1) then
       write(6,*)'espesor de la placa='
       read(5,*) espesor
       write(6,*)'numero de elementos en el espesor='
       read(5,*) nelemz
      endif
      if(icacy.eq.2) then
       write(6,*)'(Warning: the x & y coordinates are respectively'
       write(6,*)'          assumed as the radial and axial coord.)'
       write(6,*)' '
       write(6,*)'angulo circunferencial (en grados, positivo)='
       read(5,*) alfa1
       write(6,*)'numero de elementos en el arco de circunferencia='
       read(5,*) nelemz
      endif

      if(nelemz.gt.melemz) then
       write(6,*) 'nelemz gt melemz',nelemz,melemz
       stop
      endif
      allocate(lnods(nelem,nnode))
      allocate(imat(nelem))
      allocate(ln(nelem,nnode,nelemz+1))
      allocate(ipres(npoin,ndime+2))
      allocate(ifunc(npoin))
      allocate(coord(npoin,ndime))
      allocate(pres(npoin,ndime+2))
c     lee geometria

      open(85,file=filein1)
      if(ifixi.eq.2) open(89,file=filein2)
      if(igidw.eq.1.or.igidw.eq.3) then
       open(86,file=fileout1)
       open(88,file=fileout3)
      endif
      if(igidw.eq.2.or.igidw.eq.3) open(87,file=fileout2)

      if(igid.eq.1) then
       nnode8=8
       if(itfm.eq.2) read(85,*)      ! ELEMENTS
       do ielem=1,nelem
        if(nnode.le.nnode8) then
         if(itfm.eq.1)
     .    read(85,*) jelem,imat(jelem),
     .                  (lnods(jelem,inode), inode=1,nnode)
         if(itfm.eq.2)
     .    read(85,*) jelem,
     .                  (lnods(jelem,inode), inode=1,nnode)
        else
         if(itfm.eq.1)
     .    read(85,*) jelem,imat(jelem),
     .                  (lnods(jelem,inode), inode=1,nnode8)
         if(itfm.eq.2)
     .    read(85,*) jelem,
     .                  (lnods(jelem,inode), inode=1,nnode8)
         if(nnode.le.(2*nnode8)) then
          read(85,*)    (lnods(jelem,inode), inode=nnode8+1,nnode)
         else
          read(85,*)    (lnods(jelem,inode), inode=nnode8+1,2*nnode8)
          read(85,*)    (lnods(jelem,inode), inode=2*nnode8+1,nnode)
         endif
        endif
       enddo
       if(itfm.eq.2) read(85,*)      ! END_ELEMENTS
       if(itfm.eq.2) read(85,*)      ! COORDINATES
       do ipoin=1,npoin
        read(85,*) jpoin,(coord(jpoin,idime), idime=1,ndime)
        if(ifixi.eq.2) then
         if(itfm.eq.1)
     .    read(89,951) jpoin,(ipres(jpoin,idime),idime=1,ndime),
     .                    ifunc(jpoin),(pres(jpoin,idime),idime=1,ndime)
         if(itfm.eq.2) then
          if(ifront.eq.1)
     .     read(89,961) jpoin,(ipres(jpoin,idime),idime=1,ndime),
     .                         ipres(jpoin,ndime+2),
     .                               (pres(jpoin,idime),idime=1,ndime),
     .                                pres(jpoin,ndime+2)
          if(ifront.eq.2)
     .     read(89,962) jpoin,(ipres(jpoin,idime),idime=1,ndime),
     .                         ipres(jpoin,ndime+2),
     .                         ipres(jpoin,ndime+3),
     .                               (pres(jpoin,idime),idime=1,ndime),
     .                                pres(jpoin,ndime+2),
     .                                pres(jpoin,ndime+3)
         endif
        endif
       enddo
       if(itfm.eq.2) read(85,*)      ! END_COORDINATES
      endif
      if(igid.eq.2) then
       if(inoco.eq.2) then
        write(6,*) 'non-coincident mesh not implemented for GiD'
        stop
       endif

       read(85,*)            ! comentarios de GiD
       read(85,*) 

       do ipoin=1,npoin
        read(85,*) jpoin,(coord(jpoin,idime), idime=1,ndime)
       enddo

       read(85,*)            ! comentarios de GiD
       read(85,*)
       read(85,*)

       do ielem=1,nelem
        read(85,*) jelem,(lnods(jelem,inode), inode=1,nnode),imat(jelem)
       enddo
      endif
C
      nnode3d=2*nnode
      if(nnode3d.gt.mnode) then
       write(6,*) 'nnode3d gt mnode',nnode3d,mnode
       stop
      endif

      ndime3d=3
      if(ndime3d.gt.mdime) then
       write(6,*) 'ndime3d gt mdime',ndime3d,mdime
       stop
      endif

c     construye conectividades

      do i=1,nelemz+1
       do ielem=1,nelem
        nnodx=nnode
        if(inoco.eq.2) then
         if(imat(ielem).ge.imatnc) then
          nnodx=nnode/2                                  ! only quads
          ln(ielem,nnodx+1,i)=lnods(ielem,nnodx+1)       ! contact index
          ln(ielem,nnodx+2,i)=lnods(ielem,nnodx+2)       ! contact index
         endif
        endif
        do inode=1,nnodx
         ln(ielem,inode,i)=lnods(ielem,inode)+(i-1)*npoin
        enddo
       enddo
      enddo

      k=0
      do i=1,nelemz
       do ielem=1,nelem
        k=k+1
        if(itfm.eq.1) imat(k)=imat(ielem)
        nnodx=nnode
        nnode3dx=nnode3d
        if(inoco.eq.2) then
         if(imat(k).ge.imatnc) then
          nnodx=nnode/2                                  ! only quads
          nnode3dx=nnode3d/2                             ! only quads
          lnods(k,nnodx*2+1)=ln(ielem,nnodx+1,i)         ! contact index
          lnods(k,nnodx*2+2)=ln(ielem,nnodx+2,i)         ! contact index
         endif
        endif
        do inode=1,nnode3dx
         if(inode.le.nnodx) then
          lnods(k,inode)=ln(ielem,inode,i)
         else
          lnods(k,inode)=ln(ielem,inode-nnodx,i+1)
         endif
        enddo
       enddo
      enddo
      nelem3d=k                ! nelem3d should be equal to nelem*nelemz
      if(nelem3d.gt.melem) then
       write(6,*) 'nelem3d gt melem',nelem3d,melem
       stop
      endif

      if(inoco.eq.2) then  ! rotate connectivity for nc contact elements
       k=0
       do i=1,nelemz
        do ielem=1,nelem
         k=k+1
         if(imat(k).ge.imatnc) then
          iaux      =lnods(k,1)
          lnods(k,1)=lnods(k,2)
          lnods(k,2)=iaux
         endif
        enddo
       enddo
      endif

c     coordinates & boundary conditions

      k=0
      do i=1,nelemz+1
       do ipoin=1,npoin
        k=k+1
        do idime=1,ndime3d
         if(idime.le.ndime) then
          if(icacy.eq.1) then                         ! Cartesian mesh
           coord(k,idime)=coord(ipoin,idime)
          else                                        ! cylindrical mesh
           dangulo=(float(i)-1.0d0)*alfa1*twopi/360.0d0/float(nelemz)
           if(idime.eq.1) coord(k,1)=                 ! radial
     .                               coord(ipoin,1)*dcos(dangulo)
           if(idime.eq.2) coord(k,2)=coord(ipoin,2)   ! axial
          endif
          if(ifixi.eq.2) then
           ipres(k,idime)=ipres(ipoin,idime)
           pres(k,idime)=pres(ipoin,idime)
          endif
         else
          if(icacy.eq.1) then                         ! Cartesian mesh
           coord(k,idime)=(float(i)-1.0d0)*espesor/float(nelemz)
          else                                        ! cylindrical mesh
           dangulo=(float(i)-1.0)*alfa1*twopi/360.0d0/float(nelemz)
           coord(k,3)=coord(ipoin,1)*dsin(dangulo)    ! circunferential
          endif
          if(ifixi.eq.2) then
           if(itfm.eq.1) then                         ! mech. problems
            ipres(k,idime)=ipres(ipoin,idime-2)       ! assumed
            pres(k,idime)=pres(ipoin,idime-2)
           endif
           if(itfm.eq.2) then                         ! flow problems
            ipres(k,idime)=0                          ! assumed
            pres(k,idime)=0.0d0
            ipres(k,idime+1)=ipres(ipoin,idime+1)     ! pressure
            pres(k,idime+1)=pres(ipoin,idime+1)
            if(ifront.eq.2) then
             ipres(k,idime+2)=ipres(ipoin,idime+2)    ! phi
             pres(k,idime+2)=pres(ipoin,idime+2)
            endif
           endif
          endif
         endif
        enddo
        if(ifixi.eq.2) then
         if(itfm.eq.1) ifunc(k)=ifunc(ipoin)
        endif
       enddo
      enddo
      npoin3d=npoin*(nelemz+1)
      if(npoin3d.gt.mpoin) then
       write(6,*) 'npoin3d gt mpoin',npoin3d,mpoin
       stop
      endif

c     escribe .geo

      if(igidw.eq.1.or.igidw.eq.3) then
       nnode8=8
       if(itfm.eq.2) write(86,*) 'ELEMENTS'
       do ielem=1,nelem3d
        if(nnode3d.le.nnode8) then
         if(itfm.eq.1)
     .    write(86,900) ielem,imat(ielem),
     .                   (lnods(ielem,inode), inode=1,nnode3d)
         if(itfm.eq.2)
     .    write(86,902) ielem,
     .                   (lnods(ielem,inode), inode=1,nnode3d)
        else
         if(itfm.eq.1)
     .    write(86,900) ielem,imat(ielem),
     .                   (lnods(ielem,inode), inode=1,nnode8)
         if(itfm.eq.2)
     .    write(86,902) ielem,
     .                   (lnods(ielem,inode), inode=1,nnode8)
         if(nnode.le.(2*nnode8)) then
          write(86,901)  (lnods(ielem,inode), inode=nnode8+1,nnode3d)
         else
          write(86,901)  (lnods(ielem,inode), inode=nnode8+1,2*nnode8)
          write(86,901)  (lnods(ielem,inode), inode=2*nnode8+1,nnode)
         endif
        endif
       enddo
       if(itfm.eq.2) write(86,*) 'END_ELEMENTS'
       if(itfm.eq.2) write(86,*) 'COORDINATES'
       do ipoin=1,npoin3d
        write(86,920) ipoin,(coord(ipoin,idime), idime=1,ndime3d)
        if(ifixi.eq.1)
     .   write(88,*) ipoin,'  000   1  0.000000  0.000000  0.000000'
        if(ifixi.eq.2) then
         if(itfm.eq.1)
     .    write(88,952) ipoin,(ipres(ipoin,idime),idime=1,ndime3d),
     .                  ifunc(ipoin),(pres(ipoin,idime),idime=1,ndime3d)
         if(itfm.eq.2) then
          if(ifront.eq.1)
     .     write(88,962) ipoin,(ipres(ipoin,idime),idime=1,ndime3d+1),
     .                             (pres(ipoin,idime),idime=1,ndime3d+1)
          if(ifront.eq.2)
     .     write(88,963) ipoin,(ipres(ipoin,idime),idime=1,ndime3d+2),
     .                             (pres(ipoin,idime),idime=1,ndime3d+2)
         endif
        endif
       enddo
       if(itfm.eq.2) write(86,*) 'END_COORDINATES'
      endif
      if(igidw.eq.2.or.igidw.eq.3) then
       write(87,*) 'titulo: problema'
       write(87,*) 'subtitulo: problema'
       write(87,*) 'a'
       write(87,*) 'a'
       write(87,*) 'a'
       write(87,*) 'c'
       nnodx=1                                               ! hexahedra
       write(87,*) nelem3d, npoin3d, nnodx
       write(87,*) 'c'
       do ipoin=1,npoin3d
        x=coord(ipoin,1)
        y=coord(ipoin,2)
        z=coord(ipoin,3)
        write(87,920) ipoin,x,y,z      ! formatted
       end do
       write(87,*) 'c'
       do ielem=1,nelem3d                         ! nnode3d=8 is assumed
        if(inoco.eq.2) then
         if(imat(ielem).ge.imatnc) then
          write(87,903) ielem,(lnods(ielem,i),i=1,nnode3d/2),
     .                        (lnods(ielem,i),i=1,nnode3d/2),imat(ielem)
         else
          write(87,903) ielem,(lnods(ielem,i),i=1,nnode3d),imat(ielem)
         endif
        else
         write(87,903) ielem,(lnods(ielem,i),i=1,nnode3d),imat(ielem)
        endif
       enddo
      endif

  900 format(i6,i3,8(1x,i6))                            ! connectivities
  901 format(      8(1x,i6))
  902 format(i6,   8(1x,i6))                            ! flow connect.
  903 format(i6,8i6,i3)                                 ! GiD connect.

  920 format(i6,1x,3f15.6)                              ! coordinates

  951 FORMAT(I6,2x,2I1,2x,I2,2(2x,F8.6))                ! boundary cond.
  952 FORMAT(I6,2x,3I1,2x,I2,3(2x,F8.6))
 9900 format(i6,3x,3i1,3x,i2,3(e15.6,2x))               ! boundary cond.

  961 FORMAT(I6,1x,3I1,3(2x,F8.6))                      ! flow bc
  962 FORMAT(I6,2x,4I1,4(2x,F8.6))
  963 FORMAT(I6,2x,5I1,5(2x,F8.6))

      stop
      end
