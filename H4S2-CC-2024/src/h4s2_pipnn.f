      module nnparam
      implicit none
!***************************************************************************
!natom     ==>Number of atoms
!npes      ==>Number of PESs
!nmorse    ==>Number of morse-like potential
!***************************************************************************
      real*8,parameter::alpha=4.2d0,PI=3.141592653589793238d0,
     %radian=PI/180.0d0,bohr=0.5291772d0
      integer,parameter::natom=6,npes=1,nmorse=2058
      integer,parameter::nbond=natom*(natom-1)/2
      integer ninput,noutput,nscale
      integer nhid3,nlayer3,ifunc3
      integer nwe3,nodemax3
      integer,allocatable::nodes3a(:)
      real*8,allocatable::weight3a(:,:,:,:),bias3a(:,:,:),
     %pdel3a(:,:),pavg3a(:,:)
      end module nnparam
!***************************************************************************
      subroutine evvdvdx(xcart,v)
!***************************************************************************
      use nnparam
      implicit none
      integer i,j,k,ndriv
      real*8 v,vpes(npes),c(0:ninput),dvdx(1:natom*3,npes)
      real*8 xcart(1:3,1:natom),m(0:nmorse-1),xmorse(1:nbond),r(1:nbond)
      real*8 dvdr(1:nbond),x(1:natom*3),drdx(1:nbond,1:natom*3)
      real*8 txinput(1:ninput),dvdg(1:ninput,npes),p(0:ninput)
      real*8 dpdr(1:nbond,0:ninput),dvdxa(1:natom*3)
      ! ::::::::::::::::::::
      do i=1,natom
       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
      enddo

      k=0
      do i=1,natom-1
       do j=i+1,natom
        k=k+1
       r(k)=dot_product(xcart(:,i)-xcart(:,j),
     $                  xcart(:,i)-xcart(:,j))
        r(k)=dsqrt(r(k))
        xmorse(k)=dexp(-r(k)/(bohr*alpha))
       enddo
      enddo
      if(k.ne.natom*(natom-1)/2)stop "error"

      call bemsav(xmorse,p)
      do j=1,ninput
      txinput(j)=p(j)
      enddo

      call pot3a(txinput,vpes)
      v=0d0 
      do i=1,npes
       v=v+vpes(i)
      enddo      

      v=v/npes                   

      return
      end subroutine evvdvdx
!****************************************************************!
      subroutine pes_init
      use nnparam
      implicit none
      integer i,j,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile1,nfile2
      character f1*80

      nfile1=10
      nfile2=20
      open(nfile1,file='weights.txt',status='old')
      open(nfile2,file='biases.txt',status='old')

       read(nfile2,*) 
       read(nfile1,*)ninput,nhid3,noutput
       nscale=ninput+noutput
       nlayer3=nhid3+2 !additional one for input layer and one for output 
       allocate(nodes3a(nlayer3),pdel3a(nscale,npes),
     &pavg3a(nscale,npes))
       nodes3a(1)=ninput
       nodes3a(nlayer3)=noutput
       read(nfile1,*)(nodes3a(ihid),ihid=2,nhid3+1)
       nodemax3=0
       do i=1,nlayer3
        nodemax3=max(nodemax3,nodes3a(i))
       enddo
       allocate(weight3a(nodemax3,nodemax3,2:nlayer3,npes),
     &bias3a(nodemax3,2:nlayer3,npes))
       read(nfile1,*)ifunc3,nwe3

      do j=1,npes
       if(j.gt.1) then
        read(nfile1,*)
        read(nfile1,*)
        read(nfile1,*)
        read(nfile1,*)
        read(nfile2,*)
       endif
       read(nfile1,*)(pdel3a(i,j),i=1,nscale)
       read(nfile1,*)(pavg3a(i,j),i=1,nscale)
       iwe=0
       do ilay1=2,nlayer3
       ilay2=ilay1-1
       do inode1=1,nodes3a(ilay1)
       do inode2=1,nodes3a(ilay2) !
       read(nfile1,*)weight3a(inode2,inode1,ilay1,j)
       iwe=iwe+1
       enddo
       read(nfile2,*)bias3a(inode1,ilay1,j)
       iwe=iwe+1
       enddo
       enddo
       if (iwe.ne.nwe3) then
         write(*,*)'provided number of parameters ',nwe3
         write(*,*)'actual number of parameters ',iwe
         write(*,*)'nwe not equal to iwe, check input files or code'
         stop
       endif
      enddo

      close(nfile1)
      close(nfile2)
      write(*,*)'initialization done'

      end subroutine pes_init
!*************************************************************************
      subroutine pot3a(x,vpot3)
      use nnparam
      implicit none
      integer i,ndriv,inode1,inode2,ilay1,ilay2
      integer j,k,m,neu1,neu2
      real*8 x(ninput),y(nodemax3,nlayer3,npes),vpot3(npes)
      real*8,allocatable::nxw3(:,:),nxw2(:,:,:),nxw1(:,:,:)
      real*8,allocatable::ax(:,:),bx(:,:)
      real*8 vrange,girange(ninput)
      real*8 dvdg(1:ninput,npes),dvtmp(npes)
      real*8, external::tranfun

      dvdg=0d0
      dvtmp=0d0

      do m=1,npes
       do i=1,ninput
        y(i,1,m)=(x(i)-pavg3a(i,m))/pdel3a(i,m)
       enddo
!-->.....evaluate the hidden layer
       do ilay1=2,nlayer3-1
        ilay2=ilay1-1
        do inode1=1,nodes3a(ilay1)
         y(inode1,ilay1,m)=bias3a(inode1,ilay1,m)
         do inode2=1,nodes3a(ilay2)
          y(inode1,ilay1,m)=y(inode1,ilay1,m)+y(inode2,ilay2,m)
     &*weight3a(inode2,inode1,ilay1,m)
         enddo
         y(inode1,ilay1,m)=tranfun(y(inode1,ilay1,m),ifunc3)
        enddo
       enddo
!-->.....now evaluate the output
       ilay1=nlayer3
       ilay2=ilay1-1
       do inode1=1,nodes3a(ilay1)
        y(inode1,ilay1,m)=bias3a(inode1,ilay1,m)
        do inode2=1,nodes3a(ilay2)
         y(inode1,ilay1,m)=y(inode1,ilay1,m)+y(inode2,ilay2,m)
     &*weight3a(inode2,inode1,ilay1,m)
        enddo
       enddo
!-->.....the value of output layer is the fitted potential 
       vpot3(m)=y(nodes3a(nlayer3),nlayer3,m)
     &*pdel3a(nscale,m)+pavg3a(nscale,m)
      enddo

      return
      end subroutine pot3a
!**************************************************************************
        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x

        if (ifunc.eq.1) then
        tranfun=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+exp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end function tranfun
!**************************************************************************
        FUNCTION LOGSIG(X)
        REAL*8 X,LOGSIG
        LOGSIG=1.d0/(1.d0+DEXP(-X))
        RETURN
        END FUNCTION LOGSIG

