!***************************************************************************
!                         PIP-NN PES for NH4-
!***************************************************************************
!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
!-->  Modified by Kaisheng Song at 10 July 2021
      module nnparam
      implicit none
!***************************************************************************
!natom     ==>Number of atoms
!npes      ==>Number of PESs
!nmorse    ==>Number of morse-like potential
!***************************************************************************
      real*8,parameter::alpha=4.0d0,PI=3.141592653589793238d0,
     &radian=PI/180.0d0,bohr=0.5291772d0
      integer,parameter::natom=5,npes=1,nmorse=496
      integer,parameter::nbond=natom*(natom-1)/2
      integer ninput,noutput,nscale
      integer nhid3,nlayer3,ifunc3
      integer nwe3,nodemax3
      integer,allocatable::nodes3a(:)
      real*8,allocatable::weight3a(:,:,:,:),bias3a(:,:,:),
     &pdel3a(:,:),pavg3a(:,:)
      end module nnparam
!***************************************************************************
      subroutine evvdvdx(xcart,v)
!     subroutine evvdvdx(xcart,v,vpes,dvdxa,dvdx,ndriv)
!***************************************************************************
!Subroutine to calculate the average potential energy v and analytical gradient dvdxa
!Call pes_init to read files and initialize before calling evvdvdx().
!Xcart     ==>Cartesian coordinates in angstrom
!             H1: xcart(1,1), xcart(2,1), xcart(3,1)
!             H2: xcart(1,2), xcart(2,2), xcart(3,2)
!             H3: xcart(1,3), xcart(2,3), xcart(3,3)
!             H4: xcart(1,4), xcart(2,4), xcart(3,4)
!             N5: xcart(1,5), xcart(2,5), xcart(3,5)
!v         ==>Average potential energy(in eV), v=sum(vpes)/npes
!vpes      ==>Potential energy(in eV) for each PES
!dvdxa     ==>Average gradient(in eV/ang), dvdxa=sum(dvdx)/npes
!dvdx      ==>Gradient(in eV/ang) for each PES
!ndriv     ==>ndriv=1 to compute the analytical gradient, otherwise ndriv=0
!***************************************************************************
      use nnparam
      implicit none
      integer i,j,k,ndriv
      real*8 v,vpes(npes),c(0:ninput),dvdx(1:natom*3,npes)
      real*8 xcart(1:3,1:natom),m(0:nmorse-1),xbond(1:nbond),r(1:nbond)
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
        xbond(k)=dexp(-r(k)/(alpha*bohr))
       enddo
      enddo
      if(k.ne.natom*(natom-1)/2)stop "error"

      call bemsav(xbond,m,p)

      do j=1,ninput
      txinput(j)=p(j)
!      write(999,*)p(j)
      enddo

      call pot3a(txinput,vpes)
!     call pot3a(txinput,vpes,dvdg,ndriv)

      v=0d0 
      do i=1,npes
       v=v+vpes(i)
      enddo      

      v=v/npes                   
 
!      if(ndriv.eq.1)then
!       call evdvdr(r,c,m,p,dpdr!)
!       call evdrdx(xcart,r,x,drdx)
!       call evdvdx(dvdg,dpdr,drdx,dvdxa,dvdx)
!      endif
      
      return
      end subroutine evvdvdx
!****************************************************************!
!-->  read NN weights and biases from matlab output
      subroutine pes_init
      use nnparam
      implicit none
      integer i,j,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile1,nfile2
      character f1*80

      nfile1=4
      nfile2=7
      open(nfile1,file='weights.txt',status='old')
      open(nfile2,file='biases.txt',status='old')

       read(nfile1,*)
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
!     subroutine pot3a(x,vpot3,dvdg,ndriv)
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

!      if(ndriv.eq.1)then
!       neu1=nodes3a(2);neu2=nodes3a(3)
!       allocate(nxw1(1:ninput,1:neu1,npes),nxw2(1:neu1,1:neu2,npes),
!     $ax(1:neu1,npes),bx(1:neu2,npes),nxw3(1:neu2,npes))
!       do m=1,npes
!        do i=1,ninput
!         do j=1,neu1
!          nxw1(i,j,m)=weight3a(i,j,2,m)
!         enddo
!        enddo
!        do j=1,neu1
!         do k=1,neu2
!          nxw2(j,k,m)=weight3a(j,k,3,m)
!         enddo
!        enddo
!        do j=1,neu1
!         ax(j,m)=y(j,2,m)
!        enddo
!        do k=1,neu2
!         bx(k,m)=y(k,3,m)
!         nxw3(k,m)=weight3a(k,1,4,m)
!        enddo
!
!        do i=1,ninput
!         do k=1,neu2
!          dvtmp=0d0
!          do j=1,neu1
!           dvtmp(m)=dvtmp(m)+nxw2(j,k,m)*nxw1(i,j,m)*(1d0-ax(j,m)**2)
!          enddo
!          dvdg(i,m)=dvdg(i,m)+nxw3(k,m)*dvtmp(m)*(1d0-bx(k,m)**2)
!         enddo
!         dvdg(i,m)=dvdg(i,m)*pdel3a(nscale,m)/pdel3a(i,m)
!        enddo
!
!       enddo
!      endif
            
      return
      end subroutine pot3a
!**************************************************************************
        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
c    ifunc=2, transfer function is log sigmoid function, 'logsig'
c    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
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
!!**************************************************************************
!      subroutine evdvdx(dvdg,dpdr,drdx,dvdxa,dvdx)
!      use nnparam
!      implicit none
!      integer i,j,k
!      real*8 dvdx(1:natom*3,npes),xcart(1:3,1:natom)
!      real*8 drdx(1:nbond,1:natom*3),x(1:ninput)
!      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),dvdp(0:ninput)
!      real*8 dgdx(1:natom*3,1:ninput),dgdr(1:nbond,1:ninput)
!      real*8 dvdxa(1:natom*3),dvdg(1:ninput,npes),rtmp
! 
!      dgdx=0d0
!      dgdr(1:3,1:12)=dpdr(1:3,1:12)
!      do i=1,9
!       do j=1,12
!        do k=1,3
!         dgdx(i,j)=dgdx(i,j) + dgdr(k,j)*drdx(k,i)
!        enddo
!       enddo
!      enddo
!
!      dvdx=0d0;dvdxa=0
!      
!      do k=1,npes
!       do i=1,9
!        do j=1,12
!         dvdx(i,k)=dvdx(i,k) + dvdg(j,k)*dgdx(i,j)
!        enddo
!       enddo
!      enddo
!
!      do i=1,9
!       do k=1,npes
!        dvdxa(i)=dvdxa(i)+dvdx(i,k)
!       enddo
!      enddo
!
!      dvdxa=dvdxa/dble(npes)
!
!      return
!      end subroutine evdvdx
!!*********************************************************************** 
!      subroutine evdvdr(r,c,m,p,dpdr)
!      use nnparam
!      implicit none
!      integer i,j
!      real*8 r(1:nbond),dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4)
!      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),c(0:ninput)
!      real*8 m(0:nmorse-1),p(0:ninput)
!
!      dvdr(:)=0d0
!      call evdmsdr(r,dmsdr)
!      call evdmdr(m,dmsdr,dmdr)
!      call evdpdr(p,dmdr,dpdr)
!
!      return
!      end subroutine evdvdr
!!***********************************************************************
!      subroutine evdrdx(xcart,r,x,drdx)
!      use nnparam
!      implicit none
!      integer i,j,k
!      real*8 r(1:nbond),xcart(1:3,1:natom),drdx(1:nbond,1:natom*3)
!      real*8 x(1:natom*3)
!      drdx(:,:)=0d0
!
!!      do i=1,natom-1
!!       do j=i+1,natom
!!        k=k+1
!!        drdx(k,(i-1)*3+1)=(x((i-1)*3+1)-x((j-1)*3+1))/r(k)
!!        drdx(k,(i-1)*3+2)=(x((i-1)*3+2)-x((j-1)*3+2))/r(k)
!!        drdx(k,(i-1)*3+3)=(x((i-1)*3+3)-x((j-1)*3+3))/r(k)
!!        drdx(k,(j-1)*3+1)=-1d0*drdx(k,(i-1)*3+1)
!!        drdx(k,(j-1)*3+2)=-1d0*drdx(k,(i-1)*3+2)
!!        drdx(k,(j-1)*3+3)=-1d0*drdx(k,(i-1)*3+3)
!!       enddo
!!      enddo
!
!!dr1dx
!      drdx(1,1)=(x(1)-x(4))/r(1)
!      drdx(1,2)=(x(2)-x(5))/r(1)
!      drdx(1,3)=(x(3)-x(6))/r(1)
!      drdx(1,4)=-drdx(1,1)
!      drdx(1,5)=-drdx(1,2)
!      drdx(1,6)=-drdx(1,3)
!!dr2dx
!      drdx(2,1)=(x(1)-x(7))/r(2)
!      drdx(2,2)=(x(2)-x(8))/r(2)
!      drdx(2,3)=(x(3)-x(9))/r(2)
!      drdx(2,7)=-drdx(2,1)
!      drdx(2,8)=-drdx(2,2)
!      drdx(2,9)=-drdx(2,3)
!!dr3dx
!      drdx(3,4)=(x(4)-x(7))/r(3)
!      drdx(3,5)=(x(5)-x(8))/r(3)
!      drdx(3,6)=(x(6)-x(9))/r(3)
!      drdx(3,7)=-drdx(3,4)
!      drdx(3,8)=-drdx(3,5)
!      drdx(3,9)=-drdx(3,6)
! 
!      return
!      end subroutine evdrdx
!!*************************************************************************
!      subroutine evdmsdr(r,dmsdr)
!      use nnparam
!      implicit none
!      integer i,j
!      real*8 dmsdr(1:nbond,1:nbond),r(1:nbond)
!      dmsdr(:,:)=0d0
! 
!      do i=1,nbond
!       dmsdr(i,i)=-(1d0/1d0)*dexp(-r(i)/1.0d0)
!!      dmsdr(i,i)=-(1d0/alpha)*dexp(-r(i)/alpha)
!      enddo
! 
!      return
!      end subroutine evdmsdr
!!**************************************************************************
!      subroutine evdmdr(m,dmsdr,dmdr)
!      use nnparam
!      implicit none
!      integer i,j
!      real*8 dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4),m(0:nmorse-1)
! 
!      do i=1,nbond
!       dmdr(i,0)=0.D0
!       dmdr(i,1)=dmsdr(i,3)
!       dmdr(i,2)=dmsdr(i,2)
!       dmdr(i,3)=dmsdr(i,1)
!       dmdr(i,4)=dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
!      enddo
!
!      return
!      end subroutine evdmdr
!!*************************************************************************
!      subroutine evdpdr(p,dmdr,dpdr)
!      use nnparam
!      implicit none
!      integer i,j
!      real*8 dmdr(1:nbond,0:4),dpdr(1:nbond,0:ninput)
!      real*8 p(0:ninput)
! 
!      do i=1,nbond
!       dpdr(i,0)=dmdr(i,0)
!       dpdr(i,1)=dmdr(i,1)+dmdr(i,2)
!       dpdr(i,2)=dmdr(i,3)
!       dpdr(i,3)=dmdr(i,4)
!       dpdr(i,4)=dpdr(i,2)*P(1)+P(2)*dpdr(i,1)
!       dpdr(i,5)=dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
!       dpdr(i,6)=dpdr(i,2)*P(2)+P(2)*dpdr(i,2)
!       dpdr(i,7)=dpdr(i,2)*P(3)+P(2)*dpdr(i,3)
!       dpdr(i,8)=dpdr(i,3)*P(1)+P(3)*dpdr(i,1)
!       dpdr(i,9)=dpdr(i,2)*P(5)+P(2)*dpdr(i,5)
!       dpdr(i,10)=dpdr(i,2)*P(4)+P(2)*dpdr(i,4)
!       dpdr(i,11)=dpdr(i,1)*P(5)+P(1)*dpdr(i,5)-dpdr(i,8)
!       dpdr(i,12)=dpdr(i,2)*P(6)+P(2)*dpdr(i,6)
!      enddo
!
!      return
!      end subroutine evdpdr
!**********************************************************************
      
      subroutine bemsav(x,m,p)
        implicit none
        integer, parameter :: wp  = 8
        real(wp),dimension(1:10),intent(in)::x
        real(wp),dimension(0:207),intent(out)::p
        ! ::::::::::::::::::::
        real(wp),dimension(0:495)::m
        
        call evmono(x,m)
        call evpoly(m,p)
        
        return
      end subroutine bemsav
       
      subroutine evmono(x,m)
        implicit none
        integer, parameter :: wp  = 8
        real(wp),dimension(1:10),intent(in)::x
        real(wp),dimension(0:495),intent(out)::m
        !::::::::::::::::::::
        
        m(0) = 1.0D0
        m(1) = x(10)
        m(2) = x(9)
        m(3) = x(7)
        m(4) = x(4)
        m(5) = x(8)
        m(6) = x(6)
        m(7) = x(5)
        m(8) = x(3)
        m(9) = x(2)
        m(10) = x(1)
        m(11) = m(1)*m(2)
        m(12) = m(1)*m(3)
        m(13) = m(2)*m(3)
        m(14) = m(1)*m(4)
        m(15) = m(2)*m(4)
        m(16) = m(3)*m(4)
        m(17) = m(3)*m(5)
        m(18) = m(2)*m(6)
        m(19) = m(1)*m(7)
        m(20) = m(4)*m(5)
        m(21) = m(4)*m(6)
        m(22) = m(4)*m(7)
        m(23) = m(2)*m(8)
        m(24) = m(3)*m(8)
        m(25) = m(1)*m(9)
        m(26) = m(3)*m(9)
        m(27) = m(1)*m(10)
        m(28) = m(2)*m(10)
        m(29) = m(7)*m(8)
        m(30) = m(6)*m(9)
        m(31) = m(5)*m(10)
        m(32) = m(5)*m(6)
        m(33) = m(5)*m(7)
        m(34) = m(6)*m(7)
        m(35) = m(5)*m(8)
        m(36) = m(6)*m(8)
        m(37) = m(5)*m(9)
        m(38) = m(7)*m(9)
        m(39) = m(8)*m(9)
        m(40) = m(6)*m(10)
        m(41) = m(7)*m(10)
        m(42) = m(8)*m(10)
        m(43) = m(9)*m(10)
        m(44) = m(1)*m(13)
        m(45) = m(1)*m(15)
        m(46) = m(1)*m(16)
        m(47) = m(2)*m(16)
        m(48) = m(3)*m(20)
        m(49) = m(2)*m(21)
        m(50) = m(1)*m(22)
        m(51) = m(2)*m(24)
        m(52) = m(1)*m(26)
        m(53) = m(1)*m(28)
        m(54) = m(1)*m(17)
        m(55) = m(2)*m(17)
        m(56) = m(1)*m(18)
        m(57) = m(3)*m(18)
        m(58) = m(2)*m(19)
        m(59) = m(3)*m(19)
        m(60) = m(1)*m(20)
        m(61) = m(2)*m(20)
        m(62) = m(1)*m(21)
        m(63) = m(3)*m(21)
        m(64) = m(2)*m(22)
        m(65) = m(3)*m(22)
        m(66) = m(1)*m(23)
        m(67) = m(1)*m(24)
        m(68) = m(4)*m(23)
        m(69) = m(4)*m(24)
        m(70) = m(2)*m(25)
        m(71) = m(2)*m(26)
        m(72) = m(4)*m(25)
        m(73) = m(4)*m(26)
        m(74) = m(3)*m(27)
        m(75) = m(3)*m(28)
        m(76) = m(4)*m(27)
        m(77) = m(4)*m(28)
        m(78) = m(4)*m(32)
        m(79) = m(4)*m(33)
        m(80) = m(4)*m(34)
        m(81) = m(3)*m(35)
        m(82) = m(2)*m(36)
        m(83) = m(3)*m(37)
        m(84) = m(1)*m(38)
        m(85) = m(3)*m(39)
        m(86) = m(2)*m(40)
        m(87) = m(1)*m(41)
        m(88) = m(2)*m(42)
        m(89) = m(1)*m(43)
        m(90) = m(2)*m(32)
        m(91) = m(3)*m(32)
        m(92) = m(1)*m(33)
        m(93) = m(3)*m(33)
        m(94) = m(1)*m(34)
        m(95) = m(2)*m(34)
        m(96) = m(2)*m(35)
        m(97) = m(3)*m(36)
        m(98) = m(4)*m(35)
        m(99) = m(4)*m(36)
        m(100) = m(1)*m(37)
        m(101) = m(3)*m(38)
        m(102) = m(4)*m(37)
        m(103) = m(4)*m(38)
        m(104) = m(1)*m(39)
        m(105) = m(2)*m(39)
        m(106) = m(1)*m(40)
        m(107) = m(2)*m(41)
        m(108) = m(4)*m(40)
        m(109) = m(4)*m(41)
        m(110) = m(1)*m(42)
        m(111) = m(3)*m(42)
        m(112) = m(2)*m(43)
        m(113) = m(3)*m(43)
        m(114) = m(5)*m(29)
        m(115) = m(6)*m(29)
        m(116) = m(5)*m(30)
        m(117) = m(6)*m(38)
        m(118) = m(6)*m(39)
        m(119) = m(7)*m(39)
        m(120) = m(5)*m(40)
        m(121) = m(5)*m(41)
        m(122) = m(5)*m(42)
        m(123) = m(7)*m(42)
        m(124) = m(5)*m(43)
        m(125) = m(6)*m(43)
        m(126) = m(5)*m(34)
        m(127) = m(5)*m(39)
        m(128) = m(6)*m(42)
        m(129) = m(7)*m(43)
        m(130) = m(5)*m(36)
        m(131) = m(5)*m(38)
        m(132) = m(6)*m(41)
        m(133) = m(8)*m(43)
        m(134) = m(1)*m(47)
        m(135) = m(1)*m(48)
        m(136) = m(2)*m(48)
        m(137) = m(1)*m(49)
        m(138) = m(2)*m(63)
        m(139) = m(1)*m(64)
        m(140) = m(1)*m(65)
        m(141) = m(1)*m(51)
        m(142) = m(2)*m(69)
        m(143) = m(1)*m(71)
        m(144) = m(1)*m(73)
        m(145) = m(1)*m(75)
        m(146) = m(1)*m(77)
        m(147) = m(7)*m(66)
        m(148) = m(7)*m(67)
        m(149) = m(7)*m(68)
        m(150) = m(7)*m(69)
        m(151) = m(6)*m(70)
        m(152) = m(6)*m(71)
        m(153) = m(6)*m(72)
        m(154) = m(6)*m(73)
        m(155) = m(5)*m(74)
        m(156) = m(5)*m(75)
        m(157) = m(5)*m(76)
        m(158) = m(5)*m(77)
        m(159) = m(2)*m(78)
        m(160) = m(3)*m(78)
        m(161) = m(1)*m(79)
        m(162) = m(3)*m(79)
        m(163) = m(1)*m(80)
        m(164) = m(2)*m(80)
        m(165) = m(2)*m(81)
        m(166) = m(2)*m(97)
        m(167) = m(3)*m(98)
        m(168) = m(2)*m(99)
        m(169) = m(1)*m(83)
        m(170) = m(1)*m(101)
        m(171) = m(3)*m(102)
        m(172) = m(1)*m(103)
        m(173) = m(1)*m(85)
        m(174) = m(2)*m(85)
        m(175) = m(1)*m(86)
        m(176) = m(1)*m(107)
        m(177) = m(2)*m(108)
        m(178) = m(1)*m(109)
        m(179) = m(1)*m(88)
        m(180) = m(2)*m(111)
        m(181) = m(1)*m(112)
        m(182) = m(1)*m(113)
        m(183) = m(2)*m(91)
        m(184) = m(1)*m(93)
        m(185) = m(1)*m(95)
        m(186) = m(2)*m(98)
        m(187) = m(3)*m(99)
        m(188) = m(1)*m(102)
        m(189) = m(3)*m(103)
        m(190) = m(1)*m(105)
        m(191) = m(1)*m(108)
        m(192) = m(2)*m(109)
        m(193) = m(1)*m(111)
        m(194) = m(2)*m(113)
        m(195) = m(3)*m(114)
        m(196) = m(2)*m(115)
        m(197) = m(4)*m(114)
        m(198) = m(4)*m(115)
        m(199) = m(3)*m(116)
        m(200) = m(1)*m(117)
        m(201) = m(4)*m(116)
        m(202) = m(4)*m(117)
        m(203) = m(2)*m(118)
        m(204) = m(3)*m(118)
        m(205) = m(1)*m(119)
        m(206) = m(3)*m(119)
        m(207) = m(2)*m(120)
        m(208) = m(1)*m(121)
        m(209) = m(4)*m(120)
        m(210) = m(4)*m(121)
        m(211) = m(2)*m(122)
        m(212) = m(3)*m(122)
        m(213) = m(1)*m(123)
        m(214) = m(2)*m(123)
        m(215) = m(1)*m(124)
        m(216) = m(3)*m(124)
        m(217) = m(1)*m(125)
        m(218) = m(2)*m(125)
        m(219) = m(6)*m(119)
        m(220) = m(5)*m(123)
        m(221) = m(5)*m(125)
        m(222) = m(4)*m(126)
        m(223) = m(3)*m(127)
        m(224) = m(2)*m(128)
        m(225) = m(1)*m(129)
        m(226) = m(1)*m(78)
        m(227) = m(2)*m(79)
        m(228) = m(3)*m(80)
        m(229) = m(1)*m(81)
        m(230) = m(1)*m(82)
        m(231) = m(2)*m(83)
        m(232) = m(2)*m(84)
        m(233) = m(4)*m(85)
        m(234) = m(3)*m(86)
        m(235) = m(3)*m(87)
        m(236) = m(4)*m(88)
        m(237) = m(4)*m(89)
        m(238) = m(2)*m(130)
        m(239) = m(3)*m(130)
        m(240) = m(4)*m(130)
        m(241) = m(1)*m(131)
        m(242) = m(3)*m(131)
        m(243) = m(4)*m(131)
        m(244) = m(1)*m(132)
        m(245) = m(2)*m(132)
        m(246) = m(4)*m(132)
        m(247) = m(1)*m(133)
        m(248) = m(2)*m(133)
        m(249) = m(3)*m(133)
        m(250) = m(5)*m(115)
        m(251) = m(5)*m(117)
        m(252) = m(5)*m(118)
        m(253) = m(5)*m(119)
        m(254) = m(5)*m(132)
        m(255) = m(5)*m(128)
        m(256) = m(6)*m(123)
        m(257) = m(5)*m(129)
        m(258) = m(6)*m(129)
        m(259) = m(5)*m(133)
        m(260) = m(6)*m(133)
        m(261) = m(7)*m(133)
        m(262) = m(2)*m(160)
        m(263) = m(1)*m(162)
        m(264) = m(1)*m(164)
        m(265) = m(2)*m(167)
        m(266) = m(2)*m(187)
        m(267) = m(1)*m(171)
        m(268) = m(1)*m(189)
        m(269) = m(1)*m(174)
        m(270) = m(1)*m(177)
        m(271) = m(1)*m(192)
        m(272) = m(1)*m(180)
        m(273) = m(1)*m(194)
        m(274) = m(3)*m(197)
        m(275) = m(2)*m(198)
        m(276) = m(3)*m(201)
        m(277) = m(1)*m(202)
        m(278) = m(2)*m(204)
        m(279) = m(1)*m(206)
        m(280) = m(2)*m(209)
        m(281) = m(1)*m(210)
        m(282) = m(2)*m(212)
        m(283) = m(1)*m(214)
        m(284) = m(1)*m(216)
        m(285) = m(1)*m(218)
        m(286) = m(1)*m(159)
        m(287) = m(1)*m(160)
        m(288) = m(1)*m(227)
        m(289) = m(2)*m(162)
        m(290) = m(1)*m(228)
        m(291) = m(2)*m(228)
        m(292) = m(1)*m(165)
        m(293) = m(1)*m(166)
        m(294) = m(1)*m(167)
        m(295) = m(1)*m(168)
        m(296) = m(1)*m(231)
        m(297) = m(2)*m(170)
        m(298) = m(2)*m(171)
        m(299) = m(2)*m(172)
        m(300) = m(1)*m(233)
        m(301) = m(2)*m(233)
        m(302) = m(1)*m(234)
        m(303) = m(2)*m(235)
        m(304) = m(3)*m(177)
        m(305) = m(3)*m(178)
        m(306) = m(1)*m(236)
        m(307) = m(3)*m(236)
        m(308) = m(2)*m(237)
        m(309) = m(3)*m(237)
        m(310) = m(1)*m(195)
        m(311) = m(1)*m(196)
        m(312) = m(2)*m(197)
        m(313) = m(3)*m(198)
        m(314) = m(2)*m(199)
        m(315) = m(2)*m(200)
        m(316) = m(1)*m(201)
        m(317) = m(3)*m(202)
        m(318) = m(1)*m(203)
        m(319) = m(2)*m(205)
        m(320) = m(4)*m(204)
        m(321) = m(4)*m(206)
        m(322) = m(3)*m(207)
        m(323) = m(3)*m(208)
        m(324) = m(1)*m(209)
        m(325) = m(2)*m(210)
        m(326) = m(1)*m(212)
        m(327) = m(3)*m(213)
        m(328) = m(4)*m(211)
        m(329) = m(4)*m(214)
        m(330) = m(2)*m(216)
        m(331) = m(3)*m(218)
        m(332) = m(4)*m(215)
        m(333) = m(4)*m(217)
        m(334) = m(2)*m(195)
        m(335) = m(3)*m(196)
        m(336) = m(1)*m(197)
        m(337) = m(1)*m(198)
        m(338) = m(1)*m(199)
        m(339) = m(3)*m(200)
        m(340) = m(2)*m(201)
        m(341) = m(2)*m(202)
        m(342) = m(1)*m(204)
        m(343) = m(2)*m(206)
        m(344) = m(4)*m(203)
        m(345) = m(4)*m(205)
        m(346) = m(1)*m(207)
        m(347) = m(2)*m(208)
        m(348) = m(3)*m(209)
        m(349) = m(3)*m(210)
        m(350) = m(1)*m(211)
        m(351) = m(3)*m(214)
        m(352) = m(4)*m(212)
        m(353) = m(4)*m(213)
        m(354) = m(2)*m(215)
        m(355) = m(3)*m(217)
        m(356) = m(4)*m(216)
        m(357) = m(4)*m(218)
        m(358) = m(1)*m(222)
        m(359) = m(2)*m(222)
        m(360) = m(3)*m(222)
        m(361) = m(1)*m(223)
        m(362) = m(2)*m(223)
        m(363) = m(4)*m(223)
        m(364) = m(1)*m(224)
        m(365) = m(3)*m(224)
        m(366) = m(4)*m(224)
        m(367) = m(2)*m(225)
        m(368) = m(3)*m(225)
        m(369) = m(4)*m(225)
        m(370) = m(2)*m(239)
        m(371) = m(2)*m(240)
        m(372) = m(3)*m(240)
        m(373) = m(1)*m(242)
        m(374) = m(1)*m(243)
        m(375) = m(3)*m(243)
        m(376) = m(1)*m(245)
        m(377) = m(1)*m(246)
        m(378) = m(2)*m(246)
        m(379) = m(1)*m(248)
        m(380) = m(1)*m(249)
        m(381) = m(2)*m(249)
        m(382) = m(4)*m(250)
        m(383) = m(4)*m(251)
        m(384) = m(3)*m(252)
        m(385) = m(3)*m(253)
        m(386) = m(4)*m(254)
        m(387) = m(2)*m(255)
        m(388) = m(2)*m(256)
        m(389) = m(1)*m(257)
        m(390) = m(1)*m(258)
        m(391) = m(3)*m(259)
        m(392) = m(2)*m(260)
        m(393) = m(1)*m(261)
        m(394) = m(2)*m(250)
        m(395) = m(3)*m(250)
        m(396) = m(1)*m(251)
        m(397) = m(3)*m(251)
        m(398) = m(2)*m(252)
        m(399) = m(1)*m(253)
        m(400) = m(4)*m(252)
        m(401) = m(4)*m(253)
        m(402) = m(1)*m(254)
        m(403) = m(2)*m(254)
        m(404) = m(3)*m(255)
        m(405) = m(1)*m(256)
        m(406) = m(4)*m(255)
        m(407) = m(4)*m(256)
        m(408) = m(3)*m(257)
        m(409) = m(2)*m(258)
        m(410) = m(4)*m(257)
        m(411) = m(4)*m(258)
        m(412) = m(1)*m(259)
        m(413) = m(2)*m(259)
        m(414) = m(1)*m(260)
        m(415) = m(3)*m(260)
        m(416) = m(2)*m(261)
        m(417) = m(3)*m(261)
        m(418) = m(5)*m(219)
        m(419) = m(5)*m(256)
        m(420) = m(5)*m(258)
        m(421) = m(5)*m(260)
        m(422) = m(5)*m(261)
        m(423) = m(6)*m(261)
        m(424) = m(3)*m(135)
        m(425) = m(3)*m(136)
        m(426) = m(2)*m(137)
        m(427) = m(2)*m(138)
        m(428) = m(1)*m(139)
        m(429) = m(1)*m(140)
        m(430) = m(4)*m(135)
        m(431) = m(4)*m(136)
        m(432) = m(4)*m(137)
        m(433) = m(4)*m(138)
        m(434) = m(4)*m(139)
        m(435) = m(4)*m(140)
        m(436) = m(2)*m(141)
        m(437) = m(3)*m(141)
        m(438) = m(2)*m(142)
        m(439) = m(3)*m(142)
        m(440) = m(1)*m(143)
        m(441) = m(3)*m(143)
        m(442) = m(1)*m(144)
        m(443) = m(3)*m(144)
        m(444) = m(1)*m(145)
        m(445) = m(2)*m(145)
        m(446) = m(1)*m(146)
        m(447) = m(2)*m(146)
        m(448) = m(4)*m(159)
        m(449) = m(4)*m(160)
        m(450) = m(4)*m(161)
        m(451) = m(4)*m(162)
        m(452) = m(4)*m(163)
        m(453) = m(4)*m(164)
        m(454) = m(3)*m(165)
        m(455) = m(2)*m(166)
        m(456) = m(3)*m(167)
        m(457) = m(2)*m(168)
        m(458) = m(3)*m(169)
        m(459) = m(1)*m(170)
        m(460) = m(3)*m(171)
        m(461) = m(1)*m(172)
        m(462) = m(3)*m(173)
        m(463) = m(3)*m(174)
        m(464) = m(2)*m(175)
        m(465) = m(1)*m(176)
        m(466) = m(2)*m(177)
        m(467) = m(1)*m(178)
        m(468) = m(2)*m(179)
        m(469) = m(2)*m(180)
        m(470) = m(1)*m(181)
        m(471) = m(1)*m(182)
        m(472) = m(19)*m(114)
        m(473) = m(19)*m(115)
        m(474) = m(23)*m(114)
        m(475) = m(24)*m(115)
        m(476) = m(18)*m(116)
        m(477) = m(18)*m(117)
        m(478) = m(21)*m(118)
        m(479) = m(22)*m(119)
        m(480) = m(23)*m(119)
        m(481) = m(25)*m(116)
        m(482) = m(26)*m(117)
        m(483) = m(25)*m(118)
        m(484) = m(17)*m(120)
        m(485) = m(17)*m(121)
        m(486) = m(20)*m(122)
        m(487) = m(22)*m(123)
        m(488) = m(24)*m(123)
        m(489) = m(20)*m(124)
        m(490) = m(21)*m(125)
        m(491) = m(26)*m(125)
        m(492) = m(27)*m(120)
        m(493) = m(28)*m(121)
        m(494) = m(27)*m(122)
        m(495) = m(28)*m(124)

        return
      end subroutine evmono

      subroutine evpoly(m,p)
        implicit none
        integer, parameter :: wp  = 8
        real(wp),dimension(0:495),intent(in)::m
        real(wp),dimension(0:207),intent(out)::p
        !::::::::::::::::::::
        
        p(0) = m(0)
        p(1) = m(1) + m(2) + m(3) + m(4)
        p(2) = m(5) + m(6) + m(7) + m(8) 
     &         + m(9) + m(10)
        p(3) = m(11) + m(12) + m(13) + m(14) 
     &         + m(15) + m(16)
        p(4) = m(17) + m(18) + m(19) + m(20) 
     &         + m(21) + m(22) + m(23) 
     &         + m(24) + m(25) + m(26) 
     &         + m(27) + m(28)
        p(5) = m(29) + m(30) + m(31)
        p(6) = p(1)*p(2) - p(4)
        p(7) = m(32) + m(33) + m(34) + m(35) 
     &         + m(36) + m(37) + m(38) 
     &         + m(39) + m(40) + m(41) 
     &         + m(42) + m(43)
        p(8) = p(1)*p(1) - p(3) - p(3)
        p(9) = p(2)*p(2) - p(7) - p(5) - p(7) - p(5)
        p(10) = m(44) + m(45) + m(46) + m(47)
        p(11) = m(48) + m(49) + m(50) + m(51) 
     &         + m(52) + m(53)
        p(12) = m(54) + m(55) + m(56) + m(57) 
     &         + m(58) + m(59) + m(60) 
     &         + m(61) + m(62) + m(63) 
     &         + m(64) + m(65) + m(66) 
     &         + m(67) + m(68) + m(69) 
     &         + m(70) + m(71) + m(72) 
     &         + m(73) + m(74) + m(75) 
     &         + m(76) + m(77)
        p(13) = p(1)*p(5)
        p(14) = p(2)*p(3) - p(12) - p(11)
        p(15) = m(78) + m(79) + m(80) + m(81) 
     &         + m(82) + m(83) + m(84) 
     &         + m(85) + m(86) + m(87) 
     &         + m(88) + m(89)
        p(16) = m(90) + m(91) + m(92) + m(93) 
     &         + m(94) + m(95) + m(96) 
     &         + m(97) + m(98) + m(99) 
     &         + m(100) + m(101) + m(102) 
     &         + m(103) + m(104) + m(105) 
     &         + m(106) + m(107) + m(108) 
     &         + m(109) + m(110) + m(111) 
     &         + m(112) + m(113)
        p(17) = m(114) + m(115) + m(116) + m(117) 
     &         + m(118) + m(119) + m(120) 
     &         + m(121) + m(122) + m(123) 
     &         + m(124) + m(125)
        p(18) = m(126) + m(127) + m(128) + m(129)
        p(19) = p(1)*p(7) - p(16) - p(15)
        p(20) = m(130) + m(131) + m(132) + m(133)
        p(21) = p(1)*p(3) - p(10) - p(10) - p(10)
        p(22) = p(1)*p(4) - p(12) - p(11) - p(11)
        p(23) = p(2)*p(8) - p(22)
        p(24) = p(2)*p(4) - p(16) - p(15) - p(13) - p(15)
        p(25) = p(2)*p(5) - p(17)
        p(26) = p(1)*p(9) - p(24)
        p(27) = p(2)*p(7) - p(18) - p(20) - p(17) - p(18) 
     &         - p(20) - p(17) - p(18) 
     &         - p(20)
        p(28) = p(1)*p(8) - p(21)
        p(29) = p(2)*p(9) - p(27) - p(25)
        p(30) = m(134)
        p(31) = m(135) + m(136) + m(137) + m(138) 
     &         + m(139) + m(140) + m(141) 
     &         + m(142) + m(143) + m(144) 
     &         + m(145) + m(146)
        p(32) = m(147) + m(148) + m(149) + m(150) 
     &         + m(151) + m(152) + m(153) 
     &         + m(154) + m(155) + m(156) 
     &         + m(157) + m(158)
        p(33) = p(2)*p(10) - p(31)
        p(34) = p(3)*p(5) - p(32)
        p(35) = m(159) + m(160) + m(161) + m(162) 
     &         + m(163) + m(164) + m(165) 
     &         + m(166) + m(167) + m(168) 
     &         + m(169) + m(170) + m(171) 
     &         + m(172) + m(173) + m(174) 
     &         + m(175) + m(176) + m(177) 
     &         + m(178) + m(179) + m(180) 
     &         + m(181) + m(182)
        p(36) = m(183) + m(184) + m(185) + m(186) 
     &         + m(187) + m(188) + m(189) 
     &         + m(190) + m(191) + m(192) 
     &         + m(193) + m(194)
        p(37) = m(195) + m(196) + m(197) + m(198) 
     &         + m(199) + m(200) + m(201) 
     &         + m(202) + m(203) + m(204) 
     &         + m(205) + m(206) + m(207) 
     &         + m(208) + m(209) + m(210) 
     &         + m(211) + m(212) + m(213) 
     &         + m(214) + m(215) + m(216) 
     &         + m(217) + m(218)
        p(38) = m(219) + m(220) + m(221)
        p(39) = m(222) + m(223) + m(224) + m(225)
        p(40) = m(226) + m(227) + m(228) + m(229) 
     &         + m(230) + m(231) + m(232) 
     &         + m(233) + m(234) + m(235) 
     &         + m(236) + m(237)
        p(41) = p(3)*p(7) - p(36) - p(40) - p(35)
        p(42) = p(1)*p(17) - p(37)
        p(43) = p(1)*p(18) - p(39)
        p(44) = m(238) + m(239) + m(240) + m(241) 
     &         + m(242) + m(243) + m(244) 
     &         + m(245) + m(246) + m(247) 
     &         + m(248) + m(249)
        p(45) = m(250) + m(251) + m(252) + m(253) 
     &         + m(254) + m(255) + m(256) 
     &         + m(257) + m(258) + m(259) 
     &         + m(260) + m(261)
        p(46) = p(1)*p(20) - p(44)
        p(47) = p(1)*p(10) - p(30) - p(30) - p(30) - p(30)
        p(48) = p(1)*p(11) - p(31)
        p(49) = p(3)*p(4) - p(33) - p(31) - p(48) - p(31)
        p(50) = p(1)*p(12) - p(33) - p(31) - p(49) - p(33) 
     &         - p(31)
        p(51) = p(5)*p(8)
        p(52) = p(1)*p(14) - p(33)
        p(53) = p(1)*p(15) - p(40) - p(35)
        p(54) = p(1)*p(16) - p(41) - p(36) - p(35) - p(36)
        p(55) = p(1)*p(19) - p(41) - p(40)
        p(56) = p(2)*p(11) - p(35) - p(34)
        p(57) = p(4)*p(5) - p(37)
        p(58) = p(2)*p(12) - p(41) - p(36) - p(40) - p(35) 
     &         - p(32) - p(36) - p(40) 
     &         - p(32)
        p(59) = p(1)*p(25) - p(57)
        p(60) = p(2)*p(14) - p(41) - p(34)
        p(61) = p(2)*p(15) - p(39) - p(44) - p(37) - p(39) 
     &         - p(39)
        p(62) = p(4)*p(7) - p(43) - p(39) - p(44) - p(42) 
     &         - p(37) - p(61) - p(39) 
     &         - p(44) - p(39)
        p(63) = p(5)*p(7) - p(45)
        p(64) = p(2)*p(16) - p(43) - p(44) - p(42) - p(37) 
     &         - p(62) - p(43) - p(44)
        p(65) = p(2)*p(17) - p(45) - p(38) - p(63) - p(45) 
     &         - p(38) - p(38) - p(38)
        p(66) = p(2)*p(18) - p(45)
        p(67) = p(1)*p(27) - p(64) - p(62) - p(61)
        p(68) = p(2)*p(20) - p(45)
        p(69) = p(3)*p(3) - p(30) - p(47) - p(30) - p(47) 
     &         - p(30) - p(30) - p(30) 
     &         - p(30)
        p(70) = p(3)*p(8) - p(47)
        p(71) = p(1)*p(22) - p(49) - p(48)
        p(72) = p(2)*p(28) - p(71)
        p(73) = p(1)*p(24) - p(58) - p(56) - p(56)
        p(74) = p(5)*p(5) - p(38) - p(38)
        p(75) = p(8)*p(9) - p(73)
        p(76) = p(7)*p(7) - p(45) - p(38) - p(66) - p(68) 
     &         - p(65) - p(45) - p(38) 
     &         - p(66) - p(68) - p(65) 
     &         - p(45) - p(38) - p(45) 
     &         - p(38)
        p(77) = p(2)*p(24) - p(62) - p(61) - p(57)
        p(78) = p(5)*p(9) - p(65)
        p(79) = p(1)*p(29) - p(77)
        p(80) = p(7)*p(9) - p(66) - p(68) - p(63)
        p(81) = p(1)*p(28) - p(70)
        p(82) = p(2)*p(29) - p(80) - p(78)
        p(83) = p(30)*p(2)
        p(84) = p(5)*p(10)
        p(85) = m(262) + m(263) + m(264) + m(265) 
     &         + m(266) + m(267) + m(268) 
     &         + m(269) + m(270) + m(271) 
     &         + m(272) + m(273)
        p(86) = m(274) + m(275) + m(276) + m(277) 
     &         + m(278) + m(279) + m(280) 
     &         + m(281) + m(282) + m(283) 
     &         + m(284) + m(285)
        p(87) = m(286) + m(287) + m(288) + m(289) 
     &         + m(290) + m(291) + m(292) 
     &         + m(293) + m(294) + m(295) 
     &         + m(296) + m(297) + m(298) 
     &         + m(299) + m(300) + m(301) 
     &         + m(302) + m(303) + m(304) 
     &         + m(305) + m(306) + m(307) 
     &         + m(308) + m(309)
        p(88) = p(7)*p(10) - p(87) - p(85)
        p(89) = m(310) + m(311) + m(312) + m(313) 
     &         + m(314) + m(315) + m(316) 
     &         + m(317) + m(318) + m(319) 
     &         + m(320) + m(321) + m(322) 
     &         + m(323) + m(324) + m(325) 
     &         + m(326) + m(327) + m(328) 
     &         + m(329) + m(330) + m(331) 
     &         + m(332) + m(333)
        p(90) = m(334) + m(335) + m(336) + m(337) 
     &         + m(338) + m(339) + m(340) 
     &         + m(341) + m(342) + m(343) 
     &         + m(344) + m(345) + m(346) 
     &         + m(347) + m(348) + m(349) 
     &         + m(350) + m(351) + m(352) 
     &         + m(353) + m(354) + m(355) 
     &         + m(356) + m(357)
        p(91) = p(1)*p(38)
        p(92) = p(3)*p(17) - p(89) - p(90) - p(86)
        p(93) = m(358) + m(359) + m(360) + m(361) 
     &         + m(362) + m(363) + m(364) 
     &         + m(365) + m(366) + m(367) 
     &         + m(368) + m(369)
        p(94) = p(3)*p(18) - p(93)
        p(95) = m(370) + m(371) + m(372) + m(373) 
     &         + m(374) + m(375) + m(376) 
     &         + m(377) + m(378) + m(379) 
     &         + m(380) + m(381)
        p(96) = m(382) + m(383) + m(384) + m(385) 
     &         + m(386) + m(387) + m(388) 
     &         + m(389) + m(390) + m(391) 
     &         + m(392) + m(393)
        p(97) = m(394) + m(395) + m(396) + m(397) 
     &         + m(398) + m(399) + m(400) 
     &         + m(401) + m(402) + m(403) 
     &         + m(404) + m(405) + m(406) 
     &         + m(407) + m(408) + m(409) 
     &         + m(410) + m(411) + m(412) 
     &         + m(413) + m(414) + m(415) 
     &         + m(416) + m(417)
        p(98) = m(418) + m(419) + m(420) + m(421) 
     &         + m(422) + m(423)
        p(99) = p(3)*p(20) - p(95)
        p(100) = p(1)*p(45) - p(97) - p(96)
        p(101) = p(30)*p(1)
        p(102) = m(424) + m(425) + m(426) + m(427) 
     &         + m(428) + m(429) + m(430) 
     &         + m(431) + m(432) + m(433) 
     &         + m(434) + m(435) + m(436) 
     &         + m(437) + m(438) + m(439) 
     &         + m(440) + m(441) + m(442) 
     &         + m(443) + m(444) + m(445) 
     &         + m(446) + m(447)
        p(103) = p(4)*p(10) - p(83) - p(102) - p(83)
        p(104) = p(1)*p(31) - p(83) - p(102) - p(83)
        p(105) = p(1)*p(32) - p(84) - p(84)
        p(106) = p(1)*p(33) - p(83) - p(103) - p(83)
        p(107) = p(1)*p(34) - p(84)
        p(108) = m(448) + m(449) + m(450) + m(451) 
     &         + m(452) + m(453) + m(454) 
     &         + m(455) + m(456) + m(457) 
     &         + m(458) + m(459) + m(460) 
     &         + m(461) + m(462) + m(463) 
     &         + m(464) + m(465) + m(466) 
     &         + m(467) + m(468) + m(469) 
     &         + m(470) + m(471)
        p(109) = p(1)*p(35) - p(87) - p(85) - p(108) - p(85)
        p(110) = p(1)*p(36) - p(88) - p(85)
        p(111) = p(1)*p(37) - p(89) - p(90) - p(86) - p(86)
        p(112) = p(1)*p(39) - p(93)
        p(113) = p(3)*p(15) - p(87) - p(85) - p(108)
        p(114) = p(3)*p(16) - p(88) - p(87) - p(85) - p(110) 
     &         - p(109) - p(88) - p(85)
        p(115) = p(1)*p(40) - p(87) - p(113)
        p(116) = p(3)*p(19) - p(88) - p(87) - p(115)
        p(117) = p(8)*p(17) - p(111)
        p(118) = p(8)*p(18) - p(112)
        p(119) = p(1)*p(44) - p(99) - p(95) - p(95)
        p(120) = p(1)*p(46) - p(99)
        p(121) = p(5)*p(11) - p(86)
        p(122) = p(6)*p(11) - p(87) - p(109) - p(107)
        p(123) = p(5)*p(12) - p(89) - p(90)
        p(124) = p(9)*p(10) - p(122)
        p(125) = p(5)*p(14) - p(92)
        p(126) = p(7)*p(11) - p(93) - p(95) - p(90)
        p(127) = p(5)*p(15) - p(96)
        p(128) = m(472) + m(473) + m(474) + m(475) 
     &         + m(476) + m(477) + m(478) 
     &         + m(479) + m(480) + m(481) 
     &         + m(482) + m(483) + m(484) 
     &         + m(485) + m(486) + m(487) 
     &         + m(488) + m(489) + m(490) 
     &         + m(491) + m(492) + m(493) 
     &         + m(494) + m(495)
        p(129) = p(2)*p(35) - p(93) - p(95) - p(90) - p(86) 
     &         - p(126) - p(93) - p(95) 
     &         - p(86)
        p(130) = p(2)*p(36) - p(94) - p(95) - p(89)
        p(131) = p(5)*p(16) - p(97) - p(128)
        p(132) = p(2)*p(37) - p(97) - p(96) - p(91) - p(131) 
     &         - p(127) - p(96) - p(91)
        p(133) = p(2)*p(38) - p(98)
        p(134) = p(2)*p(39) - p(96)
        p(135) = p(4)*p(18) - p(97) - p(134)
        p(136) = p(5)*p(18)
        p(137) = p(2)*p(40) - p(93) - p(99) - p(89)
        p(138) = p(4)*p(19) - p(93) - p(99) - p(90) - p(118) 
     &         - p(137) - p(117) - p(99)
        p(139) = p(5)*p(19) - p(100)
        p(140) = p(7)*p(14) - p(94) - p(99) - p(90)
        p(141) = p(1)*p(65) - p(132)
        p(142) = p(1)*p(66) - p(135) - p(134)
        p(143) = p(4)*p(20) - p(100) - p(96)
        p(144) = p(5)*p(20)
        p(145) = p(2)*p(44) - p(97) - p(96) - p(143)
        p(146) = p(2)*p(45) - p(98) - p(136) - p(144) - p(98) 
     &         - p(98) - p(98)
        p(147) = p(2)*p(46) - p(100)
        p(148) = p(3)*p(10) - p(101) - p(101) - p(101)
        p(149) = p(8)*p(10) - p(101)
        p(150) = p(3)*p(11) - p(83) - p(102)
        p(151) = p(8)*p(11) - p(104)
        p(152) = p(3)*p(22) - p(103) - p(102) - p(151)
        p(153) = p(1)*p(49) - p(103) - p(102) - p(152) - p(103)
        p(154) = p(2)*p(69) - p(153) - p(150)
        p(155) = p(8)*p(12) - p(106) - p(102) - p(152)
        p(156) = p(5)*p(28)
        p(157) = p(8)*p(14) - p(103)
        p(158) = p(1)*p(53) - p(113) - p(108)
        p(159) = p(1)*p(54) - p(114) - p(110) - p(109)
        p(160) = p(1)*p(55) - p(116) - p(115)
        p(161) = p(1)*p(56) - p(122)
        p(162) = p(5)*p(22) - p(111)
        p(163) = p(3)*p(24) - p(124) - p(122) - p(161) - p(122)
        p(164) = p(1)*p(74)
        p(165) = p(1)*p(58) - p(124) - p(122) - p(163) - p(124) 
     &         - p(122)
        p(166) = p(5)*p(23) - p(117)
        p(167) = p(1)*p(60) - p(124)
        p(168) = p(1)*p(61) - p(137) - p(129) - p(126)
        p(169) = p(1)*p(62) - p(138) - p(130) - p(126)
        p(170) = p(5)*p(17) - p(98) - p(133) - p(98)
        p(171) = p(1)*p(64) - p(140) - p(130) - p(129)
        p(172) = p(1)*p(67) - p(140) - p(138) - p(137)
        p(173) = p(7)*p(15) - p(97) - p(96) - p(91) - p(134) 
     &         - p(143) - p(132) - p(96) 
     &         - p(134)
        p(174) = p(7)*p(16) - p(100) - p(97) - p(96) - p(91) 
     &         - p(142) - p(135) - p(145) 
     &         - p(143) - p(141) - p(132) 
     &         - p(100) - p(97) - p(96) 
     &         - p(91) - p(135) - p(145)
        p(175) = p(7)*p(17) - p(98) - p(146) - p(136) - p(144) 
     &         - p(133) - p(98) - p(136) 
     &         - p(144) - p(133) - p(98) 
     &         - p(98)
        p(176) = p(7)*p(18) - p(98) - p(146) - p(98)
        p(177) = p(1)*p(76) - p(174) - p(173)
        p(178) = p(7)*p(20) - p(98) - p(146) - p(98)
        p(179) = p(2)*p(56) - p(126) - p(121)
        p(180) = p(5)*p(24) - p(132)
        p(181) = p(2)*p(58) - p(138) - p(130) - p(137) - p(129) 
     &         - p(123)
        p(182) = p(1)*p(78) - p(180)
        p(183) = p(2)*p(60) - p(140) - p(125)
        p(184) = p(9)*p(15) - p(134) - p(145) - p(131)
        p(185) = p(2)*p(62) - p(135) - p(143) - p(132) - p(128) 
     &         - p(174) - p(135)
        p(186) = p(5)*p(27) - p(146) - p(175)
        p(187) = p(9)*p(16) - p(142) - p(143) - p(139) - p(127) 
     &         - p(185)
        p(188) = p(2)*p(65) - p(146) - p(133) - p(175)
        p(189) = p(9)*p(18) - p(144)
        p(190) = p(1)*p(80) - p(187) - p(185) - p(184)
        p(191) = p(9)*p(20) - p(136)
        p(192) = p(1)*p(69) - p(148)
        p(193) = p(3)*p(28) - p(149)
        p(194) = p(1)*p(71) - p(152) - p(151)
        p(195) = p(2)*p(81) - p(194)
        p(196) = p(1)*p(73) - p(163) - p(161)
        p(197) = p(9)*p(28) - p(196)
        p(198) = p(1)*p(77) - p(181) - p(179) - p(179)
        p(199) = p(2)*p(74) - p(170)
        p(200) = p(8)*p(29) - p(198)
        p(201) = p(2)*p(76) - p(176) - p(178) - p(175)
        p(202) = p(2)*p(77) - p(185) - p(184) - p(180)
        p(203) = p(5)*p(29) - p(188)
        p(204) = p(1)*p(82) - p(202)
        p(205) = p(7)*p(29) - p(189) - p(191) - p(186)
        p(206) = p(1)*p(81) - p(193)
        p(207) = p(2)*p(82) - p(205) - p(203)

        return
      end subroutine evpoly

