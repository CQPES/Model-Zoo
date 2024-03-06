!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
       integer,parameter::nbasis=83,ndim=10,natom=5
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nterm(1:nbasis),nindex(1:nbasis,1:100,1:ndim)
       integer nscale
       integer, allocatable::nodes(:)

       real*8, allocatable::weighta(:,:,:),biasa(:,:)
       real*8, allocatable::pdela(:),pavga(:)
       real*8, allocatable::weightb(:,:,:),biasb(:,:)
       real*8, allocatable::pdelb(:),pavgb(:)
       real*8, allocatable::weightc(:,:,:),biasc(:,:)
       real*8, allocatable::pdelc(:),pavgc(:)

       end module nnparam

       subroutine ch4pipNN(ct,vpes,vpesa,vpesb,vpesc) !ct in HHHHC order and angstrom
       use nnparam
       implicit none
       integer i,j,k
       real*8 rb(ndim),xbond(ndim),basis(nbasis),tmp1,txinput(nbasis-1)
       real*8 ct(3,natom),xvec(3,ndim),cx(3,natom)
       real*8 vpes,vpesa,vpesb,vpesc

       basis=0.d0

       cx=ct

       k=0
       do i=1,natom-1
        do j=i+1,natom
         k=k+1
!        write(*,*)k,i,j
         xvec(:,k)=cx(:,i)-cx(:,j)
        enddo
       enddo
       if(k.ne.10)stop 'error in bond dimension'
       do i=1,ndim
        rb(i)=dsqrt(dot_product(xvec(:,i),xvec(:,i)))
       enddo

       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav(xbond,basis)

       do j=1,nbasis-1
        txinput(j)=basis(j+1)
       enddo

       call getpota(txinput,vpesa)
       call getpotb(txinput,vpesb)
       call getpotc(txinput,vpesc)
       vpes=(vpesa+vpesb+vpesc)/3.0d0

       return

       end subroutine ch4pipNN

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling the getpot() subroutine
      subroutine pes_init
      use nnparam
      implicit none
      integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile
      character f1*80,line*80

      open(4,file='biases.txt',status='old')
      read(4,'(a80)')line
!       write(*,*)line

        nfile=7
        open(nfile,file='weights.txt')
        read(nfile,'(a80)')line
!       write(*,*)line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(nodes(nlayer),pdela(nscale),pavga(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weighta(nodemax,nodemax,2:nlayer),
     %   biasa(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdela(i),i=1,nscale)
        read(nfile,*)(pavga(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weighta(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasa(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif

        read(4,'(a80)')line
        ! write(*,*)line
        read(nfile,'(a80)')line
        ! write(*,*)line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelb(nscale),pavgb(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightb(nodemax,nodemax,2:nlayer),
     %   biasb(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelb(i),i=1,nscale)
        read(nfile,*)(pavgb(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightb(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasb(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif


        read(4,'(a80)')line
!       write(*,*)line
        read(nfile,'(a80)')line
!       write(*,*)line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelc(nscale),pavgc(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightc(nodemax,nodemax,2:nlayer),
     % biasc(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelc(i),i=1,nscale)
        read(nfile,*)(pavgc(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightc(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasc(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        close(nfile)
        close(4)

        return

        end subroutine pes_init

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)
        return
        end subroutine getpota

        subroutine getpotb(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavgb(i))/pdelb(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightb(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightb(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelb(nscale)+pavgb(nscale)
        return
        end subroutine getpotb

        subroutine getpotc(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
!       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavgc(i))/pdelc(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightc(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightc(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)
        return
        end subroutine getpotc

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
        end

      function emsav(x,c) result(v)
      implicit none
      real*8,dimension(1:10)::x
      real*8,dimension(0:82)::c
      real*8::v
      ! ::::::::::::::::::::
      real*8,dimension(0:82)::p
      call bemsav(x,p)
      v = dot_product(p,c)
      return
      end function emsav
 
      subroutine bemsav(x,p)
      implicit none
      real*8,dimension(1:10),intent(in)::x
      real*8,dimension(0:82),intent(out)::p
      ! ::::::::::::::::::::
      real*8,dimension(0:261)::m
      call evmono(x,m)
      call evpoly(m,p)
      return
      end subroutine bemsav
 
      subroutine evmono(x,m)
      implicit none
      real*8,dimension(1:10),intent(in)::x
      real*8,dimension(0:261),intent(out)::m
 
      m(0)=1.D0
      m(1)=x(10)
      m(2)=x(9)
      m(3)=x(7)
      m(4)=x(4)
      m(5)=x(8)
      m(6)=x(6)
      m(7)=x(5)
      m(8)=x(3)
      m(9)=x(2)
      m(10)=x(1)
      m(11)=m(1)*m(2)
      m(12)=m(1)*m(3)
      m(13)=m(2)*m(3)
      m(14)=m(1)*m(4)
      m(15)=m(2)*m(4)
      m(16)=m(3)*m(4)
      m(17)=m(3)*m(5)
      m(18)=m(2)*m(6)
      m(19)=m(1)*m(7)
      m(20)=m(4)*m(5)
      m(21)=m(4)*m(6)
      m(22)=m(4)*m(7)
      m(23)=m(2)*m(8)
      m(24)=m(3)*m(8)
      m(25)=m(1)*m(9)
      m(26)=m(3)*m(9)
      m(27)=m(1)*m(10)
      m(28)=m(2)*m(10)
      m(29)=m(7)*m(8)
      m(30)=m(6)*m(9)
      m(31)=m(5)*m(10)
      m(32)=m(5)*m(6)
      m(33)=m(5)*m(7)
      m(34)=m(6)*m(7)
      m(35)=m(5)*m(8)
      m(36)=m(6)*m(8)
      m(37)=m(5)*m(9)
      m(38)=m(7)*m(9)
      m(39)=m(8)*m(9)
      m(40)=m(6)*m(10)
      m(41)=m(7)*m(10)
      m(42)=m(8)*m(10)
      m(43)=m(9)*m(10)
      m(44)=m(1)*m(13)
      m(45)=m(1)*m(15)
      m(46)=m(1)*m(16)
      m(47)=m(2)*m(16)
      m(48)=m(3)*m(20)
      m(49)=m(2)*m(21)
      m(50)=m(1)*m(22)
      m(51)=m(2)*m(24)
      m(52)=m(1)*m(26)
      m(53)=m(1)*m(28)
      m(54)=m(1)*m(17)
      m(55)=m(2)*m(17)
      m(56)=m(1)*m(18)
      m(57)=m(3)*m(18)
      m(58)=m(2)*m(19)
      m(59)=m(3)*m(19)
      m(60)=m(1)*m(20)
      m(61)=m(2)*m(20)
      m(62)=m(1)*m(21)
      m(63)=m(3)*m(21)
      m(64)=m(2)*m(22)
      m(65)=m(3)*m(22)
      m(66)=m(1)*m(23)
      m(67)=m(1)*m(24)
      m(68)=m(4)*m(23)
      m(69)=m(4)*m(24)
      m(70)=m(2)*m(25)
      m(71)=m(2)*m(26)
      m(72)=m(4)*m(25)
      m(73)=m(4)*m(26)
      m(74)=m(3)*m(27)
      m(75)=m(3)*m(28)
      m(76)=m(4)*m(27)
      m(77)=m(4)*m(28)
      m(78)=m(4)*m(32)
      m(79)=m(4)*m(33)
      m(80)=m(4)*m(34)
      m(81)=m(3)*m(35)
      m(82)=m(2)*m(36)
      m(83)=m(3)*m(37)
      m(84)=m(1)*m(38)
      m(85)=m(3)*m(39)
      m(86)=m(2)*m(40)
      m(87)=m(1)*m(41)
      m(88)=m(2)*m(42)
      m(89)=m(1)*m(43)
      m(90)=m(2)*m(32)
      m(91)=m(3)*m(32)
      m(92)=m(1)*m(33)
      m(93)=m(3)*m(33)
      m(94)=m(1)*m(34)
      m(95)=m(2)*m(34)
      m(96)=m(2)*m(35)
      m(97)=m(3)*m(36)
      m(98)=m(4)*m(35)
      m(99)=m(4)*m(36)
      m(100)=m(1)*m(37)
      m(101)=m(3)*m(38)
      m(102)=m(4)*m(37)
      m(103)=m(4)*m(38)
      m(104)=m(1)*m(39)
      m(105)=m(2)*m(39)
      m(106)=m(1)*m(40)
      m(107)=m(2)*m(41)
      m(108)=m(4)*m(40)
      m(109)=m(4)*m(41)
      m(110)=m(1)*m(42)
      m(111)=m(3)*m(42)
      m(112)=m(2)*m(43)
      m(113)=m(3)*m(43)
      m(114)=m(5)*m(29)
      m(115)=m(6)*m(29)
      m(116)=m(5)*m(30)
      m(117)=m(6)*m(38)
      m(118)=m(6)*m(39)
      m(119)=m(7)*m(39)
      m(120)=m(5)*m(40)
      m(121)=m(5)*m(41)
      m(122)=m(5)*m(42)
      m(123)=m(7)*m(42)
      m(124)=m(5)*m(43)
      m(125)=m(6)*m(43)
      m(126)=m(5)*m(34)
      m(127)=m(5)*m(39)
      m(128)=m(6)*m(42)
      m(129)=m(7)*m(43)
      m(130)=m(5)*m(36)
      m(131)=m(5)*m(38)
      m(132)=m(6)*m(41)
      m(133)=m(8)*m(43)
      m(134)=m(1)*m(47)
      m(135)=m(1)*m(48)
      m(136)=m(2)*m(48)
      m(137)=m(1)*m(49)
      m(138)=m(2)*m(63)
      m(139)=m(1)*m(64)
      m(140)=m(1)*m(65)
      m(141)=m(1)*m(51)
      m(142)=m(2)*m(69)
      m(143)=m(1)*m(71)
      m(144)=m(1)*m(73)
      m(145)=m(1)*m(75)
      m(146)=m(1)*m(77)
      m(147)=m(7)*m(66)
      m(148)=m(7)*m(67)
      m(149)=m(7)*m(68)
      m(150)=m(7)*m(69)
      m(151)=m(6)*m(70)
      m(152)=m(6)*m(71)
      m(153)=m(6)*m(72)
      m(154)=m(6)*m(73)
      m(155)=m(5)*m(74)
      m(156)=m(5)*m(75)
      m(157)=m(5)*m(76)
      m(158)=m(5)*m(77)
      m(159)=m(2)*m(78)
      m(160)=m(3)*m(78)
      m(161)=m(1)*m(79)
      m(162)=m(3)*m(79)
      m(163)=m(1)*m(80)
      m(164)=m(2)*m(80)
      m(165)=m(2)*m(81)
      m(166)=m(2)*m(97)
      m(167)=m(3)*m(98)
      m(168)=m(2)*m(99)
      m(169)=m(1)*m(83)
      m(170)=m(1)*m(101)
      m(171)=m(3)*m(102)
      m(172)=m(1)*m(103)
      m(173)=m(1)*m(85)
      m(174)=m(2)*m(85)
      m(175)=m(1)*m(86)
      m(176)=m(1)*m(107)
      m(177)=m(2)*m(108)
      m(178)=m(1)*m(109)
      m(179)=m(1)*m(88)
      m(180)=m(2)*m(111)
      m(181)=m(1)*m(112)
      m(182)=m(1)*m(113)
      m(183)=m(2)*m(91)
      m(184)=m(1)*m(93)
      m(185)=m(1)*m(95)
      m(186)=m(2)*m(98)
      m(187)=m(3)*m(99)
      m(188)=m(1)*m(102)
      m(189)=m(3)*m(103)
      m(190)=m(1)*m(105)
      m(191)=m(1)*m(108)
      m(192)=m(2)*m(109)
      m(193)=m(1)*m(111)
      m(194)=m(2)*m(113)
      m(195)=m(3)*m(114)
      m(196)=m(2)*m(115)
      m(197)=m(4)*m(114)
      m(198)=m(4)*m(115)
      m(199)=m(3)*m(116)
      m(200)=m(1)*m(117)
      m(201)=m(4)*m(116)
      m(202)=m(4)*m(117)
      m(203)=m(2)*m(118)
      m(204)=m(3)*m(118)
      m(205)=m(1)*m(119)
      m(206)=m(3)*m(119)
      m(207)=m(2)*m(120)
      m(208)=m(1)*m(121)
      m(209)=m(4)*m(120)
      m(210)=m(4)*m(121)
      m(211)=m(2)*m(122)
      m(212)=m(3)*m(122)
      m(213)=m(1)*m(123)
      m(214)=m(2)*m(123)
      m(215)=m(1)*m(124)
      m(216)=m(3)*m(124)
      m(217)=m(1)*m(125)
      m(218)=m(2)*m(125)
      m(219)=m(6)*m(119)
      m(220)=m(5)*m(123)
      m(221)=m(5)*m(125)
      m(222)=m(4)*m(126)
      m(223)=m(3)*m(127)
      m(224)=m(2)*m(128)
      m(225)=m(1)*m(129)
      m(226)=m(1)*m(78)
      m(227)=m(2)*m(79)
      m(228)=m(3)*m(80)
      m(229)=m(1)*m(81)
      m(230)=m(1)*m(82)
      m(231)=m(2)*m(83)
      m(232)=m(2)*m(84)
      m(233)=m(4)*m(85)
      m(234)=m(3)*m(86)
      m(235)=m(3)*m(87)
      m(236)=m(4)*m(88)
      m(237)=m(4)*m(89)
      m(238)=m(2)*m(130)
      m(239)=m(3)*m(130)
      m(240)=m(4)*m(130)
      m(241)=m(1)*m(131)
      m(242)=m(3)*m(131)
      m(243)=m(4)*m(131)
      m(244)=m(1)*m(132)
      m(245)=m(2)*m(132)
      m(246)=m(4)*m(132)
      m(247)=m(1)*m(133)
      m(248)=m(2)*m(133)
      m(249)=m(3)*m(133)
      m(250)=m(5)*m(115)
      m(251)=m(5)*m(117)
      m(252)=m(5)*m(118)
      m(253)=m(5)*m(119)
      m(254)=m(5)*m(132)
      m(255)=m(5)*m(128)
      m(256)=m(6)*m(123)
      m(257)=m(5)*m(129)
      m(258)=m(6)*m(129)
      m(259)=m(5)*m(133)
      m(260)=m(6)*m(133)
      m(261)=m(7)*m(133)
 
      return
      end subroutine evmono
 
      subroutine evpoly(m,p)
      implicit none
      real*8,dimension(0:261),intent(in)::m
      real*8,dimension(0:82),intent(out)::p
 
      p(0)=m(0)
      p(1)=m(1)+m(2)+m(3)+m(4)
      p(2)=m(5)+m(6)+m(7)+m(8)+m(9)+m(10)
      p(3)=m(11)+m(12)+m(13)+m(14)+m(15)+m(16)
      p(4)=m(17)+m(18)+m(19)+m(20)+m(21)+m(22)+m(23)+m(24)+m(25)+m(26)+m
     &(27)+m(28)
      p(5)=m(29)+m(30)+m(31)
      p(6)=p(1)*p(2)-p(4)
      p(7)=m(32)+m(33)+m(34)+m(35)+m(36)+m(37)+m(38)+m(39)+m(40)+m(41)+m
     &(42)+m(43)
      p(8)=p(1)*p(1)-p(3)-p(3)
      p(9)=p(2)*p(2)-p(7)-p(5)-p(7)-p(5)
      p(10)=m(44)+m(45)+m(46)+m(47)
      p(11)=m(48)+m(49)+m(50)+m(51)+m(52)+m(53)
      p(12)=m(54)+m(55)+m(56)+m(57)+m(58)+m(59)+m(60)+m(61)+m(62)+m(63)+
     &m(64)+m(65)+m(66)+m(67)+m(68)+m(69)+m(70)+m(71)+m(72)+m(73)+m(74)+
     &m(75)+m(76)+m(77)
      p(13)=p(1)*p(5)
      p(14)=p(2)*p(3)-p(12)-p(11)
      p(15)=m(78)+m(79)+m(80)+m(81)+m(82)+m(83)+m(84)+m(85)+m(86)+m(87)+
     &m(88)+m(89)
      p(16)=m(90)+m(91)+m(92)+m(93)+m(94)+m(95)+m(96)+m(97)+m(98)+m(99)+
     &m(100)+m(101)+m(102)+m(103)+m(104)+m(105)+m(106)+m(107)+m(108)+m(1
     &09)+m(110)+m(111)+m(112)+m(113)
      p(17)=m(114)+m(115)+m(116)+m(117)+m(118)+m(119)+m(120)+m(121)+m(12
     &2)+m(123)+m(124)+m(125)
      p(18)=m(126)+m(127)+m(128)+m(129)
      p(19)=p(1)*p(7)-p(16)-p(15)
      p(20)=m(130)+m(131)+m(132)+m(133)
      p(21)=p(1)*p(3)-p(10)-p(10)-p(10)
      p(22)=p(1)*p(4)-p(12)-p(11)-p(11)
      p(23)=p(2)*p(8)-p(22)
      p(24)=p(2)*p(4)-p(16)-p(15)-p(13)-p(15)
      p(25)=p(2)*p(5)-p(17)
      p(26)=p(1)*p(9)-p(24)
      p(27)=p(2)*p(7)-p(18)-p(20)-p(17)-p(18)-p(20)-p(17)-p(18)-p(20)
      p(28)=p(1)*p(8)-p(21)
      p(29)=p(2)*p(9)-p(27)-p(25)
      p(30)=m(134)
      p(31)=m(135)+m(136)+m(137)+m(138)+m(139)+m(140)+m(141)+m(142)+m(14
     &3)+m(144)+m(145)+m(146)
      p(32)=m(147)+m(148)+m(149)+m(150)+m(151)+m(152)+m(153)+m(154)+m(15
     &5)+m(156)+m(157)+m(158)
      p(33)=p(2)*p(10)-p(31)
      p(34)=p(3)*p(5)-p(32)
      p(35)=m(159)+m(160)+m(161)+m(162)+m(163)+m(164)+m(165)+m(166)+m(16
     &7)+m(168)+m(169)+m(170)+m(171)+m(172)+m(173)+m(174)+m(175)+m(176)+
     &m(177)+m(178)+m(179)+m(180)+m(181)+m(182)
      p(36)=m(183)+m(184)+m(185)+m(186)+m(187)+m(188)+m(189)+m(190)+m(19
     &1)+m(192)+m(193)+m(194)
      p(37)=m(195)+m(196)+m(197)+m(198)+m(199)+m(200)+m(201)+m(202)+m(20
     &3)+m(204)+m(205)+m(206)+m(207)+m(208)+m(209)+m(210)+m(211)+m(212)+
     &m(213)+m(214)+m(215)+m(216)+m(217)+m(218)
      p(38)=m(219)+m(220)+m(221)
      p(39)=m(222)+m(223)+m(224)+m(225)
      p(40)=m(226)+m(227)+m(228)+m(229)+m(230)+m(231)+m(232)+m(233)+m(23
     &4)+m(235)+m(236)+m(237)
      p(41)=p(3)*p(7)-p(36)-p(40)-p(35)
      p(42)=p(1)*p(17)-p(37)
      p(43)=p(1)*p(18)-p(39)
      p(44)=m(238)+m(239)+m(240)+m(241)+m(242)+m(243)+m(244)+m(245)+m(24
     &6)+m(247)+m(248)+m(249)
      p(45)=m(250)+m(251)+m(252)+m(253)+m(254)+m(255)+m(256)+m(257)+m(25
     &8)+m(259)+m(260)+m(261)
      p(46)=p(1)*p(20)-p(44)
      p(47)=p(1)*p(10)-p(30)-p(30)-p(30)-p(30)
      p(48)=p(1)*p(11)-p(31)
      p(49)=p(3)*p(4)-p(33)-p(31)-p(48)-p(31)
      p(50)=p(1)*p(12)-p(33)-p(31)-p(49)-p(33)-p(31)
      p(51)=p(5)*p(8)
      p(52)=p(1)*p(14)-p(33)
      p(53)=p(1)*p(15)-p(40)-p(35)
      p(54)=p(1)*p(16)-p(41)-p(36)-p(35)-p(36)
      p(55)=p(1)*p(19)-p(41)-p(40)
      p(56)=p(2)*p(11)-p(35)-p(34)
      p(57)=p(4)*p(5)-p(37)
      p(58)=p(2)*p(12)-p(41)-p(36)-p(40)-p(35)-p(32)-p(36)-p(40)-p(32)
      p(59)=p(1)*p(25)-p(57)
      p(60)=p(2)*p(14)-p(41)-p(34)
      p(61)=p(2)*p(15)-p(39)-p(44)-p(37)-p(39)-p(39)
      p(62)=p(4)*p(7)-p(43)-p(39)-p(44)-p(42)-p(37)-p(61)-p(39)-p(44)-p(
     &39)
      p(63)=p(5)*p(7)-p(45)
      p(64)=p(2)*p(16)-p(43)-p(44)-p(42)-p(37)-p(62)-p(43)-p(44)
      p(65)=p(2)*p(17)-p(45)-p(38)-p(63)-p(45)-p(38)-p(38)-p(38)
      p(66)=p(2)*p(18)-p(45)
      p(67)=p(1)*p(27)-p(64)-p(62)-p(61)
      p(68)=p(2)*p(20)-p(45)
      p(69)=p(3)*p(3)-p(30)-p(47)-p(30)-p(47)-p(30)-p(30)-p(30)-p(30)
      p(70)=p(3)*p(8)-p(47)
      p(71)=p(1)*p(22)-p(49)-p(48)
      p(72)=p(2)*p(28)-p(71)
      p(73)=p(1)*p(24)-p(58)-p(56)-p(56)
      p(74)=p(5)*p(5)-p(38)-p(38)
      p(75)=p(8)*p(9)-p(73)
      p(76)=p(7)*p(7)-p(45)-p(38)-p(66)-p(68)-p(65)-p(45)-p(38)-p(66)-p(
     &68)-p(65)-p(45)-p(38)-p(45)-p(38)
      p(77)=p(2)*p(24)-p(62)-p(61)-p(57)
      p(78)=p(5)*p(9)-p(65)
      p(79)=p(1)*p(29)-p(77)
      p(80)=p(7)*p(9)-p(66)-p(68)-p(63)
      p(81)=p(1)*p(28)-p(70)
      p(82)=p(2)*p(29)-p(80)-p(78)
 
      return
      end subroutine evpoly
