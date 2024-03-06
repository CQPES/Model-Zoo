!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
       integer,parameter::nbasis=1332,ndim=21,natom=7
       real*8,parameter::vpesmin=-116.1036487912d0
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nscale
       real*8,allocatable:: pesnn(:),pes(:),delta1(:)
       real*8 minbs(1331),maxbs(1331)
       integer, allocatable::nodes(:)
       real*8, allocatable::weighta(:,:,:),biasa(:,:)
       real*8, allocatable::pdela(:),pavga(:)

       end module nnparam

       subroutine hch3ohpipNN(ct,vpes)
       use nnparam
       use bemsa
       implicit none
       integer i,j,k,lpes
       real*8 rb(ndim),xbond(ndim),basis(nbasis),tmp1,txinput(nbasis-1)
       real*8 ct(3,natom),xvec(3,ndim),cx(3,natom),cy(3,natom)
       real*8 vpes,vpesa,vpesb,vpesc,vpesH,vnn,threshold
       basis=0.d0
       threshold=1.0d-8
       cx=ct

       k=0
       do i=1,natom-1
        do j=i+1,natom
         k=k+1
!        write(*,*)k,i,j
         xvec(:,k)=cx(:,i)-cx(:,j)
        enddo
       enddo
       if(k.ne.ndim)stop 'error in bond dimension'
       do i=1,ndim
        rb(i)=dsqrt(dot_product(xvec(:,i),xvec(:,i)))
       enddo

       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav(xbond,basis)

       do j=1,nbasis-1
        txinput(j)=basis(j+1)
       enddo

       call getpota(txinput,vpes)
!       vnn=vpes
!       vpes=vnn

       return

       end subroutine hch3ohpipNN
!*************************************************************
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

        nfile=7
        open(nfile,file='weights.txt')
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
     $  biasa(nodemax,2:nlayer))
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

        close(nfile)
        close(4)

        return
        end subroutine pes_init 

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        integer j,k,neu1,neu2,ndriv
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
!       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*
     $  weighta(inode2,inode1,ilay1)
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
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*
     $ weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)

        return
        end subroutine getpota


        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
!c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
!c    ifunc=2, transfer function is log sigmoid function, 'logsig'
!c    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if (ifunc.eq.1) then
        tranfun=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+dexp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end

