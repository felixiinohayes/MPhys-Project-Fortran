Program Projected_band_structure
    Implicit None
    !INCLUDE 'mpif.h'
!------------------------------------------------------
    integer*4 i,j
    integer*4 IPARAM(11),IPNTR(14),iter,IDO,NCV,NEV,N,ishift,LDV,lworkl,maxiter,mode
    real*8 TOL,INFO
    real*8 testmatrix(5,5)
    real*8,allocatable::RWORK(:)
    character*1 :: bmat='I'
    character*2 :: which='LM'
    complex*16,allocatable::RESID(:),V(:,:),WORKD(:),WORKL(:)
!------------------------------------------------------
!-----Create matrix
    do i=1,5
        do j=1,5
        testmatrix(i,j)=i + j
        enddo
    enddo
    print *, testmatrix

!-----Set parameters
    maxiter=100
    mode=1
    N=size(testmatrix)
    ishift=1
    ipntr=0
    NEV=2
    NCV=4
    iter=0
    IDO=0
    TOL=0.1d0
    INFO=0
    LDV=N
    lworkl=3*(NCV**2) + 5*NCV
    iparam(1)=ishift
    iparam(3)=maxiter
    iparam(7)=mode

    allocate(RESID(N),V(N,NCV),WORKD(N*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    RESID=0d0

    WORKL=0d0
    WORKD=0d0
    RWORK=0d0

!-------Perform iterations 
    do while (IDO==0 .or. IDO==-1 .or. IDO==1)
        call znaupd(IDO,bmat,N,which,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,lworkl,RWORK,INFO)
        print *, "done"

        if(IDO==99) exit
    enddo

end Program Projected_band_structure
