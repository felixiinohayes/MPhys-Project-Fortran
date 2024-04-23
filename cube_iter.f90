module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    character*1:: bmat='I'
    character*2:: which='LM'
    real*8,parameter::ef= 4.18903772,a=1,emin=5.5,emax=6.5,bfactor=0.006,TOL=0.01
    integer*8,parameter::nblocks=3,matsize=(nblocks)**3,maxiter=1,NEV=matsize*18-2,NCV=NEV+2,ishift=1,mode=1
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    !INCLUDE 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,ie,lwork,info,ik,count,ir,ir3,ir12,nr12,r1,r2,r3,sign,il,i1,j1,i2,j2,i3,j3,xindex,yindex,rvec_data(3)
    integer*4 IPARAM(11),IPNTR(14),iter,IDO,LWORKL,LDV,LDZ,N
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), B=0.00d0
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2
    real*8,allocatable:: rvec(:,:),rvec_miller(:,:),rwork(:),ene(:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::top_Hr(:,:),triv_Hr(:,:),work(:),super_H(:,:)
    complex*16,allocatable::RESID(:),V(:,:),WORKD(:),WORKL(:),D(:),WORKEV(:),Z(:,:)
    complex*16,dimension(18,18,-6:6,-6:6,-6:6) :: interp_Hr
    complex*16 B_sigma(2,2),SIGMA
    logical:: rvecmat
    logical,allocatable:: select(:)
!------------------------------------------------------
    !call init_mpi

    write(top_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
    write(triv_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
    write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

    pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
    open(98,file=trim(adjustl(nnkp)))
111 read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,*)bvec

!------read H(R)
    open(99,file=trim(adjustl(top_file)))
    open(97,file=trim(adjustl(triv_file)))
    open(100,file='DOS.dat')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),rvec_miller(3,nr),top_Hr(nb,nb),triv_Hr(nb,nb),ndeg(nr))
    allocate(super_H(nb*matsize,nb*matsize),ene(nb*matsize))
    read(99,*)ndeg
    ! print *, blocksize

    do i=1,80
      read(97,*)
    enddo
    do ir=1,nr
        do i=1,nb
            do j=1,nb
               read(99,*)rvec_data(1),rvec_data(2),rvec_data(3),i1,i2,x1,y1
               top_Hr(i1,i2)=dcmplx(x1,y1)
               read(97,*)rvec_data(1),rvec_data(2),rvec_data(3),j1,j2,x2,y2
               triv_Hr(j1,j2)=dcmplx(x2,y2)

               interp_Hr(i1,i2,rvec_data(1),rvec_data(2),rvec_data(3))=(1-a)*triv_Hr(i1,i2) + a*top_Hr(i1,i2)
            enddo
        enddo
        rvec(:,ir) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2)
    enddo

!----- Construct supercell hamiltonian
    super_H=0d0
    do i3=0,nblocks-1
        do j3=0,nblocks-1
            r3=i3-j3
            do i2=0,nblocks-1
                do j2=0,nblocks-1
                    r2=i2-j2
                    do i1=0,nblocks-1
                        do j1=0,nblocks-1
                            r1=i1-j1
                            xindex = i3*((nblocks)**2)+i2*(nblocks)+i1
                            yindex = j3*((nblocks)**2)+j2*(nblocks)+j1
                            super_H((1+nb*xindex):(nb*(xindex+1)),(1+nb*yindex):(nb*(yindex+1))) = interp_Hr(:,:,r1,r2,r3)
                            ! print *, xindex, yindex,r1,r2,r3
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    ! call zheev('V','U',nb*blocksize,super_H,nb*blocksize,ene,work,lwork,rwork,info)
    N=nb*matsize
    allocate(RESID(N),V(N,NCV),WORKD(N*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    allocate(select(NCV),D(NEV),Z(N,NEV),WORKEV(2*NCV))
    ! print *, super_H
    ! print *, WORKD

    iparam(1)=ishift
    iparam(3)=maxiter
    iparam(7)=mode
    lworkl=3*(NCV**2) + 5*NCV

    iter=0
    IDO=0
    INFO=0
    LDV=N
    WORKL=0d0
    WORKD=0d0
    RWORK=0d0
    
    select=.true.
    rvecmat=.true.

    do while (iter<maxiter)

        iter=iter+1
        ! print *, iter
        call znaupd(IDO,bmat,N,which,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)

        if(IDO==99) exit

        if(IDO==-1 .or. IDO==1) then
            WORKD(IPNTR(2):IPNTR(2)+N-1) = matmul(super_H,WORKD(IPNTR(1):IPNTR(1)+N-1))
            ! call matmul_chunk(interp_Hr, WORKD(IPNTR(1):IPNTR(1)+N-1), WORKD(IPNTR(2):IPNTR(2)+N-1), N)
            print *, "input: ", WORKD(IPNTR(1)+2), "output", WORKD(IPNTR(2)+2)
            continue
        endif
    enddo

    if ( info .lt. 0 ) then

        print *, ' '
        print *, ' Error with _naupd, info = ', info
        print *, ' Check the documentation of _naupd'
        print *, ' '
    
    else
    
        rvecmat = .true.
    
        call zneupd (rvecmat, 'A', select, d, v, ldv, sigma,&
             workev, bmat, n, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl, &
             rwork, info)
    endif

    write(100, '(f12.6)') real(d)

    ! print *, v

    ! do i=1,nblocks

    contains

        subroutine matmul_chunk(interp_Hr,vec_in,vec_out,N)  
            integer*4,intent(in)::N
            complex*16,dimension(18,18,-6:6,-6:6,-6:6), intent(in):: interp_Hr
            complex*16,intent(in):: vec_in(N*3)
            complex*16,intent(out)::vec_out(N*3)
            complex*16::tempvec(N)

            tempvec=0d0
            do i3=0,nblocks-1
                do j3=0,nblocks-1
                    r3=i3-j3
                    do i2=0,nblocks-1
                        do j2=0,nblocks-1
                            r2=i2-j2
                            do i1=0,nblocks-1
                                do j1=0,nblocks-1
                                    r1=i1-j1
                                    xindex = i3*((nblocks)**2)+i2*(nblocks)+i1
                                    yindex = j3*((nblocks)**2)+j2*(nblocks)+j1
                                    tempvec((1+nb*yindex):(nb*(yindex+1))) = tempvec((1+nb*yindex):(nb*(yindex+1))) + matmul(interp_Hr(:,:,r1,r2,r3),vec_in((1+nb*xindex):(nb*(xindex+1))))
                                    print *, (1+nb*yindex),(nb*(yindex+1)) 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            vec_out=tempvec

        end subroutine matmul_chunk

end Program Projected_band_structure
