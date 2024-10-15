module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    character*1:: bmat='I'
    character*2:: which='SM'
    real*8,parameter::ef= 5.847,a=0,TOL=0.0001,Bx=0.0
    integer*4,parameter::nblocks=2,maxiter=100000,ishift=1,mode=1,N2=nblocks**2,N3=nblocks**3
    integer*4 nb,nloc,myid,nprocs

    complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr
    integer*4,allocatable :: npminlist(:),npmaxlist(:),nloclist(:)
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    include 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,l,nr,ie,lwork,info,ik,count,ir,ir3,ir12,nr12,r1,r2,r3,sign,il,i1,j1,i2,j2,i3,j3,xindex,yindex,rvec_data(3),index,interp_size,matsize
    integer*4 IPARAM(11),IPNTR(14),IDO,LDV,LDZ
    integer*8 LWORKL,NEV,NCV,N
    integer*4 ierr,comm,rx
    integer*8 npmin,npmax,leng,iter
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,a_spec,factor,p_l,de,dos
    real*8,allocatable:: rvec(:,:),rwork(:),real_d(:)
    integer*4,allocatable:: ndeg(:),vec_ind(:,:)
    complex*16,allocatable::top_Hr(:,:),triv_Hr(:,:),super_H(:,:),surface_vec(:),B_pt(:,:)
    complex*16,allocatable::RESID(:),V(:,:),WORKD(:),WORKL(:),D(:),WORKEV(:),Z(:,:),extrarow(:,:),extracol(:,:)
    complex*16 SIGMA,b_sigma(2,2)
    logical:: rvecmat
    logical,allocatable:: select(:)
!----Date and Time
    integer,dimension(8)::time_start
    integer,dimension(8)::time_end

    call date_and_time(VALUES=time_start)

    pi2=4.0d0*atan(1.0d0)*2.0d0

    write(top_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
    write(triv_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
    write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
    open(98,file=trim(adjustl(nnkp)))
111 read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,*)bvec
    open(99,file=trim(adjustl(top_file)))
    open(97,file=trim(adjustl(triv_file)))

    open(150, file='data/parpack_eigenvalues.dat')
    open(200, file='data/parpack_eigenvectors.dat')

!------read H(R)
    interp_size=6
    if(abs(nblocks) > interp_size) interp_size = abs(nblocks)
    matsize=(nblocks)**3

    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),top_Hr(nb,nb),triv_Hr(nb,nb),ndeg(nr))
    allocate(interp_Hr(nb,nb,-interp_size:interp_size,-interp_size:interp_size,-interp_size:interp_size))
    read(99,*)ndeg
    do i=1,80
      read(97,*)
    enddo

    allocate(B_pt(nb,nb))

    !B along X-axis
    ! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(Bx,0d0)]
    ! B_sigma(2,:) = [dcmplx(Bx,0d0) ,  dcmplx(0d0,0d0)]

    !B along Y axis
	! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-Bx)]
    ! B_sigma(2,:) = [dcmplx(0d0,Bx) ,  dcmplx(0d0,0d0)]

    !B along Z-axis
    B_sigma(1,:) = [dcmplx(Bx,0d0),  dcmplx(0d0,0d0)]
    B_sigma(2,:) = [dcmplx(0d0,0d0) ,  dcmplx(-Bx,0d0)]

	B_pt=0d0
	do i=1,nb
		do j=1,nb
			if (i==j) then
				if (i<10) then
					B_pt(i,j) = B_sigma(1,1)
				else
					B_pt(i,j) = B_sigma(2,2)
				endif
			else if (i==j+9) then
				B_pt(i,j) = B_sigma(2,1)
			else if (j==i+9) then
				B_pt(i,j) = B_sigma(1,2)
			endif
		enddo
	enddo
    do ir=1,nr
        do i=1,nb
            do j=1,nb
               read(99,*)rvec_data(1),rvec_data(2),rvec_data(3),i1,i2,x1,y1
               top_Hr(i1,i2)=dcmplx(x1,y1)
               read(97,*)rvec_data(1),rvec_data(2),rvec_data(3),j1,j2,x2,y2
               triv_Hr(j1,j2)=dcmplx(x2,y2)

               interp_Hr(i1,i2,rvec_data(1),rvec_data(2),rvec_data(3))=(1-a)*triv_Hr(i1,i2) + a*top_Hr(i1,i2) + B_pt(i1,i2)
            enddo
        enddo
        rvec(:,ir) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2)
        
        ! print *, ir, nr
    enddo
    deallocate(rvec,top_Hr,triv_Hr,ndeg)

! !------ARPACK

    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK(comm, myid, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)

    allocate(npminlist(nprocs),npmaxlist(nprocs),nloclist(nprocs))
 
    N=nb*matsize
    NEV=50
    NCV=100

    do i=1,nprocs
        nloc = matsize/nprocs
        if (mod(matsize,nprocs).ne.0) nloc= nloc +1 
        npminlist(i)=1+(i-1)*nloc
        npmaxlist(i)=min(i*nloc,matsize)
        nloclist(i) = (npmaxlist(i)-npminlist(i)+1)
    enddo

    nloc = nloclist(myid+1)

    if(myid.eq.0) print *, "nloc:", nloclist
    if(myid.eq.0) print *, "npmin:", npminlist
    if(myid.eq.0) print *, "npmax:", npmaxlist

    allocate(RESID(nloc*nb),V(nloc*nb,NCV),WORKD(nloc*nb*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    allocate(select(NCV),D(NEV),WORKEV(2*NCV))

    iparam(1)=ishift
    iparam(3)=maxiter
    iparam(7)=mode
    lworkl=3*(NCV**2) + 5*NCV
    iter=0
    IDO=0
    INFO=0
    LDV=nloc*nb
    WORKL=0d0
    WORKD=0d0
    RWORK=0d0
    rvecmat=.true.
    select=.true.


    do while (iter<maxiter)
        iter=iter+1
        if(myid.eq.0) print *, 'Iterations: ', iter
        call pznaupd(comm,IDO,bmat,nloc*nb,which,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)
        
        if(IDO==99) exit
        
        if(IDO==-1 .or. IDO==1) then
            !WORKD(IPNTR(2):IPNTR(2)+N-1) = matmul(super_H,WORKD(IPNTR(1):IPNTR(1)+N-1))
            call matmul_(WORKD(IPNTR(1):IPNTR(1)+(nloc*nb)-1), WORKD(IPNTR(2):IPNTR(2)+(nloc*nb)-1))
            ! call matmul_(interp_Hr, WORKD(IPNTR(1):IPNTR(1)+N-1), WORKD(IPNTR(2):IPNTR(2)+N-1),N,nblocks)
            ! print *, "input: ", WORKD(IPNTR(1)+2), "output", WORKD(IPNTR(2)+2)
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
        print *, "Finished iterations, calling zneupd..."
        call pzneupd (comm,rvecmat, 'A', select, d, v, ldv, sigma,&
             workev, bmat, n, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl, &
             rwork, info)
        ! print*, v(1,:)
    endif

    deallocate(RESID,WORKD,WORKL,RWORK)


    ! if(myid.eq.0) then
    !     do i=1,10
    !         real_d(i) = real(d(i))
    !     enddo
    !     do i=1,10
    !         write(150, '(1(1x,f12.6))') real_d(i)
    !         write(200, *) v(:,i)
    !     enddo
    ! endif

!     allocate(surface_vec(4*nb*(nblocks-1)),vec_ind(matsize,3))

!     print *, 'Total iterations: ', iparam(3)

!     do i=1,NEV
!         real_d(i) = real(d(i))
!     enddo
!     do i=1,NEV
!         write(150, '(1(1x,f12.10))') real_d(i)
!         write(200, *) v(:,i)
!     enddo

    call date_and_time(VALUES=time_end)

!     print *, 'Start time: ', time_start(5), ':', time_start(6), ':', time_start(7)
!     print *, 'End time: ', time_end(5), ':', time_end(6), ':', time_end(7)

    call MPI_FINALIZE(ierr)

    contains
    subroutine matmul_(vec_in, vec_out)
        use mpi
        use parameters
        implicit none
        ! Declare intent and input/output parameters
        complex*16, intent(in) :: vec_in(nloc*nb)  ! Global input vector
        complex*16, intent(out) :: vec_out(nloc*nb) ! Output vector for this process
        integer*4 :: next, prev, status(MPI_STATUS_SIZE)
        integer*4 :: irow, icol, r1, r2, r3, f1, f2, f3, rowcount, colcount, reqsend, reqrec
        complex*16 :: mv_buf(nloclist(1)*nb)
        integer :: ierr, i_start, i_end
        integer*8 npmin,npmax

        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs, ierr )
    
        ! call mpi_barrier(comm,ierr)
        call tv(myid,vec_in,vec_out)

        do i=1,nprocs-1
            next = mod(myid + i, nprocs)
            prev = mod(myid - i + nprocs, nprocs)

            ! print *, "myid=", myid, "received from prev=", prev
            ! print *, "myid=", myid, " sending to next=", next, "receiving from prev=", prev
            call mpi_sendrecv(vec_in(1), nloclist(myid+1)*nb, MPI_DOUBLE_COMPLEX, next, 0, &
                              mv_buf(1), nloclist(prev+1)*nb, MPI_DOUBLE_COMPLEX, prev, 0, comm, status, ierr)

            call tv(prev,mv_buf(1:nloclist(prev)*nb),vec_out)

        enddo
    
    end subroutine matmul_

    subroutine tv(id,vec_in,vec_out)
        use parameters
        implicit none
        ! Declare intent and input/output parameters
        integer*4,intent(in):: id
        complex*16, intent(in) :: vec_in(nloclist(id+1)*nb)
        complex*16, intent(out) :: vec_out(nloclist(myid+1)*nb)
        integer*4 :: irow, icol, r1, r2, r3, f1, f2, f3
        integer :: rowcount, colcount

        rowcount=0
        colcount=0
        do irow = npminlist(myid+1),npmaxlist(myid+1)
            f3 = (irow-1) / (N2)
            f2 = mod((irow-1) / nblocks, 3)
            f1 = mod(irow-1, 3)
            do icol = npminlist(id+1),npmaxlist(id+1)
                r3 = ((icol-1) / N2) - f3
                r2 = mod((icol-1) / nblocks, 3) - f2
                r1 = mod(icol-1, 3) - f1
                if ((abs(r1) .lt. 6) .or. (abs(r2) .lt. 6) .or. (abs(r3) .lt. 6)) then
                    vec_out(1+colcount*nb : nb*(colcount+1)) = vec_out(1+colcount*nb : nb*(colcount+1)) + &
                        matmul(interp_Hr(:,:,r1,r2,r3), vec_in(1+rowcount*nb : nb*(rowcount+1)))
                endif
                colcount = colcount + 1
            enddo
            rowcount = rowcount + 1
        enddo

    end subroutine tv

end Program Projected_band_structure
