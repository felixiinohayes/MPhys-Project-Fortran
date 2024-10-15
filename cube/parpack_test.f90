module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="../BiTeI"
    character*1:: bmat='I'
    character*2:: which='LM'
    real*8,parameter::ef= 5.847,a=0,TOL=0.1,Bx=0.0
    integer*4,parameter::nblocks=3,maxiter=100000,ishift=1,mode=1,N2=nblocks**2,N3=nblocks**3
    integer*4 nb,nloc,myid,nprocs

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
    complex*16,allocatable::RESID(:),Vloc(:,:),WORKD(:),WORKL(:),D(:),WORKEV(:),Z(:,:),extrarow(:,:),extracol(:,:)
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


! !------ARPACK

    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK(comm, myid, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)

    allocate(npminlist(nprocs),npmaxlist(nprocs),nloclist(nprocs))
 
    N=8
    NEV=N-2
    NCV=N

    do i=1,nprocs
        nloc = N/nprocs
        if (mod(N,nprocs).ne.0) nloc= nloc +1 
        npminlist(i)=1+(i-1)*nloc
        npmaxlist(i)=min(i*nloc,n)
        nloclist(i) = (npmaxlist(i)-npminlist(i)+1)
    enddo

    nloc = nloclist(myid+1)

    if(myid.eq.0) print *, "nloc:", nloclist
    if(myid.eq.0) print *, "npmin:", npminlist
    if(myid.eq.0) print *, "npmax:", npmaxlist

    allocate(RESID(N),Vloc(nloc,NCV),WORKD(N*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    allocate(select(NCV),D(NEV+1),WORKEV(2*NCV))

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
    rvecmat=.true.
    select=.true.


    do while (iter<maxiter)
        iter=iter+1
        if(myid.eq.0) print *, 'Iterations: ', iter
        call pznaupd(comm,IDO,bmat,nloc,which,NEV,TOL,RESID,NCV,Vloc,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)
        
        if(IDO==99) exit
        
        if(IDO==-1 .or. IDO==1) then
            call matmul_(comm,WORKD(IPNTR(1)), WORKD(IPNTR(2)))
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
        call pzneupd (comm,rvecmat, 'A', select, d, vloc, ldv, sigma,&
             workev, bmat, nloc, which, nev, tol, resid, ncv,&
             vloc, ldv, iparam, ipntr, workd, workl, lworkl, &
             rwork, info)
    endif

    deallocate(RESID,WORKD,WORKL,RWORK)
    
    if(myid.eq.1) then
        print *, d
    endif
    !     do i=1,NEV
    !         real_d(i) = real(d(i))
    !         print *, real_d(i)
    !     enddo
    !     do i=1,NEV
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
    subroutine matmul_(comm,vec_in, vec_out)
        use mpi
        use parameters
        implicit none
        ! Declare intent and input/output parameters
        complex*16, intent(in) :: vec_in(nloc)  ! Global input vector
        complex*16, intent(out) :: vec_out(nloc) ! Output vector for this process
        integer*4 :: comm, next, prev, status(MPI_STATUS_SIZE)
        complex*16 :: mv_buf(nloclist(1))


        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs, ierr )
    
        call tv(myid,vec_in(1),vec_out(1))

        do i=1,nprocs-1
            next = mod(myid + i, nprocs)
            prev = mod(myid - i + nprocs, nprocs)

            ! print *, "myid=", myid, "received from prev=", prev
            ! print *, "myid=", myid, " sending to next=", next, "receiving from prev=", prev
            call mpi_sendrecv(vec_in(1), nloclist(myid+1), MPI_DOUBLE_COMPLEX, next, 0, &
                              mv_buf(1), nloclist(prev+1), MPI_DOUBLE_COMPLEX, prev, 0, comm, status, ierr)

            ! if(i==2 .and. myid.eq.0 .and. iter==1) print *, vec_in(nloclist(next+1)*nb)
            ! if(i==2 .and. myid.eq.2 .and. iter==1) print *, mv_buf(nloc*nb)
            call tv(prev,mv_buf(1),vec_out(1))
        enddo
    
    end subroutine matmul_

    subroutine tv(id,vec_in,vec_out)
        use parameters
        implicit none
        ! Declare intent and input/output parameters
        integer*4,intent(in):: id
        complex*16, intent(in) :: vec_in(nloclist(id+1))
        complex*16, intent(out) :: vec_out(nloclist(myid+1))
        integer*4 :: irow, icol, r1, r2, r3, f1, f2, f3
        integer :: rowcount, colcount

        rowcount=0
        colcount=0
        do irow = npminlist(myid+1),npmaxlist(myid+1)
            do icol = npminlist(id+1),npmaxlist(id+1)
                if (irow==icol) then
                    vec_out(1+colcount) = vec_out(1+colcount) + vec_in(1+rowcount)*5
                endif
                colcount = colcount + 1
            enddo
            rowcount = rowcount + 1
        enddo
        ! if(myid.eq.0 .and. iter==1) print *, id, vec_out(1)

        ! print *, myid, vec_in(1), vec_out(1)
    end subroutine tv

end Program Projected_band_structure
