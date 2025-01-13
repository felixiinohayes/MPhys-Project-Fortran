module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    character*1:: bmat='I'
    character*2:: which='SM'
    real*8,parameter::ef_triv=4.23,ef_top=6.5,a=1,TOL=0.01,B=0.00,emin=-0.3,emax=0.3,eta=0.02
    integer*4,parameter::nxblocks=10,nyblocks=10,nzblocks=10,maxiter=100000,N3=nxblocks*nyblocks*nzblocks,Nxy=nxblocks*nyblocks
    integer*4,parameter::NEV=100,NCV=200,eres=100
    integer*4 nb,nloc,myid,nprocs
    complex*16,dimension(:,:,:,:,:),allocatable::interp_Hr
    integer*4,allocatable::npminlist(:),npmaxlist(:),nloclist(:),nloc_sum(:),nev_sum(:),nloclist_nev(:),displs(:)
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    include 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line,v_file,d_file
    character(len=5) suffix
    integer*4 i, j, k, l, nr_top, nr_triv, ie, lwork, info, ik, count, ir, ir3, ir12, nr12, r1, r2, r3, sign, il, i1, j1, i2, j2, i3, j3, xindex, yindex, rvec(3), index, interp_size
    integer*4 IPARAM(11), IPNTR(14), IDO, LDV, LDZ
    integer*8 LWORKL, N, ishift
    integer*4 ierr, comm
    integer*4 npmin, npmax, iter
    real*8 avec(3,3), bvec(3,3), pi2, x1, x2, y1, y2, a_spec, factor, p_l, de, dos, epoints(eres)
    real*8, allocatable :: rwork(:), real_d(:)
    integer*4, allocatable :: vec_ind(:,:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:)
    complex*16, allocatable :: super_H(:,:), surface_vec(:), B_pt(:,:),top_Hr_temp(:,:),triv_Hr_temp(:,:)
    complex*16, allocatable :: RESID(:), Vloc(:,:), WORKD(:), WORKL(:), D(:), WORKEV(:), vloc_flat(:), vloc_flatg(:), Zloc(:,:), v(:,:)
    complex*16 b_sigma(2,2)
    logical :: rvecmat
    logical, allocatable :: select(:)
    complex*16,dimension(:,:,:,:,:),allocatable::triv_Hr, top_Hr
!----Date and Time
    integer, dimension(8) :: time_start
    integer, dimension(8) :: time_end

    call date_and_time(VALUES=time_start)

    pi2 = 4.0d0 * atan(1.0d0) * 2.0d0

    ! write(top_file, '(a,a)') trim(adjustl(prefix)), "_hr_topological_4band.dat"
    ! write(triv_file, '(a,a)') trim(adjustl(prefix)), "_hr_trivial_4band.dat"

    write(top_file , '(a,a)') trim(adjustl(prefix)),"_hr_topological_new.dat" 
    write(triv_file, '(a,a)') trim(adjustl(prefix)),"_hr_trivial_new.dat" 

    write(nnkp, '(a,a)') trim(adjustl(prefix)), "_new.nnkp"
    open(98, file=trim(adjustl(nnkp)))
111 read(98, '(a)') line
    if (trim(adjustl(line)) .ne. "begin real_lattice") goto 111
    read(98, *) avec
    read(98, '(a)') line
    read(98, '(a)') line
    read(98, '(a)') line
    read(98, *) bvec
    open(99, file=trim(adjustl(top_file)))
    open(97, file=trim(adjustl(triv_file)))

    ! Determine the suffix based on the value of a
    if (a == 1.0d0) then
        if (B .ne. 0d0) then
            suffix = "TOP_B"
        else
            suffix = "TOP"
        endif
    else
        if (B .ne. 0d0) then
            suffix = "TRIV_B"
        else
            suffix = "TRIV"
        endif
    endif

!------read H(R) + B_pt
    interp_size=6
    ! if((nxblocks > interp_size).or.(nyblocks > interp_size).or.(nzblocks > interp_size)) interp_size = max(max(nxblocks,nyblocks),nzblocks)
    read(99,*)
    read(99,*)nb,nr_top
    read(97,*)
    read(97,*)nb,nr_triv

    allocate(top_Hr_temp(nb,nb),triv_Hr_temp(nb,nb),ndeg_top(nr_top),ndeg_triv(nr_triv))
    allocate(rvec_top(nr_top,3))
    allocate(interp_Hr(nb,nb,-6:6, -6:6, -6:6), top_Hr(nb,nb,-6:6, -6:6, -6:6), triv_Hr(nb,nb,-6:6, -6:6, -6:6))

    read(99,*)ndeg_top
    read(97,*)ndeg_triv
    allocate(B_pt(nb,nb))

    !B along Z-axis
    B_sigma(1,:) = [dcmplx(B,0d0),  dcmplx(0d0,0d0)]
    B_sigma(2,:) = [dcmplx(0d0,0d0) ,  dcmplx(-B,0d0)]

	B_pt=0d0

	do i=1,nb
		do j=1,nb
			if (i==j) then
				if (i==1 .or. i==3) then
					B_pt(i,j) = B_sigma(1,1)
				else
					B_pt(i,j) = B_sigma(2,2)
				endif
            endif
		enddo
	enddo
    ! print *, B_pt

    interp_Hr=0d0
    do ir=1,nr_top
        do i=1,nb
            do j=1,nb
               read(99,*)rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3),i1,i2,x1,y1
               top_Hr_temp(i1,i2)=dcmplx(x1,y1)
            enddo
        enddo
        top_Hr(:,:,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = top_Hr_temp(:,:)
    enddo
    do ir=1,nr_triv
        do i=1,nb
            do j=1,nb
               read(97,*)rvec(1),rvec(2),rvec(3),i1,i2,x1,y1
               triv_Hr_temp(i1,i2)=dcmplx(x1,y1)
            enddo
        enddo
        triv_Hr(:,:,rvec(1),rvec(2),rvec(3)) = triv_Hr_temp
    enddo

    ! Interpolate Hamiltonian and add magnetic field
    do ir=1,nr_top
        do i=1,nb
            do j=1,nb
                interp_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = (1-a)*triv_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) + a*top_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3))
            enddo
        enddo
    enddo
    do i=1,nb
        do j=1,nb
            interp_Hr(i,j,0,0,0) = interp_Hr(i,j,0,0,0) + B_pt(i,j)
        enddo
    enddo

    do i=1,nb
        if(a==0) then
            interp_Hr(i,i,0,0,0) = interp_Hr(i,i,0,0,0) - ef_triv
        else
            interp_Hr(i,i,0,0,0) = interp_Hr(i,i,0,0,0) - ef_top
        endif
    enddo

    deallocate(top_Hr_temp, triv_Hr_temp, top_Hr, triv_Hr, rvec_top, ndeg_top, ndeg_triv)

!--- MPI nloc arrays
    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK(comm, myid, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)

    allocate(npminlist(nprocs),npmaxlist(nprocs),nloclist(nprocs),nloc_sum(nprocs+1),nev_sum(nprocs+1),nloclist_nev(nprocs))
 
    N=nb*N3

    do i=1,nprocs
        nloc = (N3)/nprocs
        if (mod(N3,nprocs).ne.0) nloc = nloc +1 
        npminlist(i)=1+(i-1)*nloc
        npmaxlist(i)=min(i*nloc,N3)
        nloclist(i) = (npmaxlist(i)-npminlist(i)+1)*nb
        nloclist_nev(i) = (npmaxlist(i)-npminlist(i)+1)*nb*nev
    enddo

    nloc_sum(1) = 1
    nev_sum(1) = 1
    do i=1,nprocs
        nloc_sum(i+1) = nloc_sum(i) + nloclist(i)
    enddo

    do i=2,nprocs+1
        nev_sum(i) = (nloc_sum(i)-1)*nev +1
    enddo
    
    
    nloc = nloclist(myid+1) 
    if(myid.eq.0) print *, "nloc:", nloclist
    if(myid.eq.0) print *, "nloc*nev:", nloclist_nev
    if(myid.eq.0) print *, "npmin:", npminlist
    if(myid.eq.0) print *, "npmax:", npmaxlist
    if(myid.eq.0) print *, "nloc_sum:", nloc_sum
    if(myid.eq.0) print *, "nev_sum:", nev_sum

!---P_ARPACK
    allocate(RESID(nloc),Vloc(nloc,NCV),Zloc(nloc,NEV),WORKD(nloc*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    allocate(select(NCV),D(NEV),real_d(NEV),WORKEV(2*NCV))

    iparam(1)=1
    iparam(3)=maxiter
    iparam(7)=1
    lworkl=3*(NCV**2) + 5*NCV
    iter=0
    IDO=0
    INFO=0
    LDV=nloc
    WORKL=0d0
    WORKD=0d0
    RWORK=0d0
    rvecmat=.true.
    select=.true.

    ! Arnoldi iterations using pznaupd
    do while (iter < maxiter)
        iter = iter + 1
        if (myid == 0) print *, 'Iterations: ', iter

        call pznaupd(comm, IDO, 'I', nloc, 'SM', NEV, TOL, RESID, NCV, Vloc, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, WORKD, INFO)

        if (IDO == 99) exit
        if (IDO == -1 .or. IDO == 1) then
            call matmul_(comm, WORKD(IPNTR(1):IPNTR(1)+nloc-1), WORKD(IPNTR(2):IPNTR(2)+nloc-1))
            continue
        endif
    end do

    ! Check for errors
    if (INFO .ne. 0) then
        if (myid == 0) print *, 'Error with pznaupd, info =', INFO
            call MPI_FINALIZE(ierr)
        stop
    end if

    rvecmat = .true.

    call pzneupd(comm, rvecmat, 'A', select, D, Zloc, LDV, (0.0d0, 0.0d0), WORKEV, 'I', nloc, 'SM', NEV, TOL, RESID, NCV, Vloc, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, WORKD, INFO)

    if (INFO .ne. 0) then
        if (myid == 0) print *, 'Error with pzneupd, info =', INFO
            call MPI_FINALIZE(ierr)
        stop
    end if

    deallocate(RESID, Vloc, WORKD, WORKL, RWORK, select, real_d, WORKEV)

    if (myid==0) allocate(vloc_flatg(N*NEV))
    allocate(displs(nprocs))
    do i=1, nprocs
        displs(i) = nev_sum(i) - 1
    enddo

    call MPI_GATHERV(Zloc, nloc*nev, MPI_DOUBLE_COMPLEX, &
                    vloc_flatg, nloclist_nev, displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    
    deallocate(Zloc)
!---- Projections
    if (myid == 0) then
        allocate(v(N,NEV))

        open(100, file='data/cube.dx')
        do i=1,nprocs
            v(nloc_sum(i):nloc_sum(i+1)-1,1:nev) = reshape(vloc_flatg(nev_sum(i):nev_sum(i+1)-1), [nloclist(i), nev])
        enddo
    
        de = (emax-emin)/eres
        do i=1, eres
            epoints(i) = emin + de*i
        enddo

        write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nzblocks,nxblocks,nyblocks
        write(100, '(a,3(1x,f12.8))') 'origin',0d0,0d0,0d0
        write(100, '(a,3(1x,f12.8))') 'delta',0d0,0d0,1d0
        write(100, '(a,3(1x,f12.8))') 'delta',0d0,1d0,0d0
        write(100, '(a,3(1x,f12.6))') 'delta',1d0,0d0,0d0
        write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nzblocks,nxblocks,nyblocks

        print *, dot_product(v(:,1),v(:,1)) 

        print *, "Calculating DOS..."
        count = 0 
        do ie=1,eres
            count = count + 1
            print*, ie, eres

            write(100, '(a,i8,a,i8,a,i10,a)') 'object',2+count,' class array type float rank 1 shape',1,&
                                    ' item', N3, ' data follows'
            !----Spectral DOS
            do j=0,N3-1
                a_spec = 0d0
                do i=1,NEV
                    p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))
                    factor = ((epoints(ie) - d(i)))/eta
                    a_spec = a_spec + p_l * (exp(-0.5d0*factor**2)) * 1/(eta*sqrt(2*pi2))

                enddo
                write(100, '(1(1x,f12.8))') a_spec
            enddo
            write(100, '(a)') 'attribute "dep" string "positions"' 
        enddo
        do i=0,eres-1
            write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
            'object',eres+3+i,' class field', &
            'component "positions" value 1', &
            'component "connections" value 2', &
            'component "data" value ',3+i
        enddo
        write(100, '(a)') 'object "series" class series'
        do i=0,eres-1
            write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+eres+3), ' position', i
        enddo


        write(100, '(A)') 'end'

        call date_and_time(VALUES=time_end)

        print *, 'Start time: ', time_start(5), ':', time_start(6), ':', time_start(7)
        print *, 'End time: ', time_end(5), ':', time_end(6), ':', time_end(7)

    endif
    ! Finalize MPI
    call MPI_FINALIZE(ierr)

    contains

    subroutine matmul_(comm,vec_in, vec_out)
        use parameters
        implicit none
        include 'mpif.h'
        complex*16, intent(in) :: vec_in(nloclist(myid+1))
        complex*16, intent(out) :: vec_out(nloclist(myid+1))
        integer*4 :: comm, next, prev, status(MPI_STATUS_SIZE)
        integer*4 :: irow, icol, rowcount, colcount, reqsend, reqrec
        complex*16 :: mv_buf(nloclist(1))

        vec_out = 0d0

        call tv(myid, vec_in, vec_out)

       do i=1,nprocs-1
            mv_buf=0d0
            next = mod(myid + i, nprocs)
            prev = mod(myid - i + nprocs, nprocs)

            call mpi_sendrecv(vec_in, nloclist(myid+1), MPI_DOUBLE_COMPLEX, next, 0, &
                              mv_buf, nloclist(prev+1), MPI_DOUBLE_COMPLEX, prev, 0, comm, status, ierr)

            call tv(prev,mv_buf(1:nloclist(prev+1)),vec_out)
        enddo

    end subroutine matmul_

    subroutine tv(id, vec_in, vec_out)
        use parameters
        implicit none
        integer*4, intent(in) :: id
        complex*16, intent(in) :: vec_in(nloclist(id+1))
        complex*16, intent(out) :: vec_out(nloclist(myid+1))
        integer*4 :: irow, icol, r1, r2, r3, f1, f2, f3
        integer :: rowcount, colcount

        rowcount = 0
        do irow = npminlist(myid+1), npmaxlist(myid+1)
            colcount = 0
            f3 = (irow-1) / (Nxy)
            f2 = mod((irow-1) / nxblocks, nyblocks)
            f1 = mod(irow-1, nxblocks)
            do icol = npminlist(id+1), npmaxlist(id+1)
                r3 = ((icol-1) / Nxy) - f3
                r2 = mod((icol-1) / nxblocks, nyblocks) - f2
                r1 = mod(icol-1, nxblocks) - f1
                if ((abs(r1) .lt. 6) .and. (abs(r2) .lt. 6) .and. (abs(r3) .lt. 6)) then
                    vec_out(1+rowcount*nb : nb*(rowcount+1)) = vec_out(1+rowcount*nb : nb*(rowcount+1)) + matmul(interp_Hr(:,:,r1,r2,r3) + B_pt, vec_in(1+colcount*nb : nb*(colcount+1)))
                endif
                colcount = colcount + 1
            enddo
            rowcount = rowcount + 1
        enddo

    end subroutine tv

end program 
