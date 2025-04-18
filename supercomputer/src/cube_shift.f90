module parameters
    Implicit None
    character(len=80):: prefix="input/BiTeI"
    character*1:: bmat='I'
    character*2:: which='SM'
    real*8,parameter::ef_triv=4.196,ef_top=6.5,a=1,TOL=0.01,emin=-0.3,emax=0.3,eta=0.005
    integer*4,parameter::nxblocks=80,nyblocks=nxblocks,nzblocks=nxblocks,maxiter=1000000,N3=nxblocks*nyblocks*nzblocks,Nxy=nxblocks*nyblocks
    integer*4,parameter::eres=20,nshifts=eres+1,NEV=3500,NCV=2*NEV
    integer*4 nb,nloc,myid,nprocs
    complex*16,dimension(:,:,:,:,:),allocatable::interp_Hr
    integer*4,allocatable::npminlist(:),npmaxlist(:),nloclist(:),nloc_sum(:),nev_sum(:),nloclist_nev(:),displs(:)
    real*8 :: B 
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    include 'mpif.h'
!------------------------------------------------------    
    character(len=32) :: arg
    character(len=80) top_file,triv_file,nnkp,line,v_file,d_file,data_file
    character(len=5) suffix
    integer*4 i, j, k, l, nr_top, nr_triv, ie, lwork, info, ik, count, ir, ir3, ir12, nr12, r1, r2, r3, sign, il, i1, j1, i2, j2, i3, j3,is
    integer*4 IPARAM(11), IPNTR(14), IDO, LDV, LDZ, xindex, yindex, rvec(3), index, interp_size, jloc, owner
    integer*8 LWORKL, N, ishift
    integer*4 comm, ierr, num_args
    integer*4 npmin, npmax, iter
    real*8 avec(3,3), bvec(3,3), pi2, x1, x2, y1, y2, a_spec, factor, de, dos, epoints(eres+1), mem
    real*8, allocatable :: rwork(:), real_d(:)
    integer*4, allocatable :: vec_ind(:,:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:)
    complex*16, allocatable :: super_H(:,:), surface_vec(:), B_pt(:,:),top_Hr_temp(:,:),triv_Hr_temp(:,:),diag(:)
    complex*16, allocatable :: RESID(:), Vloc(:,:), WORKD(:), WORKL(:), D(:), WORKEV(:), Zloc(:,:), p_l_loc(:), p_l(:)
    complex*16 b_sigma(2,2)
    logical :: rvecmat
    logical, allocatable :: select(:)
    complex*16,dimension(:,:,:,:,:),allocatable::triv_Hr, top_Hr
!----Date and Time
    integer, dimension(8) :: time_start
    integer, dimension(8) :: time_end

!------------------------------------------------------
    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK(comm, myid, ierr)
    call MPI_COMM_SIZE(comm, nprocs, ierr)

    B = 0.0d0

    if (myid == 0) then
        num_args = command_argument_count()
        if (num_args > 0) then
            call get_command_argument(1, arg)
            read(arg, *, iostat=ierr) B
            if (ierr /= 0) then
                print *, "Error: Invalid B value. Using default B = 0.0"
                B = 0.0d0
            endif
        endif
    endif

    ! Broadcast B value to all processes
    call MPI_BCAST(B, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!------------------------------------------------------
    call date_and_time(VALUES=time_start)

    pi2 = 4.0d0 * atan(1.0d0) * 2.0d0
    write(top_file , '(a,a)') trim(adjustl(prefix)),"_hr_topological_new.dat"
    write(triv_file, '(a,a)') trim(adjustl(prefix)),"_hr_trivial_new.dat"

    write(nnkp, '(a,a)') trim(adjustl(prefix)), "_ortho.nnkp"
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

!---- Determine the suffix (and filename) based on the value of a
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
                if (mod(i,2).eq.1) then
                    B_pt(i,j) = B_sigma(1,1)
                    B_pt(i+1,j) = B_sigma(2,1)
                    B_pt(i,j+1) = B_sigma(1,2)
                else
                    B_pt(i,j) = B_sigma(2,2)
                endif
            endif
        enddo
    enddo

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

    
    print*, 'myid:', myid, 'B:', B

!----Energy points
    de = (emax-emin)/eres
    do i=1, eres+1
        epoints(i) = emin + de*(i-1)
    enddo


    if (myid == 0) then
        data_file = "data/C_"// trim(adjustl(suffix)) // ".dx"
        open(100, file=trim(adjustl(data_file)))
        write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nzblocks,nxblocks,nyblocks
        write(100, '(a,3(1x,f12.8))') 'origin',0d0,0d0,0d0
        write(100, '(a,3(1x,f12.8))') 'delta' ,0d0,0d0,1d0
        write(100, '(a,3(1x,f12.8))') 'delta' ,0d0,1d0,0d0
        write(100, '(a,3(1x,f12.8))') 'delta' ,1d0,0d0,0d0
        write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nzblocks,nxblocks,nyblocks
    endif

    allocate(npminlist(nprocs),npmaxlist(nprocs),nloclist(nprocs),nloc_sum(nprocs+1),nev_sum(nprocs+1),nloclist_nev(nprocs),diag(nb))
    allocate(select(NCV),D(NEV),real_d(NEV),WORKEV(2*NCV))

!----Shift loop
    do is = 1,nshifts

        do i=1,nb
            diag(i) = interp_Hr(i,i,0,0,0)
            interp_Hr(i,i,0,0,0) = interp_Hr(i,i,0,0,0) + epoints(is)
        enddo

        if(myid.eq.0) then
            print*, 'Shift:', is, 'Energy:', epoints(is), 'NEV:', NEV
            write(100, '(a,i8,a,i8,a,i10,a)') 'object',2+is,' class array type float rank 1 shape',1,&
                                ' item', N3, ' data follows'
        endif

!--- MPI nloc arrays
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
    ! ----Debugging
        ! if(myid.eq.0) then
        !     print *, "nloc:", nloclist
        !     print *, "nloc*nev:", nloclist_nev
        !     print *, "npmin:", npminlist
        !     print *, "npmax:", npmaxlist
        !     print *, "nloc_sum:", nloc_sum
        !     print *, "nev_sum:", nev_sum
        ! endif

    ! !---P_ARPACK

        allocate(RESID(nloc),Vloc(nloc,NCV),Zloc(nloc,NEV),WORKD(nloc*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))

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
            if ((myid == 0).and.(mod(iter,100).eq.0)) print *, 'Iterations: ', iter, 'Shift:', is,'/',nshifts

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

    !----Shift eigenvalues to correct position
        do i=1,NEV
            D(i) = D(i) + epoints(is)
        enddo

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if (INFO .ne. 0) then
            if (myid == 0) print *, 'Error with pzneupd, info =', INFO
                call MPI_FINALIZE(ierr)
            stop
        end if

        deallocate(RESID, Vloc, WORKD, WORKL, RWORK)

        count = 0
        !----Spectral DOS
        do j=0,N3-1
            count = count + 1
            owner = 0
            do i = 1, nprocs
                if (j >= (npminlist(i) - 1) .and. j < npmaxlist(i)) then
                    owner = i - 1
                    exit
                endif
            enddo

            ! Compute local contribution to p_l for each eigenvector
            allocate(p_l_loc(NEV))
            p_l_loc = 0d0  ! Array to hold local contributions for all NEV
            allocate(p_l(NEV))
            p_l = 0d0

            if (myid == owner) then
                jloc = j - (npminlist(myid+1) - 1)
                do i = 1, NEV
                    p_l_loc(i) = dot_product(Zloc(1+(jloc*nb):(jloc+1)*nb, i), &
                                            Zloc(1+(jloc*nb):(jloc+1)*nb, i))
                enddo
            endif

            call MPI_Reduce(p_l_loc, p_l, NEV, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            if (myid == 0) then
                a_spec = 0d0
                do i = 1, NEV
                    factor = ((epoints(is) - real(D(i))) / eta)
                    ! print*,i,epoints(is),real(D(i)),p_l(i),factor
                    a_spec = a_spec + real(p_l(i)) * (exp(-0.5d0*factor**2)) * 1/(eta*sqrt(2*pi2))
                enddo
                write(100, '(1(1x,f12.8))') a_spec
            endif
            deallocate(p_l_loc, p_l)
        enddo

        if (myid == 0) write(100, '(a)') 'attribute "dep" string "positions"'

        ! Restore original diagonal elements
        do i=1,nb
            interp_Hr(i,i,0,0,0) = diag(i)
        enddo
        deallocate(Zloc)
    enddo

    deallocate(npminlist,npmaxlist,nloclist,nloc_sum,nev_sum,nloclist_nev)
    deallocate(diag,select,D,real_d,WORKEV)

    if (myid == 0) then
        do i=0,eres
            write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
            'object',eres+4+i,' class field', &
            'component "positions" value 1', &
            'component "connections" value 2', &
            'component "data" value ',3+i
        enddo
        write(100, '(a)') 'object "series" class series'
        do i=0,eres
            write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+eres+3), ' position', i
        enddo

        write(100, '(A)') 'end'

        call date_and_time(VALUES=time_end)

        print *, 'Start time: ', time_start(5), ':', time_start(6), ':', time_start(7)
        print *, 'End time: ', time_end(5), ':', time_end(6), ':', time_end(7)

        mem = size(Vloc)+size(Zloc)+size(WORKD)+size(WORKL)+size(RWORK)+size(select)+size(D)+size(real_d)+size(WORKEV)
        mem = mem*16.0d0/(1024.0**3)

        print *, ''
        print *, '------------------------------------------------'
        print *, 'Memory calculations:'
        print *, 'Total:',mem*nprocs,'GB'
        print *, 'Per process:',mem,'GB'
        print *, '------------------------------------------------'
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
                    vec_out(1+rowcount*nb : nb*(rowcount+1)) = vec_out(1+rowcount*nb : nb*(rowcount+1)) + &
                                            matmul(interp_Hr(:,:,r1,r2,r3), vec_in(1+colcount*nb : nb*(colcount+1)))
                endif
                colcount = colcount + 1
            enddo
            rowcount = rowcount + 1
        enddo

    end subroutine tv

end program
