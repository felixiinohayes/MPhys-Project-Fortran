module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    character*1:: bmat='I'
    character*2:: which='SM'
    real*8,parameter::ef= 4.18903772,a=1,emin=5.5,emax=6.5,eta1=2,eta2=0.03,TOL=0.0001,Bx=0.08
    integer*8,parameter::nblocks=6,matsize=(nblocks)**3,maxiter=100000,ishift=1,mode=1,eres=200
    integer nb
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    INCLUDE 'mpif.h'
    ! INCLUDE 'debug-arpack.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,l,nr,ie,lwork,info,ik,count,ir,ir3,ir12,nr12,r1,r2,r3,sign,il,i1,j1,i2,j2,i3,j3,xindex,yindex,rvec_data(3),pindex,index,nindex,interp_size
    integer*4 IPARAM(11),IPNTR(14),iter,IDO,LDV,LDZ
    integer*4 comm,myid,nprocs,rc,ierr
    integer*8 LWORKL,NEV,NCV,N,nloc,npmin,npmax,leng
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,epoints(eres),a_spec,factor,p_l,de,dos
    real*8,allocatable:: rvec(:,:),rwork(:)
    integer*4,allocatable:: ndeg(:),vec_ind(:,:)
    complex*16,allocatable::top_Hr(:,:),triv_Hr(:,:),super_H(:,:),surface_vec(:),B_pt(:,:)
    complex*16,allocatable::RESID(:),V(:,:),WORKD(:),WORKL(:),D(:),WORKEV(:),Z(:,:)
    complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr
    complex*16 SIGMA,b_sigma(2,2)
    logical:: rvecmat
    logical,allocatable:: select(:)
!----Date and Time
    character(len=8) :: date_start, time_start
    character(len=8) :: date_end, time_end
    real*8::end_second, start_second
!------------------------------------------------------
    pi2=4.0d0*atan(1.0d0)*2.0d0
    call date_and_time(date_start, time_start)

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

    open(100,file='mpi_test.dx')
    open(200,file='ene_total.dat')
    open(300,file='ene_surface.dat')

!------read H(R)
    interp_size=6
    if(abs(nblocks) > interp_size) interp_size = abs(nblocks)

    !---- Magnetic Perturbation


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
    B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(Bx,0d0)]
    B_sigma(2,:) = [dcmplx(Bx,0d0) ,  dcmplx(0d0,0d0)]

    !B along Y axis
	! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-Bx)]
    ! B_sigma(2,:) = [dcmplx(0d0,Bx) ,  dcmplx(0d0,0d0)]

    !B along Z-axis
    ! B_sigma(1,:) = [dcmplx(Bx,0d0),  dcmplx(0d0,0d0)]
    ! B_sigma(2,:) = [dcmplx(0d0,0d0) ,  dcmplx(-Bx,0d0)]
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

!-----Energy Mesh
    de = (emax-emin)/eres

    do i=1, eres
        epoints(i) = emin + de*i
    enddo

!------ARPACK

    call MPI_INIT( ierr )
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK( comm, myid, ierr )
    call MPI_COMM_SIZE( comm, nprocs, ierr )
 
    N=nb*matsize
    NEV=N-2
    NCV=NEV+2

    nloc = (N/(nprocs*nb))
    if (mod(N,nprocs).ne.0) nloc= nloc +1 

    npmin=1+myid*nloc
    npmax=(myid+1)*nloc
    leng = (min(npmax,matsize)-npmin+1)*nb

    ! print*,'npmin:',npmin,'npmax:',npmax,'leng:',leng,'myid:',myid

    allocate(RESID(N),V(N,NCV),WORKD(leng*3),WORKL(3*NCV*NCV + 5*NCV+10),RWORK(NCV))
    allocate(select(NCV),D(NEV),Z(N,NEV),WORKEV(2*NCV))


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
        if(myid.eq.0) print *, iter
        call pznaupd(comm,IDO,bmat,leng,which,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)

        if(IDO==99) exit
        
        if(IDO==-1 .or. IDO==1) then
            !WORKD(IPNTR(2):IPNTR(2)+N-1) = matmul(super_H,WORKD(IPNTR(1):IPNTR(1)+N-1))
            ! call matmul_chunk(interp_Hr, WORKD(IPNTR(1):IPNTR(1)+nloc-1), WORKD(IPNTR(2):IPNTR(2)+nloc-1),nloc)
            call matmul_(interp_Hr, WORKD(IPNTR(1):IPNTR(1)+leng-1), WORKD(IPNTR(2):IPNTR(2)+leng-1),nloc,nblocks,N,npmin,npmax,leng)
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
        if(myid.eq.0)print *, "Finished iterations, calling zneupd..."
        call pzneupd (comm,rvecmat, 'A', select, d, v, ldv, sigma,&
             workev, bmat, leng, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl, &
             rwork, info)
        !print*, v(1,:)
    endif

    deallocate(RESID,WORKD,WORKL,RWORK)
    deallocate(Z,WORKEV)


!-------Header File
    if(myid.eq.0) then
        write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nblocks,nblocks,nblocks
        write(100, '(a,3(1x,f12.8))') 'origin',0d0,0d0,0d0
        write(100, '(a,3(1x,f12.8))') 'delta',0d0,0d0,1d0
        write(100, '(a,3(1x,f12.8))') 'delta',0d0,1d0,0d0
        write(100, '(a,3(1x,f12.6))') 'delta',1d0,0d0,0d0
        write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nblocks,nblocks,nblocks

        ! allocate(surface_vec(4*nb*(nblocks-1)),vec_ind(matsize,3))

    !------Computes total DOS for each Z layer

        allocate(surface_vec(4*nb*(nblocks-1)),vec_ind(matsize,3))

        do j=1,matsize
        !------ Indices for blocks in eigenvector
            r3 = ((j-1)/nblocks**2)
            r2 = mod((j-1)/nblocks,3)
            r1 = mod((j-1),3)

            vec_ind(j,1) = r1
            vec_ind(j,2) = r2
            vec_ind(j,3) = r3
        enddo
        ! print *, d

        count = 0 
        do ie=1,eres
            count = count + 1
            print*, ie, eres

            write(100, '(a,i8,a,i8,a,i10,a)') 'object',2+count,' class array type float rank 1 shape',1,&
                                    ' item', matsize, ' data follows'
            !----Spectral DOS
            do j=0,matsize-1
                a_spec = 0d0
                do i=1,N
                    p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))

                    factor = ((epoints(ie)- d(i)))/eta1
                    ! if(ie==1)print *, epoints(ie)-d(i)
                    ! if(ie==1)print*, p_l* (exp(-0.5d0*factor**2)) * 1/sqrt(2*pi2*eta**2)

                    a_spec = a_spec + p_l* (exp(-0.5d0*factor**2)) * 1/sqrt(2*pi2*eta1**2)
                enddo
                write(100, '(3(1x,f12.10))') a_spec
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
    endif
!
    call MPI_FINALIZE(IERR)

    ! call date_and_time(date_end, time_end)
    ! read(time_end  , '(f10.1)') end_second
    ! read(time_start, '(f10.1)') start_second
    ! ! print*,"Start Time: ", trim(time_start)
    ! print*, 'Duration: ', end_second- start_second

    contains

        subroutine matmul_chunk(interp_Hr,vec_in,vec_out,N)  
            integer*4,intent(in)::N
            complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr
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
                                    if((abs(r1).lt.6).or.(abs(r2).lt.6).or.((abs(r3).lt.6))) then
                                        tempvec((1+nb*yindex):(nb*(yindex+1))) = tempvec((1+nb*yindex):(nb*(yindex+1))) + matmul(interp_Hr(:,:,r1,r2,r3),vec_in((1+nb*xindex):(nb*(xindex+1))))
                                    endif
                                    ! tempvec(1:N) = tempvec(1:N) + 0.0001
                                    ! print *, xindex,yindex 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            vec_out=tempvec

        end subroutine matmul_chunk

        subroutine matmul_(interp_Hr,vec_in,vec_out,nloc,nblocks,N,npmin,npmax,leng)  
            integer*8,intent(in)::N,nloc,nblocks,npmin,npmax,leng
            integer*8 npmin_t,npmax_t,len_t,next,prev,send,recv,status(MPI_STATUS_SIZE)
            complex*16,dimension(18,18,-nblocks:nblocks,-nblocks:nblocks,-nblocks:nblocks), intent(in):: interp_Hr
            complex*16,intent(in):: vec_in(leng)
            complex*16,intent(out)::vec_out(leng)
            complex*16::tempvec(N)
            complex*16 mv_buf(leng),tv_out(leng)
            integer*4::icol,irow,r1,r2,r3,f1,f2,f3,N2,N3,count,count1

            call MPI_COMM_RANK( comm, myid, ierr )
            call MPI_COMM_SIZE( comm, nprocs, ierr )

            N3 = nblocks**3
            N2 = nblocks**2
            tempvec=0d0
            count1 =0

            next = myid + 1 
            prev = myid - 1

            do icol = npmin,min(npmax,N3)
                f3 = (icol-1)/(N2)
                f2 = mod((icol-1)/nblocks,3)
                f1 = mod(icol-1,3)
                do irow = 1,N3
                    r3 = ((irow-1)/N2)- f3
                    r2 = mod((irow-1)/nblocks,3) - f2
                    r1 = mod(irow-1,3) - f1
                    if((abs(r1).lt.6).or.(abs(r2).lt.6).or.((abs(r3).lt.6))) then
                        tempvec(1+(irow-1)*nb : nb*(irow)) = tempvec(1+(irow-1)*nb : nb*(irow)) + matmul(interp_Hr(:,:,r1,r2,r3), vec_in( 1+(icol-1)*nb : nb*(icol) ))
                    endif
                enddo
            enddo

            do i=1,nprocs 
                send = modulo(myid+i,nprocs)
                recv = modulo(myid-i,nprocs)


                npmin_t=1+(send)*nloc
                npmax_t=(send+1)*nloc
                len_t = (min(npmax_t,N3)-npmin_t+1)*nb

                ! print*,recv,myid,send

                call mpi_sendrecv(tempvec(1+(npmin_t-1)*nb : nb*min(npmax_t,N3)),len_t,MPI_DOUBLE_COMPLEX,&
                                    send,0, mv_buf,leng,MPI_DOUBLE_COMPLEX,recv,0,comm,status)
                
                ! print*, 'Processor:',myid,'sent to',send
                ! print*, 'Processor:',myid,'recieved from',recv
                ! call mpi_send(tempvec(1+(npmin_t-1)*nb : nb*min(npmax_t,N3)),len_t,MPI_DOUBLE_COMPLEX,&
                !                 myid +i, myid+ i ,comm, ierr)
                ! call mpi_recv(mv_buf,leng,MPI_DOUBLE_COMPLEX,&
                !                 myid -i, myid -i ,ierr)
                ! call mpi_barrier(comm)
                tv_out = tv_out + mv_buf
            enddo

            vec_out = tv_out

            call mpi_barrier(comm)



                ! do i=1,nprocs

                !     pindex = mod((i+nprocs-3),nprocs)
                !     index  = mod((i+nprocs-2),nprocs)
                !     nindex = mod((i+nprocs-1),nprocs) 

                !     print*,pindex,index,nindex

                !     npmin_t=1+(i-1)*nloc
                !     npmax_t=(i)*nloc
                !     len_t = (min(npmax_t,N3)-npmin_t+1)*nb

                !     if(myid.ne.index) then
                !         send = mod(next+(i-1),nprocs)

                !         if(myid.eq.pindex.and.i.gt.1) send = send + 1

                !         ! call mpi_send(tempvec(1+(npmin_t-1)*nb : nb*min(npmax_t,N3)),len_t,MPI_DOUBLE_COMPLEX,&
                !         !                 send,send,comm, ierr)
                !         ! print*, 'Processor:',myid,'sent to',send
                !     endif
                !     if(myid.ne.nindex) then 
                !         ! call mpi_recv(mv_buf,leng,MPI_DOUBLE_COMPLEX,&
                !         !                  (prev-(i-1)), (prev-(i-1)),comm, ierr)
                !         ! tv_out = tv_out + mv_buf
                !         ! print*, 'Processor:',myid,'recv from ',(prev-(i-1))
                !     endif
                ! enddo

                ! print*,''


        end subroutine matmul_

end Program Projected_band_structure

            ! do i=1,nprocs
                
            !     npmin_t=1+(i-1)*nloc
            !     npmax_t=(i)*nloc

            !     len_t = (min(npmax_t,N3)-npmin_t+1)*nb
            !     if(myid.ne.(i-1)) then 
            !         call mpi_send(tempvec(1+(npmin_t-1)*nb : nb*min(npmax_t,N3)),len_t,MPI_DOUBLE_COMPLEX,&
            !                         i-1, myid,comm, ierr)
            !         print*, 'Processor:',myid,'sent to',i-1
            !     endif 
            ! enddo 

            ! do i=1,nprocs
            !     if(myid.ne.(i-1)) then
            !         call mpi_recv(mv_buf,leng,MPI_DOUBLE_COMPLEX,&
            !                         i-1, myid,comm, ierr)
            !          print*, 'Processor:',myid,'recieved from',i-1
            !             ! tv_out = tv_out + mv_buf
            !     endif 
            ! enddo

            ! vec_out = tv_out







