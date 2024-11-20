module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="../BiTeI", ax = 'x'
    real*8,parameter::ef_triv=4.23,ef_top=6.5,a=1,emin=6,emax=7,bfactor=0.002, B=0.00d0
    integer,parameter::nkpath=3,np=200,eres=400,nblocks=20,nk=(nkpath-1)*np+1,nepoints=2*eres+1
    integer nb
    INTEGER IERR,MYID,NUMPROCS

end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    INCLUDE 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    character(len=5) suffix
    integer*4 i,j,k,nr,i1,i2,j1,j2,ie,lwork,info,ik,count,ir,ir3,ir12,nr12,r3,sign,il,kpool,kpmin,kpmax,ecounts,ikp,jk,kcount,sum,interp_size,nr_top,nr_triv,rvec(3),ira,irb,irc,ra
    integer*4 recv(1),nr_(3),kindex(2)
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,de,exp_factor,p_l,spectral_A,emiddle
    real*8 xk(nk),avec(3,3),bvec(3,3),rvec_data(3),kpoints(3,nkpath),dk(3),epoints(nepoints),spectral_A_comm(3,nk*nepoints),kpath(3,nk)
    real*8,allocatable:: rwork(:),k_ene(:,:),spectral_A_single(:,:)
    integer*4,allocatable:: ndeg(:),displs(:),recvcounts(:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:)
    complex*16,allocatable::Hk(:,:),Hkra(:,:,:),work(:),super_H(:,:),sH(:,:),a_p_top(:,:),a_p_bottom(:,:),B_pt(:,:),top_Hr_temp(:,:),triv_Hr_temp(:,:),extrarow(:,:)
    complex*16 B_sigma(2,2),temp1,temp2
    complex*16,dimension(4,4,-6:6,-6:6,-6:6) :: top_Hr
    complex*16,dimension(4,4,-6:6,-6:6,-6:6) :: triv_Hr
    complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr
!------------------------------------------------------
    call init_mpi

    pi2=4.0d0*atan(1.0d0)*2.0d0

    write(top_file, '(a,a)') trim(adjustl(prefix)), "_hr_topological_4band.dat"
    write(triv_file, '(a,a)') trim(adjustl(prefix)), "_hr_trivial_4band.dat"
    write(nnkp, '(a,a)') trim(adjustl(prefix)), ".nnkp"
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
    open(100,file='super_H_B12.dx')

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

!------read H(R)
    interp_size=6
    ! if((nxblocks > interp_size).or.(nyblocks > interp_size).or.(nzblocks > interp_size)) interp_size = max(max(nxblocks,nyblocks),nzblocks)
    read(99,*)
    read(99,*)nb,nr_top
    read(97,*)
    read(97,*)nb,nr_triv

    allocate(top_Hr_temp(nb,nb),triv_Hr_temp(nb,nb),ndeg_top(nr_top),ndeg_triv(nr_triv))
    allocate(rvec_top(nr_top,3))
    allocate(interp_Hr(nb,nb,-6:6, -6:6, -6:6),super_H(nb*nblocks,nb*nblocks))
    allocate(extrarow(nb,nb*nblocks))

    read(99,*)ndeg_top
    read(97,*)ndeg_triv

    lwork=max(1,2*(nb*nblocks)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*nblocks)-2)))

!-----kpath
    kpoints(:,1) = [ -0.2d0,  0.0d0,   0.0d0]  !-M
    kpoints(:,2) = [  0.0d0,  0.0d0,   0.0d0]  !Gamma
    kpoints(:,3) = [  0.2d0,  0.0d0,   0.0d0]  !M   
    
    ! ky -> -ky 
    ! kpoints(:,1) = [ 0.05d0,  -0.1d0,   0.5d0]  !H
    ! kpoints(:,2) = [ 0.0d0,   0.0d0,    0.5d0]  !A
    ! kpoints(:,3) = [ -0.05d0,   0.1d0,  0.5d0]  !-H


    ! kx -> -kx
    ! kpoints(:,1) = [ -0.5d0,   0.0d0,   0.5d0 ]  !L
    ! kpoints(:,2) = [ 0.0d0,   0.0d0,   0.5d0 ]  !A
    ! kpoints(:,3) = [ 0.5d0,  0.0d0,   0.5d0 ]  !-L



    do j = 1, nkpath-1
          sign = 1
          if(j ==1) sign = -1
        do i = 1, np+1
            ik = i + np*(j-1)
            dk = (kpoints(:,j+1)-kpoints(:,j))/np
            kpath(:, ik) = kpoints(:,(j)) + (dk*(i-1))
            xk(ik) =  sign*sqrt(dot_product(kpoints(:,2)- kpath(:, ik),kpoints(:,2) - kpath(:, ik)))
        enddo
    enddo

    emiddle = emin + (emax-emin)/2
    de = (emax-emin)/(2*eres)
    ie=0
    do i=-eres, eres
        ie=ie+1
        epoints(ie) = emiddle + de*i
    enddo

    kpool=nk/numprocs
    if (mod(nk,numprocs).ne.0) kpool=kpool+1

    kpmin=1+myid*kpool
    kpmax=(myid+1)*kpool

!----Construct magnetic perturbation
    ! call write_header()
    if(myid.eq.0) then
        write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nk,nepoints
        write(100, '(a,2(1x,f12.6))') 'origin',-0.1d0,emin
        write(100, '(a,2(1x,f12.6))') 'delta',sqrt(dot_product(dk,dk)),0d0
        write(100, '(a,2(1x,f12.6))') 'delta',0d0,de
        write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nk,nepoints
    endif
    allocate(B_pt(nb, nb))

     !B along Y axis
    ! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-B)]
    ! B_sigma(2,:) = [dcmplx(0d0,B) ,  dcmplx(0d0,0d0)]

     !B along X axis
    B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(B,0d0)]
    B_sigma(2,:) = [dcmplx(B,0d0) ,  dcmplx(0d0,0d0)]

    ! B_sigma(1,:) = [dcmplx(B,0d0),  dcmplx(0d0,0d0)]
    ! B_sigma(2,:) = [dcmplx(0d0,0d0) ,  dcmplx(-B,0d0)]
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

    recv(1)=(min(kpmax,nk)-kpmin+1)*nepoints*3

    allocate(spectral_A_single(3,recv(1)*nepoints),displs(numprocs),recvcounts(numprocs))

    displs=0
    call MPI_GATHER(recv,1,MPI_INTEGER, &
                        recvcounts,1,MPI_INTEGER, &
                        0, MPI_COMM_WORLD,IERR)
    sum=0
    do i=1,numprocs
        displs(i) = sum
        sum=sum+recvcounts(i) 
    enddo
    kcount = (min(kpmax,nk)-kpmin+1)*nepoints*3

!----- Axis selection
    if(ax == 'x') kindex = [2,3]
    if(ax == 'y') kindex = [1,3]
    if(ax == 'z') kindex = [1,2]

    allocate(Hkra(nb,nb,-6:6))

!----- Perform fourier transform
    ! nr12=nr/nr3
    do il=0,nblocks-1
        if(myid.eq.0) write(100, '(a,i8,a,i10,a)') 'object',il+3,' class array type float rank 1 shape 3 item',nk*nepoints,' data follows'
        ikp=0
        if(myid.eq.0)print *, "block", il+1, "/", nblocks
        do ik=kpmin,min(kpmax,nk)
            ikp=ikp+1
            do ira= -6,6 ! Loop over R_ vectors
                Hk=0d0    
                do irb = -6,6
                    do irc = -6,6
                        phase = 0d0

                        phase = phase + kpath(kindex(1),ik) * irb + kpath(kindex(2),ik) * irc

                        Hk=Hk+((1-a)*(triv_Hr(:,:,ira,irb,irc))+(a)*(top_Hr(:,:,ira,irb,irc)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
                    enddo
                enddo
                Hkra(:,:,ira) = Hk
            enddo

            do i=0,nblocks-1
                do j=0,nblocks-1
                    ra = i-j
                    if (ra<=6 .AND. ra>=-6) then   
                        super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1))) = Hkra(:,:,ra)
                    else
                        super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1))) = 0d0
                    endif
                enddo
            enddo
            call zheev('V','U',nb*nblocks,super_H,nb*nblocks,k_ene(:,ik),work,lwork,rwork,info)

            do ie=1,nepoints
                spectral_A = 0d0
                do i=1,nb*nblocks
                    p_l = dot_product(super_H((1+nb*il):(nb*(il+1)),i),super_H((1+nb*il):(nb*(il+1)),i))
                    exp_factor = (epoints(ie) - k_ene(i,ik))/bfactor
                    spectral_A = spectral_A + p_l * exp(-0.5d0 * (exp_factor**2))
                enddo
                spectral_A_single(1,(ikp-1)*nepoints + ie) = xk(ik)
                spectral_A_single(2,(ikp-1)*nepoints + ie) = epoints(ie)
                spectral_A_single(3,(ikp-1)*nepoints + ie) = real(spectral_A)! Top surface
            enddo
        enddo
        call MPI_GATHERV(spectral_A_single,kcount,MPI_DOUBLE_PRECISION, &
                            spectral_A_comm,recvcounts,displs,MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD,IERR)

        if(myid.eq.0) write(100, '(3(1x,f12.6))') spectral_A_comm
        if(myid.eq.0) write(100, '(a)') 'attribute "dep" string "positions"'
    enddo

    if(myid.eq.0) then
        do i=0,nblocks-1
            write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
            'object',nblocks+3+i,' class field', &
            'component "positions" value 1', &
            'component "connections" value 2', &
            'component "data" value ',3+i
        enddo
        write(100, '(a)') 'object "series" class series'
        do i=0,nblocks-1
            write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+nblocks+3), ' position', i
        enddo
        write(100, '(A)') 'end'
    endif

    call MPI_FINALIZE(IERR)
end Program Projected_band_structure

SUBROUTINE INIT_MPI
    USE PARAMETERS               ,             ONLY: IERR,MYID,NUMPROCS
    IMPLICIT NONE
    INCLUDE 'mpif.h'
        Call MPI_INIT( IERR )
        Call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
        Call MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS , IERR )
!        Write(*,*) ‘Process’, myid, ' of ’, NUMPROCS , ‘is alive.’
END SUBROUTINE INIT_MPI
