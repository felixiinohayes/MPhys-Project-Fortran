module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="../BiTeI", ax = 'x'
    real*8,parameter::ef_triv=5.2,ef_top=6.5,a=1,emin=5.5,emax=7,bfactor=0.005, B=0.00d0, passval=0.0d0
    integer,parameter::nkpath=3,np=300,eres=300,nblocks=5,nk=(nkpath-1)*np+1,nepoints=2*eres+1,N2=nblocks**2
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
    integer*4 i,j,k,l,nr,i1,i2,j1,j2,ie,il,ir,ir3,ir12,nr12,r3,ikp,jk,ira,irb,irc,ra,rb,fa,fb,n,matsize,dim,sum2,sum1
    integer*4 lwork,info,ik,count,sign,kloc,kpmin,kpmax,ecounts,kcount,interp_size,nr_top,nr_triv,rvec(3),index
    integer*4 recv(1),nr_(3),kindex(2),ai(3)
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,de,exp_factor,p_l,spectral_A,emiddle
    real*8 xk(nk),avec(3,3),bvec(3,3),rvec_data(3),kpoints(3,nkpath),dk(3),epoints(nepoints),spectral_A_comm(3,nk*nepoints),kpath(3,nk),kvec1(3),kvec2(3)
    real*8,allocatable:: rwork(:),eval(:,:),eval_flat(:),eval_flatg(:)
    integer*4,allocatable:: ndeg(:,:,:),displs1(:),recvcounts1(:),displs2(:),recvcounts2(:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:),kpminlist(:),kpmaxlist(:),kloclist(:),kloc_sum(:),buff_sum1(:),buff_sum2(:),buffsize1(:),buffsize2(:)
    complex*16,allocatable::Hk(:,:),Hkra(:,:,:,:),work(:),super_H(:,:,:),sH_flat(:),sH_flatg(:),SH(:,:,:),B_pt(:,:),top_Hr_temp(:,:),triv_Hr_temp(:,:),extrarow(:)
    complex*16 B_sigma(2,2)
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
    open(100,file='super_H_X_30_BP2D.dx')
    open(200,file='2D_EVECS.dat')
    open(300,file='2D_EVALS.dat')


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

    allocate(top_Hr_temp(nb,nb),triv_Hr_temp(nb,nb),ndeg_top(nr_top),ndeg_triv(nr_triv),ndeg(-6:6, -6:6, -6:6))
    allocate(rvec_top(nr_top,3))
    allocate(interp_Hr(nb,nb,-6:6, -6:6, -6:6))
    allocate(extrarow(nb*N2))
    allocate(Hkra(nb,nb,-6:6,-6:6))
    allocate(Hk(nb,nb))

    read(99,*)ndeg_top
    read(97,*)ndeg_triv

    lwork=max(1,2*(nb*N2+1)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*N2+1)-2)))

!-----kpath
    ! !-kx -> kx
    kpoints(:,1) = [ -0.5d0,  0.0d0,   0.0d0]  !-M
    kpoints(:,2) = [  0.0d0,  0.0d0,   0.0d0]  !Gamma
    kpoints(:,3) = [  0.5d0,  0.0d0,   0.0d0]  !M   
    
    ! ! -ky -> ky 
    ! kpoints(:,1) = [ 0.00d0, -0.5d0,  0.0d0]  !H
    ! kpoints(:,2) = [  0.0d0,  0.0d0,  0.0d0]  !A
    ! kpoints(:,3) = [ 0.00d0,  0.5d0,  0.0d0]  !-H

    ! -kz -> kz
    ! kpoints(:,1) = [ 0.00d0, 0.0d0,  -0.5d0]  !H
    ! kpoints(:,2) = [  0.0d0,  0.0d0,  0.0d0]  !A
    ! kpoints(:,3) = [0.00d0,  0.0d0,  0.5d0]  !-H


    ! Initial point in the k path
    kvec1(:)=(kpoints(1,1)-kpoints(1,2))*bvec(:,1)+(kpoints(2,1)-kpoints(2,2))*bvec(:,2)+(kpoints(3,1)-kpoints(3,2))*bvec(:,3)
    xk(1)= -sqrt(dot_product(kvec1,kvec1))

    kvec1 = 0d0
    
    ! Interpolation between points
    do i = 1, nkpath-1
        do j = 1, np
            ik = j + np*(i-1)
            dk = (kpoints(:,i+1)-kpoints(:,i))/np
            kpath(:, ik) = kpoints(:,(i)) + (dk*(j-1))
            kvec2 = kpath(1,ik)*bvec(:,1) + kpath(2,ik)*bvec(:,2) + kpath(3,ik)*bvec(:,3) 
            if(ik.gt.1) xk(ik) =  xk(ik-1) + sqrt(dot_product(kvec2-kvec1,kvec2-kvec1))
            kvec1 = kvec2
        enddo
    enddo
    ! Final point in the kpath
    kpath(:,nk) = kpoints(:,nkpath)
    kvec2=kpath(1,nk)*bvec(:,1)+kpath(2,nk)*bvec(:,2)+kpath(3,nk)*bvec(:,3)
    xk(nk)=xk(nk-1)+sqrt(dot_product(kvec2-kvec1,kvec2-kvec1))
    kpath = kpath*pi2

!----energy mesh
    emiddle = emin + (emax-emin)/2
    de = (emax-emin)/(2*eres)
    ie=0
    do i=-eres, eres
        ie=ie+1
        epoints(ie) = emiddle + de*i
    enddo


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

    ! B along Y axis
    ! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-B)]
    ! B_sigma(2,:) = [dcmplx(0d0,B) ,  dcmplx(0d0,0d0)]

    ! B along X axis
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


!----Read in Hamiltonian
    interp_Hr=0d0
    ndeg = 0d0
    do ir=1,nr_top
        do i=1,nb
            do j=1,nb
               read(99,*)rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3),i1,i2,x1,y1
               top_Hr_temp(i1,i2)=dcmplx(x1,y1)
               ndeg(rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = ndeg_top(ir)
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


    ! do i=1,nb
    !     if(a==0) then 
    !         interp_Hr(i,i,0,0,0) = interp_Hr(i,i,0,0,0) - ef_triv
    !     else 
    !         interp_Hr(i,i,0,0,0) = interp_Hr(i,i,0,0,0) - ef_top
    !     endif
    ! enddo
!-----kloc
    dim = (N2*nb)
    matsize= dim**2
    allocate(kpminlist(numprocs),kpmaxlist(numprocs),kloclist(numprocs),kloc_sum(numprocs+1), & 
             buff_sum2(numprocs+1),buffsize2(numprocs),displs2(numprocs),recvcounts2(numprocs),&
             buff_sum1(numprocs+1),buffsize1(numprocs),displs1(numprocs),recvcounts1(numprocs))
    
    sum1 = 1
    sum2 = 1
    kloc_sum(1) = 1
    buff_sum1(1)  = 1
    buff_sum2(1)  = 1

    do i=1,numprocs
        kloc=nk/numprocs
        if (mod(nk,numprocs).ne.0) kloc=kloc+1
        kpminlist(i)=1+(i-1)*kloc
        kpmaxlist(i)=min(i*kloc,nk)
        kloclist(i) = (kpmaxlist(i)-kpminlist(i)+1)

        kloc_sum(i+1) = kloc_sum(i) + kloclist(i)

        !For the MPI_GATHERV(super_H) call
        buffsize1(i) = (kpmaxlist(i)-kpminlist(i)+1)*matsize
        sum1 = sum1+ buffsize1(i)

        buff_sum1(i+1) = sum1

        displs1(i) = buff_sum1(i) - 1

        !For the MPI_GATHERV(evals) call
        buffsize2(i) = (kpmaxlist(i)-kpminlist(i)+1)*dim !For the MPI_GATHERV(eval) call
        sum2 = sum2+ buffsize2(i)

        buff_sum2(i+1) = sum2

        displs2(i) = buff_sum2(i) -1
    enddo

    recvcounts1 = buffsize1
    recvcounts2 = buffsize2
    kloc = kloclist(myid+1) 

    allocate(SH(dim,dim,nk),super_H(dim,dim,kloc),sH_flat(buffsize1(myid+1)),sH_flatg(matsize*nk)&
                                ,eval(dim,kloc),eval_flat(dim*kloc),eval_flatg(dim*nk))

    if(myid.eq.0) print *, "kloc:", kloclist
    if(myid.eq.0) print *, "kpmin:", kpminlist
    if(myid.eq.0) print *, "kpmax:", kpmaxlist
    if(myid.eq.0) print *, "kloc_sum:", kloc_sum
    if(myid.eq.0) print *, "buffsize1:", buffsize1
    if(myid.eq.0) print *, "buffsize2:", buffsize2
    if(myid.eq.0) print *, "buff_sum1:", buff_sum1, displs1
    if(myid.eq.0) print *, "buff_sum2:", buff_sum2, displs2


!----- Axis selection
    extrarow=0d0
    if(ax == 'x') index = 1
    if(ax == 'y') index = 2
    if(ax == 'z') index = 3

    ! extrarow(3:4) = +2*passval
    ! extrarow(nblocks*nb-3:nblocks*nb-2) = -3*passval

!----- Perform fourier transform
    count = 0
    ! if(myid.eq.0 .and. (il == 0 .or. il == nblocks - 1 .or. il == nblocks/2)) then
    !     write(100, '(a,i8,a,i10,a)') 'object',count+3,' class array type float rank 1 shape 3 item',nk*nepoints,' data follows'
    !     count = count + 1
    ! endif
    ikp=0
    ! if(myid.eq.0)print *, "block", il+1, "/", nblocks
    do ik=kpminlist(myid+1),kpmaxlist(myid+1)
        ikp=ikp+1
        do ira= -6,6 
            do irb = -6,6  ! Loop over R_ vectors
                Hk=0d0   
                do irc = -6,6
                    
                    ! Now 'ax' is  the only axis where periodicity is maintained
                    if(ax == 'x') ai = [irc, ira, irb]
                    if(ax == 'y') ai = [ira, irc, irb]
                    if(ax == 'z') ai = [ira, irb, irc]

                    if(ndeg(ai(1),ai(2),ai(3)).ne.0) then 
                        phase = 0d0

                        phase = phase + kpath(index,ik) * irc

                        Hk=Hk+((1-a)*(triv_Hr(:,:,ai(1),ai(2),ai(3)))+(a)*(top_Hr(:,:,ai(1),ai(2),ai(3))))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ai(1),ai(2),ai(3)))
                    endif
                enddo
                Hkra(:,:,ira,irb) = Hk
            enddo
        enddo

        do i=0,N2-1
            fb = mod((i)/nblocks,nblocks)
            fa = mod(i,nblocks)
            do j=0,N2-1
                rb = mod((j)/nblocks,nblocks) - fb
                ra = mod(j,nblocks) - fa
                ! print*,ra,rb
                if (abs(ra).lt.6 .AND. abs(rb).lt.6 ) then
                    super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1)),ikp) = Hkra(:,:,ra,rb)
                else
                    super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1)),ikp) = 0d0
                endif
            enddo
        enddo

        call zheev('V','U',dim,super_H(:,:,ikp),dim,eval(:,ikp),work,lwork,rwork,info)

        sH_flat((1+((ikp-1)*matsize)):(matsize*(ikp))) = reshape(super_H(:,:,ikp),[matsize])
        eval_flat((1+((ikp-1)*dim)):(dim*ikp)) = reshape(eval(:,ikp),[dim]) 
    enddo

    call MPI_GATHERV(sH_flat,matsize*kloc,MPI_DOUBLE_COMPLEX, &
                        sH_flatg,recvcounts1,displs1,MPI_DOUBLE_COMPLEX, &
                        0, MPI_COMM_WORLD,IERR)

    call MPI_GATHERV(eval_flat,dim*kloc,MPI_DOUBLE_PRECISION, &
                        eval_flatg,recvcounts2,displs2,MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD,IERR)

    if(myid.eq.0) then 
        write(200, *) nblocks
        write(200, *) nk

        do i =1, nk
            do j = 1, dim
                write(200, *) sH_flatg(1 + matsize*(i-1)+(j-1)*dim : dim*(j) + matsize*(i-1))
            enddo
            write(300,*) eval_flatg(1 + dim*(i-1): dim*(i))
        enddo
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