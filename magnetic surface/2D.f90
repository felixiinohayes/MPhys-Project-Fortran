module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="../BiTeI", ax = 'x'
    real*8,parameter::ef_triv=5.2,ef_top=6.5,a=1,B=0.005d0,passval=0.0d0,emin=6,emax=7,eta=0.005
    integer,parameter::nkpath=3,np=30,nblocks=15,nk=(nkpath-1)*np+1,N2=nblocks**2,eres=80,nblocks_2=nblocks/2
    integer nb
    INTEGER IERR,MYID,NUMPROCS
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    INCLUDE 'mpif.h'

!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,l,nr,i1,i2,j1,j2,ie,il,ir,ir3,ir12,nr12,r3,ikp,jk,ira,irb,irc,ra,rb,fa,fb,n,matsize,dim,sum2,sum1,ib,offset,dur
    integer*4 lwork,info,ik,count,sign,kloc,kpmin,kpmax,ecounts,kcount,interp_size,nr_top,nr_triv,rvec(3),index
    integer*4 recv(1),nr_(3),kindex(2),ai(3)
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,de,exp_factor,p_l,spectral_A,emiddle
    real*8 xk(nk),avec(3,3),bvec(3,3),rvec_data(3),kpoints(3,nkpath),dk(3),kpath(3,nk),kvec1(3),kvec2(3),epoints(eres),a_spec,factor,dos,data_g(5,3,nk*eres)
    real*8,allocatable:: rwork(:),eval(:),eval_flat(:),eval_flatg(:),data_row(:,:,:)
    integer*4,allocatable:: ndeg(:,:,:),displs1(:),recvcounts1(:),displs2(:),recvcounts2(:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:),kpminlist(:),kpmaxlist(:),kloclist(:),kloc_sum(:),buff_sum1(:),buff_sum2(:),buffsize1(:),buffsize2(:)
    complex*16,allocatable::Hk(:,:),Hkra(:,:,:,:),work(:),super_H(:,:),sH_flat(:),sH_flatg(:),SH(:,:,:),B_pt(:,:),top_Hr_temp(:,:),triv_Hr_temp(:,:),extrarow(:)
    complex*16 B_sigma(2,2)
    complex*16,dimension(4,4,-6:6,-6:6,-6:6) :: top_Hr
    complex*16,dimension(4,4,-6:6,-6:6,-6:6) :: triv_Hr
    complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr
    integer, dimension(8) :: time_start
    integer, dimension(8) :: time_prev
    integer, dimension(8) :: time_next
    integer, dimension(8) :: time_end
!------------------------------------------------------
    call init_mpi
!----Date and Time
    call date_and_time(VALUES=time_start)

    pi2=4.0d0*atan(1.0d0)*2.0d0
    
    write(triv_file, '(a,a)') trim(adjustl(prefix)), "_hr_trivial_4band.dat"
    write(top_file, '(a,a)') trim(adjustl(prefix)), "_hr_topological_4band.dat"
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
    open(100,file='block_1.dx')
    open(101,file='block_2.dx')
    open(102,file='block_3.dx')
    open(103,file='block_4.dx')
    open(104,file='block_5.dx')
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

    select case (ax)
    case ('x')
        !-kx -> kx
        kpoints(:,1) = [ -0.1d0,  0.0d0,  0.0d0]  !-M
        kpoints(:,2) = [ 0.0d0,   0.0d0,  0.0d0]  !-M
        kpoints(:,3) = [ 0.1d0,   0.0d0,  0.0d0]  !-M
    case ('y')
        ! -ky -> ky 
        kpoints(:,1) = [ 0.00d0, -0.1d0,  0.0d0]  !H
        kpoints(:,2) = [  0.0d0,  0.0d0,  0.0d0]  !A
        kpoints(:,3) = [ 0.00d0,  0.1d0,  0.0d0]  !-H
    case ('z')
        ! -kz -> kz
        kpoints(:,1) = [ 0.00d0, 0.0d0,  -0.1d0]  !H
        kpoints(:,2) = [  0.0d0,  0.0d0,  0.0d0]  !A
        kpoints(:,3) = [0.00d0,  0.0d0,  0.1d0]  !-H
    case default
        stop "Error: Select axis."
    end select
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
!---- emesh
    de = (emax-emin)/eres
    do i=1, eres
    
        epoints(i) = emin + de*i
    enddo
!----Construct magnetic perturbation
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

    interp_Hr=0d0
!----Read in Hamiltonian
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

        buffsize1(i) = kloclist(i)*3*eres
        sum1 = sum1+ buffsize1(i)

        buff_sum1(i+1) = sum1

        displs1(i) = buff_sum1(i) - 1

        ! !For the MPI_GATHERV(super_H) call
        ! buffsize1(i) = (kpmaxlist(i)-kpminlist(i)+1)*matsize
        ! sum1 = sum1+ buffsize1(i)

        ! buff_sum1(i+1) = sum1

        ! displs1(i) = buff_sum1(i) - 1

        ! !For the MPI_GATHERV(evals) call
        ! buffsize2(i) = (kpmaxlist(i)-kpminlist(i)+1)*dim !For the MPI_GATHERV(eval) call
        ! sum2 = sum2+ buffsize2(i)

        ! buff_sum2(i+1) = sum2

        ! displs2(i) = buff_sum2(i) -1
    enddo

    recvcounts1 = buffsize1
    recvcounts2 = buffsize2
    kloc = kloclist(myid+1) 

    allocate(SH(dim,dim,nk),super_H(dim,dim),sH_flat(buffsize1(myid+1)),sH_flatg(matsize*nk)&
                                ,eval(dim),eval_flat(dim*kloc),eval_flatg(dim*nk))

    if(myid.eq.0) print *, "kloc:", kloclist
    if(myid.eq.0) print *, "kpmin:", kpminlist
    if(myid.eq.0) print *, "kpmax:", kpmaxlist
    if(myid.eq.0) print *, "kloc_sum:", kloc_sum
    if(myid.eq.0) print *, "buffsize1:", buffsize1
    ! if(myid.eq.0) print *, "buffsize2:", buffsize2
    if(myid.eq.0) print *, "buff_sum1:", buff_sum1, displs1
    ! if(myid.eq.0) print *, "buff_sum2:", buff_sum2, displs2

!----DX Files
    if(myid.eq.0) then
        do i = 1,5
            write(100+(i-1), '(a,3(1x,i8))') 'object 1 class gridpositions counts',nk,eres
            write(100+(i-1), '(a,3(1x,f12.8))') 'origin',-0.1d0,emin
            write(100+(i-1), '(a,3(1x,f12.8))') 'delta',sqrt(dot_product(dk,dk)),0d0
            write(100+(i-1), '(a,3(1x,f12.6))') 'delta',0d0,de
            write(100+(i-1), '(a,3(1x,i8))') 'object 2 class gridconnections counts',nk,eres
            write(100+(i-1), '(a,i8,a,i8,a,i10,a)') 'object',3,' class array type float rank 1 shape',3,' item', nk*eres, ' data follows'
        enddo
    endif


!----- Axis selection
    extrarow=0d0
    if(ax == 'x') index = 1
    if(ax == 'y') index = 2
    if(ax == 'z') index = 3

    ! extrarow(3:4) = +2*passval
    ! extrarow(nblocks*nb-3:nblocks*nb-2) = -3*passval

!----- Perform fourier transform

    allocate(data_row(5,3,kloc*eres))

    count = 0
    ! endif
    ikp=0
    do ik=kpminlist(myid+1),kpmaxlist(myid+1)
        if(myid==0) print *,'Iter:',ik,'/',kloc
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
                    super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1))) = Hkra(:,:,ra,rb)
                else
                    super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1))) = 0d0
                endif
            enddo
        enddo

        call zheev('V','U',dim,super_H(:,:),dim,eval(:),work,lwork,rwork,info)
        
        count = 0
        do ib =1, N2
            if(ib.eq.(nblocks_2 + 1).or. ib.eq.((nblocks_2)*nblocks +1).or. ib.eq.((nblocks_2)*(nblocks+1) +1).or. ib.eq.((nblocks_2)*(nblocks+2) +1).or. ib.eq.((nblocks_2)*(2*nblocks+1)+1)) then
                ! if(myid==0) print*, ib
                do ie=1,eres
                    ! print*, ie, eres
                    a_spec = 0d0
                    do i=1,dim
                        p_l = dot_product(super_H((1+nb*(ib-1)):(nb*(ib)),i), super_H((1+nb*(ib-1)):(nb*(ib)),i))
                        factor = ((epoints(ie) - eval(i)))/eta
                        a_spec = a_spec + p_l * (exp(-0.5d0*factor**2)) * 1/(eta*sqrt(2*pi2))
                        ! print *, a_spec
                    enddo
                    data_row(count+1, 1,(ikp-1)*eres + ie) = xk(ik)
                    data_row(count+1, 2,(ikp-1)*eres + ie) = epoints(ie)
                    data_row(count+1, 3,(ikp-1)*eres + ie) = real(a_spec) ! Top surface
                    ! if(myid==0) print*, data_row(count+1,:,(ikp-1)*eres + ie), count, ik, ib
                enddo
                count = count+1
            endif
        enddo
        if(myid.eq.0) then
            call date_and_time(VALUES=time_next)
            if(ikp.gt.1) then 
                dur =((time_next(6)-time_prev(6))*60 +  (time_next(7)-time_prev(7)))
                print *, 'Duration: ', dur/60, ':', mod(dur,60)
                print *, 'Remaining: ', dur*(kloc-ikp)/60 ,':', mod(dur*(kloc-ikp),60)
                print *, 'Time: ', time_next(5), ':', time_next(6), ':', time_next(7)
                print *,''
            endif
            call date_and_time(VALUES=time_prev)
        endif
    enddo

    do i=1,5
        call MPI_GATHERV(data_row(i,:,:),3*kloc*eres,MPI_DOUBLE_PRECISION, &
                        data_g(i,:,:),recvcounts1,displs1,MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD,IERR)
    enddo

    if(myid.eq.0) then
        call date_and_time(VALUES=time_end)
        dur =((time_end(5)-time_start(5))*3600 + (time_end(6)-time_start(6))*60 +  (time_end(7)-time_start(7)))

        print* , ''
        print *, 'Start time: ', time_start(5), ':', time_start(6), ':', time_start(7)
        print *, 'End time: ', time_end(5), ':', time_end(6), ':', time_end(7)
        print *, 'Duration: ', dur/3600,':', mod(dur,3600)/60, ':', mod(dur,60)
    
        do i=1,5
            do j=1,nk*eres
                write(100+(i-1),'(3(1x,f12.6))') data_g(i,:,j)
            enddo
            write(100+(i-1),*) 'object "regular positions regular connections" class field'
            write(100+(i-1),*) 'component "positions" value 1'
            write(100+(i-1),*) 'component "connections" value 2'
            write(100+(i-1),*) 'component "data" value 3'
            write(100+(i-1),*) 'end'
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

