module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,a=1,emin=5.5,emax=6.5,bfactor=0.006
    integer,parameter::nkpath=3,np=20,nblocks=10,nr3=11,nk=(nkpath-1)*np+1,eres=10
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    INCLUDE 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,ie,lwork,info,ik,count,ir,ir3,ir12,nr12,r3,sign,il,kpool,kpmin,kpmax,ecounts,ikp
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), B=0.00d0
    real*8 phase,pi2,x1,y1,x2,y2,de,exp_factor,p_l,spectral_A
    real*8 xk(nk-1),avec(3,3),bvec(3,3),rvec_data(3),kpoints(3,nkpath),kpath(3,nk),dk(3),epoints(eres),spectral_A_comm(3,nk*eres)
    real*8,allocatable:: rvec(:,:),rvec_miller(:,:),rwork(:),k_ene(:,:),spectral_A_single(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),Hkr3(:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),super_H(:,:),sH(:,:),a_p_top(:,:),a_p_bottom(:,:),B_pt(:,:)
    complex*16 B_sigma(2,2)
!------------------------------------------------------
    call init_mpi

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
    open(100,file='super_H_B06.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),rvec_miller(3,nr),Hk(nb,nb),Hkr3(nb,nb,nr3),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr))
    allocate(super_H(nb*nblocks,nb*nblocks),sH(nb,nb*nblocks),k_ene(nb*nblocks,nk))
    allocate(a_p_top(nb*nblocks,nk),a_p_bottom(nb*nblocks,nk))
    read(99,*)ndeg

    do i=1,80
      read(97,*)
    enddo
    do ir=1,nr
        do i=1,nb
            do j=1,nb
               read(99,*)rvec_data(1),rvec_data(2),rvec_data(3),i1,i2,x1,y1
               top_Hr(i1,i2,ir)=dcmplx(x1,y1)
               read(97,*)rvec_data(1),rvec_data(2),rvec_data(3),j1,j2,x2,y2
               triv_Hr(j1,j2,ir)=dcmplx(x2,y2)
            enddo
        enddo
        rvec(:,ir) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2)
    enddo

    lwork=max(1,2*(nb*nblocks)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*nblocks)-2)))

!-----kpath
    kpoints(:,1) = [ -0.15d0,  0.0d0,   0.0d0]  !-M
    kpoints(:,2) = [  0.0d0,  0.0d0,   0.0d0]  !Gamma
    kpoints(:,3) = [  0.15d0,  0.0d0,   0.0d0]  !M    
    
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
        do i = 1, np
            ik = i + np*(j-1)
            dk = (kpoints(:,j+1)-kpoints(:,j))/np
            kpath(:, ik) = kpoints(:,(j)) + (dk*(i-1))
            xk(ik) =  sign*sqrt(dot_product(kpoints(:,2)- kpath(:, ik),kpoints(:,2) - kpath(:, ik)))
        enddo
    enddo
    de = (emax-emin)/eres
    ! if(myid.eq.0) print *, kpath(:,:)

    do i=1, eres
        epoints(i) = emin + de*i
    enddo

    kpool=nk/numprocs
    if (mod(nk,numprocs).ne.0) kpool=kpool+1

    kpmin=1+myid*kpool
    kpmax=(myid+1)*kpool
    ! print *, kpmin, kpmax, kpool, numprocs

    ecounts=kpool*eres ! to account for bands 12 and 13

!----Construct magnetic perturbation
    ! call write_header()
    if(myid.eq.0) then
        write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nk,eres
        write(100, '(a,2(1x,f12.6))') 'origin',-0.1d0,emin
        write(100, '(a,2(1x,f12.6))') 'delta',sqrt(dot_product(dk,dk)),0d0
        write(100, '(a,2(1x,f12.6))') 'delta',0d0,de
        write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nk,eres
    endif
    allocate(B_pt(nb, nb))

     !B along Y axis
    B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-B)]
    B_sigma(2,:) = [dcmplx(0d0,B) ,  dcmplx(0d0,0d0)]

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

    allocate(spectral_A_single(3,kpool*eres))
!----- Perform fourier transform
    nr12=nr/nr3
    do il=0,nblocks-1
        ! print *, myid
        if(myid.eq.0) write(100, '(a,i8,a,i10,a)') 'object',il+3,' class array type float rank 1 shape 3 item',nk*eres,' data follows'
        ikp=0
        do ik=kpmin,min(kpmax,nk)
            ikp=ikp+1
            do ir3=1,nr3 ! Loop over R3 vectors
                Hk=0d0    
                do ir12=0,nr12-1 ! Loop over (R1,R2) vectors
                    ir = ir3 + ir12*nr3 ! Calculate index of (R1,R2) vector in nr

                    phase = 0d0
                    do j = 1,2
                        phase = phase + kpath(j,ik)*rvec(j,ir)
                    enddo

                    Hk=Hk+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
                enddo
                Hk = Hk + B_pt
                Hkr3(:,:,ir3) = Hk
            enddo

            do i=0,nblocks-1
                do j=0,nblocks-1
                    r3 = i-j
                    if (r3<=5 .AND. r3>=-5) then
                        super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1))) = Hkr3(:,:,r3 + (nr3+1)/2)
                    else
                        super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1))) = 0d0
                    endif
                enddo
            enddo
            call zheev('V','U',nb*nblocks,super_H,nb*nblocks,k_ene(:,ik),work,lwork,rwork,info)

            do ie=1,eres
                spectral_A = 0d0
                do i=1,nb*nblocks
                    p_l = dot_product(super_H((1+nb*il):(nb*(il+1)),i),super_H((1+nb*il):(nb*(il+1)),i))
                    ! a_p_top(i,ik) = dot_product(super_H(1:18,i),super_H(1:18,i))
                    ! a_p_bottom(i,ik) = dot_product(super_H(nb*(nblocks-1):nb*nblocks,i),super_H(nb*(nblocks-1):nb*nblocks,i))
                    exp_factor = (epoints(ie) - k_ene(i,ik))/bfactor
                    ! print *, exp_factor
                    spectral_A = spectral_A + p_l * exp(-0.5d0 * (exp_factor**2))
                enddo
                spectral_A_single(1,(ikp-1)*eres + ie) = xk(ik)
                spectral_A_single(2,(ikp-1)*eres + ie) = epoints(ie)
                spectral_A_single(3,(ikp-1)*eres + ie) = real(spectral_A)! Top surface
                ! if (myid.eq.1) print *, (ikp-1)*eres + ie, kpool*eres
                if(myid.eq.0) print *, spectral_A_single(:,(ikp-1)*eres + ie),(ikp-1)*eres + ie
            enddo
            ! if(myid.eq.0) write(100, '(3(1x,f12.6))') spectral_A_single(:,ikp)
            
            ! if(myid.eq.0) print *, ikp, "/", nk*nblocks
        enddo
        call MPI_GATHER(spectral_A_single,ecounts,MPI_DOUBLE_PRECISION, &
                            spectral_A_comm,ecounts,MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD,IERR)
        do i=1,nk*eres
            ! if(myid.eq.0) write(100, '(3(1x,f12.6))') spectral_A_comm(:,i)
            ! if(myid.eq.0 .and. il==1) print *, spectral_A_comm(:,i)
        enddo
        ! print *, spectral_A_comm(3,:)
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
