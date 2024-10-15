module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,a=0,emin=5.5,emax=6.5,bfactor=0.006
    integer*8,parameter::nkpath=3,np=150,nr3=11,nk=(nkpath-1)*np+1,eres=5,nepoints=(2*eres)+1,nblocks=5,blocksize=(nblocks+1)**3
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
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), B=0.00d0
    real*8 phase,pi2,x1,y1,x2,y2,de,spectral_A,exp_factor,p_l,emiddle
    real*8 xk(nk),avec(3,3),bvec(3,3),kpoints(3,nkpath),kpath(3,nk),dk(3),epoints(nepoints),g_E(nepoints-1)
    real*8,allocatable:: rvec(:,:),rvec_miller(:,:),rwork(:),ene(:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),Hkr3(:,:,:),top_Hr(:,:),triv_Hr(:,:),work(:),super_H(:,:),sH(:,:),a_p_top(:,:),a_p_bottom(:,:),B_pt(:,:)
    complex*16,dimension(18,18,-6:6,-6:6,-6:6) :: interp_Hr
    complex*16 B_sigma(2,2)
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
    open(100,file='cube_eigenvalues_5_zheev.dat')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),rvec_miller(3,nr),Hk(nb,nb),Hkr3(nb,nb,nr3),top_Hr(nb,nb),triv_Hr(nb,nb),ndeg(nr))
    allocate(super_H(nb*blocksize,nb*blocksize),sH(nb,nb*nblocks),ene(nb*blocksize))
    allocate(a_p_top(nb*nblocks,nk),a_p_bottom(nb*nblocks,nk))
    read(99,*)ndeg

    do i=1,80
      read(97,*)
    enddo
    interp_Hr=0d0
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
    ! print *, interp_Hr(2,2,-1,1,1)

    lwork=max(1,2*(nb*blocksize)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*blocksize)-2)))

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

    emiddle = emin + (emax-emin)/2
    de = (emax-emin)/(2*eres)
    ie=0
    do i=-eres, eres
        ie=ie+1
        epoints(ie) = emiddle + de*i
    enddo


! !----Magnetic perturbation
     !B along Y axis
    ! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-B)]
    ! B_sigma(2,:) = [dcmplx(0d0,B) ,  dcmplx(0d0,0d0)]
    ! B_pt=0d0
    ! do i=1,nb
    !     do j=1,nb
    !         if (i==j) then
    !             if (i<10) then
    !                 B_pt(i,j) = B_sigma(1,1)
    !             else
    !                 B_pt(i,j) = B_sigma(2,2)
    !             endif
    !         else if (i==j+9) then
    !             B_pt(i,j) = B_sigma(2,1)
    !         else if (j==i+9) then
    !             B_pt(i,j) = B_sigma(1,2)
    !         endif
    !     enddo
    ! enddo

!----- Construct supercell hamiltonian

    super_H=0d0
    do i3=0,nblocks
        do j3=0,nblocks
            r3=i3-j3
            do i2=0,nblocks
                do j2=0,nblocks
                    r2=i2-j2
                    do i1=0,nblocks
                        do j1=0,nblocks
                            r1=i1-j1
                            xindex = i3*((nblocks+1)**2)+i2*(nblocks+1)+i1
                            yindex = j3*((nblocks+1)**2)+j2*(nblocks+1)+j1
                            super_H((1+nb*xindex):(nb*(xindex+1)),(1+nb*yindex):(nb*(yindex+1))) = interp_Hr(:,:,r1,r2,r3)
                            ! print *, xindex, yindex,r1,r2,r3
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    call zheev('V','U',nb*blocksize,super_H,nb*blocksize,ene,work,lwork,rwork,info)
    ! call dsaupd()

    ! g_E=0d0
    ! do i=1,nepoints-1
    !     do j=1,nb*blocksize
    !         if((ene(j) .gt. epoints(i)) .and. (ene(j) .lt. epoints(i+1))) g_E(i)=g_E(i)+1
    !     enddo
    ! enddo
    print *, ene

    write(100, '(f12.6)') ene


end Program Projected_band_structure

! SUBROUTINE INIT_MPI
!     USE PARAMETERS               ,             ONLY: IERR,MYID,NUMPROCS
!     IMPLICIT NONE
!     INCLUDE 'mpif.h'
!         Call MPI_INIT( IERR )
!         Call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
!         Call MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS , IERR )
! !        Write(*,*) ‘Process’, myid, ' of ’, NUMPROCS , ‘is alive.’
! END SUBROUTINE INIT_MPI
