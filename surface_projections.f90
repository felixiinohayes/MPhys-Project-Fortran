module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kxmax=0.03,kymax=0.03,a=0.79858
    integer,parameter::xmeshres=5,ymeshres=5,nkxpoints=(2*xmeshres+1),nkypoints=(2*ymeshres+1),nbmin=12,nbmax=13,nkp2=nkxpoints*nkypoints,nblocks=10,nr3=11
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    !INCLUDE 'mpif.h'
!------------------------------------------------------
    real*8 dx,dy,dz,da
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir,ir3,jr3,ir12,nr12,r3
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2
    real*8 avec(3,3),bvec(3,3),kpoint(2,nkp2),rvec_data(3)
    real*8,allocatable:: rvec(:,:),rwork(:)
    real*8, allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hkr3(:,:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),super_H(:,:,:)
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
    if(myid.eq.0) then
        open(100,file='super_H.dat')
    endif
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),Hk(nb,nb),Hkr3(nb,nb,nr3,nkp2),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr),super_H(nb*nblocks,nb*nblocks,nkp2))
    read(99,*)ndeg

    do i=1,80
      read(97,*)
    enddo
    do k=1,nr
		do i=1,nb
			do j=1,nb
			   read(99,*)rvec_data(1),rvec_data(2),rvec_data(3),i1,i2,x1,y1
			   top_Hr(i1,i2,k)=dcmplx(x1,y1)
			   read(97,*)rvec_data(1),rvec_data(2),rvec_data(3),j1,j2,x2,y2
			   triv_Hr(j1,j2,k)=dcmplx(x2,y2)
			enddo
		enddo
		rvec(:,k) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2) !+ rvec_data(3)*avec(:,3)
    enddo

    lwork=max(1,2*nb-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

    dx=kxmax/xmeshres
	dy=kymax/ymeshres

  !----- Create a uniform k-mesh
      ik=0
    do ikx=-xmeshres,xmeshres
      do iky=-ymeshres,ymeshres
		  ik=ik+1
		  kpoint(1,ik)=ikx*dx + kxmax
		  kpoint(2,ik)=iky*dy + kymax
      enddo
    enddo

!----- Perform fourier transform
	nr12=nr/nr3
	ikp=0
	do ik=1,nkp2
		do ir3=1,nr3	
			Hk=0d0	
			do ir12=0,nr12-1 ! Loop over (R1,R2) vectors
				ir = ir3 + ir12*11 ! Calculate index of (R1,R2) vector in nr
				phase = dot_product(kpoint(:,ik),rvec(:,ir))
				Hk=Hk+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
			enddo
			Hkr3(:,:,ir3,ik) = Hk
		enddo

		do i=0,nblocks-1
			do j=0,nblocks-1
				r3 = i - j
				do i2=1,nb
					do j2=1,nb
						super_H(18*i + i2, 18*j + j2, ik) = Hkr3(i2, j2, r3 + 6, ik)
						write(100, '(i8,2(1x,f12.6))') ik, super_H(18*i + i2, 18*j + j2, ik)
					enddo
				enddo
			enddo
		enddo
	enddo

!----- Construct supercell Hamiltonian



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
