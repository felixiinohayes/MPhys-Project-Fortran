module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kxmax=0.05,kymax=0.05,a=0
    integer,parameter::xmeshres=15,ymeshres=15,nkxpoints=(2*xmeshres+1),nkypoints=(2*ymeshres+1),nbmin=1,nbmax=18,nkp=nkxpoints*nkypoints,nblocks=10 
    integer,parameter::nr1=6,nr2=4,nr3=5
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
    integer*4 i,j,k,nr,i1,i2,j1,j2,r1,r2,r3,lwork,info,ikx,iky,ik,ir,ir1,ir2,ir3
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,sumtotal
    real*8 avec(3,3),bvec(3,3),kpoint(2,nkp)
    real*8,allocatable:: rvec(:,:,:),rwork(:),rvec_data(:,:),k_ene(:)
	integer*4,allocatable:: ndeg(:),super_H_index(:,:)
    complex*16,allocatable:: Hk(:,:),Hkr3(:,:,:),top_Hr(:,:,:,:,:),triv_Hr(:,:,:,:,:),work(:),super_H(:,:)
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
    open(200,file='top_surface_ene.dx')
    open(300,file='bottom_surface_ene.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,-nr1:nr1,-nr2:nr2),rvec_data(3,nr),Hk(nb,nb),ndeg(nr),super_H(nb*nblocks,nb*nblocks),k_ene(nb*nblocks))
    allocate(Hkr3(nb,nb,-nr3:nr3),top_Hr(nb,nb,-nr1:nr1,-nr2:nr2,-nr3:nr3),triv_Hr(nb,nb,-nr1:nr1,-nr2:nr2,-nr3:nr3))
    read(99,*)ndeg

    do i=1,80
      read(97,*)
    enddo
    do ir=1,nr
		do i=1,nb
			do j=1,nb
			   read(99,*)r1,r2,r3,i1,i2,x1,y1
               rvec_data(1,ir)=r1
               rvec_data(2,ir)=r2
               rvec_data(3,ir)=r3
               
			   top_Hr(i1,i2,r1,r2,r3)=dcmplx(x1,y1)
			   read(97,*)r1,r2,r3,j1,j2,x2,y2
			   triv_Hr(j1,j2,r1,r2,r3)=dcmplx(x2,y2)

               rvec(:,r1,r2) = r1*avec(:,1) + r2*avec(:,2)
			enddo
		enddo
    enddo

    lwork=max(1,2*(nb*nblocks)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*nblocks)-2)))

    dx=kxmax/xmeshres
	dy=kymax/ymeshres

	write(200, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkxpoints,nkypoints
    write(200, '(a,3(1x,f12.6))') 'origin',0d0,0d0
    write(200, '(a,3(1x,f12.6))') 'delta',dx,0d0
    write(200, '(a,3(1x,f12.6))') 'delta',0d0,dy
    write(200, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkxpoints,nkypoints
	write(200, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape', nbmax-nbmin+1, &
                                     ' item', nkp,' data follows'
	write(300, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkxpoints,nkypoints
    write(300, '(a,3(1x,f12.6))') 'origin',0d0,0d0
    write(300, '(a,3(1x,f12.6))') 'delta',dx,0d0
    write(300, '(a,3(1x,f12.6))') 'delta',0d0,dy
    write(300, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkxpoints,nkypoints
	write(300, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape', nbmax-nbmin+1, &
                                    ' item', nkp,' data follows'

  !----- Create a uniform k-mesh
    ik=0
    do ikx=-xmeshres,xmeshres
      do iky=-ymeshres,ymeshres
		  ik=ik+1
		  kpoint(1,ik)=ikx*dx
		  kpoint(2,ik)=iky*dy
      enddo
    enddo

!----- Perform fourier transform

	do ik=1,nkp
		do ir3= -nr3,nr3 ! Loop over R3 vectors
			Hk=0d0	
			do ir1= -nr1,nr1 ! Loop over (R1,R2) vectors
                do ir2= -nr2,nr2
                    phase = dot_product(kpoint(:,ik),rvec(:,r1,r2))
                    Hk=Hk+((1-a)*(triv_Hr(:,:,ir1,ir2,ir3))+(a)*(top_Hr(:,:,ir1,ir2,ir3)))*dcmplx(cos(phase),-sin(phase))
                enddo
			enddo
			Hkr3(:,:,ir3) = Hk
		enddo

	    do i=0,nblocks-1
			do j=0,nblocks-1
				r3 = i-j
				if (r3<6 .AND. r3>-6) then
					super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1))) = Hkr3(:,:,r3)
				else
					super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1))) = 0d0
				endif
			enddo
		enddo
		call zheev('V','U',nb*nblocks,super_H,nb*nblocks,k_ene,work,lwork,rwork,info)
		write(200, '(10(1x,f12.6))') k_ene(nbmin:nbmax) ! Top surface
		write(300, '(10(1x,f12.6))') k_ene((nb*(nblocks-1))+nbmin:(nb*(nblocks-1))+nbmax) ! Bottom surface
	enddo

    do i=1,nb*nblocks
        do j=1,nblocks
            print *, super_H(i,j)
        enddo
    enddo 

	write(200,'(A,/,A,/,A,/,A)') &
    'object "regular positions regular connections" class field', &
    'component "positions" value 1', &
    'component "connections" value 2', &
    'component "data" value 3', &
    'end'
	write(300,'(A,/,A,/,A,/,A)') &
    'object "regular positions regular connections" class field', &
    'component "positions" value 1', &
    'component "connections" value 2', &
    'component "data" value 3', &
    'end'

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
