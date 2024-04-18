module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.00015,a=0.791
    integer,parameter::meshres=35,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints,cp=meshres+1
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
	!INCLUDE 'mpif.h'
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,ikp,ir,node,pair
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkpoints,nkpoints,nkpoints),rvec_data(3),dV(3),offset(3,2,5),normal(3),v(3,3,nkpoints,nkpoints,nkpoints),v2xv3(3),total_c(3)
	real*8 dAdx(3,2),dAdy(3,2),dAdz(3,2)
	complex*8 spinor(2,2),H_2(2,2,2)
	real*8,allocatable:: rvec(:,:),rwork(:)
    real*8,allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:),eff(:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),H_k(:,:,:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eig_eff(:,:,:,:,:),eig(:,:,:,:,:),H_eff(:,:),dHdK(:,:,:)
	complex*16,allocatable::U(:,:,:,:,:),curvature(:,:,:,:,:),du(:,:,:)
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
    open(100,file='curvature.dx')
	open(200,file='connection.dx')
	open(300,file='curvature.dat')
	open(400,file='eigenvalues.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(3,nr),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr))
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
	   rvec(:,k) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2) + rvec_data(3)*avec(:,3)
    enddo

    lwork=max(1,2*nb-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))


!-----BTPs:

	node = 1
	pair = 1

	offset(:,1,1) = (/-0.017659606952654991,0.046513917396043679,0.43965460613976798/) !+ve
	offset(:,2,1) = (/ 0.017665681958398235,0.046638430945586576,0.47514974714462382/) !-ve

	offset(:,1,2) = (/ 0.04938161634772854, -0.00675093566736830,0.43883900501076140/) !+ve
	offset(:,2,2) = (/ 0.049353657200408851,0.0069912476907357081,0.43897293203235604/) !-ve
	
	offset(:,1,3) = (/ 0.01766767357940942,-0.04650968222518360,0.47482330199264761/) !+ve
	offset(:,2,3) = (/-0.01765908461242609,-0.04663940372172213,0.43932749486813732/) !-ve

	offset(:,1,4) = (/-0.04936478978846845,0.006832977105296199,0.47559337269211072/) !+ve
	offset(:,2,4) = (/-0.04933437400371396,-0.00708356750361176,0.47545289372516780/) !-ve

	offset(:,1,5) = (/-0.00879804124970561,0.04855142510428899,0.44745395357167184/)  !Dirac Point


    dk=kmax/(nkpoints-1)

!----- Create header of dx files

	write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints-2,nkpoints-2,nkpoints-2
	write(100, '(a,3(1x,f12.8))') 'origin',offset(1,node,pair)-(kmax/2),offset(2,node,pair)-(kmax/2),offset(3,node,pair)-(kmax/2)
	write(100, '(a,3(1x,f12.8))') 'delta',dk,0d0,0d0
	write(100, '(a,3(1x,f12.8))') 'delta',0d0,dk,0d0
	write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dk
	write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints-2,nkpoints-2,nkpoints-2
	write(100, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',3,&
								   ' item', (nkpoints-2)*(nkpoints-2)*(nkpoints-2), ' data follows'
	write(200, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints-1,nkpoints-1,nkpoints-1
	write(200, '(a,3(1x,f12.8))') 'origin',offset(1,node,pair)-(kmax/2),offset(2,node,pair)-(kmax/2),offset(3,node,pair)-(kmax/2)
	write(200, '(a,3(1x,f12.8))') 'delta',dk,0d0,0d0
	write(200, '(a,3(1x,f12.8))') 'delta',0d0,dk,0d0
	write(200, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dk
	write(200, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints-1,nkpoints-1,nkpoints-1
	write(200, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',3,&
									' item', (nkpoints-1)*(nkpoints-1)*(nkpoints-1), ' data follows'
	write(400, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkpoints
	write(400, '(a,3(1x,f12.8))') 'origin',offset(1,node,pair)-(kmax/2),offset(2,node,pair)-(kmax/2),offset(3,node,pair)-(kmax/2)
	write(400, '(a,3(1x,f12.8))') 'delta',dk,0d0,0d0
	write(400, '(a,3(1x,f12.8))') 'delta',0d0,dk,0d0
	write(400, '(a,3(1x,f12.8))') 'delta',0d0,0d0,dk
	write(400, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkpoints
	write(400, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',1,&
									' item', (nkpoints)*(nkpoints)*(nkpoints), ' data follows'
    
  !----- Create a uniform k-mesh
	ik=0	
	do ikx=1,nkpoints	
		do iky=1,nkpoints
			do ikz=1,nkpoints
				kpoint(1,ikx,iky,ikz) = (ikx-1)*dk + offset(1,node,pair) - kmax/2
				kpoint(2,ikx,iky,ikz) = (iky-1)*dk + offset(2,node,pair) - kmax/2
				kpoint(3,ikx,iky,ikz) = (ikz-1)*dk + offset(3,node,pair) - kmax/2
			enddo
		enddo
		! print *, kpoint(1,ikx,1,1)
	enddo

	! kpool=nkp3/numprocs
    ! if (mod(nkp3,numprocs).ne.0) kpool=kpool+1

    ! kpmin=1+myid*kpool
    ! kpmax=(myid+1)*kpool

    ! ecounts=kpool*3 !
	
	allocate(HK(nb,nb),H_k(nb,nb,ikx,iky,ikz),k_ene(nb),eig(nb,2,nkpoints,nkpoints,nkpoints))
    allocate(H_eff(2,2),eig_eff(2,2,nkpoints,nkpoints,nkpoints),eff(2),du(nb,2,3))
	allocate(U(3,2,nkpoints,nkpoints,nkpoints),curvature(3,2,nkpoints,nkpoints,nkpoints))

	print *, "Finding eigenvectors"
!----- Fourier Transform
	ikp=0
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				print*, ikp, '/', nkp3
				ikp=ikp+1
				Hk = 0d0
				do ir=1,nr
					phase = dot_product(kpoint(:,ikx,iky,ikz),rvec(:,ir))
					HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
				enddo
				H_k(:,:,ikx,iky,ikz) = HK(:,:) ! Store Hamiltonian
				call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)

!----12th and 13th bands' eigenvectors
				eig(:,1,ikx,iky,ikz) = HK(:,12)
				eig(:,2,ikx,iky,ikz) = HK(:,13)

				write(400,'(1(1x,f10.8))') k_ene(12)

!----Constructing effective Hamiltonian and diagonalizing
				! do i=1,2
				! 	do j=1,2
				! 		H_eff(i,j) = dot_product(eig(:,i,ikx,iky,ikz),matmul(H_k(:,:,ikx,iky,ikz),eig(:,j,ikx,iky,ikz)))
				! 	enddo
				! enddo
				! call zheev('V','U',2,H_eff,2,eff,work,lwork,rwork,info)

				! eig_eff(:,1,ikx,iky,ikz) = H_eff(:,1)
				! eig_eff(:,2,ikx,iky,ikz) = H_eff(:,2)
			enddo
		enddo
	enddo

!----- Berry connection
	print *, "Computing Berry connection"
	do ikx=1,nkpoints-1
		do iky=1,nkpoints-1
			do ikz=1,nkpoints-1
				! do i=1,2 !--Accounts for the 2 eigenvectors of the effective hamiltonian

				U(1,1,ikx,iky,ikz) = (dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx+1,iky,ikz)))/(abs(dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx+1,iky,ikz))))
				U(2,1,ikx,iky,ikz) = (dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx,iky+1,ikz)))/(abs(dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx,iky+1,ikz))))
				U(3,1,ikx,iky,ikz) = (dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx,iky,ikz+1)))/(abs(dot_product(eig(:,1,ikx,iky,ikz),eig(:,1,ikx,iky,ikz+1))))
				! enddo

				! total_c(:) = U(:,1,ikx,iky,ikz)*dcmplx(0d0,-1d0) 

				! write(100, '(3(1x,f20.2))') total_c
				! write(300, '(3(1x,f12.8),3(1x,f20.2))') kpoint(:,ikx,iky,ikz),total_c
			enddo
		enddo
	enddo

!----- Berry curvature
	print *, "Computing Berry curvature"
	total_c = 0d0
	do ikx=1,nkpoints-2
		do iky=1,nkpoints-2
			do ikz=1,nkpoints-2
				! do i =1,2

				curvature(1,1,ikx,iky,ikz) = log((U(2,1,ikx,iky,ikz)*U(3,1,ikx,iky+1,ikz))/(U(2,1,ikx,iky,ikz+1)*U(3,1,ikx,iky,ikz)))/(dk*dk) ! F_23
				curvature(2,1,ikx,iky,ikz) = log((U(3,1,ikx,iky,ikz)*U(1,1,ikx,iky,ikz+1))/(U(3,1,ikx+1,iky,ikz)*U(1,1,ikx,iky,ikz)))/(dk*dk) ! F_31
				curvature(3,1,ikx,iky,ikz) = log((U(1,1,ikx,iky,ikz)*U(2,1,ikx+1,iky,ikz))/(U(1,1,ikx,iky+1,ikz)*U(2,1,ikx,iky,ikz)))/(dk*dk) ! F_12

				! enddo
				! print*, curvature(:,1,ikx,iky,ikz)
				! sum(:) = curvature(:,1,ikx,iky,ikz) + curvature(:,2,ikx,iky,ikz)
				
				write(100, '(3(1x,f20.5))') aimag(curvature(:,1,ikx,iky,ikz))
				write(300, '(3(1x,f12.8),3(1x,f20.5))') kpoint(:,ikx,iky,ikz),aimag(curvature(:,1,ikx,iky,ikz))
			enddo	
		enddo
	enddo

	deallocate(H_k,Hk,eig,eig_eff)
	write(100,'(A,/,A,/,A,/,A)') &
		'object "regular positions regular connections" class field', &
		'component "positions" value 1', &
		'component "connections" value 2', &
		'component "data" value 3', &
		'end'
	write(200,'(A,/,A,/,A,/,A)') &
		'object "regular positions regular connections" class field', &
		'component "positions" value 1', &
		'component "connections" value 2', &
		'component "data" value 3', &
		'end'
	write(400,'(A,/,A,/,A,/,A)') &
		'object "regular positions regular connections" class field', &
		'component "positions" value 1', &
		'component "connections" value 2', &
		'component "data" value 3', &
		'end'
	!call MPI_FINALIZE( IERR )
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
