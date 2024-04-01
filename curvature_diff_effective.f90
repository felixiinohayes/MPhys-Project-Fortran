module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.00035,a=0.791,diff_value=kmax/10
    integer,parameter::meshres=5,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints,cp=meshres+1
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
    real*8 avec(3,3),bvec(3,3),kpoint_o(3,nkpoints,nkpoints,nkpoints),kpoint(4,4,3,nkpoints,nkpoints,nkpoints),rvec_data(3),dV(3),offset(3,2,5),normal(3),v(3,3,nkpoints,nkpoints,nkpoints),v2xv3(3),total_c(3)
	real*8 dAdx(3),dAdy(3),dAdz(3),diff(3,4),sum(3)
	complex*8 spinor(2,2),H_2(2,2,2)
	real*8,allocatable:: rvec(:,:),rwork(:)
    real*8,allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:),eff(:),curvature(:,:,:,:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),H_k(:,:,:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eig_eff(:,:,:,:,:,:,:),eig(:,:,:,:,:,:,:),H_eff(:,:),dHdK(:,:,:)
	complex*16,allocatable::connection(:,:,:,:,:,:),du(:,:,:)
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
    open(100,file='curvature20.dx')
	open(200,file='curvature20.dat')
	open(300,file='eigenvalues20.dx')
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

	node = 2
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

	write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkpoints
	write(100, '(a,3(1x,f12.8))') 'origin',offset(1,node,pair)-(kmax/2),offset(2,node,pair)-(kmax/2),offset(3,node,pair)-(kmax/2)
	write(100, '(a,3(1x,f12.8))') 'delta',dk,0d0,0d0
	write(100, '(a,3(1x,f12.8))') 'delta',0d0,dk,0d0
	write(100, '(a,3(1x,f12.8))') 'delta',0d0,0d0,dk
	write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkpoints
	write(100, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',3,&
								   ' item', nkp3, ' data follows'
	write(300, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkpoints
	write(300, '(a,3(1x,f12.8))') 'origin',offset(1,node,pair)-(kmax/2),offset(2,node,pair)-(kmax/2),offset(3,node,pair)-(kmax/2)
	write(300, '(a,3(1x,f12.8))') 'delta',dk,0d0,0d0
	write(300, '(a,3(1x,f12.8))') 'delta',0d0,dk,0d0
	write(300, '(a,3(1x,f12.8))') 'delta',0d0,0d0,dk
	write(300, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkpoints
	write(300, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',1,&
									' item', nkp3, ' data follows'
	diff(:,1) = (/0d0,0d0,0d0/)
	diff(:,2) = (/diff_value,0d0,0d0/)
	diff(:,3) = (/0d0,diff_value,0d0/)
	diff(:,4) = (/0d0,0d0,diff_value/)
    
  !----- Create a uniform k-mesh
	ik=0	
	do ikx=1,nkpoints	
		do iky=1,nkpoints
			do ikz=1,nkpoints
				! print *, (ikz-1)*dk + offset(1,node,pair) - kmax/2, ikx,iky,ikz
				kpoint(1,1,1,ikx,iky,ikz) = (ikx-1)*dk + offset(1,node,pair) - kmax/2
				kpoint(1,1,2,ikx,iky,ikz) = (iky-1)*dk + offset(2,node,pair) - kmax/2
				kpoint(1,1,3,ikx,iky,ikz) = (ikz-1)*dk + offset(3,node,pair) - kmax/2

				! Create shifted k-mesh for finite difference
				do i=2,4
					kpoint(1,i,:,ikx,iky,ikz) = kpoint(1,1,:,ikx,iky,ikz) + diff(:,i) ! 1st derivative
					do j=1,4
						kpoint(i,j,:,ikx,iky,ikz) = kpoint(1,i,:,ikx,iky,ikz) + diff(:,j) ! 2nd derivative
					enddo
				enddo
				! print *, kpoint(4,1,:,2,3,1), kpoint(1,4,:,2,3,1)
			enddo
		enddo
	enddo
	print *, kpoint(4,:,:,2,3,1)

	! kpool=nkp3/numprocs
    ! if (mod(nkp3,numprocs).ne.0) kpool=kpool+1

    ! kpmin=1+myid*kpool
    ! kpmax=(myid+1)*kpool

    ! ecounts=kpool*3 !
	
	allocate(HK(nb,nb),H_k(nb,nb,ikx,iky,ikz),k_ene(nb),eig(4,4,nb,2,nkpoints,nkpoints,nkpoints))
    allocate(H_eff(2,2),eig_eff(4,4,2,2,nkpoints,nkpoints,nkpoints),eff(2),du(2,2,3))
	allocate(connection(4,3,2,nkpoints,nkpoints,nkpoints),curvature(3,2,ikx,iky,ikz))

	print *, "Finding eigenvectors"
!----- Fourier Transform
	ikp=0
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				do k=1,4
					do ik=1,4
						print*, ikp, '/', nkp3*16
						ikp=ikp+1
						Hk = 0d0
						do ir=1,nr
							phase = dot_product(kpoint(k,ik,:,ikx,iky,ikz),rvec(:,ir))
							HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
						enddo
						H_k(:,:,ikx,iky,ikz) = HK(:,:) ! Store Hamiltonian
						call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)

		!----12th and 13th bands' eigenvectors
						eig(k,ik,:,1,ikx,iky,ikz) = HK(:,12)
						eig(k,ik,:,2,ikx,iky,ikz) = HK(:,13)

		!----Constructing effective Hamiltonian and diagonalizing
						do i=1,2
							do j=1,2
								H_eff(i,j) = dot_product(eig(k,ik,:,i,ikx,iky,ikz),matmul(H_k(:,:,ikx,iky,ikz),eig(k,ik,:,j,ikx,iky,ikz)))
							enddo
						enddo
						call zheev('V','U',2,H_eff,2,eff,work,lwork,rwork,info)

						eig_eff(k,ik,:,1,ikx,iky,ikz) = H_eff(:,1)
						eig_eff(k,ik,:,2,ikx,iky,ikz) = H_eff(:,2)
					enddo
				enddo
			enddo
		enddo
	enddo
	! print *, kpoint(4,:,:,2,3,1)

!----- Berry connection
	print *, "Computing Berry connection"
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				do ik=1,4
					do i=1,2
						du(:,i,1) = (eig_eff(ik,2,i,:,ikx,iky,ikz) - eig_eff(ik,1,i,:,ikx,iky,ikz))/diff_value
						du(:,i,2) = (eig_eff(ik,3,i,:,ikx,iky,ikz) - eig_eff(ik,1,i,:,ikx,iky,ikz))/diff_value
						du(:,i,3) = (eig_eff(ik,4,i,:,ikx,iky,ikz) - eig_eff(ik,1,i,:,ikx,iky,ikz))/diff_value

						connection(ik,1,i,ikx,iky,ikz) = dot_product(eig_eff(ik,1,i,:,ikx,iky,ikz),du(:,i,1))
						connection(ik,2,i,ikx,iky,ikz) = dot_product(eig_eff(ik,1,i,:,ikx,iky,ikz),du(:,i,2))
						connection(ik,3,i,ikx,iky,ikz) = dot_product(eig_eff(ik,1,i,:,ikx,iky,ikz),du(:,i,3))
					enddo
				enddo
				sum = connection(1,:,1,ikx,iky,ikz) + connection(1,:,2,ikx,iky,ikz)
				! write(100, '(3(1x,f20.8))') sum(:)
				! print *, eig(4,:,:,ikx,iky,ikz), ikx,iky,ikz
			enddo
		enddo
	enddo

!----- Berry curvature
	print *, "Computing Berry curvature"
	total_c = 0d0
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				do i=1,2
					dAdx(:) = (connection(2,:,i,ikx,iky,ikz) - connection(1,:,i,ikx,iky,ikz))/diff_value
					dAdy(:) = (connection(3,:,i,ikx,iky,ikz) - connection(1,:,i,ikx,iky,ikz))/diff_value
					dAdz(:) = (connection(4,:,i,ikx,iky,ikz) - connection(1,:,i,ikx,iky,ikz))/diff_value
					! print *, dAdx(2)-dAdy(1),ikx,iky,ikz

					curvature(1,i,ikx,iky,ikz) = real((dAdy(3) - dAdz(2))) ! Compute curl
					curvature(2,i,ikx,iky,ikz) = real((dAdz(1) - dAdx(3)))
					curvature(3,i,ikx,iky,ikz) = real((dAdx(2) - dAdy(1)))			
				enddo
				! print *, connection(4,:,ikx,iky,ikz)
				sum = curvature(:,1,ikx,iky,ikz) + curvature(:,2,ikx,iky,ikz)

				! print*,curvature(3,1,ikx,iky,ikz),curvature(3,2,ikx,iky,ikz),sum(3)
				
				write(100, '(3(1x,f20.8))') sum(:) * diff_value
				! write(200, '(3(1x,f12.8),3(1x,f20.2))') kpoint(1,ikx,iky,ikz),kpoint(2,ikx,iky,ikz),kpoint(3,ikx,iky,ikz),sum(:)
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
	write(300,'(A,/,A,/,A,/,A)') &
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



! write(*, '(f10.5,3I5)') real(eig(1,1,ikx,iky,ikz)), ikx, iky, ikz
! write(*, '(2f10.5,3I5)') aimag(eig(1,2,ikx,iky,ikz)), aimag(eig(2,2,ikx,iky,ikz)), ikx, iky, ikz
! write(*, '(f10.7,3I5)') real(dot_product(eig(:,2,ikx,iky,ikz),matmul(H_k(:,:,ikx,iky,ikz),eig(:,1,ikx,iky,ikz)))),ikx,iky,ikz


! print *, eig_eff(1,1,ikx,iky,ikz), ikx, iky, ikz
! write(*, '(2f10.5,3I5)') real(eig_eff(1,1,ikx,iky,ikz)),aimag(eig_eff(1,1,ikx,iky,ikz)), ikx, iky, ikz
! write(*, '(2f10.5,3I5)') aimag(eig_eff(1,1,ikx,iky,ikz)), aimag(eig_eff(1,2,ikx,iky,ikz)), ikx, iky, ikz
! if(ikx == cp.and.iky == cp.and.ikz == cp) print*, eff(1),eff(2)
! H_eff(i,j) = dot_product(eig(:,i,cp,cp,cp),matmul(H_k(:,:,ikx,iky,ikz),eig(:,j,cp,cp,cp)))
