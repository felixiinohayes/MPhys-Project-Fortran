module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.021,a=0.791
    integer,parameter::meshres=20,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir,choice
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,chern,div_F,diff_z
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkpoints,nkpoints,nkpoints),rvec_data(3),dV(3),offset(3,2),normal(3),v(3,3,nkpoints,nkpoints,nkpoints),v2xv3(3),sum(3)
	real*8 dAdx(3),dAdy(3),dAdz(3)
	complex*16 dn(18,3)
	real*8,allocatable:: rvec(:,:),rwork(:)
    real*8,allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),Hamk(:,:,:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eigenvector(:),eig_v(:,:,:,:,:),H_eff(:,:),dHdK(:,:,:)
	complex*16,allocatable::connection(:,:,:,:),curvature(:,:,:,:)
	!------------------------------------------------------

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

    dk=kmax/(nkpoints-1)

	choice = 1

	offset(:,1) = (/0.017665681958398235 -(kmax/2),0.046638430945586576-(kmax/2),0.47514974714462382-(kmax/2)/)
	offset(:,2) = (/-0.017659606952654991-(kmax/2),0.046513917396043679-(kmax/2),0.43965460613976798-(kmax/2)/)
	!offset(:,2) = (/0.0493647897886,-0.00683297710,0.43887604172172545/) 


!----- Create header of dx files

	write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints-1,nkpoints-1,nkpoints-1
	write(100, '(a,3(1x,f12.6))') 'origin',offset(1,choice),offset(2,choice),offset(3,choice)
	write(100, '(a,3(1x,f12.6))') 'delta',dk,0d0,0d0
	write(100, '(a,3(1x,f12.6))') 'delta',0d0,dk,0d0
	write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dk
	write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints-1,nkpoints-1,nkpoints-1
	write(100, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',3,&
								   ' item', (nkpoints-1)*(nkpoints-1)*(nkpoints-1), ' data follows'
    
  !----- Create a uniform k-mesh

	! Front and back face
	ik=0	
	do ikx=1,nkpoints	
		do iky=1,nkpoints
			do ikz=1,nkpoints
				kpoint(1,ikx,iky,ikz) = (ikx-1)*dk + offset(1,choice)
				kpoint(2,ikx,iky,ikz) = (iky-1)*dk + offset(2,choice)
				kpoint(3,ikx,iky,ikz) = (ikz-1)*dk + offset(3,choice)
			enddo
		enddo
	enddo

!----- Perform fourier transform
	allocate(HK(nb,nb),Hamk(nb,nb,ikx,iky,ikz))
    allocate(k_ene(nb),eig_v(nb,2,nkpoints,nkpoints,nkpoints))

	print *, "Finding eigenvectors"

	ikp=0
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				ikp=ikp+1
				print*, ikp, "/", nkp3
				Hk = 0d0
				do ir=1,nr
					phase = dot_product(kpoint(:,ikx,iky,ikz),rvec(:,ir))
					HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
				enddo
				Hamk(:,:,ikx,iky,ikz) = HK(:,:)
				call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
				eig_v(:,1,ikx,iky,ikz) = HK(:,12)
				eig_v(:,2,ikx,iky,ikz) = HK(:,13)
			enddo
		enddo
	enddo

	allocate(H_eff(2,2),dHdK(2,2,3))

	!Constructing effective Hamiltonian and finite difference in each direction
	do ikx=1,nkpoints-1
		do iky=1,nkpoints-1
			do ikz=1,nkpoints-1

				do i=1,2
					do j=1,2
						H_eff(i,j) = dot_product(eig_v(:,i,ikx,iky,ikz),matmul(Hamk(:,:,ikx,iky,ikz),eig_v(:,j,ikx,iky,ikz)))
						print*, 'Element: ',i,j,' = ', H_eff(i,j)

						dHdK(i,j,1) = (dot_product(eig_v(:,i,ikx,iky,ikz),matmul(Hamk(:,:,ikx+1,iky,   ikz  ),eig_v(:,j,ikx,iky,ikz))) - H_eff(i,j))/dk
						dHdK(i,j,2) = (dot_product(eig_v(:,i,ikx,iky,ikz),matmul(Hamk(:,:,ikx,  iky+1, ikz  ),eig_v(:,j,ikx,iky,ikz))) - H_eff(i,j))/dk
						dHdK(i,j,3) = (dot_product(eig_v(:,i,ikx,iky,ikz),matmul(Hamk(:,:,ikx,  iky,   ikz+1),eig_v(:,j,ikx,iky,ikz))) - H_eff(i,j))/dk
					enddo
				enddo

				!Constructing the V_i indices
				do k=1,3
					v(k,1,ikx,iky,ikz) =	real(dHdK(1,2,k))   
					v(k,2,ikx,iky,ikz) =  -aimag(dHdK(1,2,k))
					v(k,3,ikx,iky,ikz) =   		 dHdK(1,1,k)

					sum(k) = v(k,1,ikx,iky,ikz) + v(k,2,ikx,iky,ikz) + v(k,3,ikx,iky,ikz)
				enddo 

				
				write(100, '(3(1x,f12.2))') sum(1),sum(2),sum(3)
			enddo
		enddo
	enddo

	deallocate(Hamk,Hk,eig_v)

	! CALL MPI_GATHER( ENE   ,ECOUNTS,MPI_DOUBLE_PRECISION,   &
	!                  ENERGY,ECOUNTS,MPI_DOUBLE_PRECISION, &
	!                       0,MPI_COMM_WORLD,IERR)
	! write(100, '(2(1x,f12.6))') energy
	write(100,'(A,/,A,/,A,/,A)') &
		'object "regular positions regular connections" class field', &
		'component "positions" value 1', &
		'component "connections" value 2', &
		'component "data" value 3', &
		'end'

	! write(100, '(A)') 'end'
    
end Program Projected_band_structure

