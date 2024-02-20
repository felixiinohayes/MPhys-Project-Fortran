module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.01,a=0.77966
    integer,parameter::meshres=5,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,chern,div_F
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkpoints,nkpoints,nkpoints),rvec_data(3),dV(3),offset(3),normal(3)
	real*8 dAdx(3),dAdy(3),dAdz(3)
	complex*16 dn(18,3)
	real*8,allocatable:: rvec(:,:),rwork(:)
    real*8,allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eigenvector(:),eigenvec_data(:,:,:,:)
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
    open(100,file='berrycurvature.dx')
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

!----- Create header of dx files

	write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkpoints
	write(100, '(a,3(1x,f12.6))') 'origin',-kmax,-kmax,-kmax+0.5d0*bvec(3,3)
	write(100, '(a,3(1x,f12.6))') 'delta',dk,0d0,0d0
	write(100, '(a,3(1x,f12.6))') 'delta',0d0,dk,0d0
	write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dk
	write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkpoints
	write(100, '(a,i8,a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',3,&
								   ' item', nkp3, ' data follows'
    
  !----- Create a uniform k-mesh

	! Front and back face
	offset = (/0.013d0-(kmax/2),-0.046d0-(kmax/2),0.5d0*bvec(3,3)-(kmax/2)/)
	ik=0	
	do ikx=1,nkpoints	
		do iky=1,nkpoints
			do ikz=1,nkpoints
				kpoint(1,ikx,iky,ikz) = (ikx-1)*dk + offset(1)
				kpoint(2,ikx,iky,ikz) = (iky-1)*dk + offset(2)
				kpoint(3,ikx,iky,ikz) = (ikz-1)*dk + offset(3)
			enddo
		enddo
	enddo

!----- Perform fourier transform
	allocate(HK(nb,nb))
    allocate(k_ene(nb),eigenvec_data(nb,nkpoints,nkpoints,nkpoints))

	print *, "Finding eigenvectors"

	ikp=0
	do ikx=1,nkpoints
		do iky=1,nkpoints
			do ikz=1,nkpoints
				print *, ikp, "/", nkp3
				ikp=ikp+1
				Hk = 0d0
				do ir=1,nr
					!phase = kpoint(1,ik)*rvec(1,ir)+kpoint(2,ik)*rvec(2,ir)+kpoint(3,ik)*rvec(3,ir)
					phase = dot_product(kpoint(:,ikx,iky,ikz),rvec(:,ir))
					HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
					!HK=HK+(top_Hr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i)))
				enddo

				call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
				eigenvec_data(:,ikx,iky,ikz) = HK(:,12)
			enddo
		enddo
	enddo

	print *, "Computing Berry connection"
	allocate(connection(3,nkpoints,nkpoints,nkpoints),curvature(3,nkpoints,nkpoints,nkpoints))

!-------Compute Berry connection
	do ikx=1,nkpoints-1
		do iky=1,nkpoints-1
			do ikz=1,nkpoints-1
				eigenvector = eigenvec_data(:,ikx,iky,ikz)

				dn(:,1) = eigenvec_data(:,ikx+1,iky,ikz) - eigenvector
				dn(:,2) = eigenvec_data(:,ikx,iky+1,ikz) - eigenvector
				dn(:,3) = eigenvec_data(:,ikx,iky,ikz+1) - eigenvector

				connection(1,ikx,iky,ikz) = dot_product(eigenvector,dn(:,1))
				connection(2,ikx,iky,ikz) = dot_product(eigenvector,dn(:,2))
				connection(3,ikx,iky,ikz) = dot_product(eigenvector,dn(:,3))
				connection(:,ikx,iky,ikz) = connection(:,ikx,iky,ikz) / dk
				! print *, ikx,iky,ikz
			enddo	
		enddo
	enddo
	deallocate(eigenvec_data)

	print *, "Computing Berry curvature"

	do ikx=1,nkpoints-2
		do iky=1,nkpoints-2
			do ikz=1,nkpoints-2
				dAdx(:) = connection(:,ikx+1,iky,  ikz  ) - connection(:,ikx,iky,ikz)
				dAdy(:) = connection(:,ikx  ,iky+1,ikz  ) - connection(:,ikx,iky,ikz)
				dAdz(:) = connection(:,ikx  ,iky  ,ikz+1) - connection(:,ikx,iky,ikz)

				curvature(1,ikx,iky,ikz) = (dAdy(3) - dAdz(2)) ! Compute curl
				curvature(2,ikx,iky,ikz) = (dAdz(1) - dAdx(3))
				curvature(3,ikx,iky,ikz) = (dAdx(2) - dAdy(1))
				curvature(:,ikx,iky,ikz) = curvature(:,ikx,iky,ikz) / dk
			enddo	
		enddo
	enddo

	chern=0d0
	do ikx=1,nkpoints-3
		do iky=1,nkpoints-3
			do ikz=1,nkpoints-3
				div_F = (curvature(1,ikx+1,iky  ,ikz  )-curvature(1,ikx,iky,ikz) + &
						curvature(2,ikx  ,iky+1,ikz  )-curvature(2,ikx,iky,ikz) + &
						curvature(3,ikx  ,iky  ,ikz+1)-curvature(3,ikx,iky,ikz)) / dk

				chern = chern + (div_F)/dk
			enddo	
		enddo
	enddo

	print *, chern

	deallocate(connection,curvature)

!!-------Compute Berry curvature
!	dV = (/dk,dk,dk/)
!	chern = 0
!	do ik=1,nkp3
!		eigenvector = eigenvec_data(:,ik)	
!		if (modulo(ik,nkpoints*nkpoints*nkpoints) > (nkpoints*nkpoints*nkpoints - nkpoints*nkpoints)) then
!		   	dFdk = 0d0
!		else
!		  	dFdk(:,1) = (connection(:,ik+(nkpoints*nkpoints)) - connection(:,ik)) / dk ! dF/dx
!		endif
!		if (modulo(ik,nkpoints*nkpoints) > (nkpoints*nkpoints - nkpoints)) then 
!			dFdk = 0d0 
!		else
!		   	dFdk(:,2) = (connection(:,ik+nkpoints) - connection(:,ik)) / dk ! dF/dy
!		endif
!		if (modulo(ik,nkpoints) == 0) then
!			dFdk = 0d0
!	   	else
!		   	dFdk(:,3) = (connection(:,ik+1) - connection(:,ik)) / dk ! dF/dz
!		endif
!		curvature(1,ik) = dFdk(3,2) - dFdk(2,3) ! Compute curl
!		curvature(2,ik) = dFdk(1,3) - dFdk(3,1) 
!		curvature(3,ik) = dFdk(2,1) - dFdk(1,2) 
!		write(100, '(3(1x,f12.6))') curvature(:,ik)

!		chern = chern + dot_product(curvature(:,ik),dV) ! VOLUME INTEGRAL - wrong pretty sure
!	enddo

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

