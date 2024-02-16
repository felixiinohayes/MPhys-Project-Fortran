module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.06,a=0.7747
    integer,parameter::meshres=10,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints
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
    real*8 phase,pi2,x1,y1,x2,y2,chern
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkp3),rvec_data(3),dV(3),dAx(3),dAy(3),dAz(3),dFdx(3),dFdy(3),dFdz(3)
    real*8,allocatable:: rvec(:,:),rwork(:)
    real*8, allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eigenvector(:),eigenvec_data(:,:)
	complex*16,allocatable::connection(:,:),curvature(:,:)
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
    if(myid.eq.0) then
        open(100,file='btp_motion.dx')
    endif
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(3,nr),Hk(nb,nb),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr),eigenvector(nb))
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

    dk=kmax/meshres

!----- Create header of dx files

	! write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkzpoints
	! write(100, '(a,3(1x,f12.6))') 'origin',-kmax,-kmax,-kzmax+0.5d0*bvec(3,3)
	! write(100, '(a,3(1x,f12.6))') 'delta',dxy,0d0,0d0
	! write(100, '(a,3(1x,f12.6))') 'delta',0d0,dxy,0d0
	! write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dz
	! write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkzpoints
    
  !----- Create a uniform k-mesh
	ik=0
    do ikx=-meshres,meshres
      do iky=-meshres,meshres
        do ikz=-meshres,meshres
          ik=ik+1
          kpoint(1,ik)=ikx*dk
          kpoint(2,ik)=iky*dk
          kpoint(3,ik)=ikz*dk + 0.5d0*bvec(3,3)
        enddo
      enddo
    enddo

    ! kpool=nkp3/numprocs
    ! if (mod(nkp3,numprocs).ne.0) kpool=kpool+1

    ! kpmin=1+myid*kpool
    ! kpmax=(myid+1)*kpool

    ! ecounts=kpool*2 ! to account for bands 12 and 13


!----- Perform fourier transform
    allocate(k_ene(nb),eigenvec_data(nb,nkp3),connection(3,nkp3),curvature(3,nkp3))

	! write(100, '(a,i8,a,i8,a,i10,a)') 'object',count,' class array type float rank 1 shape',2,&
								   ! ' item', nkp3, ' data follows'
	print *, "Finding eigenvectors"

	ikp=0
	do ik=1,nkp3
		ikp=ikp+1
		Hk = 0d0
		do ir=1,nr
			!phase = kpoint(1,ik)*rvec(1,ir)+kpoint(2,ik)*rvec(2,ir)+kpoint(3,ik)*rvec(3,ir)
			phase = dot_product(kpoint(:,ik),rvec(:,ir))
			HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
			!HK=HK+(top_Hr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i)))
		enddo
		call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
		eigenvec_data(:,ik) = HK(:,12)
	enddo

	print *, "Computing Berry connection"

!-------Compute Berry connection
	do ik=1,nkp3
		eigenvector = eigenvec_data(:,ik)	
		if (modulo(ik,nkpoints*nkpoints*nkpoints) > (nkpoints*nkpoints*nkpoints - nkpoints*nkpoints)) then ! Check if last x element	
		   	dAx = 0d0
		else
		  	dAx = eigenvec_data(:,ik+(nkpoints*nkpoints)) - eigenvector
		endif
		if (modulo(ik,nkpoints*nkpoints) > (nkpoints*nkpoints - nkpoints)) then ! Check if last y element in kmesh
			dAy = 0d0 
		else
		   	dAy = eigenvec_data(:,ik+nkpoints) - eigenvector
		endif
		if (modulo(ik,nkpoints) == 0) then ! Check if last z element in kmesh
			dAz = 0d0
	   	else
		   	dAz = eigenvec_data(:,ik+1) - eigenvector
		endif
		connection(1,ik) = dot_product(eigenvector,dAx)
		connection(2,ik) = dot_product(eigenvector,dAy)
		connection(3,ik) = dot_product(eigenvector,dAz)
	enddo

	print *, "Computing Berry curvature"

!-------Compute Berry curvature
	dV = (/dk,dk,dk/)
	chern = 0
	do ik=1,nkp3
		eigenvector = eigenvec_data(:,ik)	
		if (modulo(ik,nkpoints*nkpoints*nkpoints) > (nkpoints*nkpoints*nkpoints - nkpoints*nkpoints)) then
		   	dFdx = 0d0
		else
		  	dFdx = (connection(:,ik+(nkpoints*nkpoints)) - connection(:,ik)) / dk
		endif
		if (modulo(ik,nkpoints*nkpoints) > (nkpoints*nkpoints - nkpoints)) then 
			dFdy = 0d0 
		else
		   	dFdy = (connection(:,ik+nkpoints) - connection(:,ik)) / dk
		endif
		if (modulo(ik,nkpoints) == 0) then
			dFdz = 0d0
	   	else
		   	dFdz = (connection(:,ik+1) - connection(:,ik)) / dk
		endif
		curvature(1,ik) = dFdy(3) - dFdz(2) ! Compute curl
		curvature(2,ik) = dFdz(1) - dFdx(3) 
		curvature(3,ik) = dFdx(2) - dFdy(1) 

		chern = chern + dot_product(curvature(:,ik),dV) ! VOLUME INTEGRAL - wrong pretty sure
	enddo

	print *, chern

	! CALL MPI_GATHER( ENE   ,ECOUNTS,MPI_DOUBLE_PRECISION,   &
	!                  ENERGY,ECOUNTS,MPI_DOUBLE_PRECISION, &
	!                       0,MPI_COMM_WORLD,IERR)
	write(100, '(2(1x,f12.6))') energy
	write(100, '(a)') 'attribute "dep" string "positions"'

	! write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
	! 'object',napoints+3+i,' class field', &
	! 'component "positions" value 1', &
	! 'component "connections" value 2', &
	! 'component "data" value ',3+i
	! write(100, '(a)') 'object "series" class series'
	! write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+napoints+3), ' position', i

	! write(100, '(A)') 'end'
    
end Program Projected_band_structure

