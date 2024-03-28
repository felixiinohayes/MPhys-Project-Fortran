module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.001,a=0.791
    integer,parameter::meshres=1,nkpoints=(2*meshres+1),nkp3=nkpoints*nkpoints*nkpoints
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir,node,pair
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,x1,y1,x2,y2,chern,div_F,diff_z
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkpoints,nkpoints,nkpoints),rvec_data(3),dV(3),offset(3,2,5),normal(3),v(3,3),v2xv3(3)
	real*8 dAdx(3),dAdy(3),dAdz(3)
	complex*16 dn(18,3)
	real*8,allocatable:: rvec(:,:),rwork(:)
    real*8,allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),Hamk(:,:,:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),eigenvector(:),eig_v(:,:,:,:,:),H_eff(:,:,:),dHdK(:,:,:)
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

	node = 1
	pair = 5

	offset(:,1,1) = (/-0.017659606952654991,0.046513917396043679,0.43965460613976798/) !+ve
	offset(:,2,1) = (/ 0.017665681958398235,0.046638430945586576,0.47514974714462382/) !-ve

	offset(:,2,2) = (/ 0.049353657200408851,0.0069912476907357081,0.43897293203235604/) !-ve
	
	offset(:,1,3) = (/ 0.01766767357940942,-0.04650968222518360,0.47482330199264761/) !+ve
	offset(:,2,3) = (/-0.01765908461242609,-0.04663940372172213,0.43932749486813732/) !-ve

	offset(:,1,4) = (/-0.04936478978846845,0.006832977105296199,0.47559337269211072/) !+ve
	offset(:,2,4) = (/-0.04933437400371396,-0.00708356750361176,0.47545289372516780/) !-ve

	offset(:,1,5) = (/-0.00879804124970561,0.04855142510428899,0.44745395357167184/)  !Dirac Point
    
  !----- Create a uniform k-mesh

	! Front and back face
	ik=0	
	do ikx=1,nkpoints	
		do iky=1,nkpoints
			do ikz=1,nkpoints
				kpoint(1,ikx,iky,ikz) = (ikx-1)*dk + offset(1,node,pair)
				kpoint(2,ikx,iky,ikz) = (iky-1)*dk + offset(2,node,pair)
				kpoint(3,ikx,iky,ikz) = (ikz-1)*dk + offset(3,node,pair)
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

	allocate(H_eff(2,2,4),dHdK(2,2,3))

	!Constructing effective Hamiltonian and finite difference in each direct
	do k = 1,3
		do i=1,2
			do j=1,2
				H_eff(i,j,1) = dot_product(eig_v(:,i,1,1,1),matmul(Hamk(:,:,1,1,1),eig_v(:,j,1,1,1)))
				!write(*, '(a,3(1x,i8),a,2(1x,f12.4))') 'Element: ',i,j,k,' = ', H_eff(i,j,1)

				H_eff(i,j,2) = dot_product(eig_v(:,i,1,1,1),matmul(Hamk(:,:,2,1,1),eig_v(:,j,1,1,1)))
				H_eff(i,j,3) = dot_product(eig_v(:,i,1,1,1),matmul(Hamk(:,:,1,2,1),eig_v(:,j,1,1,1)))
				H_eff(i,j,4) = dot_product(eig_v(:,i,1,1,1),matmul(Hamk(:,:,1,1,2),eig_v(:,j,1,1,1)))
				!write(*, '(a,3(1x,i8),a,2(1x,f12.4))') 'Element: ',i,j,k,' = ', H_eff(i,j,k+1)

				dHdK(i,j,1) = (H_eff(i,j,2) - H_eff(i,j,1))/dk
				dHdK(i,j,2) = (H_eff(i,j,3) - H_eff(i,j,1))/dk
				dHdK(i,j,3) = (H_eff(i,j,4) - H_eff(i,j,1))/dk
				write(*, '(a,2(1x,i8),a,2(1x,f12.4))') 'Element: ',i,j,' = ', dHdK(i,j,k)
			enddo
		enddo
		print*,''
	enddo

	!Constructing the V_i indices
	do k=1,3
		v(k,1) =	real(dHdK(1,2,k))   
		v(k,2) =  -aimag(dHdK(1,2,k))
		v(k,3) =   		 dHdK(1,1,k)
	enddo 

	!Cross product v2 x v3
	v2xv3(1) = v(2,2)*v(3,3) - v(3,2)*v(2,3)
	v2xv3(2) = v(1,2)*v(3,3) - v(3,2)*v(1,3)
	v2xv3(3) = v(1,2)*v(2,3) - v(2,2)*v(1,3)

	chern= dot_product(v(:,1),v2xv3)
	print *, 'Chirality =', chern

	deallocate(Hamk,Hk,eig_v)

	! CALL MPI_GATHER( ENE   ,ECOUNTS,MPI_DOUBLE_PRECISION,   &
	!                  ENERGY,ECOUNTS,MPI_DOUBLE_PRECISION, &
	!                       0,MPI_COMM_WORLD,IERR)
    
end Program Projected_band_structure
