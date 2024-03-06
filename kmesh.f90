module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.2,two=2.0d0,sqrt2=sqrt(two),Bx=0.06d0
    integer,parameter::meshres=100, nkpoints=(2*meshres+1),nbmin=11,nbmax=14
    integer nb
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dx, dy
    character(len=80) hamil_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,ikp
    real*8,parameter::third=1d0/3d0
    real*8 phase,pi2,a,b
    real*8 avec(3,3),bvec(3,3),rvec(3),kpoint(3),B_sigma(2,2)
    real*8,allocatable:: rvec_data(:,:),ene(:),rwork(:),k_ene(:),kpoints(:,:), sam(:,:), oam(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:),B_pt(:,:)
    complex*8, parameter:: one = complex(1.d0,0.d0),im = complex(0.d0,1.d0), zero = complex(0.d0,0.d0)
!------------------------------------------------------
    write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
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
    open(99,file=trim(adjustl(hamil_file)))
    open(100,file='energy.dx')
    open(200,file='sam.dx')
    open(300,file='oam.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),Hk(nb,nb),Hamr(nb,nb,nr),ndeg(nr),ene(nb))
    read(99,*)ndeg
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             hamr(i1,i2,k)=dcmplx(a,b)
          enddo
       enddo
    enddo

   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

!----- Create K-mesh
    dx = kmax / meshres
    dy = kmax / meshres
    
!----- Create header of dx files
    write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(200, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(300, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(100, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(200, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(300, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(100, '(a,2(1x,f12.6))') 'delta',dx,0d0
    write(100, '(a,2(1x,f12.6))') 'delta',0d0,dy
    write(200, '(a,2(1x,f12.6))') 'delta',dx,0d0
    write(200, '(a,2(1x,f12.6))') 'delta',0d0,dy
    write(300, '(a,2(1x,f12.6))') 'delta',dx,0d0
    write(300, '(a,2(1x,f12.6))') 'delta',0d0,dy
    write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(200, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(300, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',nbmax-nbmin+1,&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(200, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(300, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'

                            
!----- Perform Fourier transform
    allocate(sam(3,nbmin:nbmax), oam(3,nbmin:nbmax), k_ene(nb))

	allocate(B_pt(nb,nb))
	data B_sigma /0d0,Bx,Bx,0d0/
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

    ikp = 1
    do ikx=-meshres,meshres
        do iky=-meshres,meshres
            print*, ikp,'/',nkpoints**2 
            ikp = ikp + 1
            kpoint(1)= ikx*dx
            kpoint(2)= iky*dy
            kpoint(3)= 0.5d0*bvec(3,3)

            HK=(0d0,0d0)

            do i=1,nr
                rvec = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)

                phase = dot_product(kpoint,rvec)

                HK=HK+Hamr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
            enddo
			HK=HK+B_pt
            call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
            call projections(HK,sam,oam)
            write(100, '(4(1x,f12.6))') k_ene(nbmin), k_ene(nbmin+1),k_ene(nbmin+2),k_ene(nbmin+3)


            write(200, '(6(1x,f12.6))') sam(:,nbmin), sam(:,nbmax)
            write(300, '(6(1x,f12.6))') oam(:,nbmin), oam(:,nbmax)
        enddo
    enddo
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
    write(300,'(A,/,A,/,A,/,A)') &
    'object "regular positions regular connections" class field', &
    'component "positions" value 1', &
    'component "connections" value 2', &
    'component "data" value 3', &
    'end'
end program Projected_band_structure

subroutine projections(H,sam,oam)
    use parameters
    Implicit None
    complex*16 H(nb,nb),chi(2,1),phi(3)
    real*8 sam(3,nbmin:nbmax),oam(3,nbmin:nbmax),sx(1,1),sy(1,1),sz(1,1),lx(1,1),ly(1,1),lz(1,1)
    complex*8 pauli_x(2, 2), pauli_y(2, 2), pauli_z(2, 2), Lhat_x(3,3), Lhat_y(3,3), Lhat_z(3,3), Y_lm(3,1)
    integer ib,jb,orbital_index
!-----Spin projection
   !-Define Pauli matrices
   
    data pauli_x / (0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
    data pauli_y / (0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
    data pauli_z / (1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/
    data Lhat_x /  (0d0,0d0),(1d0,0d0),(0d0, 0d0),( 1d0, 0d0),(0d0,0d0),(1d0,0d0),(0d0,0d0),(1d0, 0d0),( 0d0,0d0)/     
    data Lhat_y /  (0d0,0d0),(0d0,1d0),(0d0, 0d0),( 0d0,-1d0),(0d0,0d0),(0d0,1d0),(0d0,0d0),(0d0,-1d0),( 0d0,0d0)/     
    data Lhat_z /  (1d0,0d0),(0d0,0d0),(0d0, 0d0),( 0d0, 0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0,0d0)/     

	Lhat_x = Lhat_x * 1/sqrt2
	Lhat_y = Lhat_y * 1/sqrt2

    sam=0d0 
    oam=0d0

    do ib=nbmin,nbmax
            do jb=1,nb/2
                chi(1,1) = H(jb     ,ib) 
                chi(2,1) = H(jb+nb/2,ib)

                sx = matmul(conjg(transpose(chi)),matmul(pauli_x, chi))
                sy = matmul(conjg(transpose(chi)),matmul(pauli_y, chi))
                sz = matmul(conjg(transpose(chi)),matmul(pauli_z, chi))
                sam(1,ib)=sam(1,ib)+sx(1,1)
                sam(2,ib)=sam(2,ib)+sy(1,1)
                sam(3,ib)=sam(3,ib)+sz(1,1)
           
            enddo
            do jb=1,nb,3
                phi(1) = H(jb  ,ib) 
                phi(2) = H(jb+1,ib)
                phi(3) = H(jb+2,ib)                               
                Y_lm(3,1) = (-1/sqrt2) * (phi(1) + dcmplx(0,1)*phi(2))
                Y_lm(1,1) = ( 1/sqrt2) * (phi(1) - dcmplx(0,1)*phi(2))
                Y_lm(2,1) = phi(3)
                
                lx = matmul(conjg(transpose(Y_lm)),matmul(Lhat_x, Y_lm))
                ly = matmul(conjg(transpose(Y_lm)),matmul(Lhat_y, Y_lm))
                lz = matmul(conjg(transpose(Y_lm)),matmul(Lhat_z, Y_lm))
                oam(1,ib)=oam(1,ib)+lx(1,1)
                oam(2,ib)=oam(2,ib)+ly(1,1)
                oam(3,ib)=oam(3,ib)+lz(1,1)
                                                                        
            enddo
    enddo
end subroutine projections
