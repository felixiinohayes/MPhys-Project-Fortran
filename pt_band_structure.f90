module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::kmax=0.07, alpha = 0.77966, ef = 5.5632247424157093
    integer,parameter::meshres=200, nkpoints=(2*meshres+1),nbmin=13,nbmax=13
    integer nb
    
end module parameters
Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dx, dy
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,a,b,x1,y1
    real*8 avec(3,3),bvec(3,3),rvec(3),kpoint(3)
    real*8,allocatable:: rvec_data(:,:),rvec_data_t(:,:),ene(:),rwork(:),k_ene(:),kpoints(:,:), sam(:,:), oam(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:)
    complex*8, parameter:: one = complex(1.d0,0.d0),im = complex(0.d0,1.d0), zero = complex(0.d0,0.d0)
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
    open(100,file='pt_band_structure.dx')
    open(200,file='sam_pt.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),rvec_data_t(3,nr),Hk(nb,nb),Top_hr(nb,nb,nr),Triv_hr(nb,nb,nr),ndeg(nr),ene(nb))
    read(99,*)ndeg
    do i = 1, 80
        read(97, *)! Read and discard 80 lines
    end do
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             top_hr(i1,i2,k)=dcmplx(a,b)
             read(97,*)rvec_data_t(1,k),rvec_data_t(2,k),rvec_data_t(3,k),j1,j2,x1,y1
             triv_hr(j1,j2,k)=dcmplx(x1,y1)
          enddo
       enddo
    enddo
   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
!----- Create K-mesh
    dx = kmax / meshres
    dy = kmax / meshres
!----- Fermi Energy
!----- Create header of dx files
    write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(100, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(100, '(a,2(1x,f12.6))') 'delta',dx,0d0
    write(100, '(a,2(1x,f12.6))') 'delta',0d0,dy
    write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',nbmax-nbmin+1,&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(200, '(a,2(1x,i8))') 'object 1 class gridpoints counts',nkpoints,nkpoints
    write(200, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(200, '(a,2(1x,f12.6))') 'delta',dx,0d0
    write(200, '(a,2(1x,f12.6))') 'delta',0d0,dy
    write(200, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(200, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',nbmax-nbmin+1,&
                                    ' item', 3*nkpoints*nkpoints,' data follows'

!----- Perform fourier transform
    allocate(sam(3,nbmin:nbmax), oam(3,nbmin:nbmax), k_ene(nb))
    do ikx=-meshres,meshres
        do iky=-meshres,meshres
            kpoint(1)= ikx*dx
            kpoint(2)= iky*dy
            kpoint(3)= 0.5d0*bvec(3,3)
            HK=(0d0,0d0)
            do i=1,nr
                rvec = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)
                phase = dot_product(kpoint,rvec)
                HK=HK+((1-alpha)*(triv_hr(:,:,i))+alpha*(top_hr(:,:,i)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
            enddo
            call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
            write(100, '(2(1x,f12.6))') k_ene(nbmin) - ef, k_ene(nbmax) - ef
            call projections(HK,sam,oam)
            write(200, '(6(1x,f12.6))') sam(:,nbmin), sam(:,nbmax)
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
end program Projected_band_structure

subroutine projections(H,sam,oam)
    use parameters
    Implicit None
    complex*16 H(nb,nb),chi(2,1)
    real*8 sam(3),oam(3),sx(1,1),sy(1,1),sz(1,1)
    complex*8 pauli_x(2, 2), pauli_y(2, 2), pauli_z(2, 2), Lx(3,3), Ly(3,3), Lz(3,3), Y_lm(3,1)
    integer ib,jb
!-----Spin projection
   !-Define Pauli matrices
   
    data pauli_x / (0d0,0d0),(1d0, 0d0),(1d0,0d0),( 0d0,0d0)/
    data pauli_y / (0d0,0d0),(0d0,-1d0),(0d0,1d0),( 0d0,0d0)/
    data pauli_z / (1d0,0d0),(0d0, 0d0),(0d0,0d0),(-1d0,0d0)/

    sam=0d0 
    oam=0d0

    do ib=nbmin,nbmax
            do jb=1,nb/2
                chi(1,1) = H(jb     ,ib) 
                chi(2,1) = H(jb+nb/2,ib)

              
                sx = matmul(conjg(transpose(chi)),matmul(pauli_x, chi))
                sy = matmul(conjg(transpose(chi)),matmul(pauli_y, chi))
                sz = matmul(conjg(transpose(chi)),matmul(pauli_z, chi))
                sam(1)=sam(1)+sx(1,1)
                sam(2)=sam(2)+sy(1,1)
                sam(3)=sam(3)+sz(1,1)
           
            enddo
    enddo
end subroutine projections
