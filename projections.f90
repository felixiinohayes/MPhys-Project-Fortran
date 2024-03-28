module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.2,two=2.0d0,sqrt2=sqrt(two),Bx=0.0d0,beta0=1d0,alpha=0
    integer,parameter::meshres=20, nkpoints=(2*meshres+1),nbmin=11,nbmax=14,nkp2=nkpoints*nkpoints
    integer nb
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikp,ir,ik
    real*8,parameter::third=1d0/3d0,pi_8=4*atan(1.0_8),cos60=cos(pi_8/3),tan30=tan(pi_8/6),sin30=sin(pi_8/6)
    real*8 phase,pi2,a,b,x1,y1,theta_r,beta,ratio
    real*8 avec(3,3),bvec(3,3),r_frac(3,3)
    real*8,allocatable:: rvec_data(:,:),ene(:),rwork(:),k_ene(:), sam(:,:), oam(:,:),rs(:,:),kpoint(:,:),rvec(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:),B_pt(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),B_sigma(:,:)
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
    open(100,file='energy.dx')
    open(200,file='sam.dx')
    open(300,file='oam.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),Hk(nb,nb),Top_hr(nb,nb,nr),Triv_hr(nb,nb,nr),ndeg(nr),ene(nb),rvec(3,nr),rs(3,nr),Hamr(nb,nb,nr))
    read(99,*)ndeg
    do i = 1, 80
        read(97, *)! Read and discard 80 lines
    end do
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             top_hr(i1,i2,k)=dcmplx(a,b)
             read(97,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),j1,j2,x1,y1
             triv_hr(j1,j2,k)=dcmplx(x1,y1)
          enddo
       enddo
       rvec(:,k) = rvec_data(1,k)*avec(:,1) + rvec_data(2,k)*avec(:,2) + rvec_data(3,k)*avec(:,3)
    enddo
   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

   dk=kmax/meshres
!----- Create header of dx files
    write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(200, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(300, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(100, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(200, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(300, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(100, '(a,2(1x,f12.6))') 'delta',dk,0d0
    write(100, '(a,2(1x,f12.6))') 'delta',0d0,dk
    write(200, '(a,2(1x,f12.6))') 'delta',dk,0d0
    write(200, '(a,2(1x,f12.6))') 'delta',0d0,dk
    write(300, '(a,2(1x,f12.6))') 'delta',dk,0d0
    write(300, '(a,2(1x,f12.6))') 'delta',0d0,dk
    write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(200, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(300, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',4,&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(200, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(300, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'

!----- Create a uniform k-mesh
    ik=0
    do ikx=-meshres,meshres
      do iky=-meshres,meshres
          ik=ik+1
          kpoint(1,ik)=ikx*dk
          kpoint(2,ik)=iky*dk 
          kpoint(3,ik)= 0.5d0*bvec(3,3)
      enddo
    enddo

!----- Perform Fourier transform
    ikp=0                
    do ikx=-meshres,meshres
        do iky=-meshres,meshres
			ikp = ikp+1		
			print *, ikp, "/", nkp2
            HK=(0d0,0d0)
            do i=1,nr
                phase = dot_product(kpoint(:,ikp),rvec(:,i))
                HK=HK+((1-alpha)*(triv_hr(:,:,i))+alpha*(top_hr(:,:,i)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
            enddo
			HK=HK+B_pt
            call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
            call projections(HK,sam,oam)
            write(100, '(4(1x,f12.6))') k_ene(nbmin), k_ene(nbmin+1),k_ene(nbmin+2),k_ene(nbmin+3)

            write(200, '(6(1x,f12.6))') sam(:,nbmin), sam(:,nbmax)
        enddo
    enddo
    print*,'beta(0)=',2-(1-beta0*abs(cos(0d0)))
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

