module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.4,two=2.0d0,sqrt2=sqrt(two),Bx=0.0d0,beta0=0.2d0,alpha=0
    integer,parameter::meshres=20, nkpoints=(2*meshres+1),nbmin=11,nbmax=14,nkp2=nkpoints*nkpoints
    integer nb
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dk
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikp,ir,ik,ii,jj,l,m
    real*8,parameter::third=1d0/3d0,pi_8=4*atan(1.0_8),cos60=cos(pi_8/3),tan30=tan(pi_8/6),sin30=sin(pi_8/6)
    real*8 phase,pi2,a,b,x1,y1,theta_r,beta,ratio
    real*8 avec(3,3),bvec(3,3),r_frac(3,3),rvec_stress(2)
    real*8,allocatable:: rvec_data(:,:),ene(:),rwork(:),k_ene(:),rs(:),kpoint(:,:),rvec(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:),B_pt(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),B_sigma(:,:)
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

    r_frac(1,:) = (/third  ,third*2,0d0/) !I
	r_frac(2,:) = (/third*2,third  ,0d0/) !Te
	r_frac(3,:) = (/0d0    ,0d0    ,0d0/) !Bi

!------read H(R)
    open(99,file=trim(adjustl(top_file)))
    open(97,file=trim(adjustl(triv_file)))
    open(100,file='energy.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),Hk(nb,nb),Top_hr(nb,nb,nr),Triv_hr(nb,nb,nr),ndeg(nr),ene(nb),rvec(3,nr),rs(3),Hamr(nb,nb,nr))
    read(99,*)ndeg
    do i = 1, 80
        read(97, *)! Read and discard 80 lines
    end do
    do k=1,nr
        do i=1,nb
            do j=1,nb
                !-----Topological
                read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
                
                ii = mod((i1-1)/3,3)+1
                jj = mod((i2-1)/3,3)+1

                rvec_stress(:) = rvec_data(:,k) + r_frac(jj,:) - r_frac(ii,:)
                rs(:) = rvec_stress(1)*avec(:,1) + rvec_stress(2)*avec(:,2) 

                top_hr(i1,i2,k)=dcmplx(a,b)*(1d0 + beta0*abs(rs(2))/max(1d-8,sqrt(dot_product(rs,rs))))

                !------Trivial
                read(97,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),j1,j2,x1,y1

                ii = mod((j1-1)/3,3)+1
                jj = mod((j2-1)/3,3)+1

                rvec_stress(:) = rvec_data(:,k) + r_frac(jj,:) - r_frac(ii,:)
                rs(:) = rvec_stress(1)*avec(:,1) + rvec_stress(2)*avec(:,2) 

                triv_hr(j1,j2,k)=dcmplx(x1,y1)*(1d0 + beta0* abs(rs(2))/max(1d-8,sqrt(dot_product(rs,rs))))

            enddo
        enddo
       
        ! avec(2,:) = avec(2,:) * 1.000
        rvec(:,k) = rvec_data(1,k)*avec(:,1) + rvec_data(2,k)*avec(:,2) + rvec_data(3,k)*avec(:,3)
    enddo
   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

   dk=kmax/meshres
!----- Create header of dx files
    write(100, '(a,2(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints
    write(100, '(a,2(1x,f12.6))') 'origin',-kmax,-kmax
    write(100, '(a,2(1x,f12.6))') 'delta',dk,0d0
    write(100, '(a,2(1x,f12.6))') 'delta',0d0,dk
    write(100, '(a,2(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints
    write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',4,&
                                    ' item', nkpoints*nkpoints,' data follows'

!----Magnetic Perturbation
    allocate(k_ene(nb),kpoint(3,nkp2),B_pt(nb,nb),B_sigma(2,2))

    !B along X-axis
    B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(Bx,0d0)]
    B_sigma(2,:) = [dcmplx(Bx,0d0) ,  dcmplx(0d0,0d0)]

    ! !B along Y axis
	! B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-Bx)]
    ! B_sigma(2,:) = [dcmplx(0d0,Bx) ,  dcmplx(0d0,0d0)]
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
            write(100, '(4(1x,f12.6))') k_ene(nbmin), k_ene(nbmin+1),k_ene(nbmin+2),k_ene(nbmin+3)
        enddo
    enddo
    write(*,'(a,(1x,f10.3))') 'beta(0)=',2-(1-beta0*abs(cos(0d0)))
    write(100,'(A,/,A,/,A,/,A)') &
    'object "regular positions regular connections" class field', &
    'component "positions" value 1', &
    'component "connections" value 2', &
    'component "data" value 3', &
    'end' 
end program Projected_band_structure

