	module parameters
		Implicit None
	!--------to be modified by the user
		real*8,parameter::ef= 4.18903772,alpha=1,kmax=0.1,beta0=1d0,Bx=0.00d0
		INTEGER IERR,MYID,NUMPROCS
		
	end module parameters

     Program Projected_band_structure
	 use parameters
     Implicit None
!-------to be midified by the usere
     character(len=80):: prefix="BiTeI"
     integer,parameter::nkpath=3,np=600,nbmin=12,nbmax=13,meshres=20
!-----------------------------------------------------
     integer*4 ik,ikmax, skip,sign
     real*8 kz 
     real*8 :: Te_sum, Bi_sum, I_sum
     character(len=30)::klabel(nkpath)
     character(len=80) nnkp,line,top_file,triv_file
     integer*4,parameter::nkpoints=(2*meshres+1),nkp2=nkpoints*nkpoints
     integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2,nb,l,ikp,ir
     real*8,parameter::third=1d0/3d0,pi_8=4*atan(1.0_8),cos60=cos(pi_8/3),tan30=tan(pi_8/6),sin30=sin(pi_8/6),sqrt2=sqrt(2d0)
     real*8 phase,pi2,jk,a,b,x1,y1,dx,dy,theta_r,beta
     real*8 xk(nkp2),bvec(3,3),avec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),kpoints(3,nkpath),kpath(3,nkp2),dk(3),kpoint(3),rvec_3(3)
     real*8,allocatable:: rvec(:,:),rvec_data(:,:),rvec_cart(:,:),ene(:,:),rwork(:),od(:,:,:)
     integer*4,allocatable:: ndeg(:)
     complex*16,allocatable:: Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:),H_col(:),B_pt(:,:)
     complex*16 temp1,temp2,B_sigma(2,2)
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

!-----read H(R) 
     open(99,file=trim(adjustl(top_file)))
     open(97,file=trim(adjustl(triv_file)))
     open(100,file='energy.dx')
     read(99,*)
     read(99,*)nb,nr
     allocate(rvec(3,nr),rvec_data(3,nr),rvec_cart(3,nr),Hk(nb,nb),triv_hr(nb,nb,nr),Top_hr(nb,nb,nr),ndeg(nr),ene(nb,nkp2))
     read(99,*)ndeg
     do i = 1, 80
            read(97, *)! Read and discard 80 lines
     enddo
     do k=1,nr
        do i=1,nb
           do j=1,nb
              read(99,*)rvec(1,k),rvec(2,k),rvec(3,k),i1,i2,a,b
              top_hr(i1,i2,k)=dcmplx(a,b)
              read(97,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),j1,j2,x1,y1
              triv_hr(j1,j2,k)=dcmplx(x1,y1)
           enddo
        enddo
     enddo

   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
	! do ir=1,nr
    	! do i=1,6
        	! do j=1,6
	! 			rvec_cart(:,ir) = rvec_data(1,ir)*avec(:,1) + rvec_data(2,ir)*avec(:,2)

	! 			if(modulo(i,3)==modulo(j,3)) then ! Bi-I
	! 				continue
	! 			else if((modulo(i,3)==0 .and. modulo(j,3)==2) .or. (modulo(i,3)==2 .and. modulo(j,3)==0)) then ! 
	! 				rvec_cart(1,ir) = rvec_cart(1,ir) - avec(2,2)*0.5d0*tan30/sin30
	! 				rvec_cart(2,ir) = rvec_cart(2,ir) + avec(2,2)*0.5d0
	! 			else if((modulo(i,3)==0 .and. modulo(j,3)==1) .or. (modulo(i,3)==1 .and. modulo(j,3)==0)) then ! Te-I
	! 				rvec_cart(1,ir) = rvec_cart(1,ir) + avec(2,2)*0.5d0*tan30/sin30
	! 			else if((modulo(i,3)==1 .and. modulo(j,3)==2) .or. (modulo(i,3)==2 .and. modulo(j,3)==1)) then ! Te-I
	! 				rvec_cart(1,ir) = rvec_cart(1,ir) + avec(2,2)*0.5d0*tan30
	! 				rvec_cart(2,ir) = rvec_cart(2,ir) + avec(2,2)*0.5d0
	! 			endif

	! 			if(rvec_cart(2,ir)==0) then ! To stop divide by zero error
	! 				theta_R = 0
	! 			else
	! 				theta_R = atan(rvec_cart(1,ir)/rvec_cart(2,ir))
	! 			endif
	! 			beta=beta0*abs(cos(theta_R))

	! 			do i1=1,3
	! 				do i2=1,3
	! 					top_hr(3*(i-1)+i1, 3*(j-1)+i2, ir) = top_hr(3*(i-1)+i1, 3*(j-1)+i2, ir) * beta
	! 					triv_hr(3*(i-1)+i1, 3*(j-1)+i2, ir) = triv_hr(3*(i-1)+i1, 3*(j-1)+i2, ir) * beta
	! 				enddo
	! 			enddo
	! 			print *, "test"
	! 		enddo
	! 	enddo
	! enddo
	print *, "test"


!----- Create K-mesh
    dx = kmax / meshres
    dy = kmax / meshres
      !----- Create a uniform k-mesh
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
    write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',4,&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(200, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'
    write(300, '(a,a,i10,a)') 'object 3 class array type float rank 1 shape 6',&
                                    ' item', nkpoints*nkpoints,' data follows'

                            
!----- Perform Fourier transform

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
	print *, "test"

    ikp = 0
    do ikx=-meshres,meshres
        do iky=-meshres,meshres
            ikp = ikp + 1
			print *, ikp, "/", nkp2
            kpoint(1)= ikx*dx
            kpoint(2)= iky*dy
            kpoint(3)= 0.5d0*bvec(3,3)
            HK=(0d0,0d0)
            do i=1,nr
                rvec_3 = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)
                phase = dot_product(kpoint,rvec_3)
                ! HK=HK+Hamr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
	            HK=HK+((1-alpha)*(triv_hr(:,:,j))+alpha*(top_hr(:,:,j)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))
            enddo
			HK=HK+B_pt
            call zheev('V','U',nb,HK,nb,ene(:,ikp),work,lwork,rwork,info)
            write(100, '(4(1x,f12.6))') ene(nbmin,ikp), ene(nbmin+1,ikp),ene(nbmin+2,ikp),ene(nbmin+3,ikp)
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

