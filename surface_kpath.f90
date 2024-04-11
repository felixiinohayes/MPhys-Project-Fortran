module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,a=0.7
    integer,parameter::nkpath=3,np=20,nblocks=20,nr3=11,nk=(nkpath-1)*np+1
	integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    !INCLUDE 'mpif.h'
!------------------------------------------------------
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ik,count,ir,ir3,ir12,nr12,r3,sign
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), B=0.0d0
    real*8 phase,pi2,x1,y1,x2,y2
    real*8 xk(nk),avec(3,3),bvec(3,3),rvec_data(3),kpoints(3,nkpath),kpath(3,nk),dk(3)
    real*8,allocatable:: rvec(:,:),rvec_miller(:,:),rwork(:),k_ene(:,:)
	integer*4,allocatable:: ndeg(:)
    complex*16,allocatable::Hk(:,:),Hkr3(:,:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),super_H(:,:),sH(:,:),a_p_top(:,:),a_p_bottom(:,:),B_pt(:,:)
	complex*16 B_sigma(2,2)
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
    open(100,file='super_H_top.dat')
    open(200,file='super_H_bottom.dat')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(2,nr),rvec_miller(3,nr),Hk(nb,nb),Hkr3(nb,nb,nr3),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr))
	allocate(super_H(nb*nblocks,nb*nblocks),sH(nb,nb*nblocks),k_ene(nb*nblocks,nk))
	allocate(a_p_top(nb*nblocks,nk),a_p_bottom(nb*nblocks,nk))
    read(99,*)ndeg

    do i=1,80
      read(97,*)
    enddo
    do ir=1,nr
		do i=1,nb
			do j=1,nb
			   read(99,*)rvec_data(1),rvec_data(2),rvec_data(3),i1,i2,x1,y1
			   top_Hr(i1,i2,ir)=dcmplx(x1,y1)
			   read(97,*)rvec_data(1),rvec_data(2),rvec_data(3),j1,j2,x2,y2
			   triv_Hr(j1,j2,ir)=dcmplx(x2,y2)
			enddo
		enddo
		rvec(:,ir) = rvec_data(1)*avec(:,1) + rvec_data(2)*avec(:,2)
    enddo

    lwork=max(1,2*(nb*nblocks)-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*(nb*nblocks)-2)))

!-----kpath
    ! ky -> -ky 
	kpoints(:,1) = [ 0.125d0,  -0.25d0,   0.5d0]  !H
	kpoints(:,2) = [ 0.0d0,   0.0d0,    0.5d0]  !A
	kpoints(:,3) = [ -0.125d0,   0.25d0,  0.5d0]  !-H


	! kx -> -kx
	! kpoints(:,1) = [ -0.5d0,   0.0d0,   0.5d0 ]  !L
	! kpoints(:,2) = [ 0.0d0,   0.0d0,   0.5d0 ]  !A
	! kpoints(:,3) = [ 0.5d0,  0.0d0,   0.5d0 ]  !-L

	do j = 1, nkpath-1
		  sign = 1
		  if(j ==1) sign = -1
		do i = 1, np
			ik = i + np*(j-1)
			dk = (kpoints(:,j+1)-kpoints(:,j))/np
			kpath(:, ik) = kpoints(:,(j)) + (dk*(i-1))
			xk(ik) =  sign*sqrt(dot_product(kpoints(:,2)- kpath(:, ik),kpoints(:,2) - kpath(:, ik)))
		enddo
	enddo

!----Construct magnetic perturbation
	allocate(B_pt(nb, nb))

     !B along Y axis
	B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-B)]
    B_sigma(2,:) = [dcmplx(0d0,B) ,  dcmplx(0d0,0d0)]

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
!----- Perform fourier transform
	nr12=nr/nr3
	count = 0
	do ik=1,nk
		count = count + 1
		do ir3=1,nr3 ! Loop over R3 vectors
			Hk=0d0	
			do ir12=0,nr12-1 ! Loop over (R1,R2) vectors
				ir = ir3 + ir12*nr3 ! Calculate index of (R1,R2) vector in nr

				phase = 0d0
				do j = 1,2
					phase = phase + kpath(j,ik)*rvec(j,ir)
				enddo

				Hk=Hk+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
			enddo
			Hk = Hk + B_pt
			Hkr3(:,:,ir3) = Hk
		enddo

		do i=0,nblocks-1
			do j=0,nblocks-1
				r3 = i-j
				if (r3<=5 .AND. r3>=-5) then
					super_H((1+nb*i):(nb*(i+1)),(1+nb*j):(nb*(j+1))) = Hkr3(:,:,r3 + (nr3+1)/2)
				else
					super_H((1+nb*i):(nb*(i+1)), (1+nb*j):(nb*(j+1))) = 0d0
				endif
			enddo
		enddo
		call zheev('V','U',nb*nblocks,super_H,nb*nblocks,k_ene(:,ik),work,lwork,rwork,info)

		do i=1,nb*nblocks
			a_p_top(i,ik) = dot_product(super_H(1:18,i),super_H(1:18,i))
			a_p_bottom(i,ik) = dot_product(super_H(nb*(nblocks-1):nb*nblocks,i),super_H(nb*(nblocks-1):nb*nblocks,i))
		enddo
		
		do j=1, nb*nblocks
			do i=1 ,nk-1
			if(k_ene(j,i).gt.0)write(100, '(f12.6,a,f12.6,a,f12.6)') xk(i),",", k_ene(j,i),",", real(a_p_top(j,i))! Top surface
			if(k_ene(j,i).gt.0)write(200, '(f12.6,a,f12.6,a,f12.6)') xk(i),",", k_ene(j,i),",", real(a_p_bottom(j,i))! Bottom surface
			enddo
			write(100,*)
			write(100,*)
			write(200,*)
			write(200,*)
		enddo
		print *, count
	enddo
	call write_plt()
	stop
end Program Projected_band_structure

subroutine write_plt()
	open(99,file='band.plt')
	write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [ -0.09 : 0.09]'
	write(99,'(a)') &
		 'set terminal pngcairo enhanced font "DejaVu,30" fontscale 1 size 2000, 2000'
	write(99,'(a,f4.2,a)')'set output "band.png"'
	write(99,'(13(a,/),a)') &
		 'set border',&
		 'unset xtics',&
		 'unset ytics',&
		 'set encoding iso_8859_1',&
		 'set size ratio 0 1.0,1.0',&
		 'set yrange [5.5: 6.5]',&
		 'unset key',&
		 'set mytics 2',&
		 'set parametric',&
		 'set palette defined (0 "#deebf7", 1 "#c6dbef", 2 "#9ecae1", 3 "#6baed6", 4 "#4292c6", 5 "#2171b5", 6 "#084594")',&
		 'set trange [-10:10]',&
		 'set multiplot',&
		 'plot "super_H_top.dat" u 1:2:3 with l lt 1 lw 1 linecolor palette notitle',&
		 'unset multiplot'

   end subroutine write_plt
! SUBROUTINE INIT_MPI
!     USE PARAMETERS               ,             ONLY: IERR,MYID,NUMPROCS
!     IMPLICIT NONE
!     INCLUDE 'mpif.h'
!         Call MPI_INIT( IERR )
!         Call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
!         Call MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS , IERR )
! !        Write(*,*) ‘Process’, myid, ' of ’, NUMPROCS , ‘is alive.’
! END SUBROUTINE INIT_MPI
