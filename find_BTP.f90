module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,a=0.791
    integer,parameter::meshres=3,nkp=(2*meshres+1),nkp3=nkp*nkp*nkp
    integer nb,nr
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dx,dy,dz,da
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), kmax=0.008d0
    real*8 phase,pi2,x1,y1,x2,y2,bandgap,minbandgap
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkp3),rvec_data(3),minkpoint(3),kmiddle_initial(3)
    real*8,allocatable:: rvec(:,:),rwork(:)
    real*8, allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:)
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
    allocate(rvec(3,nr),Hk(nb,nb),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr))
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

	kmiddle_initial = (/-0.017,-0.05,0.435/)
	lwork=max(1,2*nb-1)
	allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
	allocate(k_ene(nb))

	minbandgap=0.5d0
	call find_minbandgap(kmax,kmiddle_initial,triv_Hr,top_Hr,rvec,ndeg,minbandgap)


contains

	RECURSIVE SUBROUTINE find_minbandgap(kmax,kmiddle, triv_Hr, top_Hr, rvec,ndeg, minbandgap)
		use parameters
		Implicit None
		real*8 kmiddle(3), kpoint(3,nkp3), rvec(3,nr)
		real*8 kmax, dk, bandgap, minbandgap, phase
		integer*4 ndeg(nr),ik,ir,ikx,iky,ikz
		complex*16 top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr)
	  !----- Create a uniform k-mesh
		if (kmax < 0.0000000001) then
			write(*,'(3(1x,f24.17))') kmiddle
			return
		endif

		dk=kmax/meshres

		ik=0
		do ikx=-meshres,meshres
		  do iky=-meshres,meshres
			do ikz=-meshres,meshres
			  ik=ik+1
			  kpoint(1,ik)=ikx*dk + kmiddle(1)
			  kpoint(2,ik)=iky*dk + kmiddle(2)
			  kpoint(3,ik)=ikz*dk + kmiddle(3)
			enddo
		  enddo
		enddo

	!----- Perform fourier transform

		ikp=0
		do ik=1, nkp3
			print *, ik, "/", nkp3, kmiddle, minbandgap
			ikp=ikp+1
			Hk = 0d0
			do ir=1,nr
				phase = dot_product(kpoint(:,ik),rvec(:,ir))
				HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
			enddo
			call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)

			bandgap = abs(k_ene(13)-k_ene(12))	
			if (bandgap < minbandgap) then
				minbandgap = bandgap
				kmiddle = kpoint(:,ik)	
			endif
		enddo
		call find_minbandgap(kmax/1.5,kmiddle,triv_Hr,top_Hr,rvec,ndeg,minbandgap)


	END SUBROUTINE find_minbandgap

end Program Projected_band_structure
!1.7665681958398235E-002   4.6638430945586576E-002  0.47514974714462382  