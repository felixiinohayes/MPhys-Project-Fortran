module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kxmax=0.06,kymax=0.06,kzmax=0.06,amax=0.01892,acritical=0.791
    integer,parameter::meshres=20,xmeshres=meshres,ymeshres=meshres,zmeshres=meshres,ares=0,nkxpoints=(2*xmeshres+1),nkypoints=(2*ymeshres+1),nkzpoints=(2*zmeshres+1),napoints=(2*ares+1),nbmin=12,nbmax=13,nkp3=nkxpoints*nkypoints*nkzpoints
    integer nb
    INTEGER IERR,MYID,NUMPROCS
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    INCLUDE 'mpif.h'
!------------------------------------------------------
    real*8 dx,dy,dz,da
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,j1,j2,lwork,info,ikx,iky,ikz,ia,ik,count,kpool,kpmin,kpmax,ecounts,ikp,ir
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two), B = 0.06d0
    real*8 phase,pi2,x1,y1,x2,y2,a,minbandgap
    real*8 avec(3,3),bvec(3,3),kpoint(3,nkp3),rvec_data(3)
    real*8,allocatable:: rvec(:,:),rwork(:)
    real*8, allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:),energy(:,:),ene(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:),B_pt(:,:),B_sigma(:,:)
!------------------------------------------------------
    call init_mpi

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
        open(100,file='fermi_surface.dx')
    endif
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

    lwork=max(1,2*nb-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

    dx=kxmax/xmeshres
    dy=kymax/ymeshres
    dz=kzmax/zmeshres
    !da=amax/ares
    da = 0

!----- Create header of dx files

    if(myid.eq.0) then
        write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkxpoints,nkypoints,nkzpoints
        write(100, '(a,3(1x,f12.6))') 'origin',-kxmax,-kymax,-kzmax+0.5d0*bvec(3,3)
        write(100, '(a,3(1x,f12.6))') 'delta',dx,0d0,0d0
        write(100, '(a,3(1x,f12.6))') 'delta',0d0,dy,0d0
        write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dz
        write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkxpoints,nkypoints,nkzpoints
    endif
    
  !----- Create a uniform k-mesh
      ik=0
    do ikx=-xmeshres,xmeshres
      do iky=-ymeshres,ymeshres
        do ikz=-zmeshres,zmeshres
          ik=ik+1
		  kpoint(1,ik)=ikx*dx! + 0.017665681958398235
          kpoint(2,ik)=iky*dy! + 0.046638430945586576
          kpoint(3,ik)=ikz*dz + 0.5*bvec(3,3)! + 0.47514974714462382
		!   print *, kpoint(:,ikx,iky,ikz)
        enddo
      enddo
    enddo

    kpool=nkp3/numprocs
    if (mod(nkp3,numprocs).ne.0) kpool=kpool+1

    kpmin=1+myid*kpool
    kpmax=(myid+1)*kpool

    ecounts=kpool*2 ! to account for bands 12 and 13

!------- Construct B perturbation
    allocate(B_pt(nb, nb),B_sigma(2,2))
    ! B_sigma(1,:) = [0d0,  B]
    ! B_sigma(2,:) = [B,  0d0]

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
    allocate(sam(3,nbmin:nbmax), oam(3,nbmin:nbmax), ene(2,kpool),k_ene(nb),energy(2,nkp3))

    count=3
    !do ia=-ares,ares
        print *,'interpolation =',ares,'processor =',myid
        !a=ia*da + acritical
        a = acritical

        if(myid.eq.0) then
            write(100, '(a,i8,a,i8,a,i10,a)') 'object',count,' class array type float rank 1 shape',2,&
                                        ' item', nkp3, ' data follows'
        endif
        count=count+1
        ikp=0
        minbandgap=0d0
        do ik=kpmin,min(kpmax,nkp3)
            ikp=ikp+1
            if(myid.eq.0) then
				print*, ikp, "/", nkp3/NUMPROCS
			endif
            Hk = 0d0
            do ir=1,nr
            phase = dot_product(kpoint(:,ik),rvec(:,ir))
            HK=HK+((1-a)*(triv_Hr(:,:,ir))+(a)*(top_Hr(:,:,ir)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ir))
            enddo
            HK = HK+B_pt
            call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
            ene(:,ikp)=k_ene(12:13)
        enddo

        CALL MPI_GATHER( ENE   ,ECOUNTS,MPI_DOUBLE_PRECISION,   &
                        ENERGY,ECOUNTS,MPI_DOUBLE_PRECISION, &
                            0,MPI_COMM_WORLD,IERR)
        if(myid.eq.0)  write(100, '(2(1x,f12.6))') energy
        if(myid.eq.0)  write(100, '(a)') 'attribute "dep" string "positions"'
    !enddo

    if(myid.eq.0) then
        do i=0,napoints-1
            write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
            'object',napoints+3+i,' class field', &
            'component "positions" value 1', &
            'component "connections" value 2', &
            'component "data" value ',3+i
        enddo
        write(100, '(a)') 'object "series" class series'
        do i=0,napoints-1
            write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+napoints+3), ' position', i
        enddo

        write(100, '(A)') 'end'
    endif
    
    deallocate(energy)
    call MPI_FINALIZE( IERR )
end Program Projected_band_structure

SUBROUTINE INIT_MPI
    USE PARAMETERS               ,             ONLY: IERR,MYID,NUMPROCS
    IMPLICIT NONE
    INCLUDE 'mpif.h'
        Call MPI_INIT( IERR )
        Call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
        Call MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS , IERR )
!        Write(*,*) ‘Process’, myid, ' of ’, NUMPROCS , ‘is alive.’
END SUBROUTINE INIT_MPI
