     Program Wannier_band_structure
     Implicit None
!-------to be midified by the usere
     character(len=80):: prefix="BiTeI"
     integer,parameter::nkpath=3,np=600
!-----------------------------------------------------
     integer*4 ik,ikmax, skip,sign
     real*8 kz,ef 
     real*8 :: Te_sum, Bi_sum, I_sum
     character(len=30)::klabel(nkpath)
     character(len=80) nnkp,line,top_file,triv_file
     integer*4,parameter::nk=(nkpath-1)*np+1
     integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2,nb,l
     real*8,parameter::third=1d0/3d0, alpha = 1, Bx = 0.06d0
     real*8 phase,pi2,jk,a,b,x1,y1
     real*8 xk(nk),bvec(3,3),avec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),kpoints(3,nkpath),kpath(3,nk),dk(3)
     real*8,allocatable:: rvec(:,:),rvec_data(:,:),ene(:,:),rwork(:),od(:,:,:)
     integer*4,allocatable:: ndeg(:)
     complex*16,allocatable:: Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:),H_col(:),B_pt(:,:)
     complex*16 temp1,temp2,B_sigma(2,2)
!-----------------------------------------------------
     write(top_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
     write(triv_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
     write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
     pi2=4.0d0*atan(1.0d0)*2.0d0
!--------------  reciprocal vectors
     open(98,file=trim(adjustl(nnkp)))
     111 read(98,'(a)')line
         if(trim(adjustl(line)).ne."begin real_lattice") goto 111
         read(98,*)avec
         read(98,'(a)')line
         read(98,'(a)')line
         read(98,'(a)')line
         read(98,*)bvec
!--------------kpath

     ! ky -> -ky 
     kpoints(:,1) = [ 0.25d0,  -0.5d0,   0.5d0 ]  !H
     kpoints(:,2) = [ 0.0d0,   0.0d0,   0.5d0 ]  !A
     kpoints(:,3) = [ -0.25d0,   0.5d0,   0.5d0 ]  !-H


     ! kx -> -kx
     ! kpoints(:,1) = [ -0.5d0,   0.0d0,   0.5d0 ]  !L
     ! kpoints(:,2) = [ 0.0d0,   0.0d0,   0.5d0 ]  !A
     ! kpoints(:,3) = [ 0.5d0,  0.0d0,   0.5d0 ]  !-L

     data klabel     /'L','A','H'/
  
     do j = 1, nkpath-1
          sign = 1
          if(j ==1) sign = -1
          do i = 1, np
               ik = i + np*(j-1)
               dk = (kpoints(:,j+1)-kpoints(:,j))/np
               kpath(:, ik) = kpoints(:,(j)) + (dk*(i-1))

               kpath(:,ik) = kpath(1,ik)*bvec(:,1) + kpath(2,ik)*bvec(:,2) + kpath(3,ik)*bvec(:,3) 
               xk(ik) =  sign*sqrt(dot_product(kpoints(:,2)*bvec(:,3) - kpath(:, ik),kpoints(:,2)*bvec(:,3) - kpath(:, ik)))

               if(ik==2*np) then
                    kpath(:,nk) = kpoints(1,nkpath)*bvec(:,1) + kpoints(2,nkpath)*bvec(:,2) + kpoints(3,nkpath)*bvec(:,3) 
                    xk(nk) = xk(nk-1) + sqrt(dot_product(kpoints(:,2)*bvec(:,3) - kpath(:, nk),kpoints(:,2)*bvec(:,3) - kpath(:, nk)))
               endif
         enddo
     enddo

!-----read H(R) 
     open(99,file=trim(adjustl(top_file)))
     open(97,file=trim(adjustl(triv_file)))
     open(100,file='band.dat')
     read(99,*)
     read(99,*)nb,nr
     allocate(rvec(3,nr),rvec_data(3,nr),Hk(nb,nb),triv_hr(nb,nb,nr),Top_hr(nb,nb,nr),ndeg(nr),ene(nb,nk),od(nb,nk,3),H_col(nb))
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

!------- Construct B perturbation
	allocate(B_pt(nb, nb))

     !B along Y axis
	B_sigma(1,:) = [dcmplx(0d0,0d0),  dcmplx(0d0,-Bx)]
     B_sigma(2,:) = [dcmplx(0d0,Bx) ,  dcmplx(0d0,0d0)]

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

     ! do i = 1,nb
     !      write(*, '(18(1x,f12.4))') real(B_pt(i,:))
     ! enddo
!---- Fourier transform H(R) to H(k)
    ene=0d0
    do k=1,nk
       HK=(0d0,0d0)
       do j=1,nr
          phase=0.0d0

          rvec(:,j) = rvec_data(1,j)*avec(:,1) + rvec_data(2,j)*avec(:,2) + rvec_data(3,j)*avec(:,3)

          do i=1,3
             phase=phase+kpath(i,k)*rvec(i,j)  
          enddo
          HK=HK+((1-alpha)*(triv_hr(:,:,j))+alpha*(top_hr(:,:,j)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))
       enddo
	   HK=HK+B_pt
       call zheev('V','U',nb,HK,nb,ene(:,k),work,lwork,rwork,info) 
    enddo

!---Fermi level:
    ef = (MAXVAL(ene(12, :)) + MINVAL(ene(13, :)))/2.0d0
    print * ,"E_F: ", ef

!---Orbital probability:


    deallocate(HK,work)
    
    do i=1,nb
       do k=1,nk-1
         write(100,'(5(x,f12.6))') xk(k), ene(i,k) 
       enddo
         write(100,*)
         write(100,*)
    enddo
    call write_plt(nkpath,xkl,klabel,ef,Bx)
    stop
333 write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
    stop
444 write(*,'(3a)')'ERROR: input file "',trim(adjustl(triv_file)),' not found'
    stop

    end

    subroutine write_plt(nkp,xkl,kl,ef,Bx)
          implicit none
          integer nkp,i
          real*8 xkl(nkp),ef,Bx
          character(len=30)kl(nkp)
          
          open(99,file='band.plt')
          write(99,'(a,f12.8)')'ef=',ef
          write(99,'(a,f12.6,a,f12.6,a)') '#set xrange [ -0.12 : 0.12]'
          write(99,'(a)') &
               'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 0.5 size 5.00in, 5.00in'
          write(99,'(a,f4.2,a)')'set output "band.pdf"'
          write(99, '(a,f12.3,a)') 'set title "B =',Bx,'T (along y-axis)"'
          write(99,'(17(a,/),a)') &
               'set border',&
               'set xtics',&
               'set ytics',&
               'set encoding iso_8859_1',&
               'set xlabel "k_y"',&
               'set ylabel "E"',&
               'set size ratio 0 1.0,1.0',&
               'set yrange [-0.4: 0.4 ]',&
               'set xrange [-0.12: 0.12 ]',&
               'unset key',&
               'set mytics 2',&
               'set parametric',&
               'set trange [-10:10]',&
               'set multiplot',&
			   'set arrow from 0,-0.4 to 0, 0.4 nohead lc rgb "red"',&
                  'set arrow from -0.12,0 to 0.12, 0  nohead lc rgb "red"',&
               'plot "band.dat" u 1:($2-ef) with l lt 1 lw 1 lc rgb "black"',&
               'unset multiplot'
     
         end subroutine write_plt
