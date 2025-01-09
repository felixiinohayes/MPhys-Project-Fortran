      Program Wannier_band_structure
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="../BiTeI"
      integer,parameter::nkpath=10,np=600
!------------------------------------------------------
      integer*4 ik,ikmax, skip,sign,interp_size,nr_top,nr_triv,index,ir,ix,iy,iz,rvec(3),rvec_data(3)
      real*8 kz,ef,eg
      real*8 :: Te_sum, Bi_sum, I_sum
      character(len=30)::klabel(nkpath)
      character(len=80) nnkp,line,top_file,triv_file
      integer*4,parameter::nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2,nb,l
      real*8,parameter::third=1d0/3d0, a = 1
      real*8 phase,pi2,jk,b,x1,y1
      real*8 xk(nk),bvec(3,3),avec(3,3),kvec1(3),kvec2(3),xkl(nkpath),kpoints(3,nkpath),kpath(3,nk),dk(3)
      real*8,allocatable:: ene(:,:),rwork(:),od(:,:,:)
      integer*4,allocatable:: ndeg(:,:,:),ndeg_top(:),ndeg_triv(:),rvec_top(:,:)
      complex*16,allocatable:: Hk(:,:),work(:),H_col(:),top_Hr_temp(:,:),triv_Hr_temp(:,:)
      complex*16 temp1,temp2
      complex*16,dimension(:,:,:,:,:),allocatable :: interp_Hr,top_Hr,triv_Hr
!------------------------------------------------------
      write(top_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological_new.dat"
      write(triv_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial_new.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),"_ortho.nnkp"
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
      open(99,file=trim(adjustl(top_file)))
      open(97,file=trim(adjustl(triv_file)))
      open(100,file='band.dat')
!---------------kpath
    !   data kpoints(:,1) /     0.5d0,      0.0d0,    0.5d0/  !L
    !   data kpoints(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
    !   data kpoints(:,3) /     third,      third,    0.5d0/  !H

    ! data kpoints(:,1) /     0.5d0,      0.0d0,    0.5d0/  !L
    ! data kpoints(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
    ! data kpoints(:,3) /     0.0d0,      0.5d0,    0.5d0/  !H

    kpoints(:,1) = [ 0.0d0,    0.0d0,   0.0d0]  !Gamma  
    kpoints(:,2) = [ 0.5d0,    0.0d0,   0.0d0]  !X
    kpoints(:,3) = [ 0.5d0,    0.5d0,   0.0d0]  !S
    kpoints(:,4) = [ 0.0d0,    0.5d0,   0.0d0]  !Y
    kpoints(:,5) = [ 0.0d0,    0.0d0,   0.0d0]  !Gamma
    kpoints(:,6) = [ 0.0d0,    0.0d0,   0.5d0]  !Z
    kpoints(:,7) = [ 0.5d0,    0.0d0,   0.5d0]  !U
    kpoints(:,8) = [ 0.5d0,    0.5d0,   0.5d0]  !R
    kpoints(:,9) = [ 0.0d0,    0.5d0,   0.5d0]  !T
    kpoints(:,10) = [ 0.0d0,    0.0d0,   0.5d0]  !Z
    
    !   data klabel     /'L','A','H'/

      kvec1(:)=(kpoints(1,1)-kpoints(1,2))*bvec(:,1)+(kpoints(2,1)-kpoints(2,2))*bvec(:,2)+(kpoints(3,1)-kpoints(3,2))*bvec(:,3)
      xk(1)= -sqrt(dot_product(kvec1,kvec1))
  
      kvec1 = 0d0
      do i = 1, nkpath-1
          do j = 1, np
              ik = j + np*(i-1)
              dk = (kpoints(:,i+1)-kpoints(:,i))/np
              kpath(:, ik) = kpoints(:,(i)) + (dk*(j-1))
              kvec2 = kpath(1,ik)*bvec(:,1) + kpath(2,ik)*bvec(:,2) + kpath(3,ik)*bvec(:,3) 
              if(ik.gt.1) xk(ik) =  xk(ik-1) + sqrt(dot_product(kvec2-kvec1,kvec2-kvec1))
              kvec1 = kvec2
          enddo
      enddo
      ! Final point in the kpath
      kpath(:,nk) = kpoints(:,nkpath)
      kvec2=kpath(1,nk)*bvec(:,1)+kpath(2,nk)*bvec(:,2)+kpath(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(kvec2-kvec1,kvec2-kvec1))
      kpath = kpath*pi2

!------read H(R) 



      read(99,*)
      read(99,*)nb,nr_top
      read(97,*)
      read(97,*)nb,nr_triv

      allocate(Hk(nb,nb),ene(nb,nk),od(nb,nk,3),H_col(nb))
      allocate(top_Hr_temp(nb,nb),triv_Hr_temp(nb,nb),ndeg_top(nr_top),ndeg_triv(nr_triv),ndeg(-6:6, -6:6, -6:6))
      allocate(interp_Hr(nb,nb,-6:6, -6:6, -6:6),top_Hr(nb,nb,-6:6, -6:6, -6:6),triv_Hr(nb,nb,-6:6, -6:6, -6:6))
      allocate(rvec_top(nr_top,3))
     
      read(99,*)ndeg_top
      read(97,*)ndeg_triv

      interp_Hr=0d0
      !----Read in Hamiltonian
          ndeg = 0d0
          
          do ir=1,nr_top
              do i=1,nb
                  do j=1,nb
                     read(99,*)rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3),i1,i2,x1,y1
                     top_Hr_temp(i1,i2)=dcmplx(x1,y1)
                     ndeg(rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = ndeg_top(ir)
                  enddo
              enddo
              top_Hr(:,:,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = top_Hr_temp(:,:)
          enddo
          do ir=1,nr_triv
              do i=1,nb
                  do j=1,nb
                     read(97,*)rvec(1),rvec(2),rvec(3),i1,i2,x1,y1
                     triv_Hr_temp(i1,i2)=dcmplx(x1,y1)
                  enddo
              enddo
              triv_Hr(:,:,rvec(1),rvec(2),rvec(3)) = triv_Hr_temp(:,:)
          enddo
      
          ! Interpolate Hamiltonian and add magnetic field
          do ir=1,nr_top
              do i=1,nb
                  do j=1,nb
                      interp_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) = (1-a)*triv_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3)) + a*top_Hr(i,j,rvec_top(ir,1),rvec_top(ir,2),rvec_top(ir,3))
                  enddo
              enddo
          enddo

     lwork=max(1,2*nb-1)
        allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
        print*,nr_top
        print*,'interp_Hr(5,4,-3,-2,-1)', interp_Hr(5,4,-3,-2,-1),'   -3   -2   -1    5    4 -0.00275907 -0.00075338'  
        !                           1    1  -2   -2   -1    

!---- Fourier transform H(R) to H(k)
      ene=0d0
      do k=1,nk
        HK=(0d0,0d0)
        do ix=-6,6
            do iy=-6,6
                do iz=-6,6
                    if(ndeg(ix,iy,iz).ne.0) then
                        phase=0.0d0
                        phase = phase + kpath(1,k) * ix + kpath(2,k) * iy +kpath(3,k) * iz
                       
                        Hk=Hk+((1-a)*(triv_Hr(:,:,ix,iy,iz))+(a)*(top_Hr(:,:,ix,iy,iz)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(ix,iy,iz))
                    endif
                enddo
            enddo
        enddo
      
         call zheev('V','U',nb,HK,nb,ene(:,k),work,lwork,rwork,info) 

!------Orbital Probability Calculation:
      !    do l=1, nb
      !       H_col = HK(:,l)
      !       Te_sum= 0.0d0
      !       Bi_sum = 0.0d0
      !       I_sum = 0.0d0
      !       do i=1, 2
      !          skip = 0
      !          skip = (i-1)*9;
      !          do j=1, 3
      !             Te_sum = Te_sum + real(conjg(H_col(j+skip)) * (H_col(j+skip)))
      !             Bi_sum = Bi_sum + real(conjg(H_col(j+3+skip)) * (H_col(j+3+skip)))
      !             I_sum =  I_sum + real((conjg(H_col(j+6+skip)) * (H_col(j+6+skip))))
      !          enddo
      !       enddo
      !       od(l,k,:) = [Te_sum, Bi_sum, I_sum]
      !    enddo
      enddo

!-----Fermi level:
      ! ef = (MAXVAL(ene(12, :)) + MINVAL(ene(13, :)))/2.0d0
      ! eg = MINVAL(ene(13, :)) - MAXVAL(ene(12, :))

      ! print * , ef,eg

      deallocate(HK,work)
      
      do i=1,nb
         do k=1,nk-1
           write(100,'(5(x,f10.4))') xk(k), ene(i,k) 
         enddo
           write(100,*)
           write(100,*)
      enddo
      call write_plt(nkpath,xkl,klabel,ef)
      stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(triv_file)),' not found'
      stop

      end

      subroutine write_plt(nkp,xkl,kl,ef)
            implicit none
            integer nkp,i
            real*8 xkl(nkp),ef
            character(len=30)kl(nkp)
            
            open(99,file='band.plt')
            write(99,'(a,f12.8)')'ef=',ef
            write(99,'(a,f12.6,a,f12.6,a)') '#set xrange [ -5 : 5]'
            write(99,'(a)') &
                 'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 5.00in'
            write(99,'(a,f4.2,a)')'set output "band.pdf"'
            write(99,'(14(a,/),a)') &
                 'set border',&
                 'set xtics',&
                 'set ytics',&
                 'set encoding iso_8859_1',&
                 'set size ratio 0 1.0,1.0',&
                 'set yrange [-1: 8]',&
                 'unset key',&
                 'set mytics 2',&
                 'set parametric',&
                 'set trange [-10:10]',&
                 'set multiplot',&
                 '#plot "band.dat" every 4 u 1:($2-ef):(column(3)*4) with points pt 7 ps variable lc rgb "royalblue"',&
                 '#plot "band.dat" every 4 u 1:($2-ef):(column(4)*4) with points pt 7 ps variable lc rgb "light-red"',&
                 '#plot "band.dat" every 4 u 1:($2-ef):(column(5)*4) with points pt 7 ps variable lc rgb "forest-green"',&
                 'plot "band.dat" u 1:($2-ef) with l lt 1 lw 3.5 lc rgb "black"',&
                 'unset multiplot'
       
           end subroutine write_plt
