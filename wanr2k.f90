      Program Wannier_band_structure
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=600

!------------------------------------------------------
      integer*4 ik,ikmax, skip
      real*8 kz,ef 
      real*8 :: Te_sum, Bi_sum, I_sum
      character(len=30)::klabel(nkpath)
      character(len=80) nnkp,line,top_file,triv_file
      integer*4,parameter::nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2,nb,l
      real*8,parameter::third=1d0/3d0, alpha = 0.77966
      real*8 phase,pi2,jk,a,b,x1,y1
      real*8 klist(3,1:nk),xk(nk),kpath(3,np),bvec(3,3),avec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath)
      real*8,allocatable:: rvec(:,:),rvec_data(:,:),ene(:,:),rwork(:),od(:,:,:)
      integer*4,allocatable:: ndeg(:)
      complex*16,allocatable:: Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:),H_col(:)
      complex*16 temp1,temp2
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
!---------------kpath
      data kpath(:,1) /     0.5d0,      0.0d0,    0.5d0/  !L
      data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
      data kpath(:,3) /     third,      third,    0.5d0/  !H

!      open(77,file='tmp')
 !     read(77,*)ik,ikmax
 !     kz=float(ik)*0.5d0/float(ikmax)
 !     kpath(3,:)=kz

      data klabel     /'L','A','H'/

      ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

!     xk(1)= 0d0 !-sqrt(dot_product(ktemp1,ktemp1))
      xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
      xkl(1)=xk(1)
      

      k=0
      ktemp1=0d0
      do i=1,nkpath-1
       do j=1,np
        k=k+1
        jk=dfloat(j-1)/dfloat(np)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
        if(j.eq.1) xkl(i)=xk(k)
        ktemp1=ktemp2
       enddo
      enddo
      klist(:,nk)=kpath(:,nkpath)
      ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      xkl(nkpath)=xk(nk)
!      write(*,*)klist
      klist=klist*pi2

!------read H(R)
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

!---- Fourier transform H(R) to H(k)
      ene=0d0
      do k=1,nk
         HK=(0d0,0d0)
         do j=1,nr

            phase=0.0d0
            do i=1,3
               phase=phase+klist(i,k)*rvec(i,j)
            enddo

            HK=HK+((1-alpha)*(triv_hr(:,:,j))+alpha*(top_hr(:,:,j)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

         enddo
         call zheev('V','U',nb,HK,nb,ene(:,k),work,lwork,rwork,info) 

!------Orbital Probability Calculation:
         do l=1, nb
            H_col = HK(:,l)
            Te_sum= 0.0d0
            Bi_sum = 0.0d0
            I_sum = 0.0d0
            do i=1, 2
               skip = 0
               skip = (i-1)*9;
               do j=1, 3
                  Te_sum = Te_sum + real(conjg(H_col(j+skip)) * (H_col(j+skip)))
                  Bi_sum = Bi_sum + real(conjg(H_col(j+3+skip)) * (H_col(j+3+skip)))
                  I_sum =  I_sum + real((conjg(H_col(j+6+skip)) * (H_col(j+6+skip))))
               enddo
            enddo
            od(l,k,:) = [Te_sum, Bi_sum, I_sum]
         enddo
      enddo

!-----Fermi level:
      ef = (MAXVAL(ene(12, :)) + MINVAL(ene(13, :)))/2.0d0
      print * , ef

!-----Orbital probability:


      deallocate(HK,work)
      
      do i=1,nb
         do k=1,nk
           write(100,'(5(x,f12.6))') xk(k),ene(i,k),od(i,k,:)
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
            write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [ -0.12 : 0.12]'
            write(99,'(a)') &
                 'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 5.00in'
            write(99,'(a,f4.2,a)')'set output "band.pdf"'
            write(99,'(14(a,/),a)') &
                 'set border',&
                 'unset xtics',&
                 'unset ytics',&
                 'set encoding iso_8859_1',&
                 'set size ratio 0 1.0,1.0',&
                 'set yrange [-0.4: 0.4 ]',&
                 'unset key',&
                 'set mytics 2',&
                 'set parametric',&
                 'set trange [-10:10]',&
                 'set multiplot',&
                 'plot "band.dat" every 4 u 1:($2-ef):(column(3)*4) with points pt 7 ps variable lc rgb "royalblue"',&
                 'plot "band.dat" every 4 u 1:($2-ef):(column(4)*4) with points pt 7 ps variable lc rgb "light-red"',&
                 'plot "band.dat" every 4 u 1:($2-ef):(column(5)*4) with points pt 7 ps variable lc rgb "forest-green"',&
                 'plot "band.dat" u 1:($2-ef) with l lt 1 lw 3.5 lc rgb "yellow"',&
                 'unset multiplot'
       
           end subroutine write_plt
