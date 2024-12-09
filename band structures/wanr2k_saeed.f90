      Program Wannier_band_structure
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="../BiTeI"
      integer,parameter::nkpath=3,np=100
!------------------------------------------------------
      integer*4 ik,ikmax
      real*8 kz
      character(len=30)::klabel(nkpath)
      character(len=80) hamil_file,nnkp,line
      integer*4,parameter::nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,i1,i2,nb,lwork,info
      real*8,parameter::third=1d0/3d0!,kz=0d0
      real*8 phase,pi2,jk,a,b,ef,eg
      real*8 klist(3,1:nk),xk(nk),kpath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath)
      real*8,allocatable:: rvec(:,:),ene(:,:),rwork(:)
      integer*4,allocatable:: ndeg(:)
      complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:)
      complex*16 temp1,temp2
!------------------------------------------------------
      write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial_4band.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

      pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
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

!      xk(1)= 0d0 !-sqrt(dot_product(ktemp1,ktemp1))
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
      open(99,file=trim(adjustl(hamil_file)),err=444)
      open(100,file='band.dat')
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec(3,nr),Hk(nb,nb),Hamr(nb,nb,nr),ndeg(nr),ene(nb,nk))
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec(1,k),rvec(2,k),rvec(3,k),i1,i2,a,b
               hamr(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo

     lwork=max(1,2*nb-1)
     allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

!---- Fourrier transform H(R) to H(k)
      ene=0d0
      do k=1,nk
         HK=(0d0,0d0)
         do j=1,nr

            phase=0.0d0
            do i=1,3
               phase=phase+klist(i,k)*rvec(i,j)
            enddo

            HK=HK+Hamr(:,:,j)*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

         enddo

         call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
         
      enddo

      !-----Fermi level:
      ef = (MAXVAL(ene(2, :)) + MINVAL(ene(3, :)))/2.0d0
      eg = MINVAL(ene(3, :)) - MAXVAL(ene(2, :))

      print*, ef,eg

      deallocate(HK,work)
      
      do i=1,nb
         do k=1,nk
           write(100,'(2(x,f12.6))') xk(k),ene(i,k)
         enddo
           write(100,*)
           write(100,*)
      enddo
      call write_plt(nkpath,xkl,klabel,ef)
      stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file)),' not found'
      stop

      end

     subroutine write_plt(nkp,xkl,kl,ef)
     implicit none
     integer nkp,i
     real*8 xkl(nkp),ef
     character(len=30)kl(nkp)
     
     open(99,file='band.plt')
     write(99,'(a,f12.8)')'ef=',ef
     write(99,'(a)') 'set xtics ( \'
     do i=1,nkp
        if(trim(adjustl(kl(i))).eq.'g'.or.trim(adjustl(kl(i))).eq.'G')kl(i)="{/Symbol \107}"
        if(i.ne.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i),", \"
        if(i.eq.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i)," )"
     enddo
     write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [',xkl(1),':',xkl(nkp),']'
     write(99,'(a)') &
          'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in'
     write(99,'(a,f4.2,a)')'set output "band.pdf"'
     write(99,'(9(a,/),a)') &
          'set encoding iso_8859_1',&
          'set size ratio 0 1.0,1.0',&
          'set ylabel "E-E_{CBM} (eV)"',&
          'set yrange [ -2 : 2.0 ]',&
          'unset key',&
          'set ytics 1.0 scale 1 nomirror out',&
          'set mytics 2',&
          'set parametric',&
          'set trange [-10:10]',&
          'plot "band.dat" u 1:($2-ef) with l lt 1 lw 3,\'
    do i=2,nkp-1
      write(99,'(f12.6,a)') xkl(i),',t with l lt 2  lc -1,\'
    enddo
    write(99,'(a)') 't,0 with l lt 2  lc -1'
    end subroutine write_plt
