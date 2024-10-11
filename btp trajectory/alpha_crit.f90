Program Wannier_band_structure
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=200
!------------------------------------------------------
      character(len=30)::klabel(nkpath)
      character(len=80) top_file,triv_file,nnkp,line
      integer*4,parameter::nk=(nkpath-1)*np+1,alpha_res = 100
      integer*4 i,j,k,l,nr,i1,i2,j1,j2,nb,lwork,info,count, min_index
      real*8,parameter::third=1d0/3d0,alpha_max=0.05, alpha_crit = 0.8, z_plane = 0.5d0
      real*8 phase,pi2,jk,a,b,x1,y1,alpha,dalpha, min_eg, min_alpha
      real*8 avec(3,3),bvec(3,3)
      real*8 klist(3,1:nk),xk(nk),kpath(3,np),ktemp1(3),ktemp2(3),xkl(nkpath)
      real*8,allocatable:: rvec(:,:),rvec_data(:,:),rvec_data_t(:,:),ene(:,:),rwork(:), e_g_data(:),alphas(:)
      integer*4,allocatable:: ndeg(:)
      complex*16,allocatable::Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:)
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
      data kpath(:,1) /     0.5d0,      0.0d0,    z_plane/  !L
      data kpath(:,2) /     0.0d0,      0.0d0,    z_plane/  !A
      data kpath(:,3) /     third,      third,    z_plane/  !H
      data klabel     /'L','A','H'/

      ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

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
      klist=klist*pi2

!------read H(R)
      open(99,file=trim(adjustl(top_file)))
      open(97,file=trim(adjustl(triv_file)))
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec_data(3,nr),rvec_data_t(3,nr),Hk(nb,nb),top_hr(nb,nb,nr),triv_hr(nb,nb,nr),ndeg(nr),ene(nb,nk))
      read(99,*)ndeg
      do i = 1, 80
          read(97, *)! Read and discard 80 lines
      end do
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
               top_hr(i1,i2,k)=dcmplx(a,b)
               read(97,*)rvec_data_t(1,k),rvec_data_t(2,k),rvec_data_t(3,k),j1,j2,x1,y1
               triv_hr(j1,j2,k)=dcmplx(x1,y1)
            enddo
         enddo
      enddo
     lwork=max(1,2*nb-1)
     allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)),e_g_data((2*alpha_res)+1),alphas((2*alpha_res)+1))
!---- Fourier transform H(R) to H(k)
      dalpha = alpha_max/alpha_res
      ene=0d0
      count =1
      do l = -alpha_res, alpha_res
         alpha = (l*dalpha) + alpha_crit
         do k=1,nk
            HK=(0d0,0d0)
            do j=1,nr

               phase=0.0d0
               do i=1,3
                  phase=phase+klist(i,k)*rvec_data_t(i,j)
               enddo
               HK=HK+((1-alpha)*(triv_hr(:,:,j))+alpha*(top_hr(:,:,j)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))
            enddo
            call zheev('V','U',nb,HK,nb,ene(:,k),work,lwork,rwork,info) 
         enddo
         e_g_data(count) = MINVAL(ene(13,:))-MAXVAL(ene(12,:))
         print *, e_g_data(count) , alpha
         alphas(count) = alpha
         count = count +1
      enddo

      min_eg = MINVAL(e_g_data)
      min_index = MINLOC(e_g_data, 1)
      min_alpha = alphas(min_index)
      
      print * , "Min E_G for ", alpha_crit - alpha_max," < alpha < ",alpha_crit + alpha_max,": ", min_eg
      print * , "Alpha for this value: ", min_alpha
end program  
            
