module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="cube"
    real*8,parameter::emin=-0.3,emax=0.3,eta=0.02
    integer*8,parameter::nb=4,eres=500
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    integer*4 i,j,k,l,nr,ie,count,i1,i2,i3,j1,j2,j3,r1,r2,r3,xindex,yindex,NEV,nxblocks,nyblocks,nzblocks,matsize
    character(len=80) v_file,d_file,out_file
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,epoints(eres),a_spec,factor,p_l,de,dos
    real*8,allocatable:: rvec(:,:),rwork(:),real_d(:)
    complex*16,allocatable:: V(:,:)
    real*8,allocatable:: d(:)

    pi2=4.0d0*atan(1.0d0)*2.0d0

    write(v_file,'(a,a)') "data/C30_TOP_B_EVECS.dat"
    write(d_file,'(a,a)') "data/C30_TOP_B_EVALS.dat"
    write(out_file,'(a,a)') "data/C30_TOP_B_BULK.dat"
    open(98,file=trim(adjustl(d_file)))
    open(99,file=trim(adjustl(v_file)))
    open(100,file=trim(adjustl(out_file)))

    read(98, *) NEV, nxblocks, nyblocks, nzblocks
    matsize=nxblocks*nyblocks*nzblocks
    print *, NEV, nxblocks, nyblocks, nzblocks

    allocate(d(NEV), V(nb*matsize,NEV))

    read(98, *) d
    do i=1, NEV
        read(99, *) v(:,i)
    enddo

    ! print *, v(:,1)
    ! print *, d

    ! print *, v(:,1)

    de = (emax-emin)/eres
    do i=1, eres
        epoints(i) = emin + de*i
    enddo

    ! i=500

    ! do j=0,matsize-1
    !     p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))
    !     write(100, '(f12.10)') p_l
    ! enddo

    print *, dot_product(v(:,1),v(:,1))

    ! write(100, '(a)') 'attribute "dep" string "positions"' 
    do i3=0,nzblocks-1
        do i2=0,nyblocks-1
                do i1=0,nxblocks-1
                    xindex = i3*(nxblocks*nyblocks)+i2*(nxblocks)+i1
                    if(i3==nzblocks/2 .and. i2 ==nyblocks/2 .and. i1==nxblocks/2) then
                        j = xindex
                    endif
            enddo
        enddo
    enddo

    print *, j

    print *, "Calculating DOS..."
    count = 0 
    do ie=1,eres
        count = count + 1
        ! print*, ie, eres
        !----Spectral DOS
        a_spec = 0d0
        do i=1,NEV
            p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))
            factor = ((epoints(ie) - d(i)))/eta
            a_spec = a_spec + p_l * (exp(-0.5d0*factor**2)) * 1/(eta * sqrt(2*pi2))
            ! if(ie==1) print *, d(i), p_l, factor, a_spec
            ! if (ie==50) print *, epoints(ie), d(i), factor
        enddo
        write(100, '(2(1x,f12.8))') a_spec, epoints(ie)
        ! write(100, '(1(1x,f12.8))') p_l
    enddo

end program Projected_band_structure
