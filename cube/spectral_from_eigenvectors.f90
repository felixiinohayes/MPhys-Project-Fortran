module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="cube"
    real*8,parameter::emin=-0.1,emax=0.1,eta=0.03
    integer*8,parameter::nb=18,eres=100
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    integer*4 i,j,k,l,nr,ie,count,ierr,NEV,nxblocks,nyblocks,nzblocks,matsize
    character(len=80) v_file,d_file, line
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,epoints(eres),a_spec,factor,p_l,de,dos
    real*8,allocatable:: rvec(:,:),rwork(:),real_d(:)
    complex*16,allocatable:: V(:,:)
    real*8,allocatable:: d(:)

    pi2=4.0d0*atan(1.0d0)*2.0d0

    open(98,file='data/slab_10_top_eigenvalues.dat')
    open(99,file='data/slab_10_top_eigenvectors.dat')
    open(100,file='data/slab_10.dx')

    read(98, *) NEV, nxblocks, nyblocks, nzblocks
    matsize=nxblocks*nyblocks*nzblocks
    print *, NEV, nxblocks, nyblocks, nzblocks

    allocate(d(NEV), V(nb*matsize,NEV))

    read(98, *) d
    do i=1, NEV
        read(99, *) v(:,i)
    enddo
    ! print *, v(:,1)
    ! print *, d(1)

    ! print *, v(:,1)

    de = (emax-emin)/eres
    do i=1, eres
        epoints(i) = emin + de*i
    enddo

    write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nzblocks,nxblocks,nyblocks
    write(100, '(a,3(1x,f12.8))') 'origin',0d0,0d0,0d0
    write(100, '(a,3(1x,f12.8))') 'delta',0d0,0d0,1d0
    write(100, '(a,3(1x,f12.8))') 'delta',0d0,1d0,0d0
    write(100, '(a,3(1x,f12.6))') 'delta',1d0,0d0,0d0
    write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nzblocks,nxblocks,nyblocks

    ! write(100, '(a,i10,a)') 'object 3 class array type float rank 1 shape 1 item', matsize, ' data follows'

    ! i=500

    ! do j=0,matsize-1
    !     p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))
    !     write(100, '(f12.10)') p_l
    ! enddo

    print *, dot_product(v(:,1),v(:,1))

    ! write(100, '(a)') 'attribute "dep" string "positions"' 

    print *, "Calculating DOS..."
    count = 0 
    do ie=1,eres
        count = count + 1
        print*, ie, eres

        write(100, '(a,i8,a,i8,a,i10,a)') 'object',2+count,' class array type float rank 1 shape',1,&
                                ' item', matsize, ' data follows'
        !----Spectral DOS
        do j=0,matsize-1
            a_spec = 0d0
            do i=1,NEV
                p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, i), v( 1+(j*nb) : (j+1)*nb, i))
                factor = ((epoints(ie) - d(i)))/eta
                a_spec = a_spec + p_l * (exp(-0.5d0*factor**2)) * 1/(eta*sqrt(2*pi2))
                ! if(ie==1) print *, d(i), p_l, factor, a_spec
                ! if(j==1) print *, a_spec
            enddo
            write(100, '(1(1x,f12.8))') a_spec
        enddo
        write(100, '(a)') 'attribute "dep" string "positions"' 
    enddo
    do i=0,eres-1
        write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
        'object',eres+3+i,' class field', &
        'component "positions" value 1', &
        'component "connections" value 2', &
        'component "data" value ',3+i
    enddo
    write(100, '(a)') 'object "series" class series'
    do i=0,eres-1
        write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+eres+3), ' position', i
    enddo

    ! write(100,'(A,/,A,/,A,/,A,/)') &
    ! 'object 4 class field', &
    ! 'component "positions" value 1', &
    ! 'component "connections" value 2', &
    ! 'component "data" value 3'
    ! write(100, '(a)') 'object "series" class series'
    ! write(100, '(a)') 'member 0 value 4 position 0'

    write(100, '(A)') 'end'

end program Projected_band_structure
