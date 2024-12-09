module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="cube"
    real*8,parameter::emin=5,emax=8,eta=0.005
    integer*8,parameter::nb=4,eres=100
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    integer*4 i,j,k,l,nr,ie,count,ierr,NEV,nblocks,matsize,nk,ik,ib
    character(len=80) v_file,d_file, line
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,epoints(eres),a_spec,factor,p_l,de,dos,dk(3)
    real*8,allocatable:: rvec(:,:),rwork(:),real_d(:)
    complex*16,allocatable:: V(:,:)
    real*8,allocatable:: d(:,:),xk(:)

    pi2=4.0d0*atan(1.0d0)*2.0d0

    open(98,file='2D_EVALS.dat')
    open(99,file='2D_EVECS.dat')
    open(100,file='2D.dx')
    open(200,file='2D_xk.dat')

    read(98, *) nblocks, nk, dk
    NEV = nblocks * nblocks * nb
    matsize = nblocks * nblocks

    allocate(d(nk,NEV), V(NEV,NEV))

    allocate(xk(nk))
    read(200, *) xk

    do i=1, nk
        read(98, *) d(i,:)
    enddo

    de = (emax-emin)/eres
    do i=1, eres
        epoints(i) = emin + de*i
    enddo

    write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nk,eres
    write(100, '(a,3(1x,f12.8))') 'origin',-0.1d0,emin
    write(100, '(a,3(1x,f12.8))') 'delta',sqrt(dot_product(dk,dk)),0d0
    write(100, '(a,3(1x,f12.6))') 'delta',0d0,de
    write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nk,eres


    print *, dot_product(v(:,1),v(:,1))

    print *, "Calculating DOS..."
    count = 0 
    do ib = 1, matsize
        count = count + 1
        write(100, '(a,i8,a,i8,a,i10,a)') 'object',2+count,' class array type float rank 1 shape',3,' item', nk*eres, ' data follows'
        print *, ib

        do ik = 1, nk
            do j=1,NEV
                read(99, *) v(:,j)
            enddo
            do ie=1,eres
                print*, ie, eres

                a_spec = 0d0
                do i=1,NEV
                    p_l = dot_product(v((1+nb*ib):(nb*(ib+1)),i), v((1+nb*ib):(nb*(ib+1)),i))
                    factor = ((epoints(ie) - d(ik,i)))/eta
                    a_spec = a_spec + p_l * (exp(-0.5d0*factor**2)) * 1/(eta*sqrt(2*pi2))
                    ! print *, a_spec
                enddo

                write(100, '(3(1x,f20.12))') xk(ik), epoints(ie), real(a_spec)
            enddo
        enddo
        write(100, '(a)') 'attribute "dep" string "positions"' 
    enddo
    do i=0,matsize-1
        write(100,'(A,i8,A,/,A,/,A,/,A,i8,/)') &
        'object',matsize+3+i,' class field', &
        'component "positions" value 1', &
        'component "connections" value 2', &
        'component "data" value ',3+i
    enddo
    write(100, '(a)') 'object "series" class series'
    do i=0,matsize-1
        write(100, '(a,i8,a,i8,a,i8)') 'member', i, ' value', (i+matsize+3), ' position', i
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
