module parameters
    Implicit None
!--------to be modified by the user
    character(len=80):: prefix="cube"
    real*8,parameter::ef= 4.18903772,emin=3,emax=5,eta1=0.05
    integer*8,parameter::nb=18,nblocks=4,matsize=(nblocks)**3,eres=150,NEV=1
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
    integer*4 i,j,k,l,nr,ie,count
    character(len=80) v_file,d_file
    real*8 avec(3,3),bvec(3,3),pi2,x1,x2,y1,y2,epoints(eres),a_spec,factor,p_l,de,dos
    real*8,allocatable:: rvec(:,:),rwork(:),real_d(:)
    complex*16 V(matsize*nb,NEV)
    real*8 d(NEV)

    pi2=4.0d0*atan(1.0d0)*2.0d0

    write(v_file,'(a,a)')trim(adjustl(prefix)),"_eigenvectors_shift.dat"
    write(d_file,'(a,a)')trim(adjustl(prefix)),"_eigenvalues_shift.dat"
    open(98,file=trim(adjustl(v_file)))
    open(99,file=trim(adjustl(d_file)))

    open(100,file='DOS_cube_single.dx')

    read(98,*) v
    read(99,*) d

    write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nblocks,nblocks,nblocks
    write(100, '(a,3(1x,f12.8))') 'origin',0d0,0d0,0d0
    write(100, '(a,3(1x,f12.8))') 'delta',0d0,0d0,1d0
    write(100, '(a,3(1x,f12.8))') 'delta',0d0,1d0,0d0
    write(100, '(a,3(1x,f12.6))') 'delta',1d0,0d0,0d0
    write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nblocks,nblocks,nblocks

    write(100, '(a,i10,a)') 'object 3 class array type float rank 1 shape 1 item', matsize, ' data follows'

    ! i=500

    do j=0,matsize-1
        p_l = dot_product( v( 1+(j*nb) : (j+1)*nb, 1), v( 1+(j*nb) : (j+1)*nb, 1))
        write(100, '(f12.10)') p_l
    enddo

    print *, dot_product(v(:,1),v(:,1))

    write(100, '(a)') 'attribute "dep" string "positions"' 


    write(100,'(A,/,A,/,A,/,A,/)') &
    'object 4 class field', &
    'component "positions" value 1', &
    'component "connections" value 2', &
    'component "data" value 3'
    write(100, '(a)') 'object "series" class series'
    write(100, '(a)') 'member 0 value 4 position 0'

    write(100, '(A)') 'end'

end program Projected_band_structure