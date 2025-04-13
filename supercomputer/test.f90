program test
    use iso_c_binding
    implicit none
    integer*4,parameter::n=10
    complex*16,allocatable::a(:,:),b(:,:),c(:,:)
    integer*4::i,j,numprocs,myid,ierr,x,y
    integer(c_long) :: global_maxrss_sum = 0
    integer(c_long) :: global_maxrss_max = 0
    real*8::memory_sum_gb,memory_max_gb,mem
    include 'mpif.h'
    interface
    function get_maxrss() bind(c)
        use iso_c_binding
        integer(c_long) :: get_maxrss
    end function get_maxrss
    end interface
    integer(c_long) :: local_maxrss = 0
    
    ! Allocate array before MPI_Init
    i = sqrt(4.0d0 * 10.0d0**9/16.0d0 )
    allocate(a(i,i), stat=ierr)
    if (ierr /= 0) then
        print *, "Error allocating array a before MPI_Init"
        stop
    endif
    local_maxrss = max(local_maxrss, get_maxrss())
    print *, "Memory after first allocation (before MPI_Init):", real(local_maxrss) / (1024.0 ** 3), "GB"
    
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
    
    ! Allocate another array after MPI_Init
    allocate(b(i,i), stat=ierr)
    if (ierr /= 0) then
        print *, "Error allocating array b after MPI_Init"
        stop
    endif
    local_maxrss = max(local_maxrss, get_maxrss())
    print *,"myid:",myid,"Memory after second allocation (after MPI_Init):", real(local_maxrss) / (1024.0 ** 3), "GB"

    do x = 1, size(a,1)
        do y = 1, size(a,2)
            a(x,y) = cmplx(1.0d0, 1.0d0)
            b(x,y) = cmplx(1.0d0, 1.0d0)
        end do
        local_maxrss = max(local_maxrss, get_maxrss())
    end do
    print*, "Size of a: ", size(a),' =',size(a)*16.0d0/10**9,"GB, i:", i

    print *,"myid:",myid, "Memory after loop :", real(local_maxrss) / (1024.0 ** 3), "GB"
    local_maxrss = max(local_maxrss, get_maxrss())
    
    call MPI_Reduce(local_maxrss, global_maxrss_sum, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(local_maxrss, global_maxrss_max, 1, MPI_INTEGER8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if (myid == 0) then
        memory_sum_gb = real(global_maxrss_sum) / (1024.0 ** 3)
        memory_max_gb = real(global_maxrss_max) / (1024.0 ** 3)
        print *, 'Final memory usage:'
        print *, 'Total across all processes: ', memory_sum_gb, ' GB'
        print *, 'Maximum per process: ', memory_max_gb, ' GB'
    end if

    deallocate(a, b, stat=ierr)

    call MPI_Finalize(ierr)
end program test