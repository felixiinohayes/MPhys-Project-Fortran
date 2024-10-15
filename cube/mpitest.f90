program hello
    use mpi
    implicit none
    integer :: ierr, myid, numprocs
    integer :: local_value, global_sum

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    ! Each process has a local value equal to its rank (myid)
    local_value = myid

    ! Perform an all-reduce operation to sum the local values across all processes
    call MPI_ALLREDUCE(local_value, global_sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Each process prints its own rank and the global sum
    print *, "Hello from process ", myid, " of ", numprocs, " | Global sum of ranks = ", global_sum

    call MPI_FINALIZE(ierr)
end program hello