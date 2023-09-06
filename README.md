# poisson2D

brief: Poisson solver which uses a 2D processor decomposition  
Used for parallel solving of the discretized Poisson equation. The Jacobi iterative method is employed as the solver, and MPI is utilized for parallel processing to enhance computational efficiency.  

Implement four versions of the 2D ghost exchange routine:  
sendrecv.c: A version which uses sendrecv()  
nonblocking.c: A version which uses non-blocking sends and receives  
rma_fence.c: A version which uses MPI_Win_fence based synchronization along with MPI_Put and MPI_Get  
rma_pscw.c: A version which uses general active target synchronization: MPI_Win_start, MPI_Win_complete, MPI_Win_post and
MPI_Win_wait along with MPI_Put and MPI_Get  

load MPI module on chuck: module load cports openmpi  
then compile with command: make  
to run rma_fence.c, input the command: mpirun -np 4 ./rma_fence 31  
to run rma_pscw.c, input the command: mpirun -np 4 ./rma_pscw 31  
to run sendrecv.c, input the command: mpirun -np 4 ./sendrecv 31  
to run nonblocking.c, input the command: mpirun -np 4 ./nonblocking 31  

then you will get the outcome file from the first two c file,namely"grid_outcomes_fence" and "grid_outcomes_pscw"  
and then you can see it mathes the outcome of non-RMA version, which named"grid_outcomes_2d".  
  
