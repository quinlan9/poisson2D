#include<stdio.h>
#include<mpi.h>


int MPE_Decomp1d(int n, int nprocs, int myid, int* s, int* e) {
    // Decomposes a 1D array of size n among size processes, assigning a 
    // contiguous chunk of elements to each process
 
        int nlocal, leftover, offset;

        // Calculate the size of the local chunk of elements for each process
        nlocal = n / nprocs;

        // Calculate the number of leftover elements after dividing the work evenly among the processes
        leftover = n % nprocs;

        // Calculate the starting index of the local chunk of elements for the current process
        if (myid < leftover) {
            // If the current process is one of the first few processes with leftover elements,
            // it gets one extra element to distribute among its local chunk of elements
            nlocal++;
            offset = myid * nlocal;
        }
        else {
            // Otherwise, the starting index is calculated as usual
            offset = myid * nlocal + leftover;
        }

        // Set the starting and ending indices for the local chunk of elements
        *s = offset + 1;
        *e = offset + nlocal;

        // If the ending index goes beyond the end of the array, adjust it to the maximum index
        if (*e > n) {
            *e = n;
        }
	printf("myid:%d\n",myid);
        return MPI_SUCCESS;
 }

int MPE_Decomp2d(int nx, int dims[2],int coords[2], int myid, int *sx,int *ex, int *sy,int *ey) {
	int xnum,xleftover,xoffset;
	xnum=nx/dims[0];
	xleftover=nx%dims[0];
	
	if(coords[0]<xleftover){
		xnum++;
		xoffset=coords[0]*xnum;
	}else{
		xoffset=coords[0]*xnum+xleftover;	
	}
	*sx=xoffset+1;
	*ex=xoffset+xnum;	

	int ynum,yleftover,yoffset;
	ynum=nx/dims[1];
	yleftover=nx%dims[1];

	if(coords[1]<yleftover){
		ynum++;
		yoffset=coords[1]*ynum;
	}else{
		yoffset=coords[1]*ynum+yleftover;
	}
	*sy=yoffset+1;
	*ey=yoffset+ynum;
    return MPI_SUCCESS;
}


