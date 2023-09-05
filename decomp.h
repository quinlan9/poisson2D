


/* #include "gfunc.h" */
int MPE_Decomp1d(int n, int nprocs,int myid, int *s, int *e);

//int MPE_Decomp2d(int n, int dims[2], int myid, int* s0, int* e0, int* s1, int* e1);
int MPE_Decomp2d(int nx, int dims[2],int coords[2], int myid, int *sx,int *ex, int *sy,int *ey);
//int MPE_Decomp2d(int nx, int ny, int* dims, int rank, int* coords, int* xs, int* xe, int* ys, int* ye);
