#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <mpi.h>

#include "poisson.h"
#include "jacobi.h"

#define maxit 20000
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

#include "decomp.h"
void init_full_grid(double g[][maxn]);
void init_full_grids(double a[][maxn], double b[][maxn], double f[][maxn]);

void twodinit_basic(double a[][maxn], double b[][maxn], double f[][maxn],
    int nx, int ny, int sx, int ex, int sy, int ey);

void print_full_grid(double x[][maxn]);
void print_in_order(double x[][maxn], MPI_Comm comm);
void print_grid_to_file(char* fname, double x[][maxn], int nx, int ny);
void write_grid(char* fname, double a[][maxn], int s, int e, int nx, MPI_Comm comm);
void GatherGrid2D(double a[][maxn], int nx, int nprocs, int dims[2], int coords[2], int myid, int sx, int ex, int sy, int ey, MPI_Comm comm);


int main(int argc, char** argv) {
    double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
    int nx, ny;

    int myid, nprocs;
    /* MPI_Status status; */
    int it;
    int sx, ex, sy, ey;
    double glob_diff;
    double ldiff;
    double t1, t2;
    double tol = 1.0E-11;


    /*initialize MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //check the number of command line arguments and set the nx accordingly
    if (myid == 0) {
        /* set the size of the problem */
        if (argc > 2) {
            fprintf(stderr, "---->Usage: mpirun -np <nproc> %s <nx>\n", argv[0]);
            fprintf(stderr, "---->(for this code nx=ny)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (argc == 2) {
            nx = atoi(argv[1]);
        }
        if (argc == 1) {
            nx = 15;
        }

        if (nx > maxn - 2) {
            fprintf(stderr, "grid size too large\n");
            exit(1);
        }
    }

    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcast nx from process 0 to all processes
    printf("(myid: %d) nx = %d\n", myid, nx);
    ny = nx;

    init_full_grids(a, b, f); // a b f

    MPI_Comm comm_cart;
    int ndims = 2; //2D
    int dims[2];
    int periodic[2] = { 0,0 };
    int reorder = 0;
    int coords[2];
    int nbrs[4]; //store up down left right rank

    MPI_Dims_create(nprocs, 2, dims);
    if (myid == 0) {
        printf("dims[0]: %d, dims[1]: %d\n", dims[0], dims[1]);
    }

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, &comm_cart);//create a Cartesian communicator

    MPI_Cart_coords(comm_cart, myid, 2, coords);//shift the rank to coorinate
    printf("myid:%d,coords[0]=%d,coords[1]=%d\n", myid, coords[0], coords[1]);

    MPI_Cart_shift(comm_cart, 0, 1, &nbrs[LEFT], &nbrs[RIGHT]); //left and right
    MPI_Cart_shift(comm_cart, 1, 1, &nbrs[UP], &nbrs[DOWN]); //get the rank of up and down of the current location


    MPE_Decomp2d(nx, dims, coords, myid, &sx, &ex, &sy, &ey);//divied the area into sub block, each block assigned to specific myid
    printf("nx:%d (myid: %d) coordinate:(%d,%d)  xs: %d,xe: %d; ys:%d, ye: %d; nbup: %d; nbdown: %d; nbrleft: %d; nbrright: %d\n", nx, myid, coords[0], coords[1], sx, ex, sy, ey,
        nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);

    twodinit_basic(a, b, f, nx, ny, sx, ex, sy, ey);
    //print_in_order(a, MPI_COMM_WORLD);
    t1 = MPI_Wtime(); //get current time

    
    MPI_Win wina, winb;
    MPI_Win_create(&a[sx - 1][0], (maxn) * (ex - sx + 3) * sizeof(double),sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &wina);
    MPI_Win_create(&b[sx - 1][0], (maxn) * (ex - sx + 3) * sizeof(double), sizeof(double),MPI_INFO_NULL, MPI_COMM_WORLD, &winb);
  

    MPI_Group ranks[4];
    MPI_Group group;

    MPI_Comm_group(MPI_COMM_WORLD, &group);
    
    if(coords[0]==0){
        MPI_Group_incl(group, 1,&(nbrs[RIGHT]),&ranks[0]);
    }
    if (coords[1]!=0) {
        MPI_Group_incl(group, 1,&(nbrs[UP]), &ranks[1]);
    }
    if (coords[0] != 0) {
        MPI_Group_incl(group, 1, &(nbrs[LEFT]), &ranks[2]);
    }
    if (coords[1]==0) { 
        MPI_Group_incl(group, 1,&(nbrs[DOWN]), &ranks[3]);
    }
    glob_diff = 1000; 
    for (it = 0; it < maxit; it++) { 

        exchange_rma_pscw_2d(a, wina, ranks,coords,nx, sx, ex, sy, ey, nbrs[LEFT], nbrs[RIGHT], nbrs[UP], nbrs[DOWN]);
        sweep2d(a, f, nx, sx, ex, sy, ey, b);


        exchange_rma_pscw_2d(b, winb, ranks,coords,nx, sx, ex, sy, ey, nbrs[LEFT], nbrs[RIGHT], nbrs[UP], nbrs[DOWN]);
        sweep2d(b, f, nx, sx, ex, sy, ey, a);

        ldiff = griddiff_2d(a, b, nx, sx, ex, sy, ey);//Calculate the difference between the subgrids for this process

        MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// Reduce the local difference to a global difference using MPI_Allreduce
        if (myid == 0 && it % 10 == 0) {
            // Print the current global and local differences every 10 iterations (only for process 0)
            printf("(myid %d) locdiff: %lf; glob_diff: %lf\n", myid, ldiff, glob_diff);
        }
        /* if( it%5==0 ){ */
        /*   print_in_order(a, MPI_COMM_WORLD); */
        /* } */

        //Break the loop if the global difference is below the tolerance
        if (glob_diff < tol) {
            if (myid == 0) {
                printf("iterative solve converged\n");
            }
            break;
        }
    }
    t2 = MPI_Wtime();

    printf("DONE! (it: %d)\n", it);

    if (myid == 0) {
        if (it == maxit) {
            fprintf(stderr, "Failed to converge\n");
        }
        printf("Run took %lf s\n", t2 - t1);
    }

    print_in_order(a, MPI_COMM_WORLD);
    
    if (nprocs == 1) {
        print_grid_to_file("grid", a, nx, ny);
        print_full_grid(a);
    }

    //Gather_grid_2d(a, nx, s, e, comm_cart);
    GatherGrid2D(a, nx, nprocs, dims, coords, myid, sx, ex, sy, ey, comm_cart);
    if (myid == 0) {
        write_grid("grid_outcomes_pscw.txt", a, 0, nx + 1, nx, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

void twodinit_basic(double a[][maxn], double b[][maxn], double f[][maxn],
    int nx, int ny, int sx, int ex, int sy, int ey)
{
    int i, j;

    /* set everything to 0 first */
    for (i = sx - 1; i <= ex + 1; i++) {
        for (j = sy - 1; j <= ey + 1; j++) {
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

    /* deal with boundaries */


    if (sx == 1) {
        for (j = sy - 1; j <= ey + 1; j++) {
            double y = 1.0 * j / ((double)(ny + 1));
            a[0][j] = y / (1.0 + y * y);
            b[0][j] = y / (1.0 + y * y);
        }
    }

    if (ex == nx) {
        for (j = sy - 1; j <= ey + 1; j++) {
            double y = 1.0 * j / ((double)(1 + nx));
            a[nx + 1][j] = y / (y * y + 4.0);
            b[nx + 1][j] = y / (y * y + 4.0);
        }

    }

    if (sy == 1) {
        for (i = sx - 1; i <= ex + 1; i++) {
            a[i][0] = 0;
            b[i][0] = 0;
        }
    }


    if (ey == nx) {
        for (i = sx - 1; i <= ex + 1; i++) {
            double x = 1.0 * i / ((double)(nx + 1));
            a[i][nx + 1] = 1 / (1 + (1 + x) * (1 + x));
            b[i][nx + 1] = 1 / (1 + (1 + x) * (1 + x));
        }

    }
}

void init_full_grid(double g[][maxn])
{
    int i, j;
    const double junkval = -5;

    for (i = 0; i < maxn; i++) {
        for (j = 0; j < maxn; j++) {
            g[i][j] = junkval;
        }
    }
}

/* set global a,b,f to initial arbitrarily chosen junk value */
void init_full_grids(double a[][maxn], double b[][maxn], double f[][maxn])
{
    int i, j;
    const double junkval = -5;

    for (i = 0; i < maxn; i++) {
        for (j = 0; j < maxn; j++) {
            a[i][j] = junkval;
            b[i][j] = junkval;
            f[i][j] = junkval;
        }
    }

}

/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn])
{
    int i, j;
    for (j = maxn - 1; j >= 0; j--) {
        for (i = 0; i < maxn; i++) {
            if (x[i][j] < 10000.0) {
                printf("|%2.6lf", x[i][j]);
            }
            else {
                printf("%9.2lf ", x[i][j]);
            }
        }
        printf("\n");
    }

}

void print_in_order(double x[][maxn], MPI_Comm comm)
{
    int myid, size;
    int i;

    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    MPI_Barrier(comm);
    printf("Attempting to print in order\n");
    sleep(1);
    MPI_Barrier(comm);

    for (i = 0; i < size; i++) {
        if (i == myid) {
            printf("proc %d\n", myid);
            print_full_grid(x);
        }
        fflush(stdout);
        usleep(500);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void print_grid_to_file(char* fname, double x[][maxn], int nx, int ny)
{
    FILE* fp;
    int i, j;

    fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "Error: can't open file %s\n", fname);
        exit(4);
    }

    for (j = ny + 1; j >= 0; j--) {
        for (i = 0; i < nx + 2; i++) {
            fprintf(fp, "%lf ", x[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/*
allows each processor to write
it  s part of the grid to file or stdout. The output should be in
  mesh/grid   format
*/
void write_grid(char* fname, double a[][maxn], int s, int e, int nx, MPI_Comm comm)
{
    FILE* fp;
    int i, j, n;
    int myid, size;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    //file or stdout
    fp = stdout;
    if (fp != NULL) {
        fp = fopen(fname, "w");
    }

    for (n = 0; n < size; n++) {
        if (myid == n) {
            for (j = maxn - 1; j >= 0; j--) {
                for (i = s; i <= e; i++) {
                    fprintf(fp, "%lf ", a[i][j]);
                }
                fprintf(fp, "\n");
            }
        }
    }

    fclose(fp);
}



/*
gathers the parallel solution which
is distributed across multiple processors onto one processor (rank0)
*/
void GatherGrid2D(double a[][maxn], int nx, int nprocs, int dims[2], int coords[2], int myid, int sx, int ex, int sy, int ey, MPI_Comm comm) {

    if (myid != 0) {
        for (int i = 0; i < ex - sx + 1; i++) {
            MPI_Send(&a[sx + i][sy], ey - sy + 1, MPI_DOUBLE, 0, myid, comm);

        }
    }
    if (myid == 0) {
        /*initialize matrix a */
        for (int i = ex + 1; i < nx + 2; i++) {
            for (int j = 0; j < nx + 2; j++) {
                a[i][j] = 0.0;
            }
        }
        for (int i = 0; i < nx + 2; i++) {
            for (int j = ex + 1; j < nx + 2; j++) {
                a[i][j] = 0.0;
            }
        }
        //receive fron other processors
        for (int k = 1; k < nprocs; k++) {
            MPI_Cart_coords(comm, k, 2, coords);
            int xleftover = nx % dims[0];
            int yleftover = nx % dims[1];
            int xnum = (coords[0] < xleftover) ? (nx / dims[0] + 1) : (nx / dims[0]);
            int ynum = (coords[1] < yleftover) ? (nx / dims[1] + 1) : (nx / dims[1]);
            int xstart = (coords[0] < xleftover) ? (coords[0] * xnum + 1) : (coords[0] * xnum + xleftover + 1);
            int ystart = (coords[1] < yleftover) ? (coords[1] * ynum + 1) : (coords[1] * ynum + yleftover + 1);
            printf("k=%d,xnum:%d,ynum:%d,xstart=%d,ystart=%d\n", k, xnum, ynum, xstart, ystart);
            for (int i = 0; i < xnum; i++) {
                MPI_Recv(&a[xstart + i][ystart], ynum, MPI_DOUBLE, k, k, comm, MPI_STATUS_IGNORE);
            }
        }

        //set the boundary
        //right
        for (int j = 1; j < nx + 1; j++) {
            double y = 1.0 * j / ((double)(1 + nx));
            a[nx + 1][j] = y / (y * y + 4.0);
        }
        //left
        for (int j = 1; j < nx + 1; j++) {
            double y = 1.0 * j / ((double)(nx + 1));
            a[0][j] = y / (1.0 + y * y);
        }
        //top
        for (int i = 1; i < nx + 1; i++) {
            double x = 1.0 * i / ((double)(nx + 1));
            a[i][nx + 1] = 1 / (1 + (1 + x) * (1 + x));
        }

    }
}
