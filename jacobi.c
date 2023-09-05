#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "poisson.h"
#include "jacobi.h"



void sweep1d(double a[][maxn], double f[][maxn], int nx,
    int s, int e, double b[][maxn])
{
    double h;
    int i, j;

    h = 1.0 / ((double)(nx + 1));

    for (i = s; i <= e; i++) {
        for (j = 1; j < nx + 1; j++) {
            b[i][j] = 0.25 * (a[i - 1][j] + a[i + 1][j] + a[i][j + 1] + a[i][j - 1] - h * h * f[i][j]);
        }
    }
}

void sweep2d(double a[][maxn], double f[][maxn], int nx,
    int sx, int ex, int sy, int ey, double b[][maxn])
{
    double h;
    int i, j;

    h = 1.0 / ((double)(nx + 1));

    for (i = sx; i <= ex; i++) {
        for (j = sy; j <= ey; j++) {
            b[i][j] = 0.25 * (a[i - 1][j] + a[i + 1][j] + a[i][j + 1] + a[i][j - 1] - h * h * f[i][j]);
        }
    }
}

void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
    int nbrleft, int nbrright)
{

    MPI_Sendrecv(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, &x[s - 1][1], nx, MPI_DOUBLE, nbrleft,
        0, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, &x[e + 1][1], nx, MPI_DOUBLE, nbrright,
        1, comm, MPI_STATUS_IGNORE);

}

void exchange3_2d(double x[][maxn], int nx, int sx, int ex, int sy, int ey, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown)
{
    
    MPI_Sendrecv(&x[ex][sy],ey-sy+1 , MPI_DOUBLE, nbrright, 0, 
		&x[sx - 1][sy], ey - sy + 1, MPI_DOUBLE,nbrleft, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[sx][sy], ey - sy + 1, MPI_DOUBLE, nbrleft, 1, 
		&x[ex + 1][sy], ey - sy + 1, MPI_DOUBLE,nbrright, 1, comm, MPI_STATUS_IGNORE);

    MPI_Datatype doubletype;
    MPI_Type_vector(ex-sx+1,1 , nx+2, MPI_DOUBLE, &doubletype);
    MPI_Type_commit(&doubletype);
    MPI_Sendrecv(&x[sx][sy], 1, doubletype, nbrup, 3, &x[sx][ey + 1], 1, doubletype,
        nbrdown, 3, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&x[sx][ey], 1, doubletype, nbrdown, 2, &x[sx][sy - 1], 1, doubletype,
        nbrup, 2, comm, MPI_STATUS_IGNORE);
}

void exchange_nb_2d(double x[][maxn], int nx, int sx, int ex, int sy, int ey, MPI_Comm comm,int nbrleft, int nbrright, int nbrup,int nbrdown) {

    MPI_Request requests[8];
    
    MPI_Datatype doubletype;
    MPI_Type_vector(ex - sx + 1, 1, nx + 2, MPI_DOUBLE, &doubletype);
    MPI_Type_commit(&doubletype);

    MPI_Isend(&x[ex][sy], ey - sy + 1, MPI_DOUBLE, nbrright, 0, comm, &requests[0]);
    MPI_Isend(&x[sx][sy], ey - sy + 1, MPI_DOUBLE, nbrleft, 1, comm, &requests[1]);
    MPI_Isend(&x[sx][sy], 1, doubletype, nbrup, 3, comm, &requests[2]);
    MPI_Isend(&x[sx][ey], 1, doubletype, nbrdown, 2, comm, &requests[3]);


    MPI_Irecv(&x[sx - 1][sy], ey - sy + 1, MPI_DOUBLE, nbrleft, 0, comm, &requests[4]);
    MPI_Irecv(&x[ex + 1][sy], ey - sy + 1, MPI_DOUBLE, nbrright, 1, comm, &requests[5]);
    MPI_Irecv(&x[sx][ey + 1], 1, doubletype, nbrdown, 3, comm, &requests[6]);
    MPI_Irecv(&x[sx][sy - 1], 1, doubletype, nbrup, 2, comm, &requests[7]);
    MPI_Waitall(8, requests, MPI_STATUS_IGNORE);

  

}

void exchange_rma_fence_2d(double a[][maxn], MPI_Win win,int nx, int sx, int ex, int sy, int ey, int nbrleft, int nbrright, int nbrup, int nbrdown) {
    MPI_Datatype doubletype;
    MPI_Type_vector(ex - sx + 1, 1, nx + 2, MPI_DOUBLE, &doubletype);
    MPI_Type_commit(&doubletype);

    MPI_Aint offset;
    MPI_Win_fence(0, win);

    offset = maxn + sy;
    MPI_Get(&a[ex + 1][sy], ey - sy + 1, MPI_DOUBLE, nbrright, offset, ey - sy + 1, MPI_DOUBLE, win);
    offset = sy;
    MPI_Put(&a[ex][sy], ey - sy + 1, MPI_DOUBLE, nbrright, offset, ey - sy + 1, MPI_DOUBLE, win);

    offset = maxn + 1 + ey;
    MPI_Get(&a[sx][ey + 1], 1, doubletype, nbrdown, offset, 1, doubletype, win);
    offset = maxn + ey;
    MPI_Put(&a[sx][ey], 1, doubletype, nbrdown, offset, 1, doubletype, win);
    MPI_Win_fence(0, win);
}


void exchange_rma_pscw_2d(double a[][maxn], MPI_Win win, MPI_Group ranks[4],int coords[2],int nx,int sx,int ex,int sy,int ey,int nbrleft,int nbrright,int nbrup,int nbrdown) {
    MPI_Datatype doubletype;
    MPI_Type_vector(ex - sx + 1, 1, nx + 2, MPI_DOUBLE, &doubletype);
    MPI_Type_commit(&doubletype);

    MPI_Aint offset;
    if (coords[0] != 0) {//it has left neighbor
        MPI_Win_post(ranks[2], 0, win);//ready to receive data from left neighbor
    }
    if (coords[0] == 0) {
        MPI_Win_start(ranks[0], 0, win);//ready to start RMA operations
    }

    if (coords[0] == 0) {
        offset = maxn + sy;
        MPI_Get(&a[ex + 1][sy], ey - sy + 1, MPI_DOUBLE, nbrright, offset, ey - sy + 1, MPI_DOUBLE, win);
        offset = sy;
        MPI_Put(&a[ex][sy], ey - sy + 1, MPI_DOUBLE, nbrright, offset, ey - sy + 1, MPI_DOUBLE, win);
    }

    if (coords[0] != 0) {
        MPI_Win_wait(win);//blocks until all incoming RMA operations are completed
    }
    if (coords[0] == 0) {
        MPI_Win_complete(win);// finished its RMA operations
    }

    if (coords[1] == 0) {
        MPI_Win_post(ranks[3], 0, win);
    }
    if (coords[1] != 0) {
        MPI_Win_start(ranks[1], 0, win);
    }

    if (coords[1] != 0) {
        offset = maxn + 1 + sy;
        MPI_Get(&a[sx][ey + 1], 1, doubletype, nbrdown, offset, 1, doubletype, win);
        offset = maxn + ey;
        MPI_Put(&a[sx][ey], 1, doubletype, nbrdown, offset, 1, doubletype, win);
    }

    if (coords[1] == 0) {
        MPI_Win_wait(win);
    }
    if (coords[1] != 0) {
        MPI_Win_complete(win);
    }
}
               

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e)
{
    double sum;
    double tmp;
    int i, j;

    sum = 0.0;

    for (i = s; i <= e; i++) {
        for (j = 1; j < nx + 1; j++) {
            tmp = (a[i][j] - b[i][j]);
            sum = sum + tmp * tmp;
        }
    }

    return sum;
}


double griddiff_2d(double a[][maxn], double b[][maxn], int nx, int sx, int ex,
    int sy,int ey)
{
    double sum;
    double tmp;
    int i, j;

    sum = 0.0;

    for (i = sx; i <= ex; i++) {
        for (j = sy; j <= ey; j++) {
            tmp = (a[i][j] - b[i][j]);
            sum = sum + tmp * tmp;
        }
    }

    return sum;

}
