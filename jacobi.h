void sweep1d(double a[][maxn], double f[][maxn], int nx,int s, int e, double b[][maxn]);
void sweep2d(double a[][maxn], double f[][maxn], int nx,int sx, int ex, int sy, int ey, double b[][maxn]);
void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm,int nbrleft, int nbrright);
void exchange3_2d(double x[][maxn], int nx, int sx, int ex, int sy, int ey, MPI_Comm comm, int nbrleft, int nbrright, int nbrup, int nbrdown);
void exchange_nb_2d(double x[][maxn], int nx, int sx, int ex, int sy, int ey, MPI_Comm comm,int nbrleft, int nbrright, int nbrup,int nbrdown);
void exchange_rma_fence_2d(double a[][maxn], MPI_Win win,int nx, int sx, int ex, int sy, int ey, int nbrleft, int nbrright, int nbrup, int nbrdown);
void exchange_rma_pscw_2d(double a[][maxn], MPI_Win win, MPI_Group ranks[4],int coords[2],int nx,int sx,int ex,int sy,int ey,int nbrleft,int nbrright,int nbrup,int nbrdown);
double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e);
double griddiff_2d(double a[][maxn], double b[][maxn], int nx, int sx, int ex,int sy,int ey);

