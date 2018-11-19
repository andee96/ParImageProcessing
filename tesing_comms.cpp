#include <iostream>
#include "comms_lib.h"
#include "mpi.h"
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"

using namespace std; 

RealNumber boundaryval(int i, int m);


int main(int argv, char **argc)
{
    int rank, size;
    int M, N; 
    int coords[2] = {0,0};
    char filename[100]; 
    sprintf(filename, "%s", argc[1]); 

    int M_local, N_local; 
    initialize(rank, size, coords); 

    cout << rank << " " << size << " " << coords[0] << " " << coords[1] << endl;

    RealNumber** buf = read(filename, rank, M_local, N_local, M, N); 

    char destname[100];
    sprintf(destname, "testcomms_%d.pgm", rank);
    double delta_max = 0.1; 
    RealNumber** old = create2darray<RealNumber>(M_local+2, N_local+2, 255);
    RealNumber** new_arr = create2darray<RealNumber>(M_local+2, N_local+2, 255);
    // Set sawtooth values for old 
    double j_local, val;
    for (int j = 1; j <N_local+1; j++)
    {
        j_local =  int(j + N_local * coords[1]);
        val = boundaryval(j_local, N);
        old[0][j] = int(255*(1.0-val));
        old[M_local+1][j] = int(255*val);
    }
    reconstruct(delta_max, old, new_arr, buf, M_local, N_local, rank);
    write(destname, rank, new_arr, M_local, N_local); 
    delete[] buf;

    
    MPI_Finalize();
}

RealNumber boundaryval(int i, int m)
{
    // Sets the boundary value 
    RealNumber val; 

    val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
    if (i >= m/2+1) {val = 2.0 - val;};
    return val;
}
