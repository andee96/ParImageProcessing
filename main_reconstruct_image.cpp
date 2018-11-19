#include <iostream>
#include "comms_lib.h"
#include "mpi.h"
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"

using namespace std; 

RealNumber boundaryval(int i, int m);
void set_sawtooth_values(int N_local, int M_local, int N, RealNumber** arr, int* coords);

int main(int argc, char **argv)
{
    int rank, size;
    int M, N; 
    int coords[2] = {0,0};
    // Define source and destination names 
    char filename[100]; 
    sprintf(filename, "%s", argv[1]); 
    char destname[100];
    double delta_max = atof(argv[3]);
    int M_local, N_local; 

    // Initialize 
    initialize(rank, size, coords); 
    if (rank == 0 ){cout << "Reconstructig original image from edge file: " << filename << endl;}
    sprintf(destname, "%s_%d_%2f.pgm", argv[2], size, delta_max);
    // Read the local buffer
    if (rank == 0 ){cout << "Reading file..." << endl;}
    RealNumber** edge_local = read(filename, rank, M_local, N_local, M, N); 
    if (rank == 0 ){cout << "Done." << endl;}
    RealNumber** old = create2darray<RealNumber>(M_local+2, N_local+2, 255);
    RealNumber** new_arr = create2darray<RealNumber>(M_local+2, N_local+2, 255);
    // Set sawtooth values for old 
    set_sawtooth_values(N_local, M_local, N, old, coords);
    if (rank == 0 ){cout << "Reconstructing..." << endl;}
    reconstruct(delta_max, old, new_arr, edge_local, M_local, N_local, rank);
    if (rank == 0 ){cout << "Writing file to: " << destname << endl;}
    write(destname, rank, new_arr, M_local, N_local); 
    if (rank == 0 ){cout << "Done." << endl;}
    delete[] edge_local;
    finalize();
}

RealNumber boundaryval(int i, int m)
{
    // Sets the boundary value 
    RealNumber val; 
    val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
    if (i >= m/2+1) {val = 2.0 - val;};
    return val;
}

void set_sawtooth_values(int N_local, int M_local, int N, RealNumber** arr, int* coords)
{
    // Set sawtooth values for arr 
    double j_local, val;
    for (int j = 1; j <N_local+1; j++)
    {
        j_local =  int(j + N_local * coords[1]);
        val = boundaryval(j_local, N);
        arr[0][j] = int(255*(1.0-val));
        arr[M_local+1][j] = int(255*val);
    }
}