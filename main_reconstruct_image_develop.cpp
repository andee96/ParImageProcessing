#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include "comms_lib.h"
#include "mpi.h"
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"

using namespace std; 

RealNumber boundaryval(int i, int m);
void set_sawtooth_values(int N_local, int M_local, int N, RealNumber** arr, int* coords);
void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN);
void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN, int &ni, int &nj);

int main(int argc, char **argv)
{
    int rank, size;
    int M, N; 
    double pM, pN, delta_max; 
    char filename[100]; 
    char destname[100];
    int coords[2] = {0,0};
    int auto_flag = atoi(argv[2]);
    int ni=0, nj=0;
    // Read parameters according to auto flag
    if (auto_flag == 0)
    {
        read_params(argv[1], filename, destname, delta_max, pM, pN);
    }
    else if (auto_flag == 1)
    {
        read_params(argv[1], filename, destname, delta_max, pM, pN, ni, nj);
    }
    int M_local, N_local; 

    // Initialize according to the value of the flag 
    if (auto_flag == 0){ initialize(rank, size, coords);}
    else if (auto_flag == 1){ initialize(rank, size, coords, ni, nj);}
    if (rank == 0 ){cout << "Reconstructig original image from edge file: " << filename << endl;}
    sprintf(destname, "%s_%d_%2f.pgm", destname, size, delta_max);
    // Read the local buffer
    if (rank == 0 ){cout << "Reading file..." << endl;}
    RealNumber** edge_local = read(filename, rank, M_local, N_local, M, N, pM, pN); 
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

void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN)
{
    // Reads the parameters from the specified file 
    string line; 
    ifstream file(filename);

    // Read parameters line by line 
    getline(file, line);
    strcpy(file_in, line.c_str());
    getline(file, line);
    strcpy(file_out, line.c_str());
    getline(file, line);
    delta_max = atof(line.c_str());
    getline(file, line);
    pM = atof(line.c_str());
    getline(file, line);
    pN = atof(line.c_str());
}

void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN, int &ni, int &nj)
{
    // Reads the parameters from the specified file 
    string line; 
    ifstream file(filename);

    // Read parameters line by line 
    getline(file, line);
    strcpy(file_in, line.c_str());
    getline(file, line);
    strcpy(file_out, line.c_str());
    getline(file, line);
    delta_max = atof(line.c_str());
    getline(file, line);
    pM = atof(line.c_str());
    getline(file, line);
    pN = atof(line.c_str());
    getline(file, line);
    ni = atoi(line.c_str());
    getline(file, line);
    nj = atoi(line.c_str()); 
}

