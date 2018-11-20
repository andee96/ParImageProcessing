#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
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
    double t1, t2, ttotal = 0;
    int n_runs;
    if (argc == 4){n_runs = atoi(argv[3]);}
    else{n_runs = 1;}    

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
    int** proc_dims = create2darray<int>(size, 2, 0);
    if (rank == 0 ){cout << "Reconstructig original image from edge file: " << filename << endl;}
    sprintf(destname, "%s_%d_%2f.pgm", destname, size, delta_max);
    // Read the local buffer
    for (int n = 0; n < n_runs; n++)
    {
        if (rank == 0 ){cout << "Run " << n << endl;}
        // add function for timing in comms library
        RealNumber** edge_local = read(filename, rank, M_local, N_local, M, N, pM, pN, size, proc_dims); 
        // cout <<"Rank " << rank << ": Local sizes: " << M_local << ", " << N_local << endl;
        // if (rank == 0 ){cout << "Done." << endl;}
        RealNumber** old = create2darray<RealNumber>(M_local+2, N_local+2, 255);
        RealNumber** new_arr = create2darray<RealNumber>(M_local+2, N_local+2, 255);
        // Set sawtooth values for old 
        set_sawtooth_values(N_local, M_local, N, old, coords);
        // if (rank == 0 ){cout << "Reconstructing..." << endl;}
        t1 = record_time();
        reconstruct(delta_max, old, new_arr, edge_local, M_local, N_local, rank);
        t2 = record_time();
        ttotal = ttotal + t2 - t1; 
        // if (rank == 0 ){cout << "Writing file to: " << destname << endl;}
        delete[] edge_local; delete[] old; 
        write(destname, rank, new_arr, proc_dims, size); 
    }
    if (rank == 0 )
    {
        cout << "Done." << endl;
        cout << "Total time spent reconstructing: " << ttotal << "s." << endl;
        cout <<"Time per reconstruction: " << ttotal/n_runs << "s." << endl; 
    }
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

    file.close();
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

    file.close();
}

// void write_log_file(char* filename, char* file_in, char* file_out, double ttotal, int n_runs)
// {
//     return 0;
// }

