#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
// #include "comms_lib.h"
#include "mpi.h"
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"
#include "ParImageProcessor.h"

using namespace std; 

RealNumber boundaryval(int i, int m);
void set_sawtooth_values(int N_local, int M_local, int N, RealNumber** arr, int start_index);
void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN);
void read_params(char* filename, char* file_in, char* file_out, double &delta_max, double &pM, double &pN, int &ni, int &nj);
void write_log_file(char* filename, char* file_in_name, char* file_out_name, double ttotal, int n_runs, int delta_max, ParImageProcessor &processor);


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
    // Timing calls 
    double t1_reconstruct, t2_reconstruct, ttotal_reconstruct=0;
    double t1_read, t2_read, ttotal_read=0;
    double t1_write, t2_write, ttotal_write=0;
    double ttotal = 0;
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
    ParImageProcessor processor;
    if (auto_flag == 0){processor.initialize(pM, pN);}
    else if (auto_flag == 1){processor.initialize(pM, pN, ni, nj);}
    // Assign coordinates and rank
    coords[0] = processor.coords[0];
    coords[1] = processor.coords[1];
    rank = processor.rank;
    size = processor.size;
    int** proc_dims = create2darray<int>(size, 2, 0);
    sprintf(destname, "%s_%d_%2f.pgm", destname, size, delta_max);
    // Read the local buffer
    for (int n = 0; n < n_runs; n++)
    {
        //if (processor.rank == 0 ){cout << "Run " << n << endl;}
        // add function for timing in comms library
        t1_read = processor.record_time();
        RealNumber** edge_local = processor.read(filename); 
        t2_read = processor.record_time();
        // Assign local and global sizes
        M_local = processor.M_local; 
        N_local = processor.N_local;
        M = processor.M;
        N = processor.N; 
        RealNumber** old = create2darray<RealNumber>(M_local+2, N_local+2, 255);
        RealNumber** new_arr = create2darray<RealNumber>(M_local+2, N_local+2, 255);
        // Set sawtooth values for old 
        set_sawtooth_values(N_local, M_local, N, old, processor.proc_indices[rank][1]);
        t1_reconstruct = processor.record_time();
        processor.reconstruct(delta_max, old, new_arr, edge_local);
        t2_reconstruct = processor.record_time();
        delete[] edge_local; delete[] old; 
        t1_write = processor.record_time();
        processor.write(destname, new_arr);
        t2_write = processor.record_time();
        // Record times 
        ttotal_read =  ttotal_read + t2_read - t1_read;
        ttotal_reconstruct = ttotal_reconstruct + t2_reconstruct - t1_reconstruct;
        ttotal_write = ttotal_write + t2_write - t1_write;
    }
	ttotal = ttotal_read + ttotal_reconstruct + ttotal_write; 
    if (rank == 0 )
    {
        cout << "Done." << endl;
        cout << "Total time spent reconstructing: " << ttotal_reconstruct << "s." << endl;
        cout <<"\tTime per reconstruction: " << ttotal_reconstruct/n_runs << "s." << endl; 
        cout << "Total time spent reading: " << ttotal_read << "s." << endl;
        cout << "\tTime per read: " << ttotal_read/n_runs << "s." << endl;
        cout << "Total time spent writing: " << ttotal_write << "s." << endl;
        cout << "\tTime per write: " << ttotal_write/n_runs << "s." << endl;
        cout << "Total runtime: " << ttotal << "s." << endl;
        cout <<"\tTotal time per run: " << ttotal/n_runs << "s." << endl;
         // Write logfile
        char log_name[100];
        sprintf(log_name, "./test_log.txt");
        write_log_file(log_name, filename, destname, ttotal, n_runs, delta_max, processor);
    }
    processor.finalize();
}

RealNumber boundaryval(int i, int m)
{
    // Sets the boundary value 
    RealNumber val; 
    val = 2.0*((RealNumber)(i-1))/((RealNumber)(m-1));
    if (i >= m/2+1) {val = 2.0 - val;};
    return val;
}

void set_sawtooth_values(int N_local, int M_local, int N, RealNumber** arr, int start_index)
{
    // Set sawtooth values for arr 
    double j_local, val;
    for (int j = 1; j <N_local+1; j++)
    {
        j_local =  int(j + start_index);
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

void write_log_file(char* filename, char* file_in_name, char* file_out_name, double ttotal, int n_runs, int delta_max, ParImageProcessor &processor)
{
    // Open file 
    ofstream file;
    file.open(filename, ios::out);
    file.open(filename, ios::out|ios::in|ios::binary);

    // Get date 
    file << "Date: " << "\n";
    file << endl;
    // Input file details
    file << "Input file: " << file_in_name << endl;
    file <<"\tM: " << processor.M << endl;
    file <<"\tN: " << processor.N << endl;
    file <<"Output File: " << file_out_name << endl;
    file << endl;
    file << "Parameters:\n";
    file << "\tNumber of cores: " << processor.size << endl; 
    file << "\tMaximum delta: " << delta_max << endl;
    file <<"\tCartesian dimensions: " << endl;
    file << endl;
    file << "Individual process parameters:\n";
    for (int r = 0; r<processor.size; r++)
    {
        file << "\tProcess " << r <<":\n";
        file << "\t\tCoordinates: " ;
        file <<"\t\tM_local: " << processor.proc_dims[r][0] << endl;
        file << "\t\tN_local: " << processor.proc_dims[r][1] << endl;
        file << endl;
    }

    file <<"Run details:\n";
    file << "\tTotal runtime: " << ttotal <<"s.\n";
    file << "\tTime per individual run: " << ttotal/n_runs <<"s.\n";

    file.close();
}


