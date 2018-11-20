#ifndef ParImageProcessor_H
#define ParImageProcessor_H

#include "mpi.h"

class ParImageProcessor
{
    // Private attributes
    MPI_Comm comm, comm_old; 

    // Private member functions
    void loop(RealNumber** old, RealNumber** new_array, RealNumber** edge, int i_init, int i_end, int j_init, int j_end);
    void create_topology();
    void create_topology(int ni, int nj);
    void scatter(RealNumber **local_buf, RealNumber **buf);
    void gather(RealNumber** local_buf, RealNumber** buf);
    void decompose();
    void decompose_single_dim(int &A_local, int A, double pA, int coord, int P);
    void get_proc_indices();
 
public:
    // Public attributes 
    int size, rank; 
    int dim[2] = {0,0};
    int M, N; 
    int M_local, N_local; 
    int coords[2];
    int ranks[4]; 
    int** proc_dims;
    int** proc_indices;
    double pM, pN; 


    // Public methods
    void initialize(int PM, int PN);
    void initialize(int PM, int PN, int ni, int nj);
    RealNumber** read(char* filename);
    void reconstruct(double delta_max, RealNumber** old, RealNumber** new_array, RealNumber** edge);
    void write(char* filename, RealNumber** local_buf);
    void finalize();
    double record_time();

    // Class constructors
    ParImageProcessor();
};



#endif