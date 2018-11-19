#include "mpi.h"
#include "precision.h"

void initialize(int &rank, int &size, int* coords);
void initialize(int &rank, int &size, int* coords, int ni, int nj);
void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int &rank, int* coords);
void scatter(int rank, RealNumber **buf, RealNumber **local_buf, int M_local, int N_local);
RealNumber** read(char* filename, int rank, int &M_local, int &N_local, int &M_master, int &N_master, double pM, double pN);
void gather(int rank,  RealNumber** local_buf, RealNumber** buf, int M_local, int N_local);
void write(char* filename, int rank, RealNumber** local_buf, int M_local, int N_local);
void reconstruct(double delta_max,  RealNumber** old, RealNumber** new_array, RealNumber** edge, int M_local, int N_local, int rank);
void finalize();





