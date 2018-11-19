/* Communication library for parallel image reconstruction 
Date: 17/11/2018 
Author: Andreas Malekos 
Written for course: Message Passing Programming 
*/ 

#include "mpi.h"
#include <iostream> 
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"

using namespace std;

// Declare global variables 
MPI_Comm comm; 
int ranks[4]; 
int dim[2] = {0,0};
int M, N; 

// Helper function declarations 
void loop(double** old, double** new_array, double** edge, int i_init, int i_end, int j_init, int j_end);
void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int &rank, int* coords);
void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int &rank, int* coords, int ni, int nj);
void scatter(int rank, RealNumber **buf, RealNumber **local_buf, int M_local, int N_local);
void gather(int rank,  RealNumber** local_buf, RealNumber** buf, int M_local, int N_local);
void decompose(int &M_local, int &N_local, double pM, double pN, int rank);
void decompose_single_dim(int &A_local, int A, double pA, int coord, int P);
void define_datatypes(MPI_Datatype &MPI_row, MPI_Datatype &MPI_col, int M, int N);

// Main functions 

void initialize(int &rank, int &size, int* coords)
{
    MPI_Comm comm_old, comm_cart; 
    // Initializes MPI, returns rank, size and coordinates 
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    comm_old = MPI_COMM_WORLD;
    create_topology(size, comm, comm_old, rank, coords);
}

void initialize(int &rank, int &size, int* coords, int ni, int nj)
{
    // Same as above, but the length along each dimension is supplied beforehand
    MPI_Comm comm_old, comm_cart; 
    // Initializes MPI, returns rank, size and coordinates 
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    comm_old = MPI_COMM_WORLD;
    create_topology(size, comm, comm_old, rank, coords, ni, nj);
}

RealNumber** read(char* filename, int rank, int &M_local, int &N_local, int &M_master, int &N_master, double pM, double pN)
{
    // Returns the local buffer and the dimensions of the "unhalo-ed" local file
    RealNumber **masterbuf; 
    pgmsize(filename, &M, &N);
    M_master = M; 
    N_master = N;
    decompose(M_local, N_local, pM, pN, rank); 
    cout << M_local << " " << N_local << " " << rank << endl;
    RealNumber** buf = create2darray<RealNumber>(M_local+2, N_local+2, 255.0);
    // Rank 0 reads and scatters the process
    if (rank == 0)
    {
        masterbuf = create2darray<RealNumber>(M, N, 0.0);
        pgmread(filename, &masterbuf[0][0], M, N);
    }
    scatter(rank, masterbuf, buf, M_local, N_local); 
    if (rank == 0){delete[] masterbuf;}
    return buf; 
}

void write(char* filename, int rank, RealNumber** local_buf, int M_local, int N_local)
{
    RealNumber **buf = create2darray<RealNumber>(M, N, 0.0);
    gather(rank, local_buf, buf, M_local, N_local);

    // Rank 0 writes the file to filename
    if (rank == 0)
    {
        pgmwrite(filename, &buf[0][0], M, N);
    }
    delete[] buf; 
}

void reconstruct(double delta_max,  RealNumber** old, RealNumber** new_array, RealNumber** edge, int M_local, int N_local, int rank)
{
    MPI_Request request;
    MPI_Status status; 
    // Get the types for halo swaping 
    MPI_Datatype MPI_REAL_ROW; 
    MPI_Datatype MPI_REAL_COL; 
    define_datatypes(MPI_REAL_ROW, MPI_REAL_COL, M_local, N_local); 

    int rank_up = ranks[0];
    int rank_down = ranks[1];
    int rank_right = ranks[2];
    int rank_left = ranks[3];

    int M1 = (int)((M_local+2)/2);
    int N1 = (int)((N_local+2)/4);
    int N2 = 2*N1; 
    
    double** delta = create2darray<double>(M_local, N_local, 0.0); 
    double** max_delta = create2darray<double>(M_local, N_local, 0.0); 
    max_delta[0][0] = 100;
    int n = 0; 
    int n_comp_delta = 100; // compute delta at every n_comp_delta iterations; 

    //Loop reconstruct the image in buffer by using the old image
    while (max_delta[0][0] > delta_max)
    {
        // Receive row 0 from rank_up and send row M_local to rank_down
        MPI_Issend( &old[M_local][0], 1, MPI_REAL_ROW, rank_down, 0, comm, &request);
        MPI_Recv( &old[0][0], 1, MPI_REAL_ROW, rank_up, 0, comm, &status);
        // Loop over area 1
        loop(old, new_array, edge, 1, M1, N1, N2);
        MPI_Wait(&request, &status);
        // Receive row M_local+1 from rank_down and send row 1 to rank_up
        MPI_Issend( &old[1][0], 1, MPI_REAL_ROW, rank_up, 0, comm, &request);
        MPI_Recv( &old[M_local+1][0], 1, MPI_REAL_ROW, rank_down, 0, comm, &status);
        // Loop over area 2
        loop(old, new_array, edge, M1, M_local+1, N1, N2);
        MPI_Wait(&request, &status);
        // Send column 1 to rank_left and receive column N_local+1 to rank_right
        MPI_Issend( &old[0][1], 1, MPI_REAL_COL, rank_left, 0, comm, &request);
        MPI_Recv ( &old[0][N_local+1], 1, MPI_REAL_COL, rank_right, 0, comm, &status); 
        // Loop over area 3
        loop(old, new_array, edge, 1, M_local+1, N2, N_local+1);
        MPI_Wait(&request, &status);
        // Receive column 0 from rank_left and send column N_local to rank_right
        MPI_Issend( &old[0][N_local], 1, MPI_REAL_COL, rank_right, 0, comm, &request);
        MPI_Recv( &old[0][0], 1, MPI_REAL_COL, rank_left, 0, comm, &status);
        // Loop over area 4
        loop(old, new_array, edge, 1, M_local+1, 1, N1);
        MPI_Wait(&request, &status);
        
        // At the end of each iteration, set the old array equal to the new array
        // old = new_array; is there a better way to do this??? 
        // cout << "Rank " << rank << ": Run " << n << " finished." << endl;
        for (int i = 1 ; i < M_local+1; i++)
        {
            for (int j = 1; j < N_local+1; j++)
            {
                // Calculate absolute difference between old and new 
                if (n%n_comp_delta == 0)
                {

                    delta[i-1][j-1] = abs(old[i][j] - new_array[i][j]);
                }
            old[i][j] = new_array[i][j];
            }
        }
        if (n%n_comp_delta == 0)
        {
            // Reduce max 
            MPI_Allreduce(&delta[0][0], &max_delta[0][0], M_local*N_local, MPI_DOUBLE, MPI_MAX, comm);
            // cout << max_delta[0][0] << endl;
         }
        n++;
    }
    delete[] max_delta; 
    delete[] delta; 
    if (rank == 0)
    {
        cout << "Reconstruction done, number of iterations: " << n << endl;
    }
}

void finalize()
{
    MPI_Finalize();
}
// Helper functions 

void loop(double** old, double** new_array, double** edge, int i_init, int i_end, int j_init, int j_end)
{
    // Perform image reconstruction by calculating the new array elements for i = [i_init, i_end)
    // j = [j_init, j_end).
    for (int i = i_init; i<i_end; i++)
    {
        for (int j = j_init; j<j_end; j++)
        {
            new_array[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
        }
    } 
}

void scatter(int rank, RealNumber **buf, RealNumber **local_buf, int M_local, int N_local)
{
    // Scatters buf among all processes 
    int ni = dim[0]; int nj = dim[1]; 
    int send_rank_coords[2] = {0,0}; 
    int send_rank; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_MASTER;
    MPI_Type_vector(M_local, N_local, N, MPI_REALNUMBER, &SUB_ARRAY_MASTER);  
    MPI_Type_commit(&SUB_ARRAY_MASTER);

    MPI_Datatype SUB_ARRAY_LOCAL; 
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_REALNUMBER, &SUB_ARRAY_LOCAL);
    MPI_Type_commit(&SUB_ARRAY_LOCAL);
    if (rank == 0)
    {
        // Asyncrhonous send and receive to itself 
        MPI_Issend(&buf[0][0], 1, SUB_ARRAY_MASTER, 0, 0, comm, &request);
        MPI_Recv( &local_buf[1][1], 1, SUB_ARRAY_LOCAL, 0, 0, comm, &status);
        MPI_Wait(&request, &status);
       for (int n = 1; n < ni*nj; n++)
       {
           send_rank = n; 
           MPI_Cart_coords(comm, send_rank, 2, send_rank_coords); 
           MPI_Ssend(&buf[send_rank_coords[0] * M_local][send_rank_coords[1]*N_local],
                    1, SUB_ARRAY_MASTER, send_rank, 0, comm);
       }
    }
    else
    {
        MPI_Recv( &local_buf[1][1], 1, SUB_ARRAY_LOCAL, 0, 0, comm, &status);
    }
}

void gather(int rank,  RealNumber** local_buf, RealNumber** buf, int M_local, int N_local)
{
    // Sends all local bufs back to master buf 
    int ni = dim[0]; int nj = dim[1]; 
    int rec_coords[2] = {0,0}; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_MASTER;
    MPI_Type_vector(M_local, N_local, N, MPI_REALNUMBER, &SUB_ARRAY_MASTER);  
    MPI_Type_commit(&SUB_ARRAY_MASTER);

    MPI_Datatype SUB_ARRAY_LOCAL; 
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_REALNUMBER, &SUB_ARRAY_LOCAL);
    MPI_Type_commit(&SUB_ARRAY_LOCAL); 

    if (rank == 0)
    {
        // Non blocking send from rank 0 to itself
        MPI_Issend(&local_buf[1][1], 1, SUB_ARRAY_LOCAL, 0, 0, comm, &request);
        MPI_Recv( &buf[0][0], 1, SUB_ARRAY_MASTER, 0, 0, comm, &status);
        MPI_Wait(&request, &status); 
        // Blocking receive from each rank 
        for (int n_rank = 1; n_rank<ni*nj; n_rank++ )
        {
            MPI_Cart_coords(comm, n_rank, 2, rec_coords);
            MPI_Recv(&buf[rec_coords[0]*M_local][rec_coords[1]*N_local],
                    1, SUB_ARRAY_MASTER, n_rank, 0, comm, &status);
        }
    }
    else
    {
        // Non send data from all ranks to rank 0
        MPI_Ssend(&local_buf[1][1], 1, SUB_ARRAY_LOCAL, 0, 0, comm);
    }
}

void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int &rank, int* coords)
{
    // Given the number of processes and the size of the problem, create appropriate 2D (or 1D) topology
    // We want a longer dimension in the j direction, as that is how images are stored in memory. This means 
    // more processes along the i direction. 
    // Need to split n into two factors (out of which one can be 1). Biggest factor will be along j direction 
    int period[2] = {0, 1};
    int reorder = 1; 
    int rank_up, rank_down, rank_left, rank_right; 
    MPI_Dims_create(n, 2, dim);
    MPI_Cart_create(comm_old, 2, dim, period, reorder, &comm_new);
    MPI_Comm_rank(comm_new, &rank); 
    MPI_Cart_coords(comm_new, rank, 2, coords);
    MPI_Cart_shift(comm_new, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm_new, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left; 
}

void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int &rank, int* coords, int ni, int nj)
{
    // Same as create topology above, but MPI_Dims_create is not used, and the nubmer of dimensions is suppled by 
    // the user. It is the user's responsability to make sure that these are defined correctly.
    int period[2] = {0, 1};
    int reorder = 1; 
    int rank_up, rank_down, rank_left, rank_right; 
    dim[0] = ni; 
    dim[1] = nj; 
    MPI_Cart_create(comm_old, 2, dim, period, reorder, &comm_new);
    MPI_Comm_rank(comm_new, &rank); 
    MPI_Cart_coords(comm_new, rank, 2, coords);
    MPI_Cart_shift(comm_new, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm_new, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left;  

}


void decompose(int &M_local, int &N_local, double pM, double pN, int rank)
{
    // Decomposes the image. For now assume the number of processes divides N
    // and M evenly 
    // M_local = M/dim[0];
    // N_local = N/dim[1];

    int coords[2] = {0,0};
    MPI_Cart_coords(comm, rank, 2, coords); 
    decompose_single_dim(M_local, M, pM, coords[0], dim[0]);
    decompose_single_dim(N_local, N, pN, coords[1], dim[1]);

}

void decompose_single_dim(int &A_local, int A, double pA, int coord, int P)
{
    // Decomposes A into A_local, given the edge/inside fraction for that 
    // dimension pA. 
    int rem = 0;
    if (P > 2)
    {
        int M_edge = int(pA * A/(2*pA + P - 2));
        int M_inner = int(A - 2*M_edge)/(P -2);

        rem = (2*M_edge + (P-2)*M_inner);
        // Check to see if process is in the corners 
        if (coord == 0 || coord == P-1){A_local = M_edge;}
        else{A_local = M_inner;}
        if (coord < rem){A_local++;}
    }
    else
    {
        // Case where we have 2 or 1 corners, i.e. no middle processes 
        A_local = int(A/P);
        rem = A%(P);

        if (coord < rem)
        {
            A_local++;
        } 
    }
    
}

void define_datatypes(MPI_Datatype &MPI_row, MPI_Datatype &MPI_col, int M, int N)
{
    // Create the derived datatypes 
    MPI_Type_vector(1, N+2, 1, MPI_REALNUMBER, &MPI_row);
    MPI_Type_commit(&MPI_row);

    MPI_Type_vector(M+2, 1, N+2, MPI_REALNUMBER, &MPI_col);
    MPI_Type_commit(&MPI_col);
}