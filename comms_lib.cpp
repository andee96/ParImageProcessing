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
void scatter(int rank, RealNumber **buf, RealNumber **local_buf, int** proc_dims, int size);
void gather(int rank,  RealNumber** local_buf, RealNumber** buf, int** proc_dims, int size);
void decompose(int &M_local, int &N_local, double pM, double pN, int rank, int size, int** proc_dims);
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

RealNumber** read(char* filename, int rank, int &M_local, int &N_local, int &M_master, int &N_master, double pM, double pN, int size, int** proc_dims)
{
    // Returns the local buffer and the dimensions of the "unhalo-ed" local file
    RealNumber **masterbuf; 
    pgmsize(filename, &M, &N);
    M_master = M; 
    N_master = N;
    decompose(M_local, N_local, pM, pN, rank, size, proc_dims); 
    RealNumber** buf = create2darray<RealNumber>(M_local+2, N_local+2, 255.0);
    // Rank 0 reads and scatters the process
    if (rank == 0)
    {
        masterbuf = create2darray<RealNumber>(M, N, 0.0);
        pgmread(filename, &masterbuf[0][0], M, N);
    }
    scatter(rank, masterbuf, buf, proc_dims, size); 
    if (rank == 0){delete[] masterbuf;}
    return buf; 
}

void write(char* filename, int rank, RealNumber** local_buf, int** proc_dims,int size)
{
    RealNumber **buf = create2darray<RealNumber>(M, N, 0.0);
    gather(rank, local_buf, buf, proc_dims, size);

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
    
    RealNumber** delta = create2darray<RealNumber>(M_local, N_local, 0.0); 
    RealNumber max_delta_local = 100;
    RealNumber max_delta_global = 100;

    int n = 0; 
    int n_comp_delta = 100; // compute delta at every n_comp_delta iterations; 

    //Loop reconstruct the image in buffer by using the old image
    while (max_delta_global > delta_max)
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
            max_delta_local = max_elem(delta, M_local, N_local);
            MPI_Allreduce(&max_delta_local, &max_delta_global, 1, MPI_REALNUMBER, MPI_MAX, comm);
         }
        n++;
    }
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

void scatter(int rank, RealNumber **buf, RealNumber **local_buf, int** proc_dims, int size)
{
    // Scatters buf among all processes 
    int my_localM = proc_dims[rank][0];
    int my_localN = proc_dims[rank][1]; 
    int curr_coords[2] = {0,0};
    int currM = my_localM, currN = my_localN;
    int index_i = 0, index_j = 0; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_SEND;
    MPI_Datatype SUB_ARRAY_RECV;
    MPI_Type_vector(my_localM, my_localN, my_localN+2, MPI_REALNUMBER, &SUB_ARRAY_RECV); 
    MPI_Type_commit(&SUB_ARRAY_RECV);

    if (rank == 0 )
    {
        // Define datatype 
        MPI_Type_vector(proc_dims[0][0], proc_dims[0][1], N, MPI_REALNUMBER, &SUB_ARRAY_SEND);
        MPI_Type_commit(&SUB_ARRAY_SEND);
        // Asyncrhonous send to itself 
        MPI_Issend(&buf[0][0], 1, SUB_ARRAY_SEND, 0, 0, comm, &request);
        MPI_Recv( &local_buf[1][1], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
        MPI_Type_free(&SUB_ARRAY_SEND);
        for (int r = 1; r < size; r++)
        {
            // Update the indices
            if ((r)%dim[1] != 0)
            {
                index_j = index_j + currN; 
            }
            else
            {
                index_j = 0; 
                index_i = index_i + currM;
            }
            // Send the right array to every element 
            currM = proc_dims[r][0];
            currN = proc_dims[r][1];
            MPI_Type_vector(currM, currN, N, MPI_REALNUMBER, &SUB_ARRAY_SEND);
            MPI_Type_commit(&SUB_ARRAY_SEND);
            MPI_Wait(&request, &status);
            MPI_Issend(&buf[index_i][index_j],1, SUB_ARRAY_SEND, r, 0, comm, &request);
            MPI_Type_free(&SUB_ARRAY_SEND);
        }
        MPI_Wait(&request, &status);
    }
    else
    {
        MPI_Recv(&local_buf[1][1], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
    }
}

void gather(int rank,  RealNumber** local_buf, RealNumber** buf, int** proc_dims, int size)
{   
    // Sends all local bufs back to master buf 
    int ni = dim[0]; int nj = dim[1]; 
    int my_localM = proc_dims[rank][0];
    int my_localN = proc_dims[rank][1];
    int currM = my_localM,  currN=my_localN;
    int index_i = 0, index_j=0; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_SEND, SUB_ARRAY_RECV;
    // Custom datatype for sending arrays
    MPI_Type_vector(my_localM, my_localN, my_localN+2, MPI_REALNUMBER, &SUB_ARRAY_SEND);  
    MPI_Type_commit(&SUB_ARRAY_SEND);

    if (rank == 0)
    {
        MPI_Type_vector(my_localM, my_localN, N, MPI_REALNUMBER, &SUB_ARRAY_RECV);
        MPI_Type_commit(&SUB_ARRAY_RECV);
        // Non blocking send from rank 0 to itself
        MPI_Issend(&local_buf[1][1], 1, SUB_ARRAY_SEND, 0, 0, comm, &request);
        MPI_Recv( &buf[0][0], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
        MPI_Type_free(&SUB_ARRAY_RECV);
        // Blocking receive from each rank 
        for (int r = 1; r<size; r++ )
        {
            // Update indices
            if ((r)%dim[1] != 0)
            {
                cout << index_j << endl;
                index_j = index_j + currN; 
                cout << index_j << endl;
            }
            else
            {
                cout << index_i << endl;
                index_j = 0; 
                index_i = index_i + currM;
                cout << index_i << endl;
            }
            currM = proc_dims[r][0];
            currN = proc_dims[r][1];

            MPI_Type_vector(currM, currN, N, MPI_REALNUMBER, &SUB_ARRAY_RECV);
            MPI_Type_commit(&SUB_ARRAY_RECV);
            MPI_Wait(&request, &status);
            MPI_Irecv(&buf[index_i][index_j], 1, SUB_ARRAY_RECV, r, 0, comm, &request);
            MPI_Type_free(&SUB_ARRAY_RECV);
        }
        MPI_Wait(&request, &status);
    }
    else
    {
        // Non send data from all ranks to rank 0
        cout << "Rank " << rank << "getting ready to send data to rank 0." << endl;
        MPI_Ssend(&local_buf[1][1], 1, SUB_ARRAY_SEND, 0, 0, comm);
        cout <<"Rank " << rank <<": Done." << endl;
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


void decompose(int &M_local, int &N_local, double pM, double pN, int rank, int size, int** proc_dims)
{
    // Decomposes the image. For now assume the number of processes divides N
    // and M evenly 
    // M_local = M/dim[0];
    // N_local = N/dim[1];

    int coords[2] = {0,0};
    MPI_Cart_coords(comm, rank, 2, coords); 
    decompose_single_dim(M_local, M, pM, coords[0], dim[0]);
    decompose_single_dim(N_local, N, pN, coords[1], dim[1]);

    int local_dims[2];
    // rank 0 waits to receive from every type 
    MPI_Status status; 
    MPI_Request request; 
    if (rank == 0)
    {
        proc_dims[0][0] = M_local;
        proc_dims[0][1] = N_local; 
        for (int i = 1; i < size ; i++)
        {
            MPI_Recv(&proc_dims[i][0], 2, MPI_REALNUMBER, i, 0, comm, &status);
        }
    }
    else
    {
        local_dims[0] = M_local;
        local_dims[1] = N_local; 
        MPI_Ssend(&local_dims[0], 2, MPI_REALNUMBER, 0, 0, comm);
    }
    // Probably temporary, process 0 then sends the array to everyone 
    MPI_Bcast(&proc_dims[0][0], size*2, MPI_INT, 0, comm);
}

void decompose_single_dim(int &A_local, int A, double pA, int coord, int P)
{
    // Decomposes A into A_local, given the edge/inside fraction for that 
    // dimension pA. 
    int rem = 0;
    if (P > 2)
    {
        int A_edge = int(pA * A/(2*pA + P - 2));
        int A_inner = int(A - 2*A_edge)/(P -2);
        cout << A_edge << " " << A_inner << endl;
        rem = A%(2*A_edge + (P-2)*A_inner);
        // Check to see if process is in the corners 
        if (coord == 0 || coord == P-1){A_local = A_edge;}
        else{A_local = A_inner;}
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