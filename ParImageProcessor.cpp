/*
Date: 26/11/2018
Author: Andreas Malekos
Files that contains the implementation of the ParImagProcessor 
methods declared in ParImageProcessor.h
*/

#include "mpi.h"
#include <iostream> 
#include "precision.h"
#include "2darray.h"
#include "pgmio.h"
#include <math.h>
#include "ParImageProcessor.h"

// Constructor definitions
ParImageProcessor::ParImageProcessor()
{
    /* Initializes MPI and assigns the size and initial rank to the process.
    The proc_dims and prod_indices arrays are also initialized. */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    comm_old = MPI_COMM_WORLD;
    // Initialize arrays with indices and dimensions of each process
    proc_dims = create2darray<int>(size, 2, 0);
    proc_indices = create2darray<int>(size, 2, 0);
    dim[0] = 0; dim[1] = 0;
}

// Public member function defintions 
void ParImageProcessor::initialize(double PM, double PN)
{ /* Assigns values to pM and pN and creates the cartesian topology
    using the build int MPI_Dims_create function
    Inputs:
        double pM: ratio of edge to inner column size.
        double pN: ratio of edge to inner row size.*/
    pM = PM;
    pN = PN;
    create_topology();
}

void ParImageProcessor::initialize(double PM, double PN, int ni, int nj)
{/* Assigns values to pM and pN and creates the cartesian topology
    using manual dimension creation
    Inputs:
        double pM: ratio of edge to inner column size.
        double pN: ratio of edge to inner row size.*/
    pM = PM;
    pN = PN;
    create_topology(ni, nj);
}

RealNumber** ParImageProcessor::read(char* filename)
{
    /* Returns the local buffer with the halo and its dimensions without the 
    halo. In the process, the image is decomposed and the proc_indices and 
    proc_dims arrays are popualted.
    Inputs:
        char* filename: Name of the file containing the edge image.
    Returns:
        RealNumber** buf: local buffer with haloes. 
    */
    RealNumber **buf; 
    // Get size of problem space
    pgmsize(filename, &M, &N);
    // Decompose the problem space 
    decompose(); 
    // Get the starting indices on the masterbuf array for each process.
    get_proc_indices();
    RealNumber** local_buf = create2darray<RealNumber>(M_local+2, N_local+2, 255.0);
    // Rank 0 reads the array
    if (rank == 0)
    {
        buf = create2darray<RealNumber>(M, N, 0.0);
        pgmread(filename, &buf[0][0], M, N);
    }
    // Scatter master buffer to each local process
    scatter(local_buf, buf); 
    if (rank == 0){delete[] buf;}
    return local_buf;    
}

void ParImageProcessor::write(char* filename, RealNumber** local_buf)
{/* Gather local buffer back to master process, and writes the master buffer 
    to specified result file.
    Inputs:
        char* filename: Name of result file. 
        RealNumber** local_buf: local buffer with haloes.
*/
    // Create empty array to store master buffer. 
    RealNumber **buf = create2darray<RealNumber>(M, N, 0.0);
    // Rank 0 gathers all the data to buf 
    gather(local_buf, buf);
    // Rank 0 writes the file to filename
    if (rank == 0)
    {
        pgmwrite(filename, &buf[0][0], M, N);
    }
    delete[] buf; 
}

RealNumber** ParImageProcessor::reconstruct(double delta_max, int n_comp_delta, RealNumber** old, RealNumber** edge)
{/* Performs image reconstruction until the value of max delta is reached. 
    Inputs:
        double delta_max: maximum value of delta
        int n_comp_delta: nubmer of iterations after which to calculate the current value of delta.
        RealNumber** old: old array containing sawtooth values.
        Realnumber** edge: local edge array with haloes.
    Returns:
        RealNumber** new_array: reconstructed local array with haloes.
    */
    RealNumber** new_array = create2darray<RealNumber>(M_local+2, N_local+2, 0.0);
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    // Define vector types for halo swapping
    MPI_Datatype MPI_REAL_ROW; 
    MPI_Datatype MPI_REAL_COL; 
    MPI_Type_vector(1, N_local+2, 1, MPI_REALNUMBER, &MPI_REAL_ROW);
    MPI_Type_commit(&MPI_REAL_ROW);

    MPI_Type_vector(M_local+2, 1, N_local+2, MPI_REALNUMBER, &MPI_REAL_COL);
    MPI_Type_commit(&MPI_REAL_COL);
    // Get the raks of the neighboring processes
    int rank_up = ranks[0];
    int rank_down = ranks[1];
    int rank_right = ranks[2];
    int rank_left = ranks[3];

    // Get the lengths of the local loop sections. 
    int M1 = (int)((M_local+2)/2);
    int N1 = (int)((N_local+2)/4);
    int N2 = 2*N1; 

    // Initialize delta and average 
    RealNumber** delta = create2darray<RealNumber>(M_local, N_local, 0.0); 
    RealNumber max_delta_local = 100;
    RealNumber max_delta_global = 100;
    RealNumber average_local = 0; 
    RealNumber average_global = 0; 

    // Initialize iteration counter 
    int n = 0; 
   
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
        MPI_Recv( &old[0][N_local+1], 1, MPI_REAL_COL, rank_right, 0, comm, &status); 
        // Loop over area 3
        loop(old, new_array, edge, 1, M_local+1, N2, N_local+1);
        MPI_Wait(&request, &status);
        // MPI_Issend column 0 from rank_left and send column N_local to rank_right
        MPI_Issend( &old[0][N_local], 1, MPI_REAL_COL, rank_right, 0, comm, &request);
        MPI_Recv( &old[0][0], 1, MPI_REAL_COL, rank_left, 0, comm, &status);
        // Loop over area 4
        loop(old, new_array, edge, 1, M_local+1, 1, N1);
        MPI_Wait(&request, &status);
        
        // At the end of each iteration, set the old array equal to the new array
        for (int i = 1 ; i < M_local+1; i++)
        {
            for (int j = 1; j < N_local+1; j++)
            {
                // Calculate absolute difference between old and new if n divideds n_comp_delta
                if (n%n_comp_delta == 0)
                {
                    delta[i-1][j-1] = fabs(old[i][j] - new_array[i][j]);
                }
            old[i][j] = new_array[i][j];
            }
        }
        if (n%n_comp_delta == 0)
        {
            // Calculate max and average of all processes using MPI_Allreduce
            max_delta_local = max_elem(delta, M_local, N_local);
            average_local = average(new_array, M_local, N_local, M*N);
            MPI_Allreduce(&max_delta_local, &max_delta_global, 1, MPI_REALNUMBER, MPI_MAX, comm);
            MPI_Allreduce(&average_local, &average_global, 1, MPI_REALNUMBER, MPI_SUM, comm);
            // Rank 0 prints the average out:
            if (rank == 0){std::cout << "The average is :" << average_global << std::endl;}
         }
        n++;
    }
    delete[] delta; 
    return new_array;
}

void ParImageProcessor::finalize()
{// Encapsulator for MPI Finalize
    MPI_Finalize();
}

double ParImageProcessor::record_time()
{
    // Encapsulator for MPI_Wtime command. 
    MPI_Barrier(comm);
    double t = MPI_Wtime();
    return t;
}

// Helper member functions defintions 
void ParImageProcessor::loop(RealNumber** old, RealNumber** new_array, RealNumber** edge, int i_init, int i_end, int j_init, int j_end)
{/* Performs image reconstruction by calculating the new array elements for i = [i_init, i_end)
    j = [j_init, j_end).
    Inputs:
        RealNumber** old: old array
        RealNumber** new_array: new_array
        RealNumber** edge: edge array 
        int i_init: starting row index
        int i_end: end row index
        int j_init: starting row index
        int j_end: ending row index
    */ 
    for (int i = i_init; i<i_end; i++)
    {
        for (int j = j_init; j<j_end; j++)
        {
            new_array[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
        }
    } 
}

void ParImageProcessor::create_topology()
{/* Creates cartesian topology and assigns new rank using automatic MPI_Dims_create. The ranks of the neighboring 
    processes are also calculated. 
    */
    int period[2] = {0, 1};
    int reorder = 1; 
    int rank_up, rank_down, rank_left, rank_right; 
    MPI_Dims_create(size, 2, dim);
    MPI_Cart_create(comm_old, 2, dim, period, reorder, &comm);
    MPI_Comm_rank(comm, &rank); 
    MPI_Cart_coords(comm, rank, 2, coords);
    MPI_Cart_shift(comm, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left; 
}

void ParImageProcessor::create_topology(int ni, int nj)
{/* Creates cartesian topology and assigns new rank by creating a cartesian topology
    with ni, nj processes along the first and second dimensions respectively. It 
    is up to the user to make sure that these are consistent with the number of processes.
    */
    int period[2] = {0, 1};
    int reorder = 1; 
    int rank_up, rank_down, rank_left, rank_right; 
    dim[0] = ni; 
    dim[1] = nj; 
    MPI_Cart_create(comm_old, 2, dim, period, reorder, &comm);
    MPI_Comm_rank(comm, &rank); 
    if (rank == 0){std::cout << dim[0] << " " << dim[1] << std::endl;}
    MPI_Cart_coords(comm, rank, 2, coords);
    MPI_Cart_shift(comm, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left;  
}

void ParImageProcessor::scatter(RealNumber** local_buf, RealNumber** buf)
{   /* Scatters buf from the master process to the local buffer of every 
    other process. 
    Inputs: 
        RealNumber** local_buf: local buffer where array will be scattered.
        RealNumber** buf: master buffer. Empty for every process apart rank 0.
    */
    // Declare and define current local size and indices 
    int currM = M_local, currN = N_local;
    int index_i = 0, index_j = 0; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_SEND;
    MPI_Datatype SUB_ARRAY_RECV;
    // Delcare custom vector type for receiving sub aray
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_REALNUMBER, &SUB_ARRAY_RECV); 
    MPI_Type_commit(&SUB_ARRAY_RECV);

    // Rank 0 sends the right sub-arrays to the other processes.
    if (rank == 0 )
    {
        // Define custom vector type for sending
        MPI_Type_vector(proc_dims[0][0], proc_dims[0][1], N, MPI_REALNUMBER, &SUB_ARRAY_SEND);
        MPI_Type_commit(&SUB_ARRAY_SEND);
        // Asyncrhonous send to itself 
        MPI_Issend(&buf[0][0], 1, SUB_ARRAY_SEND, 0, 0, comm, &request);
        MPI_Recv( &local_buf[1][1], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
        MPI_Type_free(&SUB_ARRAY_SEND);
        for (int r = 1; r < size; r++)
        {
            // Update current starting indices on buf for rank r. 
            index_i = proc_indices[r][0];
            index_j = proc_indices[r][1];
            // Update current local sizes for rank r. 
            currM = proc_dims[r][0];
            currN = proc_dims[r][1];
            // Create vector datatype for sending.
            MPI_Type_vector(currM, currN, N, MPI_REALNUMBER, &SUB_ARRAY_SEND);
            MPI_Type_commit(&SUB_ARRAY_SEND);
            // Wawit for previous send request to complete and send to rank r.
            MPI_Wait(&request, &status);
            MPI_Issend(&buf[index_i][index_j],1, SUB_ARRAY_SEND, r, 0, comm, &request);
            MPI_Type_free(&SUB_ARRAY_SEND);
        }
        MPI_Wait(&request, &status);
    }
    // Other ranks wait to receive the local buffer from rank 0.
    else
    {
        MPI_Recv(&local_buf[1][1], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
    }
}

void ParImageProcessor::gather(RealNumber** local_buf, RealNumber** buf)
{/* Gathers local buffers from all the ranks to rank 0. 
    Inputs:
        RealNumber** local_buf: local buffer with halos
        RealNumber** buf: master buffer with no halos.
*/
    int ni = dim[0]; int nj = dim[1]; 
    int currM, currN; 
    int index_i = 0, index_j=0; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_SEND, SUB_ARRAY_RECV;
    // Custom datatype for sending arrays
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_REALNUMBER, &SUB_ARRAY_SEND);  
    MPI_Type_commit(&SUB_ARRAY_SEND);

    // Rank 0 receives all the local_buf arrays and stores them in buf.
    if (rank == 0)
    {
        // Define vector type for receiving
        MPI_Type_vector(M_local, N_local, N, MPI_REALNUMBER, &SUB_ARRAY_RECV);
        MPI_Type_commit(&SUB_ARRAY_RECV);
        // Non blocking send from rank 0 to itself
        MPI_Issend(&local_buf[1][1], 1, SUB_ARRAY_SEND, 0, 0, comm, &request);
        MPI_Recv( &buf[0][0], 1, SUB_ARRAY_RECV, 0, 0, comm, &status);
        MPI_Type_free(&SUB_ARRAY_RECV);
        // Blocking receive from each rank 
        for (int r = 1; r<size; r++ )
        {
            // Update starting indices on buf and local dimensions for r. 
            index_i = proc_indices[r][0];
            index_j = proc_indices[r][1];
            currM = proc_dims[r][0];
            currN = proc_dims[r][1];

            // Define vector type for receibing sub-array from rank r. 
            MPI_Type_vector(currM, currN, N, MPI_REALNUMBER, &SUB_ARRAY_RECV);
            MPI_Type_commit(&SUB_ARRAY_RECV);
            // Wait for previous receive to complete and receive the next sub-array from r.
            MPI_Wait(&request, &status);
            MPI_Irecv(&buf[index_i][index_j], 1, SUB_ARRAY_RECV, r, 0, comm, &request);
            MPI_Type_free(&SUB_ARRAY_RECV);
        }
        MPI_Wait(&request, &status);
    }
    // Othe ranks wait to send their local buffers to 0. 
    else
    {
        MPI_Ssend(&local_buf[1][1], 1, SUB_ARRAY_SEND, 0, 0, comm);
    }   
}

void ParImageProcessor::decompose()
{/*
    Decomppse problem space. Each rank calculates its own M_local and N_local. All 
    are then sent to rank 0, which gathers them and broadcats the whole array back.
    This is done in order to link the indexes of the int_proc array to the ranks. 
*/ 
    int coords[2] = {0,0};
    MPI_Cart_coords(comm, rank, 2, coords); 
    // Decompose each dimenion separetely according to pM and rank coordinates
    decompose_single_dim(M_local, M, pM, coords[0], dim[0]);
    decompose_single_dim(N_local, N, pN, coords[1], dim[1]);

    //  Create list containing the local dimensions for each rank 
    int local_dims[2];
    local_dims[0] = M_local; local_dims[1] = N_local;
    // All gather to populate proc_dims array. 
    MPI_Allgather(&local_dims[0], 2, MPI_INT, &proc_dims[0][0], 2, MPI_INT, comm);
}

void ParImageProcessor::decompose_single_dim(int &A_local, int A, double pA, int coord, int P)
{
    /* Decomposes A into A_local, given the edge/inner fraction for that 
     dimension, pA. 
     Inputs:
        int A_local: local size along dimension to be defined.
        int A: global size of dimension
        double pA: edge/inner fractional size for given dimension
        int coord: coordinate of process along given dimension
        int P: number of processes along given dimension
     */
     // Initialize remainder
    int rem = 0;
    // Case where there are more than two process along given dimension
    if (P > 2)
    {
        // Calculate the integer parts of number of edge and inner processes
        int A_edge = int(pA * A/(2*pA + P - 2));
        int A_inner = int(A - 2*A_edge)/(P -2);
        rem = A%(2*A_edge + (P-2)*A_inner);
        // Check to see if process is in the corners 
        if (coord == 0 || coord == P-1){A_local = A_edge;}
        else{A_local = A_inner;}
        // If the remainder is non-zero, assign naively.
        if (coord < rem){A_local++;}
    }
    // Case where only 1 or 2 processes along given dimension, i.e. no middle processes.
    else
    {
        // Case where we have 2 or 1 corners, i.e. no middle processes 
        A_local = int(A/P);
        rem = A%(P);
        // If remainder is non-zero assign naively.
        if (coord < rem){A_local++;} 
    }
}

void ParImageProcessor::get_proc_indices()
{
    /* Returns an sizex2 array where arr[i][0] and arr[i][1] are the starting 
    indices in the masterbuf array for the ith process in each dimension. 
    The first elements are always 0.*/
    int currM = M_local; int currN = N_local;
    int index_i = 0, index_j=0;
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
        // Calculates new values for currM and currN
        currM = proc_dims[r][0];
        currN = proc_dims[r][1];
        // Calculate new indices 
        proc_indices[r][0] = index_i;
        proc_indices[r][1] = index_j;
    }
}
