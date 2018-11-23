/*
File containing class that will handle all communicaitons
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
    // Initializes MPI, returns rank, size and coordinates 
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
{
    pM = PM;
    pN = PN;
    create_topology();
}

void ParImageProcessor::initialize(double PM, double PN, int ni, int nj)
{
    pM = PM;
    pN = PN;
    create_topology(ni, nj);
}

RealNumber** ParImageProcessor::read(char* filename)
{
    // Returns the local buffer and the dimensions of the "unhalo-ed" local file
    RealNumber **buf; 
    pgmsize(filename, &M, &N);
    decompose(); 
    // Get the starting indices for each process. This will be important later 
    get_proc_indices();
    RealNumber** local_buf = create2darray<RealNumber>(M_local+2, N_local+2, 255.0);
    // Rank 0 reads and scatters the process
    if (rank == 0)
    {
        buf = create2darray<RealNumber>(M, N, 0.0);
        pgmread(filename, &buf[0][0], M, N);
    }
    // Scatter master buf to each local process
    scatter(local_buf, buf); 
    if (rank == 0){delete[] buf;}
    return local_buf;    
}

void ParImageProcessor::write(char* filename, RealNumber** local_buf)
{
    // Gathers all data to one buffer and write it in filename
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

void ParImageProcessor::reconstruct(double delta_max, int n_comp_delta, RealNumber** old, RealNumber** new_array, RealNumber** edge)
{
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    // Get the types for halo swaping 
    MPI_Datatype MPI_REAL_ROW; 
    MPI_Datatype MPI_REAL_COL; 
    MPI_Type_vector(1, N_local+2, 1, MPI_REALNUMBER, &MPI_REAL_ROW);
    MPI_Type_commit(&MPI_REAL_ROW);

    MPI_Type_vector(M_local+2, 1, N_local+2, MPI_REALNUMBER, &MPI_REAL_COL);
    MPI_Type_commit(&MPI_REAL_COL);

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
    RealNumber average_local = 0; 
    RealNumber average_global = 0; 

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
        // old = new_array; is there a better way to do this??? 
        for (int i = 1 ; i < M_local+1; i++)
        {
            for (int j = 1; j < N_local+1; j++)
            {
                // Calculate absolute difference between old and new 
                if (n%n_comp_delta == 0)
                {
                    delta[i-1][j-1] = fabs(old[i][j] - new_array[i][j]);
                }
            old[i][j] = new_array[i][j];
            }
        }
        if (n%n_comp_delta == 0)
        {
            // Calculate max and average 
            max_delta_local = max_elem(delta, M_local, N_local);
            average_local = average(new_array, M_local, N_local, M*N);
            MPI_Allreduce(&max_delta_local, &max_delta_global, 1, MPI_REALNUMBER, MPI_MAX, comm);
            MPI_Allreduce(&average_local, &average_global, 1, MPI_REALNUMBER, MPI_SUM, comm);
         }
        n++;
    }
    delete[] delta; 
}

void ParImageProcessor::finalize()
{
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
void ParImageProcessor::loop(double** old, double** new_array, double** edge, int i_init, int i_end, int j_init, int j_end)
{
    // Perform image reconstruction by calculating the new array elements for i = [i_init, i_end)
    // j = [j_init, j_end).
    // std::cout << "RanK " << rank << " starting new loop from i: " << i_init << ", " << i_end <<", j: " << j_init <<", " <<j_end <<std::endl;
    for (int i = i_init; i<i_end; i++)
    {
        for (int j = j_init; j<j_end; j++)
        {
            new_array[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
        }
    } 
}

void ParImageProcessor::create_topology()
{
    // Given the number of processes and the size of the problem, create appropriate 2D (or 1D) topology
    // We want a longer dimension in the j direction, as that is how images are stored in memory. This means 
    // more processes along the i direction. 
    // Need to split n into two factors (out of which one can be 1). Biggest factor will be along j direction 
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
{
    // Same as create topology above, but MPI_Dims_create is not used, and the nubmer of dimensions is suppled by 
    // the user. It is the user's responsability to make sure that these are defined correctly.
    int period[2] = {0, 1};
    int reorder = 1; 
    int rank_up, rank_down, rank_left, rank_right; 
    dim[0] = ni; 
    dim[1] = nj; 
    MPI_Cart_create(comm_old, 2, dim, period, reorder, &comm);
    MPI_Comm_rank(comm, &rank); 
    MPI_Cart_coords(comm, rank, 2, coords);
    MPI_Cart_shift(comm, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left;  
}

void ParImageProcessor::scatter(RealNumber** local_buf, RealNumber** buf)
{
  // Scatters buf among all processes 
    int currM = M_local, currN = N_local;
    int index_i = 0, index_j = 0; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_SEND;
    MPI_Datatype SUB_ARRAY_RECV;
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_REALNUMBER, &SUB_ARRAY_RECV); 
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
            index_i = proc_indices[r][0];
            index_j = proc_indices[r][1];
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

void ParImageProcessor::gather(RealNumber** local_buf, RealNumber** buf)
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
            index_i = proc_indices[r][0];
            index_j = proc_indices[r][1];
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
        MPI_Ssend(&local_buf[1][1], 1, SUB_ARRAY_SEND, 0, 0, comm);
    }   
}

void ParImageProcessor::decompose()
{
    // Decomppse image. Each rank calculates its own M_local and N_local. All 
    // are then sent to rank 0, which gathers them and broadcats the whole array back.
    // This is done in order to link the indexes of the int_proc array to the ranks.  
    int coords[2] = {0,0};
    MPI_Cart_coords(comm, rank, 2, coords); 
    decompose_single_dim(M_local, M, pM, coords[0], dim[0]);
    decompose_single_dim(N_local, N, pN, coords[1], dim[1]);

    //  Create list containing the local dimensions for each rank 
    int local_dims[2];
    local_dims[0] = M_local; local_dims[1] = N_local;
    MPI_Allgather(&local_dims[0], 2, MPI_INT, &proc_dims[0][0], 2, MPI_INT, comm);
}

void ParImageProcessor::decompose_single_dim(int &A_local, int A, double pA, int coord, int P)
{

    // Decomposes A into A_local, given the edge/inside fraction for that 
    // dimension pA. 
    int rem = 0;
    if (P > 2)
    {
        int A_edge = int(pA * A/(2*pA + P - 2));
        int A_inner = int(A - 2*A_edge)/(P -2);
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

void ParImageProcessor::get_proc_indices()
{
    // Returns an sizex2 array where arr[i][0] and arr[i][1] are the starting 
    // indices for the ith process in each dimension. The first elements are always
    // 0.
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
