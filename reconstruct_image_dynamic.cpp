/*
Parallel code for image reconstructio
*/ 

#include "pgmio.h"
#include "2darray.h"
#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std; 

void define_datatypes(MPI_Datatype &double_row, MPI_Datatype &double_col, int M, int N);
void loop(double** old, double** new_array, double** edge, int i_init, int i_end, int j_init, int j_end, int N);
double boundaryval(int i, int m);

template <typename T>
T** create2darray(int n_rows, int n_cols, T init_val)
{
    // Creates dynamically allocated 2D contiguous array 
    T **arr = new T*[n_rows];
    int size = n_rows * n_cols;
    // First element stores the whole array 
    arr[0] = new T[size];

    // Assign the memory 
    for (int i = 0; i < n_rows; i++)
    {
        arr[i] = &arr[0][i*n_cols];
    }

    // Initialize values of array to given value 
    for (int i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            arr[i][j] = init_val;
        }
    }

    return arr; 
}

template <typename T>
void delete_array(T** t)
{
    // Deletes memory for t
    delete[] t; 
}

void reconstruct(double delta_max,  double** old, double** new_array, double** edge, int M, int N, int* ranks, MPI_Comm comm, int rank)
{
    MPI_Request request;
    MPI_Status status; 
    // Get the types for halo swaping 
    MPI_Datatype MPI_DOUBLE_ROW; 
    MPI_Datatype MPI_DOUBLE_COL; 
    define_datatypes(MPI_DOUBLE_ROW, MPI_DOUBLE_COL, M, N); 

    int rank_up = ranks[0];
    int rank_down = ranks[1];
    int rank_right = ranks[2];
    int rank_left = ranks[3];

    int M1 = (int)((M+2)/2);
    int N1 = (int)((N+2)/4);
    int N2 = 2*N1; 
    
    double** delta = create2darray<double>(M, N, 0.0); 
    double** max_delta = create2darray<double>(M, N, 0.0); 
    max_delta[0][0] = 100;
    int n = 0; 
    int n_comp_delta = 100; // compute delta at every n_comp_delta iterations; 

    //Loop reconstruct the image in buffer by using the old image
    while (max_delta[0][0] > delta_max)
    {
        // Receive row 0 from rank_up and send row M to rank_down
        MPI_Issend( &old[M][0], 1, MPI_DOUBLE_ROW, rank_down, 0, comm, &request);
        MPI_Recv( &old[0][0], 1, MPI_DOUBLE_ROW, rank_up, 0, comm, &status);
        // Loop over area 1
        loop(old, new_array, edge, 1, M1, N1, N2, N);
        MPI_Wait(&request, &status);
        // Receive row M+1 from rank_down and send row 1 to rank_up
        MPI_Issend( &old[1][0], 1, MPI_DOUBLE_ROW, rank_up, 0, comm, &request);
        MPI_Recv( &old[M+1][0], 1, MPI_DOUBLE_ROW, rank_down, 0, comm, &status);
        // Loop over area 2
        loop(old, new_array, edge, M1, M+1, N1, N2, N);
        MPI_Wait(&request, &status);
        // Send column 1 to rank_left and receive column N+1 to rank_right
        MPI_Issend( &old[0][1], 1, MPI_DOUBLE_COL, rank_left, 0, comm, &request);
        MPI_Recv ( &old[0][N+1], 1, MPI_DOUBLE_COL, rank_right, 0, comm, &status); 
        // Loop over area 3
        loop(old, new_array, edge, 1, M+1, N2, N+1, N);
        MPI_Wait(&request, &status);
        // Receive column 0 from rank_left and send column N to rank_right
        MPI_Issend( &old[0][N], 1, MPI_DOUBLE_COL, rank_right, 0, comm, &request);
        MPI_Recv( &old[0][0], 1, MPI_DOUBLE_COL, rank_left, 0, comm, &status);
        // Loop over area 4
        loop(old, new_array, edge, 1, M+1, 1, N1, N);
        MPI_Wait(&request, &status);
        
        // At the end of each iteration, set the old array equal to the new array
        // old = new_array; is there a better way to do this??? 
        // cout << "Rank " << rank << ": Run " << n << " finished." << endl;
        for (int i = 1 ; i < M+1; i++)
        {
            for (int j = 1; j < N+1; j++)
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
            MPI_Allreduce(&delta[0][0], &max_delta[0][0], M*N, MPI_DOUBLE, MPI_MAX, comm);
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

void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int* dim, int &rank, int* ranks, int* coords)
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
     

void scatter(MPI_Comm &comm, int rank, double **buf, double **local_buf, int* dim, int M, int N, int M_local, int N_local)
{
    // Scatters buf among all processes 
    int ni = dim[0]; int nj = dim[1]; 
    int send_rank_coords[2] = {0,0}; 
    int send_rank; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_MASTER;
    MPI_Type_vector(M_local, N_local, N, MPI_DOUBLE, &SUB_ARRAY_MASTER);  
    MPI_Type_commit(&SUB_ARRAY_MASTER);

    MPI_Datatype SUB_ARRAY_LOCAL; 
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_DOUBLE, &SUB_ARRAY_LOCAL);
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

void gather(MPI_Comm &comm, int rank, double** buf, double** local_buf, int* dim, int M, int N, int M_local, int N_local)
{
    // Sends all local bufs back to master buf 
    int ni = dim[0]; int nj = dim[1]; 
    int rec_coords[2] = {0,0}; 
    MPI_Status status; MPI_Request request; 
    MPI_Datatype SUB_ARRAY_MASTER;
    MPI_Type_vector(M_local, N_local, N, MPI_DOUBLE, &SUB_ARRAY_MASTER);  
    MPI_Type_commit(&SUB_ARRAY_MASTER);

    MPI_Datatype SUB_ARRAY_LOCAL; 
    MPI_Type_vector(M_local, N_local, N_local+2, MPI_DOUBLE, &SUB_ARRAY_LOCAL);
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


int main(int argc, char **argv)
{
    int rank, size; 
    int ranks[4]; 
    double t1, t2, dt;
    double val;
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int M, N; 
    char dest_name[100]; 
    char filename[100]; 
    double delta_max = atof(argv[3]);
    int n_runs = 1;
    if (argc == 5)
    {
        n_runs = atoi(argv[4]);
    }
    sprintf(dest_name, "%s_%d_%2f.pgm", argv[2], size, delta_max);
    sprintf(filename, "%s", argv[1] );

    // Create cartesian topology and assign rank 
    MPI_Comm cart_comm; 
    MPI_Comm comm_old = MPI_COMM_WORLD; 
    int dim[2] = {0,0};  
    int coords[2] = {0,0};
    create_topology(size, cart_comm, comm_old, dim, rank, ranks, coords);

    // Process 0 gets size and reads the image
    pgmsize(filename, &M, &N);
    double **buf_out = create2darray<double>(M, N, 0.0);
    double **buf_in = create2darray<double>(M, N, 0.0);
    if (rank == 0)
    {
        // buf_in = create2darray<double>(M, N, 0.0);
        // buf_out = create2darray<double>(M, N, 0.0);
        pgmread(filename, &buf_in[0][0], M, N);
    }

    // For now assume image can be decomposed evenly among all the processes. 
    // i.e. size divides M.
    int M_local = int(M/dim[0]);
    int N_local = int(N/dim[1]); 
    // decompose(dim, M, N, rank, &M_local, &N_local); 

    // Now each process processes the local array 
    double **old = create2darray<double>(M_local+2, N_local+2, 0.0);
    double **new_array = create2darray<double>(M_local+2, N_local+2, 0.0);
    double **edge_local = create2darray<double>(M_local+2, N_local+2, 0.0);
    int j_local; 
    MPI_Barrier(cart_comm);
    t1 = MPI_Wtime(); 
    // Loop for the given number of iterations
    for (int i = 0; i<n_runs ; i++)
        {
        for (int i = 0; i< (M_local+2) ; i++)
        {
            for (int j = 0; j< (N_local+2) ; j++)
            {
                edge_local[i][j] = 255; 
                old[i][j] = 255; 
                new_array[i][j] = 255;
            }
        }
        // Set sawtooth values for old 
        for (int j = 1; j <N_local+1; j++)
        {
            j_local =  int(j + N_local * coords[1]);
            val = boundaryval(j_local, N);
            old[0][j] = int(255*(1.0-val));
            old[M_local+1][j] = int(255*val);
        }
        // Scatter the data 
        scatter(cart_comm, rank, buf_in, edge_local, dim, M, N, M_local, N_local);
        // Reconstruct 
        reconstruct(delta_max, old, new_array, edge_local, M_local, N_local, ranks, cart_comm, rank); 
        // Gather the data back 
        gather(cart_comm, rank, buf_out, new_array, dim, M, N, M_local, N_local);
    }
    delete_array<double>(buf_in);
    if (rank != 0)
    {
        delete_array<double>(buf_out);
    }
    delete_array<double>(edge_local);
    delete_array<double>(new_array);
    delete_array<double>(old);

    MPI_Barrier(cart_comm);
    t2 = MPI_Wtime();
    dt = t2 - t1; 
    if (rank == 0)
    {
        cout << "Run finished, total time  ellapsed: " << dt <<"s." << endl; 
        cout << "Average runtime per reconstruction: " << dt/n_runs << endl;
    }
    // Master rank saves the image 
    if (rank == 0){pgmwrite(dest_name, &buf_out[0][0], M, N); delete_array<double>(buf_out);}
    MPI_Finalize(); 
}

double boundaryval(int i, int m)
{
    // Sets the boundary value 
    double val; 

    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) {val = 2.0 - val;};
    return val;
}

void define_datatypes(MPI_Datatype &double_row, MPI_Datatype &double_col, int M, int N)
{
    // Create the derived datatypes 
    MPI_Type_vector(1, N+2, 1, MPI_DOUBLE, &double_row);
    MPI_Type_commit(&double_row);

    MPI_Type_vector(M+2, 1, N+2, MPI_DOUBLE, &double_col);
    MPI_Type_commit(&double_col);
}

void loop(double** old, double** new_array, double** edge, int i_init, int i_end, int j_init, int j_end, int N)
{
    // Perform image reconstruction by calculating the new array elements for i = [i_init, i_end)
    // j = [j_init, j_end).
    for (int i = i_init; i<i_end; i++)
    {
        for (int j = j_init; j<j_end; j++)
        {
            new_array[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
            // *(new_array + i*(N+2) + j) = 0.25*(*(old + (i-1)*(N+2) + j) + *(old + (i+1)*(N+2) + j) + *(old + i*(N+2) + j-1) + *(old + i*(N+2) + j+1) - *(edge + i*(N+2) + j) );
        }
    } 
}
