/*
Parallel code for image reconstructio
*/ 

#include "pgmio.h"
#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std; 

void define_datatypes(MPI_Datatype &double_row, MPI_Datatype &double_col, int M, int N);

void loop(double* old, double* new_array, double* edge, int i_init, int i_end, int j_init, int j_end, int N);

void reconstruct(double delta_max, double* old, double* new_array, double* edge, int M, int N, int* ranks, MPI_Comm comm, int rank)
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
    
    double delta[M][N]; 
    double max_delta[M][N]; 
    max_delta[0][0] = 100;
    int n = 0; 
    int n_comp_delta = 50; // compute delta at every n_comp_delta iterations; 

    // temp testing 
    // int temp_coord[2] = {0,0};
    // int temp_coord_down[2] = {0,0};
    // int temp_coord_up[2] = {0,0};
    // int temp_coord_left[2] = {0,0};
    // int temp_coord_right[2] = {0,0};

    // MPI_Cart_coords(comm, rank, 2, temp_coord);
    // if (rank_down != MPI_PROC_NULL)
    // {
    //     MPI_Cart_coords(comm, rank_down, 2, temp_coord_down);
    // }
    // else
    // {
    //     temp_coord_down[0] = -15; 
    //     temp_coord_down[1] = -15; 
    // }

    // if (rank_right != MPI_PROC_NULL)
    // {
    //     MPI_Cart_coords(comm, rank_right, 2, temp_coord_right);
    // }
    // else
    // {
    //     temp_coord_right[0] = -15; 
    //     temp_coord_right[1] = -15; 
    // }

    // if (rank_up != MPI_PROC_NULL)
    // {
    //     MPI_Cart_coords(comm, rank_up, 2, temp_coord_up);
    // }
    // else
    // {
    //     temp_coord_up[0] = -15; 
    //     temp_coord_up[1] = -15; 
    // }

    // if (rank_left != MPI_PROC_NULL)
    // {
    //     MPI_Cart_coords(comm, rank_left, 2, temp_coord_left);
    // }
    // else
    // {
    //     temp_coord_left[0] = -15; 
    //     temp_coord_left[1] = -15; 
    // }

    // cout << "Rank " << temp_coord[0] << ", " << temp_coord[1] << ": Rank to my left is " << temp_coord_left[0] << ", " << temp_coord_left[1] <<"\nRank to my right is " << temp_coord_right[0] << ", " << temp_coord_right[1]
    // << "\nRank to my down is " << temp_coord_down[0] << ", " << temp_coord_down[1] << "\nRank to my up is " << temp_coord_up[0] << ", " << temp_coord_up[1] << endl;
    //Loop reconstruct the image in buffer by using the old image
    while (max_delta[0][0] > delta_max)
    {
        // Receive row 0 from rank_up and send row M to rank_down
        MPI_Issend( (old + M*(N+2)), 1, MPI_DOUBLE_ROW, rank_down, 0, comm, &request);
        MPI_Recv( (old), 1, MPI_DOUBLE_ROW, rank_up, 0, comm, &status);
        // Loop over area 1
        loop(old, new_array, edge, 1, M1, N1, N2, N);
        MPI_Wait(&request, &status);
        // Receive row M+1 from rank_down and send row 1 to rank_up
        MPI_Issend( (old + (N+2)), 1, MPI_DOUBLE_ROW, rank_up, 0, comm, &request);
        MPI_Recv( (old + (M+1)*(N+2)), 1, MPI_DOUBLE_ROW, rank_down, 0, comm, &status);
        // Loop over area 2
        loop(old, new_array, edge, M1, M+1, N1, N2, N);
        MPI_Wait(&request, &status);
        // Receive N+1 column from rank_right and send column 1 to rank_left 
        // if (temp_coord[1] == 1)
        // {

        MPI_Issend( (old + 1), 1, MPI_DOUBLE_COL, rank_left, 0, comm, &request);
        // }
        // if (temp_coord[1] == 0)
        // {
        MPI_Recv ( (old + (N+1)), 1, MPI_DOUBLE_COL, rank_right, 0, comm, &status); 
        // }
        // Loop over area 3
        loop(old, new_array, edge, 1, M+1, N2, N+1, N);
        MPI_Wait(&request, &status);
        // Receive column 0 from rank_left and send column N to rank_right
        // cout << "Rank " << temp_coord[0] << ", " << temp_coord[1] << "preparing to send column N to " << temp_coord_right[0] << ", " << temp_coord_right[1] << endl; 
        // if (temp_coord[1] == 0)
        // {

        MPI_Issend( (old + N), 1, MPI_DOUBLE_COL, rank_right, 0, comm, &request);
        // }
        // if (temp_coord[1] == 1)
        // {
        MPI_Recv( (old), 1, MPI_DOUBLE_COL, rank_left, 0, comm, &status);
        // }
        // Loop over area 4
        loop(old, new_array, edge, 1, M+1, 1, N1, N);
        MPI_Wait(&request, &status);
        
        // At the end of each iteration, set the old array equal to the new array
        // old = new_array; is there a better way to do this??? 
        for (int i = 1 ; i < M+1; i++)
        {
            for (int j = 1; j < N+1; j++)
            {
                // Calculate absolute difference between old and new 
                if (n%n_comp_delta == 0)
                {
                    *((double*)delta + (i-1)*(N) + j-1) = abs(*(old + i*(N+2) + j) - *(new_array + i*(N+2) + j));
                }
                *(old + i*(N+2) + j ) = *(new_array + i*(N+2) + j);
            }
        }
        if (n%n_comp_delta == 0)
        {
            // Reduce max 
            MPI_Allreduce(&delta[0][0], &max_delta, M*N, MPI_DOUBLE, MPI_MAX, comm);
            // cout << max_delta[0][0] << endl;
         }
        n++;
    }
    if (rank == 0)
    {
        cout << "Reconstruction done, number of iterations: " << n << endl;
    }
}

void create_topology(int n, MPI_Comm &comm_new, MPI_Comm &comm_old, int* dim, int &rank, int* ranks)
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
    MPI_Cart_shift(comm_new, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm_new, 1, 1, &rank_left, &rank_right); 
    ranks[0] = rank_up; ranks[1] = rank_down;
    ranks[2] = rank_right; ranks[3] = rank_left; 
}

// void decompose(int dim, int M, int N, int rank, int &M_local, int &N_local, MPI_Comm &comm)
// {
//     // Decompose process space given the number of processes in each dimension 
//     int n = dim[0] * dim[1]; // total number of processes 
//     int coords[2]; // process coordinates 
//     MPI_Cart_coords(comm, rank, 2, &coords); 
//     M_local = int(M/dim[0]);
//     N_local = int(N/dim[1]); 
//     int M_remainder = int(M%dim[0]); 
//     int N_remainder = int(N%dim[1]);
//     n_edges = 2*(dim[0] + dim[1] ) - 4 // number of processes on outside of image
//     n_inside = n - n_edges // number of processes on the inside of the image

//     // For now, simply distribute the extra work around the edges 
//     if (coords[0] < M_remainder)
//     {
//         M_local = M_local + 1; 
//     }
//     if (coords[1] < N_remainder)
//     {
//         N_local = N_local + 1; 
//     }
     

void scatter(MPI_Comm &comm, int rank, double* buf, double* local_buf, int* dim, int M, int N, int M_local, int N_local)
{
    // Scatters buf among all processes 
    int ni = dim[0]; int nj = dim[1]; 
    int coords[2]; 
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
        MPI_Issend(buf, 1, SUB_ARRAY_MASTER, 0, 0, comm, &request);
        MPI_Recv( (local_buf + (N_local+2) + 1), 1, SUB_ARRAY_LOCAL, 0, 0, comm, &status);
        MPI_Wait(&request, &status);
       for (int n = 1; n < ni*nj; n++)
       {
           send_rank = n; 
           MPI_Cart_coords(comm, send_rank, 2, send_rank_coords); 
           MPI_Ssend( (buf + send_rank_coords[0]*(M_local)*(N) + send_rank_coords[1]*(N_local) ),
                    1, SUB_ARRAY_MASTER, send_rank, 0, comm);
       }
    }
    else
    {
        MPI_Recv( (local_buf + (N_local+2) + 1), 1, SUB_ARRAY_LOCAL, 0, 0, comm, &status);
    }
}

void gather(MPI_Comm &comm, int rank, double* buf, double* local_buf, int* dim, int M, int N, int M_local, int N_local)
{
    // Sends all local bufs back to master buf 
    int ni = dim[0]; int nj = dim[1]; 
    int coords[2]; 
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
        MPI_Issend( (local_buf + (N_local+2) + 1), 1, SUB_ARRAY_LOCAL, 0, 0, comm, &request);
        MPI_Recv( buf, 1, SUB_ARRAY_MASTER, 0, 0, comm, &status);
        MPI_Wait(&request, &status); 
        // Blocking receive from each rank 
        for (int n_rank = 1; n_rank<ni*nj; n_rank++ )
        {
            MPI_Cart_coords(comm, n_rank, 2, rec_coords);
            MPI_Recv( (buf + rec_coords[0]*M_local*N + rec_coords[1]*N_local), 1, 
                        SUB_ARRAY_MASTER, n_rank, 0, comm, &status);

        }
    }
    else
    {
        // Non send data from all ranks to rank 0
        MPI_Ssend( (local_buf + (N_local+2) + 1), 1, SUB_ARRAY_LOCAL, 0, 0, comm);
    }
}




int main(int argc, char **argv)
{
    int rank, size; 
    int ranks[4]; 
    double t1, t2, dt; 
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
    create_topology(size, cart_comm, comm_old, dim, rank, ranks);

    // Process 0 gets size and reads the image
    pgmsize(filename, &M, &N);
    double buf_in[M][N];
    double buf_out[M][N];
    if (rank == 0)
    {
        cout << filename;
        pgmread(filename, buf_in, M, N);

    }
    else
    {
        for (int i = 0; i < M; i++)
        {
            for (int j =0; j<N; j++)
            {
                buf_in[i][j] = 0;
            }
        }
    }

    // For now assume image can be decomposed evenly among all the processes. 
    // i.e. size divides M.
    int M_local = int(M/dim[0]);
    int N_local = int(N/dim[1]); 
    // decompose(dim, M, N, rank, &M_local, &N_local); 
    
    // Define new datatype that corresponds to an M_local by N subsection of the bufer, that will then be passed 
    // directly to edge_local. 

    // Now each process processes the local array 
    double old[M_local+2][N_local+2]; 
    double new_array[M_local+2][N_local+2] ;
    double edge_local[M_local+2][N_local+2];


    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime(); 
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
    // Scatter the data 
    scatter(cart_comm, rank, (double*)buf_in, (double*)edge_local, dim, M, N, M_local, N_local);
    // Initialize arrays to 0 

    // fill_edge_array((double*)edge_local, (double*)local_buf, M_local, N, rank_right, rank_left, comm);  
    reconstruct(delta_max, (double*)old, (double*)new_array, (double*)edge_local, M_local, N_local, ranks, cart_comm, rank); 

    // Send local buffers back to global buffer
    // MPI_Gather( &new_array[1][1], 1, SUB_ARRAY_LOCAL, 
    //             &buf, 1, SUB_ARRAY_MASTER, 
    //             0, cart_comm );
    gather(cart_comm, rank, (double *)buf_out, (double *)new_array, dim, M, N, M_local, N_local);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();
    dt = t2 - t1; 
    if (rank == 0)
    {
        cout << "Run finished, total time  ellapsed: " << dt <<"s." << endl; 
        cout << "Average runtime per reconstruction: " << dt/n_runs << endl;
    }


    // Master rank saves the image 
    if (rank == 0){pgmwrite(dest_name, buf_out, M, N);}
    MPI_Finalize(); 
}

void define_datatypes(MPI_Datatype &double_row, MPI_Datatype &double_col, int M, int N)
{
    // Create the derived datatypes 
    MPI_Type_vector(1, N+2, 1, MPI_DOUBLE, &double_row);
    MPI_Type_commit(&double_row);

    MPI_Type_vector(M+2, 1, N+2, MPI_DOUBLE, &double_col);
    MPI_Type_commit(&double_col);
}

void loop(double* old, double* new_array, double* edge, int i_init, int i_end, int j_init, int j_end, int N)
{
    // Perform image reconstruction by calculating the new array elements for i = [i_init, i_end)
    // j = [j_init, j_end).
    for (int i = i_init; i<i_end; i++)
    {
        for (int j = j_init; j<j_end; j++)
        {
            *(new_array + i*(N+2) + j) = 0.25*(*(old + (i-1)*(N+2) + j) + *(old + (i+1)*(N+2) + j) + *(old + i*(N+2) + j-1) + *(old + i*(N+2) + j+1) - *(edge + i*(N+2) + j) );
        }
    } 
}
