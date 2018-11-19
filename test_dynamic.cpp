/* 
Testing dynamic allocation memory scattering and gathering in c++
*/ 

#include "mpi.h"
// #include "2darray.h"
#include <iostream>

using namespace std; 

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

void process_array(double** arr, int M, int N, int rank)
{
    cout << &arr[0][0] << endl;
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            // cout << "testing " << endl;
            // cout <<  i << " " << j << " " << rank << endl;  
            arr[i][j] = rank; 
        }
    }
}


int main(int argv, char **argc)
{
    int rank, size; 
    MPI_Init(&argv, &argc);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int M = 10; 
    int N = 10; 
    int M_local = M/2; 
    int N_local = N/2;

    double **masterbuf = create2darray<double>(M, N, 0.0);
    double **localbuf = create2darray<double>(M_local, N, 0.0); 
    // Scatter the 2D array to local processes 
    MPI_Datatype SUB_ARRAY; 
    MPI_Type_vector(M_local, N, N, MPI_DOUBLE, &SUB_ARRAY);
    MPI_Type_commit(&SUB_ARRAY); 

    MPI_Datatype MPI_ROW, MPI_COL;
    MPI_Type_vector(1, N, 1, MPI_DOUBLE, &MPI_ROW);
    MPI_Type_commit(&MPI_ROW);
    MPI_Type_vector(M_local, 1, N, MPI_DOUBLE, &MPI_COL);
    MPI_Type_commit(&MPI_COL);

    MPI_Scatter(&masterbuf[0][0], 1, SUB_ARRAY, &localbuf[0][0], 1, SUB_ARRAY, 0, MPI_COMM_WORLD);
    // delete_array<double>(masterbuf); 

    cout << "Rank " << rank <<" has recceived the data." << endl; 
    process_array(localbuf, M_local, N, rank);
    cout << &localbuf[0][0] << endl;
    cout << "testing" << endl;

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    MPI_Request request; 

    // Send a row from proces 0 to process 1 
    if (rank == 0)
    {
        MPI_Ssend(&localbuf[0][0], 1, MPI_ROW, 1, 0, comm);
    }
    else{
        MPI_Recv(&localbuf[0][0], 1 , MPI_ROW, 0, 0, comm, &status);
        for (int j = 0; j < N ; j++)
        {
            cout << localbuf[0][j] <<  " "<<endl; 
        }
    }

    // Send a column from process 0 to process 1 
    if (rank == 0 )
    {
        MPI_Issend(&localbuf[0][N-1], 1, MPI_COL, 1, 0, comm, &request);
        MPI_Recv(&localbuf[0][N-1], 1, MPI_COL, 1, 0, comm, &status);
        MPI_Wait(&request, &status);
    }
    else{
        MPI_Issend(&localbuf[0][N-1], 1, MPI_COL, 0, 0, comm, &request);
        MPI_Recv(&localbuf[0][N-1], 1, MPI_COL, 0, 0, comm, &status);
        MPI_Wait(&request, &status);
        for (int i = 0; i < M_local ; i++)
        {
            cout << localbuf[i][N-1] << " " << endl;
        }
    }



    // double **masterbuf = create2darray<double>(M, N, 0.0);

    MPI_Gather(&localbuf[0][0], 1, SUB_ARRAY, &masterbuf[0][0], 1, SUB_ARRAY, 0, MPI_COMM_WORLD);
    delete_array<double>(localbuf);
    if (rank !=0 )
    {
        delete_array<double>(masterbuf);
    }
    if (rank == 0 )
    {
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<M; j++)
        {
            // cout << i << " " << j << endl;
            cout << masterbuf[i][j] << " ";
        }
        cout << endl; 
        
    }   
    delete_array<double>(masterbuf);
    } 

    MPI_Finalize();
}


    
  
