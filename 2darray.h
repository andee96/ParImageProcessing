/* 
Date: 26/11/2018
Author: Andreas Malekos
File that contains declaration and implementation 
functions to create a custom 2D array, calcuate its maximum value and local 
average.  
*/

#ifndef __2DARRAY__H__
#define __2DARRAY__H__

// Declara functions
template <typename T>
T** create2darray(int n_rows, int n_cols, T init_val);
template <typename T>
T max_elem(T** t, int M, int N);
template <typename T>
T average(T** arr, int M, int N, int  size);

// Function implementaiton
template <typename T>
T** create2darray(int n_rows, int n_cols, T init_val)
{/* Creates dynamically allocated 2D contiguous array along dimension j.
    Inputs:
        int n_rows: number of rows
        int n_cols: number of colums
        T init_val: initial value of every array element.
    Returns:
        T**: 2D array of type T
*/    
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
T max_elem(T** arr, int M, int N)
{/* Finds the maximum element of array. 
    Inputs:
        T** arr: 2D array of type T
        int M: number of array rows
        int N: number of array columns
    Returns:
        T m: maximum element of array
    */
    T m = 0; 
    for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (arr[i][j] > m){m = arr[i][j];}
            }
        }
    return m; 
}

template <typename T>
T average(T** arr, int M, int N, int  size)
{/* Calculates the local average of the an array with haloes.
    Inputs:
        T** arr: 2D array of type T
        int M: number of rows of arr without haloes.
        int N: number of columns of arr without haloes. 
        int size: total size of array 
    Returns:
        T avg: local average
    */
    T sum = 0; 
    for (int i = 1; i <M+1; i++)
    {
        for (int j = 1; j <N+1; j++)
        {
            sum = sum + arr[i][j];
        }
    }
    T avg = sum/size;
    return avg;
}
  
#endif