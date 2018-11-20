template <typename T>
T** create2darray(int n_rows, int n_cols, T init_val);
template <typename T>
void delete_array(T **t);
template <typename T>
T max_elem(T** t, int M, int N);template <typename T>
T max_elem(T** arr, int M, int N)
{
    // Finds max element of array 
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


