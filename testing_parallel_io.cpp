//testing parallel io 
#include "mpi.h"
#include "pgmio.h"
#include <iostream>

using namespace std;

int main(int argv, char **argc)
{
    char filename[100];
    sprintf(filename, "%s", argc[1]);

    int rank, size; 
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_Offset filesize; 
    MPI_File_get_size(fh, &filesize);
    int bufsize = filesize/size; 
    int ndoubles = bufsize/sizeof(double);

    cout << sizeof(int) << " " << sizeof(double) << " " <<  sizeof(float) << endl ; 

    MPI_File_seek(fh, rank * bufsize, MPI_SEEK_SET);
    double buf[96][92]; 
    MPI_Status status; 
    // MPI_File_read(fh, &buf[0][0], ndoubles, MPI_DOUBLE, &status);

    cout << ndoubles << endl; 
    cout << bufsize << endl;

    MPI_Finalize();

}