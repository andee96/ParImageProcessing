#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char* fname = argv[1];
    string line; 
    ifstream myfile(fname);
    double a, b;
    if (myfile.is_open())
    {
        // Get line 1 
        getline(myfile, line);
        a = atof(line.c_str());
        getline(myfile, line);
        b = atof(line.c_str());
    }
    cout << a << " " << b << " " << rank << endl;
    myfile.close();

    MPI_Finalize();

}