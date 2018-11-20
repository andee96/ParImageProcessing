#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    int rank;
    int dim[2] = {2,2};
    int reorder = 1; 
    int period[2] = {0,0};
    MPI_Comm comm_new; 
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm_new);
    MPI_Comm_rank(comm_new, &rank);
    int coord[2];
    int rank_up, rank_down, rank_left, rank_right; 
    MPI_Cart_shift(comm_new, 0, 1, &rank_up, &rank_down); 
    MPI_Cart_shift(comm_new, 1, 1, &rank_left, &rank_right); 
    MPI_Cart_coords(comm_new, rank, 2, coord);
    cout << rank << " " << coord[0] << " " << coord[1]<< endl;
    cout << "Rank " << rank <<": To my up is " << rank_up << endl;
    cout << "Rank " << rank <<": To my left is " << rank_left << endl;
    cout << int(1/2) << endl;
    MPI_Finalize();

}