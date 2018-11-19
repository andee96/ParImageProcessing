/* Simple program that reads in the smallest pgm image, stores 
it in a 2D array and then writes it in the same folder again. */ 

#include "pgmio.h"
#include <iostream>

using namespace std; 


int main()
{

int M, N; 
char* filename = (char*) "edgenew768x768.pgm";
char* dest_name = (char*) "test.pgm";
// Get the image sizes 
pgmsize(filename, &M, &N);

// Declare the image 
double buf[M][N];

// Read pgm image 
pgmread(filename, buf, M, N); 

for (int i = 0; i<N; i++)
{
    for (int j = 0; j<M; j++)
    {
        cout << buf[i][j] << endl; 
    }
}

// Write the file again an a pgm file 
pgmwrite(dest_name, buf, M, N);

} 