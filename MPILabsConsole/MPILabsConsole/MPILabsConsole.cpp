// MPILabsConsole.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <mpi.h>

int main()
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the rank of the process
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // Print the message
    printf("Hello World! My rank is %d\n", my_rank);
    // Finalize the MPI environment.
    MPI_Finalize();
}
