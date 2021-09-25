
#include <iostream>
#include <mpi.h>

using namespace std;

const int n = 10;
const int m = 10;

void CalculateResult(int matrix[n][m], int vector[m]) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);


}

int main()
{
    int threadsCount;
    cout << "Enter count of threads" << endl;
    cin >> threadsCount;

    int matrix[n][m];
    int vector[m];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j] = rand() % 100;

    for (int i = 0; i < m; i++)
        vector[i] = rand() % 100;

    MPI_Init(NULL, NULL);

    CalculateResult(matrix, vector);

    MPI_Finalize();
}
