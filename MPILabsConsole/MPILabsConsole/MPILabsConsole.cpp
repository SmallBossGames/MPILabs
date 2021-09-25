
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>

using namespace std;

void CalculateResult(int** matrix, int* vector) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);


}

void InitRandomMatrix(double_t** matrix, size_t m, size_t n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            matrix[i][j] = rand() % 100;
        }
    }   
}

void InitRandomVector(double_t* vector, size_t m)
{
    for (int i = 0; i < m; i++)
    {
        vector[i] = rand() % 100;
    }
}

int main()
{
    MPI_Init(NULL, NULL);

    int rank;

    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t m = 0, n = 0;

    

    if (rank == 0)
    {
        cout << "Write m:" << endl;
        cin >> m;

        cout << "Write n:" << endl;
        cin >> n;
    }

    MPI_Bcast(&m, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    double_t* vector = new double[m];

    if (rank == 0)
    {
        InitRandomVector(vector, m);
    }

    MPI_Bcast(vector, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*string sb;

    for (size_t i = 0; i < m; i++)
    {
        sb.append(to_string(vector[i])).append(", ");
    }

    cout << sb << endl;*/


    if (rank == 0)
    {
        cout << "Threads: " << size << endl;

        double_t** matrix = new double*[m];

        for (size_t i = 0; i < m; i++)
        {
            matrix[i] = new double[n];
        }

        InitRandomMatrix(matrix, m, n);
        

        for (size_t i = 0; i < m; i++)
        {
            delete[] matrix[i];
        }

        delete[] matrix;
    } 
    else
    {
        cout << "I'm fucking slave" << endl;
    }

    

    cout << "Rank: " << rank << endl;
    cout << "M: " << m << endl;
    cout << "N: " << n << endl;


    delete[] vector;

    MPI_Finalize();
}
