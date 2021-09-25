
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

void InitRandomArray(double_t* vector, size_t m)
{
    for (int i = 0; i < m; i++)
        vector[i] = rand() % 100;
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
        InitRandomArray(vector, m);

    MPI_Bcast(vector, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double_t* matrix = new double[m * n];

    if (rank == 0)
    {
        cout << "Threads: " << size << endl;
        InitRandomArray(matrix, m * n);
    }
    else
        cout << "I'm fucking slave #" << rank << endl;

    cout << "Rank: " << rank << endl;
    cout << "M: " << m << endl;
    cout << "N: " << n << endl;

    int* sendNumber;
    int* sendIndex;
    double_t* receivedPartialMatrix;
    int notSendRows = n;

    sendNumber = new int[size];
    sendIndex = new int[size];
    receivedPartialMatrix = new double_t[size];

    int rowNumber = n / size;
    sendNumber[0] = rowNumber * n;
    sendIndex[0] = 0;
    receivedPartialMatrix[0] = rowNumber * n;

    for (int i = 1; i < size; i++) {
        notSendRows -= rowNumber;
        rowNumber = notSendRows / (size - i);
        sendNumber[i] = rowNumber * n;
        sendIndex[i] = sendIndex[i - 1] + sendNumber[i - 1];
    }

    MPI_Scatterv(matrix, sendNumber, sendIndex, MPI_DOUBLE, receivedPartialMatrix, sendNumber[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*for (size_t i = 0; i < m; i++)
    {
    }*/

    string sb;

    for (int i = 1; i < rowNumber * n; i++)
        sb.append(to_string(receivedPartialMatrix[i])).append(", ");

    cout << "Rank #" << rank << " received matrix: " << sb << endl;

    //cout << sb << endl;

    delete[] sendNumber;
    delete[] sendIndex;

    delete[] vector;

    if (rank == 0)
        delete[] matrix;


    MPI_Finalize();
}
