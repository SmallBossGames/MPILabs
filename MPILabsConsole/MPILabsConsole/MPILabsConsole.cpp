#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>

using namespace std;

void GatherResult(double_t* resultMatrix, double_t* resultPartialMatrix, int n, int rowNumber, int rank, int size) {
    int* receiveNumber;
    int* receiveIndex;
    int notGatheredRows = n;

    receiveNumber = new int[size];
    receiveIndex = new int[size];

    receiveNumber[0] = n / size;
    receiveIndex[0] = 0;

    for (int i = 1; i < size; i++) {
        notGatheredRows -= receiveNumber[i - 1];
        receiveNumber[i] = notGatheredRows / (size - i);
        receiveIndex[i] = receiveIndex[i - 1] + receiveNumber[i - 1];
    }

    MPI_Allgatherv(resultPartialMatrix, receiveNumber[rank], MPI_DOUBLE, resultMatrix, receiveNumber, receiveIndex, MPI_DOUBLE, MPI_COMM_WORLD);

    delete[] receiveNumber;
    delete[] receiveIndex;
}

void CalculateResult(double_t* resultPartialMatrix, double_t* sourceMatrix, double_t* sourceVector, int n, int rowNumber) {
    for (int i = 0; i < rowNumber; i++)
    {
        resultPartialMatrix[i] = 0;

        for (int j = 0; j < n; j++)
            resultPartialMatrix[i] += sourceMatrix[i * n + j] * sourceVector[j];
    }
}

void InitRandomArray(double_t* vector, size_t m)
{
    srand((unsigned)time(0));
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

    double_t* sourceVector = new double[m];

    if (rank == 0)
        InitRandomArray(sourceVector, m);

    MPI_Bcast(sourceVector, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double_t* sourceMatrix = new double[m * n];
    double_t* resultVector = new double[n];

    if (rank == 0)
    {
        cout << "Threads: " << size << endl;
        InitRandomArray(sourceMatrix, m * n);
    }
    else
        cout << "I'm fucking slave #" << rank << endl;

    cout << "Rank: " << rank << endl;
    cout << "M: " << m << endl;
    cout << "N: " << n << endl;

    int* sendNumber;
    int* sendIndex;
    double_t* receivedPartialMatrix;
    double_t* calculatedPartialMatrix;
    int notSendRows = n;

    sendNumber = new int[size];
    sendIndex = new int[size];
    receivedPartialMatrix = new double_t[size];
    calculatedPartialMatrix = new double_t[size];

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

    MPI_Scatterv(sourceMatrix, sendNumber, sendIndex, MPI_DOUBLE, receivedPartialMatrix, sendNumber[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*for (size_t i = 0; i < m; i++)
    {
    }*/

    string sb;

    for (int i = 0; i < rowNumber * n; i++)
        sb.append(to_string(receivedPartialMatrix[i])).append(", ");

    cout << "Rank #" << rank << " received matrix: " << sb << endl;

    //cout << sb << endl;

    delete[] sendNumber;
    delete[] sendIndex;

    CalculateResult(calculatedPartialMatrix, receivedPartialMatrix, sourceVector, n, rowNumber);

    GatherResult(resultVector, calculatedPartialMatrix, n, rowNumber, rank, size);

    if (rank == 0) {
        string resultSb;

        for (int i = 0; i < n; i++)
            resultSb.append(to_string(resultVector[i])).append(", ");

        cout << "Result is " << resultSb << endl;
    }

    delete[] sourceVector;

    if (rank == 0)
        delete[] sourceMatrix;

    MPI_Finalize();
}
