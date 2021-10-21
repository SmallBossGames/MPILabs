#pragma once

#include<mpi.h>
#include<stdio.h>
#include<cmath>

namespace mpi_labs::algorythms::fox
{
    using namespace std;

    // ��������� 7.1
// �������� ����� ��������� ������ � ������� ������������� ������
//   ������� ���������� ���������: ��� ������� ����������, 
//   ������ ������ � �� ���������� �� ����������� � ���������
//   ���������, �������� �������� ���������� �������
    int ProcNum = 0;   	// ���������� ��������� ��������� 
    int ProcRank = 0;  	// ���� �������� ��������
    int GridSize;      	// ������ ����������� ������� ���������
    int GridCoords[2];  	// ���������� �������� �������� � ����������
    // �������
    MPI_Comm GridComm;  	// ������������ � ���� ���������� �������
    MPI_Comm ColComm;   	// ������������ � ������� �������
    MPI_Comm RowComm;   	// ������������ � ������ �������
    MPI_Datatype MPI_BLOCK;

    // �������� ������������� � ���� ��������� ���������� ������� 
    // � �������������� ��� ������ ������ � ������� ������� �������
    void CreateGridCommunicators() {
        int DimSize[2];  	// ���������� ��������� � ������ ��������� �������
        int Periodic[2];	// =1 ��� ������� ���������, ����������� ������������� 
        int Subdims[2];  	// =1 ��� ������� ���������, ������������ � ����������
        DimSize[0] = GridSize;
        DimSize[1] = GridSize;
        Periodic[0] = 0;
        Periodic[1] = 0;

        // �������� ������������� � ���� ���������� ������� 
        MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);

        // ����������� ��������� �������� � ������� 
        MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);

        // �������� �������������� ��� ����� ���������� �������
        Subdims[0] = 0;  // �������� ���������
        Subdims[1] = 1;  // ������� ������� ��������� � ����������
        MPI_Cart_sub(GridComm, Subdims, &RowComm);

        // �������� �������������� ��� �������� ���������� �������
        Subdims[0] = 1;
        Subdims[1] = 0;
        MPI_Cart_sub(GridComm, Subdims, &ColComm);
    }

    // ������� ��� ��������� ������ � ������������� �������� ������
    void ProcessInitialization(double*& pAMatrix, double*& pBMatrix, double*& pCMatrix, double*& pAblock, double*& pBblock, double*& pCblock, double*& pTemporaryAblock, int& Size, int& BlockSize)
    {
        if (ProcRank == 0) {
            do {
                printf("\nEnter matrix size: ");
                scanf("%d", &Size);

                if (Size % GridSize != 0) {
                    printf("������ ������ ������ ���� ������ ������� �����! \n");
                }
            } while (Size % GridSize != 0);
        }
        MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        BlockSize = Size / GridSize;

        pAblock = new double[BlockSize * BlockSize];
        pBblock = new double[BlockSize * BlockSize];
        pCblock = new double[BlockSize * BlockSize];
        pTemporaryAblock = new double[BlockSize * BlockSize];

        for (int i = 0; i < BlockSize * BlockSize; i++)
        {
            pCblock[i] = 0;
        }
        if (ProcRank == 0)
        {
            pAMatrix = new double[Size * Size];
            pBMatrix = new double[Size * Size];
            pCMatrix = new double[Size * Size];

            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++)
                {
                    pAMatrix[i * Size + j] = (double)rand() / RAND_MAX * 9;
                    // ������� ������� ��������� �������
                    // ��� ���������� ��������� �����: i + j/10.;
                    pBMatrix[i * Size + j] = -((double)rand() / RAND_MAX * 9);
                    // ������� ������� ��������� �������
                    // ��� ���������� ��������� �����: -(i + j/10.);
                }
            }
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++) {
                    // ���� ������� ������� �������, �� ������ ���
                }
            }
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++) {
                    // ���� ������� ������� �������, �� ������ ���
                }
            }
        }
    }

    void DataDistribution(double*& pAMatrix, double*& pBMatrix, double*& pTemporaryAblock, double*& pBblock, int& Size, int& BlockSize)
    {
        MPI_Type_vector(BlockSize, BlockSize, Size, MPI_DOUBLE, &MPI_BLOCK);
        MPI_Type_commit(&MPI_BLOCK);
        if (ProcRank == 0)
        {
            for (int r = 0; r < ProcNum; r++)
            {
                int c[2];
                MPI_Cart_coords(GridComm, r, 2, c);
                MPI_Send(pAMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1,
                    MPI_BLOCK, r, 0, GridComm);
                MPI_Send(pBMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1,
                    MPI_BLOCK, r, 0, GridComm);
            }
        }
        MPI_Status s;
        MPI_Recv(pAMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, GridComm, &s);
        MPI_Recv(pBMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, GridComm, &s);
        for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {}
        }
        for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {}
        }
    }

    // �������� ������ ������� � �� ������� ������� ��������� 
    void ABlockCommunication(int iter, double* pAblock,
        double* pMatrixAblock, int BlockSize) {

        // ����������� �������� �������� � ������ ���������� ������� 
        int Pivot = (GridCoords[0] + iter) % GridSize;

        // ����������� ������������� ����� � ��������� ����� ������
        if (GridCoords[1] == Pivot) {
            for (int i = 0; i < BlockSize * BlockSize; i++)
                pAblock[i] = pMatrixAblock[i];
        }

        // �������� �����
        MPI_Bcast(pAblock, BlockSize * BlockSize, MPI_DOUBLE, Pivot,
            RowComm);
    }

    // ��������� ��������� ������
    void BlockMultiplication(double* pAblock, double* pBblock,
        double* pCblock, int BlockSize) {
        // ���������� ������������ ��������� ������
        for (int i = 0; i < BlockSize; i++) {
            for (int j = 0; j < BlockSize; j++) {
                double temp = 0;
                for (int k = 0; k < BlockSize; k++)
                    temp += pAblock[i * BlockSize + k] * pBblock[k * BlockSize + j];
                pCblock[i * BlockSize + j] += temp;
            }
        }
    }

    // ����������� ����� ������ ������� � ����� ������� ���������� 
    // ������� 
    void BblockCommunication(double* pBblock, int BlockSize) {
        MPI_Status Status;
        int NextProc = GridCoords[0] + 1;
        if (GridCoords[0] == GridSize - 1) NextProc = 0;
        int PrevProc = GridCoords[0] - 1;
        if (GridCoords[0] == 0) PrevProc = GridSize - 1;
        MPI_Sendrecv_replace(pBblock, BlockSize * BlockSize, MPI_DOUBLE,
            NextProc, 0, PrevProc, 0, ColComm, &Status);
    }

    // ������� ��� ������������� ��������� ������
    void ParallelResultCalculation(double* pAblock, double* pMatrixAblock,
        double* pBblock, double* pCblock, int BlockSize) {
        for (int iter = 0; iter < GridSize; iter++) {
            // �������� ������ ������� A �� ������� ���������� �������
            ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
            // ��������� ������
            BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
            // ����������� ����� ������ ������� B � �������� ���������� 
            // �������
            BblockCommunication(pBblock, BlockSize);
        }
    }

    void ResultCollection(double*& pCMatrix, double*& pCblock, int& Size, int& BlockSize)
    {
        for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {}
        }
        MPI_Send(pCMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, GridComm);
        if (ProcRank == 0)
        {
            MPI_Status s;
            for (int r = 0; r < ProcNum; r++)
            {
                int c[2];
                MPI_Cart_coords(GridComm, r, 2, c);
                MPI_Recv(pCblock + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, GridComm, &s);
            }
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++) {}
            }
        }
    }


    void demo_function(int argc, char* argv[]) {
        double* pAMatrix; 	// ������ �������� ���������� ���������
        double* pBMatrix; 	// ������ �������� ���������� ���������
        double* pCMatrix; 	// �������������� �������
        int Size;        	// ������ ������
        int BlockSize;   	// ������ ��������� ������, ������������� 
                          // �� ���������
        double* pAblock;  	// ���� ������� � �� ��������
        double* pBblock;  	// ���� ������� � �� ��������
        double* pCblock;  	// ���� �������������� ������� � �� ��������
        double* pMatrixAblock;
        double Start, Finish, Duration;

        setvbuf(stdout, 0, _IONBF, 0);

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
        MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

        GridSize = sqrt((double)ProcNum);
        if (ProcNum != GridSize * GridSize)
        {
            if (ProcRank == 0)
                printf("Number of processes must be a perfect square \n");
        }
        else
        {
            if (ProcRank == 0)
                printf("Parallel matrix multiplication program\n");

            // �������� ����������� ������� ��������� � �������������� 
            // ����� � ��������
            CreateGridCommunicators();

            // ��������� ������ � ������������� ��������� ������
            ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock, pCblock, pMatrixAblock, Size, BlockSize);
            // ������� ������������� ������ ����� ����������
            DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size, BlockSize);

            // ���������� ������������� ������ �����
            ParallelResultCalculation(pAblock, pMatrixAblock, pBblock, pCblock, BlockSize);

            // ���� �������������� ������� �� ������� ��������
            ResultCollection(pCMatrix, pCblock, Size, BlockSize);

            // ���������� �������� ����������
            //ProcessTermination(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock, pCblock, pMatrixAblock);
        }

        MPI_Finalize();
    }

}