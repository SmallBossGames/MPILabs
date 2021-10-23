#pragma once

#include<mpi.h>
#include<iostream>
#include<cmath>

namespace mpi_labs::algorythms::fox
{
    using namespace std;

    void print_debug_message(const char* msg, int rank) 
    {
        if (rank == 0)
        {
            cout << msg << endl;
        }
    }

    // ��������� 7.1
// �������� ����� ��������� ������ � ������� ������������� ������
//   ������� ���������� ���������: ��� ������� ����������, 
//   ������ ������ � �� ���������� �� ����������� � ���������
//   ���������, �������� �������� ���������� �������
    int ProcNum;   	// ���������� ��������� ��������� 
    int ProcRank;  	// ���� �������� ��������
    int GridSize;      	// ������ ����������� ������� ���������
    int GridCoords[2];  	// ���������� �������� �������� � ���������� �������
    MPI_Comm Communicator;  	// ������������ � ���� ���������� �������
    MPI_Comm ColComm;   	// ������������ � ������� �������
    MPI_Comm RowComm;   	// ������������ � ������ �������
    MPI_Datatype MPI_BLOCK;

    // �������� ������������� � ���� ��������� ���������� ������� 
    // � �������������� ��� ������ ������ � ������� ������� �������
    void CreateGridCommunicators() {
        int DimSize[2] = { GridSize , GridSize };  	// ���������� ��������� � ������ ��������� �������
        int Periodic[2] = 
        { 
            0, // =1 ��� ������� ���������, ����������� ������������� 
            0, // =1 ��� ������� ���������, ������������ � ����������
        };

        // �������� ������������� � ���� ���������� ������� 
        MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 0, &Communicator);

        // ����������� ��������� �������� � ������� 
        MPI_Cart_coords(Communicator, ProcRank, 2, GridCoords);

        // �������� �������������� ��� ����� ���������� �������
        const int SubdimsRows[2] = { 0, 1 };

        MPI_Cart_sub(Communicator, SubdimsRows, &RowComm);

        // �������� �������������� ��� �������� ���������� �������
        const int SubdimsColumns[2] = { 1, 0 };

        MPI_Cart_sub(Communicator, SubdimsColumns, &ColComm);

        print_debug_message("Creation of comunicators completed", ProcRank);
    }

    // ������� ��� ��������� ������ � ������������� �������� ������
    void ProcessInitialization(double*& aMatrix, double*& bMatrix, double*& cMatrix, double*& aBlock, double*& bBlock, double*& cBlock, double*& pTemporaryAblock, int& Size, int& BlockSize)
    {
        if (ProcRank == 0) {
            do {
                cout << endl << "\nEnter matrix size: ";
                cin >> Size;

                if (Size % GridSize != 0) {
                    printf("������ ������ ������ ���� ������ ������� �����! \n");
                }
            } while (Size % GridSize != 0);
        }
        MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        BlockSize = Size / GridSize;

        aBlock = new double[BlockSize * BlockSize];
        bBlock = new double[BlockSize * BlockSize];
        cBlock = new double[BlockSize * BlockSize];
        pTemporaryAblock = new double[BlockSize * BlockSize];

        for (int i = 0; i < BlockSize * BlockSize; i++)
        {
            cBlock[i] = 0;
        }
        aMatrix = new double[Size * Size];
        bMatrix = new double[Size * Size];
        cMatrix = new double[Size * Size];

        if (ProcRank == 0)
        {

            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++)
                {
                    aMatrix[i * Size + j] = (double)rand() / RAND_MAX * 9;
                    // ������� ������� ��������� �������
                    // ��� ���������� ��������� �����: i + j/10.;
                    bMatrix[i * Size + j] = -((double)rand() / RAND_MAX * 9);
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

        print_debug_message("Process init completed", ProcRank);
    }

    void DataDistribution(double*& aMatrix, double*& bMatrix, double*& pTemporaryAblock, double*& pBblock, int& Size, int& BlockSize)
    {
        MPI_Type_vector(BlockSize, BlockSize, Size, MPI_DOUBLE, &MPI_BLOCK);

        print_debug_message("DataDistribution MPI_Type_vector completed", ProcRank);

        MPI_Type_commit(&MPI_BLOCK);

        print_debug_message("DataDistribution MPI_Type_commit completed", ProcRank);

        auto temp_send_counts = make_unique<int[]>(ProcNum);
        auto temp_send_shifts = make_unique<int[]>(ProcNum);

        for (int rank = 0; rank < ProcNum; rank++) {
            int c[2];
            MPI_Cart_coords(Communicator, rank, 2, c);
            temp_send_counts[rank] = 1;
            temp_send_shifts[rank] = c[0] * Size * BlockSize + c[1] * BlockSize;
        }

        MPI_Scatterv(aMatrix, temp_send_counts.get(), temp_send_shifts.get(), MPI_BLOCK, aMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, Communicator);
        MPI_Scatterv(bMatrix, temp_send_counts.get(), temp_send_shifts.get(), MPI_BLOCK, bMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, Communicator);

        //if (ProcRank == 0)
        //{
        //    //for (int r = 1; r < ProcNum; r++)
        //    //{
        //    //    MPI_Status s;
        //    //    int c[2];
        //    //    MPI_Cart_coords(GridComm, r, 2, c);
        //    //    MPI_Scatterv()
        //    //    /*MPI_Send(pAMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, GridComm);
        //    //    MPI_Send(pBMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, GridComm);*/
        //    //    MPI_Sendrecv(pAMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, pAMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, GridComm, &s);
        //    //    MPI_Sendrecv(pBMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, pBMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, GridComm, &s);
        //    //}

        //    MPI_Status s;
                //    const auto r = 0;
                //    int c[2];
//    MPI_Cart_coords(Communicator, r, 2, c);
        //    MPI_Sendrecv(pAMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, pAMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, Communicator, &s);
        //    MPI_Sendrecv(pBMatrix + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, r, 0, pBMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, Communicator, &s);
        //}

        /*if (ProcRank != 0)
        {
            MPI_Status s;
            MPI_Recv(aMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, Communicator, &s);
            MPI_Recv(bMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, Communicator, &s);
            cout << "rank: " << ProcRank << endl;
        }*/

        /*for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {
                cout << aMatrix[i + j];
            }
        }
        for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {
                cout << bMatrix[i + j] << endl;
            }
        }*/

        print_debug_message("DataDistribution completed", ProcRank);
    }

    // �������� ������ ������� � �� ������� ������� ��������� 
    void ABlockCommunication(int iter, double* pAblock, double* pMatrixAblock, int BlockSize) {
        // ����������� �������� �������� � ������ ���������� ������� 
        int Pivot = (GridCoords[0] + iter) % GridSize;

        // ����������� ������������� ����� � ��������� ����� ������
        if (GridCoords[1] == Pivot) {
            for (int i = 0; i < BlockSize * BlockSize; i++)
                pAblock[i] = pMatrixAblock[i];
        }
        // �������� �����
        MPI_Bcast(pAblock, BlockSize * BlockSize, MPI_DOUBLE, Pivot, RowComm);
    }

    // ��������� ��������� ������
    void BlockMultiplication(double* pAblock, double* pBblock, double* pCblock, int BlockSize) {
        // ���������� ������������ ��������� ������
        for (int i = 0; i < BlockSize; i++) 
        {
            for (int j = 0; j < BlockSize; j++) 
            {
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
        
        if (GridCoords[0] == GridSize - 1) 
            NextProc = 0;
        
        int PrevProc = GridCoords[0] - 1;

        if (GridCoords[0] == 0) 
            PrevProc = GridSize - 1;
        
        cout << "send to " << NextProc << " from " << PrevProc << endl;
        MPI_Sendrecv_replace(pBblock, BlockSize * BlockSize, MPI_DOUBLE, NextProc, 0, PrevProc, 0, ColComm, &Status);
        print_debug_message("shift", ProcRank);
    }

    // ������� ��� ������������� ��������� ������
    void ParallelResultCalculation(double* pAblock, double* pMatrixAblock, double* pBblock, double* pCblock, int BlockSize) {
        for (int iter = 0; iter < GridSize; iter++) {
            // �������� ������ ������� A �� ������� ���������� �������
            ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
            // ��������� ������
            BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
            // ����������� ����� ������ ������� B � �������� ���������� 
            // �������
            BblockCommunication(pBblock, BlockSize);
            print_debug_message("calculated", ProcRank);
        }
    }

    void ResultCollection(double*& cMatrix, double*& cBlock, int& Size, int& BlockSize)
    {
        /*for (int i = 0; i < BlockSize; i++)
        {
            for (int j = 0; j < BlockSize; j++) {}
        }*/
        MPI_Send(cMatrix, BlockSize * BlockSize, MPI_DOUBLE, 0, 0, Communicator);
        if (ProcRank == 0)
        {
            MPI_Status s;
            for (int rank = 1; rank < ProcNum; rank++)
            {
                int c[2];
                MPI_Cart_coords(Communicator, rank, 2, c);
                MPI_Recv(cBlock + c[0] * Size * BlockSize + c[1] * BlockSize, 1, MPI_BLOCK, rank, 0, Communicator, &s);
            }
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++) {
                    cout << cBlock[i + j] << endl;
                }
            }
        }
    }


    void demo_function(int argc, char* argv[]) {
        double* pAMatrix; 	// ������ �������� ���������� ���������
        double* pBMatrix; 	// ������ �������� ���������� ���������
        double* pCMatrix; 	// �������������� �������
        int Size;        	// ������ ������
        int BlockSize;   	// ������ ��������� ������, ������������� �� ���������
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

            // �������� ����������� ������� ��������� � �������������� ����� � ��������
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