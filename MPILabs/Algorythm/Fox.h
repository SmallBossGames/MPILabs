#pragma once

namespace mpi_labs::algorythms::fox
{
    void CreateGridCommunicators();

    void ProcessInitialization(double*& pAMatrix, double*& pBMatrix, double*& pCMatrix, double*& pAblock, double*& pBblock, double*& pCblock, double*& pTemporaryAblock, int& Size, int& BlockSize);

    void DataDistribution(double*& pAMatrix, double*& pBMatrix, double*& pTemporaryAblock, double*& pBblock, int& Size, int& BlockSize);

    void ABlockCommunication(int iter, double* pAblock, double* pMatrixAblock, int BlockSize);

    void BlockMultiplication(double* pAblock, double* pBblock, double* pCblock, int BlockSize);
 
    void BblockCommunication(double* pBblock, int BlockSize);

    void ParallelResultCalculation(double* pAblock, double* pMatrixAblock, double* pBblock, double* pCblock, int BlockSize);

    void ResultCollection(double*& pCMatrix, double*& pCblock, int& Size, int& BlockSize);

    void demo_function(int argc, char* argv[]);
}