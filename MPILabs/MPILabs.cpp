#pragma once

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>
#include <valarray>
#include <Windows.h>

using namespace std;
// ЛЕНТОЧНЫЙ АЛГОРИТМ
// 
//void LogMatrix(const unique_ptr<double_t[]> &matrix, const size_t m, const size_t n)
//{
//    string sb;
//
//    for (size_t i = 0; i < m; i++)
//    {
//        for (size_t j = 0; j < n; j++)
//        {
//            sb.append(to_string(matrix[i * n + j])).append(" ");
//        }
//
//        sb.append("\n");
//    }
//
//    sb.append("\n");
//
//    cout << sb;
//}
//
//
//void LogVector(const unique_ptr<double_t[]> &vector, const size_t m)
//{
//    string sb;
//
//    for (size_t i = 0; i < m; i++)
//    {
//        sb.append(to_string(vector[i])).append(" ");
//    }
//
//    sb.append("\n\n");
//
//    cout << sb;
//}
//
//void LogVector(const unique_ptr<int[]>& vector, const size_t m)
//{
//    string sb;
//
//    for (size_t i = 0; i < m; i++)
//    {
//        sb.append(to_string(vector[i])).append(" ");
//    }
//
//    sb.append("\n\n");
//
//    cout << sb;
//}
//
//tuple<size_t, size_t> read_matrix_sizes(const int rank) 
//{
//    size_t m = 0, n = 0;
//
//    if (rank == 0)
//    {
//        cout << "Write m:" << endl;
//        cin >> m;
//
//        cout << endl;
//
//        cout << "Write n:" << endl;
//        cin >> n;
//
//        cout << endl;
//    }
//
//    MPI_Bcast(&m, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
//
//    MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
//
//    return make_tuple(m, n);
//}
//
//unique_ptr<double[]> create_random_vector(const size_t m)
//{
//    auto vector = make_unique<double[]>(m);
//
//    for (int i = 0; i < m; i++)
//    {
//        vector[i] = rand() % 100;
//    }
//
//    return vector;
//}
//
//unique_ptr<double[]> init_matrix(const size_t m, const size_t n, const int rank)
//{
//    if (rank == 0)
//    {
//        return create_random_vector(m * n);
//    }
//
//    return unique_ptr<double[]>(nullptr);
//}
//
//void broadcast_matrix(const size_t m, const size_t n, unique_ptr<double[]>& matrix, int rank)
//{
//    auto full_length = m * n;
//
//    if (rank != 0)
//    {
//        matrix = move(make_unique<double[]>(full_length));
//    }
//
//    MPI_Bcast(matrix.get(), full_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//}
//
//void get_send_counts(const size_t m, const size_t n, const size_t processes_count, unique_ptr<int[]>& send_counts, unique_ptr<int[]>& displacements)
//{
//    const auto max_rows_for_process = m / processes_count;
//
//    auto temp_send_counts = make_unique<int[]>(processes_count);
//
//    auto temp_displacements = make_unique<int[]>(processes_count);
//
//    auto sentItems = 0;
//
//    for (size_t i = 0; i < processes_count; i++)
//    {
//        auto offset = i + 1 == processes_count ? m - sentItems : max_rows_for_process;
//
//        temp_send_counts[i] = offset * n;
//
//        temp_displacements[i] = sentItems * n;
//
//        sentItems += offset;
//    }
//
//    send_counts = move(temp_send_counts);
//
//    displacements = move(temp_displacements);
//}
//
//int main()
//{
//    MPI_Init(NULL, NULL);
//
//    int rank;
//
//    int size;
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    auto [m, n] = read_matrix_sizes(rank);
//
//    auto first_matrix = init_matrix(m, n, rank);
//
//    auto second_matrix = init_matrix(m, n, rank);
//
//    if (rank == 0)
//    {
//        cout << "First matrix: " << endl;
//        LogMatrix(first_matrix, m, n);
//
//        cout << "Second matrix: " << endl;
//        LogMatrix(second_matrix, m, n);
//    }
//
//    auto start_execution_timestamp = MPI_Wtime();
//
//    broadcast_matrix(m, n, second_matrix, rank);
//
//    unique_ptr<int[]> send_counts, displacements;
//    get_send_counts(m, n, size, send_counts, displacements);
//
//    auto matrix_part_buffer = make_unique<double_t[]>(send_counts[rank]);
//
//    MPI_Scatterv(first_matrix.get(), send_counts.get(), displacements.get(), MPI_DOUBLE, matrix_part_buffer.get(), send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    auto process_result_vector = make_unique<double_t[]>(send_counts[rank]);
//
//    for (size_t i = 0; i < send_counts[rank]; i++)
//    {
//        process_result_vector[i] = 0.0;
//    }
//
//    for (size_t i = 0; i < send_counts[rank]/n; i++)
//    {
//        for (size_t j = 0; j < m; j++)
//        {
//            auto sum = 0.0;
//
//            for (size_t k = 0; k < n; k++)
//            {
//                sum += matrix_part_buffer[i * n + k] * second_matrix[k * m + j];
//            }
//
//            process_result_vector[i * n + j] = sum;
//        }
//    }
//
//    unique_ptr<double_t[]> result_matrix;
//
//    if (rank == 0)
//    {
//        result_matrix = move(make_unique<double_t[]>(m * n));
//    }
//
//    MPI_Gatherv(process_result_vector.get(), send_counts[rank], MPI_DOUBLE, result_matrix.get(), send_counts.get(), displacements.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    auto end_execution_timestamp = MPI_Wtime();
//
//    if (rank == 0)
//    {
//        cout << "Result vector: " << endl;
//
//        LogMatrix(result_matrix, m, n);
//
//        cout << "Processes: " << size << endl;
//        cout << "Sizes: " << m << '*' << n << endl;
//        cout << "Execution time: " << (end_execution_timestamp - start_execution_timestamp) << endl;
//    }
//
//    MPI_Finalize();
//}

// Программа 7.1
// Алгоритм Фокса умножения матриц – блочное представление данных
//   Условия выполнения программы: все матрицы квадратные, 
//   размер блоков и их количество по горизонтали и вертикали
//   одинаково, процессы образуют квадратную решетку
int ProcNum = 0;   	// Количество доступных процессов 
int ProcRank = 0;  	// Ранг текущего процесса
int GridSize;      	// Размер виртуальной решетки процессов
int GridCoords[2];  	// Координаты текущего процесса в процессной
// решетке
MPI_Comm GridComm;  	// Коммуникатор в виде квадратной решетки
MPI_Comm ColComm;   	// коммуникатор – столбец решетки
MPI_Comm RowComm;   	// коммуникатор – строка решетки
MPI_Datatype MPI_BLOCK;

// Создание коммуникатора в виде двумерной квадратной решетки 
// и коммуникаторов для каждой строки и каждого столбца решетки
void CreateGridCommunicators() {
    int DimSize[2];  	// Количество процессов в каждом измерении решетки
    int Periodic[2];	// =1 для каждого измерения, являющегося периодическим 
    int Subdims[2];  	// =1 для каждого измерения, оставляемого в подрешетке
    DimSize[0] = GridSize;
    DimSize[1] = GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;

    // Создание коммуникатора в виде квадратной решетки 
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);

    // Определение координат процесса в решетке 
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);

    // Создание коммуникаторов для строк процессной решетки
    Subdims[0] = 0;  // Фиксация измерения
    Subdims[1] = 1;  // Наличие данного измерения в подрешетке
    MPI_Cart_sub(GridComm, Subdims, &RowComm);

    // Создание коммуникаторов для столбцов процессной решетки
    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

// Функция для выделения памяти и инициализации исходных данных
void ProcessInitialization(double*& pAMatrix, double*& pBMatrix, double*& pCMatrix, double*& pAblock, double*& pBblock, double*& pCblock, double*& pTemporaryAblock, int& Size, int& BlockSize) 
{
    if (ProcRank == 0) {
        do {
            printf("\nEnter matrix size: ");
            scanf("%d", &Size);

            if (Size % GridSize != 0) {
                printf("Размер матриц должен быть кратен размеру сетки! \n");
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
                // вариант задания элементов матрицы
                // без генератора случайных чисел: i + j/10.;
                pBMatrix[i * Size + j] = -((double)rand() / RAND_MAX * 9);
                // вариант задания элементов матрицы
                // без генератора случайных чисел: -(i + j/10.);
            }
        }
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++) {
                // если хочется вывести матрицу, то делать тут
            }
        }
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++) {
                // если хочется вывести матрицу, то делать тут
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

// Рассылка блоков матрицы А по строкам решетки процессов 
void ABlockCommunication(int iter, double* pAblock,
    double* pMatrixAblock, int BlockSize) {

    // Определение ведущего процесса в строке процессной решетки 
    int Pivot = (GridCoords[0] + iter) % GridSize;

    // Копирование передаваемого блока в отдельный буфер памяти
    if (GridCoords[1] == Pivot) {
        for (int i = 0; i < BlockSize * BlockSize; i++)
            pAblock[i] = pMatrixAblock[i];
    }

    // Рассылка блока
    MPI_Bcast(pAblock, BlockSize * BlockSize, MPI_DOUBLE, Pivot,
        RowComm);
}

// Умножение матричных блоков
void BlockMultiplication(double* pAblock, double* pBblock,
    double* pCblock, int BlockSize) {
    // Вычисление произведения матричных блоков
    for (int i = 0; i < BlockSize; i++) {
        for (int j = 0; j < BlockSize; j++) {
            double temp = 0;
            for (int k = 0; k < BlockSize; k++)
                temp += pAblock[i * BlockSize + k] * pBblock[k * BlockSize + j];
            pCblock[i * BlockSize + j] += temp;
        }
    }
}

// Циклический сдвиг блоков матрицы В вдоль столбца процессной 
// решетки 
void BblockCommunication(double* pBblock, int BlockSize) {
    MPI_Status Status;
    int NextProc = GridCoords[0] + 1;
    if (GridCoords[0] == GridSize - 1) NextProc = 0;
    int PrevProc = GridCoords[0] - 1;
    if (GridCoords[0] == 0) PrevProc = GridSize - 1;
    MPI_Sendrecv_replace(pBblock, BlockSize * BlockSize, MPI_DOUBLE,
        NextProc, 0, PrevProc, 0, ColComm, &Status);
}

// Функция для параллельного умножения матриц
void ParallelResultCalculation(double* pAblock, double* pMatrixAblock,
    double* pBblock, double* pCblock, int BlockSize) {
    for (int iter = 0; iter < GridSize; iter++) {
        // Рассылка блоков матрицы A по строкам процессной решетки
        ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
        // Умножение блоков
        BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
        // Циклический сдвиг блоков матрицы B в столбцах процессной 
        // решетки
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


void main(int argc, char* argv[]) {
    double* pAMatrix; 	// Первый аргумент матричного умножения
    double* pBMatrix; 	// Второй аргумент матричного умножения
    double* pCMatrix; 	// Результирующая матрица
    int Size;        	// Размер матриц
    int BlockSize;   	// Размер матричных блоков, расположенных 
                      // на процессах
    double* pAblock;  	// Блок матрицы А на процессе
    double* pBblock;  	// Блок матрицы В на процессе
    double* pCblock;  	// Блок результирующей матрицы С на процессе
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

        // Создание виртуальной решетки процессов и коммуникаторов 
        // строк и столбцов
        CreateGridCommunicators();

        // Выделение памяти и инициализация элементов матриц
        ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock, pCblock, pMatrixAblock, Size, BlockSize);
        // Блочное распределение матриц между процессами
        DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size, BlockSize);

        // Выполнение параллельного метода Фокса
        ParallelResultCalculation(pAblock, pMatrixAblock, pBblock, pCblock, BlockSize);

        // Сбор результирующей матрицы на ведущем процессе
        ResultCollection(pCMatrix, pCblock, Size, BlockSize);

        // Завершение процесса вычислений
        //ProcessTermination(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock, pCblock, pMatrixAblock);
    }

    MPI_Finalize();
}
