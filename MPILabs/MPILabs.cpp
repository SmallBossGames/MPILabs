#pragma once

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

using namespace std;

void LogMatrix(const unique_ptr<double_t[]> &matrix, const size_t m, const size_t n)
{
    string sb;

    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            sb.append(to_string(matrix[i * n + j])).append(" ");
        }

        sb.append("\n");
    }

    sb.append("\n");

    cout << sb;
}

void LogVector(const unique_ptr<double_t[]> &vector, const size_t m)
{
    string sb;

    for (size_t i = 0; i < m; i++)
    {
        sb.append(to_string(vector[i])).append(" ");
    }

    sb.append("\n\n");

    cout << sb;
}

tuple<size_t, size_t> read_matrix_sizes(const int rank) 
{
    size_t m = 0, n = 0;

    if (rank == 0)
    {
        cout << "Write m:" << endl;
        cin >> m;

        cout << endl;

        cout << "Write n:" << endl;
        cin >> n;

        cout << endl;
    }

    MPI_Bcast(&m, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    return make_tuple(m, n);
}

unique_ptr<double[]> create_random_vector(const size_t m)
{
    auto vector = make_unique<double[]>(m);

    for (int i = 0; i < m; i++)
    {
        vector[i] = rand() % 100;
    }

    return vector;
}

unique_ptr<double[]> init_matrix(const size_t m, const size_t n, const int rank)
{
    if (rank == 0)
    {
        return create_random_vector(m * n);
    }

    return unique_ptr<double[]>(nullptr);
}

unique_ptr<double[]> init_vector(const size_t n, const int rank)
{
    auto vector = rank == 0 ? create_random_vector(n) : make_unique<double[]>(n);

    MPI_Bcast(vector.get(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return vector;
}

int main()
{
    MPI_Init(NULL, NULL);

    int rank;

    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto [m, n] = read_matrix_sizes(rank);

    auto matrix = init_matrix(m, n, rank);

    auto vector = init_vector(m, rank);

    if (rank == 0)
    {
        cout << "Input vector: " << endl;
        LogVector(vector, m);
    }

    if (rank == 0)
    {
        cout << "Input matrix: " << endl;

        LogMatrix(matrix, m, n);
    }

    auto start_execution_timestamp = MPI_Wtime();

    auto val_per_process = make_unique<int[]>(size);

    auto displacepents = make_unique<int[]>(size);

    const auto max_rows_for_process = m / size;

    auto sentItems = 0;

    for (size_t i = 0; i < size; i++)
    {
        auto offset = min(max_rows_for_process, m - sentItems);

        val_per_process[i] = offset * n;

        displacepents[i] = sentItems * n;

        sentItems += offset;
    }

    auto matrix_part_buffer = make_unique<double_t[]>(val_per_process[rank]);

    MPI_Scatterv(matrix.get(), val_per_process.get(), displacepents.get(), MPI_DOUBLE, matrix_part_buffer.get(), val_per_process[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    unique_ptr<double_t[]> result_vector(nullptr);

    if (rank == 0)
    {
        result_vector.reset(new double_t[m]);
    }

    const size_t process_result_vector_size = val_per_process[rank] / n;

    auto process_result_vector = make_unique<double_t[]>(process_result_vector_size);

    for (size_t i = 0; i < process_result_vector_size; i++)
    {
        process_result_vector[i] = 0.0;

        for (size_t j = 0; j < n; j++)
        {
            process_result_vector[i] += vector[j] * matrix_part_buffer[i * n + j];
        }
    }

    auto result_val_per_process = make_unique<int[]>(size);

    auto result_dispacepents = make_unique<int[]>(size);

    auto resultRowsLeft = m;

    for (size_t i = 0; i < size; i++)
    {
        auto rows_for_process = i + 1 == size ? resultRowsLeft : max_rows_for_process;

        result_val_per_process[i] = rows_for_process;

        result_dispacepents[i] = m - resultRowsLeft;

        resultRowsLeft -= rows_for_process;
    }

    MPI_Gatherv(process_result_vector.get(), process_result_vector_size, MPI_DOUBLE, result_vector.get(), result_val_per_process.get(), result_dispacepents.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto end_execution_timestamp = MPI_Wtime();

    if (rank == 0)
    {
        cout << "Result vector: " << endl;

        LogVector(result_vector, m);

        cout << "Processes: " << size << endl;
        cout << "Sizes: " << m << '*' << n << endl;
        cout << "Execution time: " << (end_execution_timestamp - start_execution_timestamp) << endl;
    }

    MPI_Finalize();
}