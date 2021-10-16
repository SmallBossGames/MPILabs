#pragma once

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>
#include <valarray>

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

void LogVector(const unique_ptr<int[]>& vector, const size_t m)
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

void broadcast_matrix(const size_t m, const size_t n, unique_ptr<double[]>& matrix, int rank)
{
    auto full_length = m * n;

    if (rank != 0)
    {
        matrix = move(make_unique<double[]>(full_length));
    }

    MPI_Bcast(matrix.get(), full_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void get_send_counts(const size_t m, const size_t n, const size_t processes_count, unique_ptr<int[]>& send_counts, unique_ptr<int[]>& displacements)
{
    const auto max_rows_for_process = m / processes_count;

    auto temp_send_counts = make_unique<int[]>(processes_count);

    auto temp_displacements = make_unique<int[]>(processes_count);

    auto sentItems = 0;

    for (size_t i = 0; i < processes_count; i++)
    {
        auto offset = i + 1 == processes_count ? m - sentItems : max_rows_for_process;

        temp_send_counts[i] = offset * n;

        temp_displacements[i] = sentItems * n;

        sentItems += offset;
    }

    send_counts = move(temp_send_counts);

    displacements = move(temp_displacements);
}

int main()
{
    MPI_Init(NULL, NULL);

    int rank;

    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto [m, n] = read_matrix_sizes(rank);

    auto first_matrix = init_matrix(m, n, rank);

    auto second_matrix = init_matrix(m, n, rank);

    if (rank == 0)
    {
        cout << "First matrix: " << endl;
        LogMatrix(first_matrix, m, n);

        cout << "Second matrix: " << endl;
        LogMatrix(second_matrix, m, n);
    }

    auto start_execution_timestamp = MPI_Wtime();

    broadcast_matrix(m, n, second_matrix, rank);

    unique_ptr<int[]> send_counts, displacements;
    get_send_counts(m, n, size, send_counts, displacements);

    auto matrix_part_buffer = make_unique<double_t[]>(send_counts[rank]);

    MPI_Scatterv(first_matrix.get(), send_counts.get(), displacements.get(), MPI_DOUBLE, matrix_part_buffer.get(), send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto process_result_vector = make_unique<double_t[]>(send_counts[rank]);

    for (size_t i = 0; i < send_counts[rank]; i++)
    {
        process_result_vector[i] = 0.0;
    }

    for (size_t i = 0; i < send_counts[rank]/n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            auto sum = 0.0;

            for (size_t k = 0; k < n; k++)
            {
                sum += matrix_part_buffer[i * n + k] * second_matrix[k * m + j];
            }

            process_result_vector[i * n + j] = sum;
        }
    }

    unique_ptr<double_t[]> result_matrix;

    if (rank == 0)
    {
        result_matrix = move(make_unique<double_t[]>(m * n));
    }

    MPI_Gatherv(process_result_vector.get(), send_counts[rank], MPI_DOUBLE, result_matrix.get(), send_counts.get(), displacements.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto end_execution_timestamp = MPI_Wtime();

    if (rank == 0)
    {
        cout << "Result vector: " << endl;

        LogMatrix(result_matrix, m, n);

        cout << "Processes: " << size << endl;
        cout << "Sizes: " << m << '*' << n << endl;
        cout << "Execution time: " << (end_execution_timestamp - start_execution_timestamp) << endl;
    }

    MPI_Finalize();
}