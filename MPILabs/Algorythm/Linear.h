#pragma once

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

namespace mpi_labs::algorythms::linear
{
    using namespace std;

    void LogMatrix(const unique_ptr<double_t[]>& matrix, const size_t m, const size_t n);

    void LogVector(const unique_ptr<double_t[]>& vector, const size_t m);

    void LogVector(const unique_ptr<int[]>& vector, const size_t m);

    tuple<size_t, size_t> read_matrix_sizes(const int rank);

    unique_ptr<double[]> create_random_vector(const size_t m);

    unique_ptr<double[]> init_matrix(const size_t m, const size_t n, const int rank);

    void broadcast_matrix(const size_t m, const size_t n, unique_ptr<double[]>& matrix, int rank);

    void get_send_counts(const size_t m, const size_t n, const size_t processes_count, unique_ptr<int[]>& send_counts, unique_ptr<int[]>& displacements);

    int demo_function();
}