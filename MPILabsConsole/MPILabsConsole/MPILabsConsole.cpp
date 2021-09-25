
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>

using namespace std;

void InitRandomVector(double_t* vector, size_t m)
{
    for (int i = 0; i < m; i++)
    {
        vector[i] = rand() % 100;
    }
}

void LogMatrix(double_t* matrix, size_t m, size_t n)
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

    cout << sb;
}

void LogVector(double_t* vector, size_t m)
{
    string sb;

    for (size_t i = 0; i < m; i++)
    {
        sb.append(to_string(vector[i])).append(" ");
    }

    sb.append("\n");

    cout << sb;
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

    unique_ptr<double_t[]> matrix(nullptr);

    unique_ptr<double_t[]> vector(new double[m]);

    if (rank == 0)
    {
        InitRandomVector(vector.get(), m);
    }

    if (rank == 0) 
    {
        cout << "Input vector: " << endl;
        LogVector(vector.get(), m);
    }

    MPI_Bcast(vector.get(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    unique_ptr<int[]> val_per_process(new int[size]);

    unique_ptr<int[]> dispacepents(new int[size]);

    if (rank == 0)
    {
        matrix.reset(new double[m * n]);

        InitRandomVector(matrix.get(), m * n);

        cout << "Input matrix: " << endl;

        LogMatrix(matrix.get(), m, n);
    }

    const auto max_rows_for_process = m / size;

    auto rowsLeft = m;

    for (size_t i = 0; i < size; i++)
    {
        auto rows_for_process = i + 1 == size ? rowsLeft : max_rows_for_process;

        val_per_process[i] = rows_for_process * n;

        dispacepents[i] = (m - rowsLeft) * n;

        rowsLeft -= rows_for_process;
    }

    unique_ptr<double_t[]> matrix_part_buffer(new double_t[val_per_process[rank]]);

    MPI_Scatterv(matrix.get(), val_per_process.get(), dispacepents.get(), MPI_DOUBLE, matrix_part_buffer.get(), val_per_process[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    unique_ptr<double_t[]> result_vector(nullptr);

    if (rank == 0)
    {
        result_vector.reset(new double_t[m]);
    }

    const size_t process_result_vector_size = val_per_process[rank] / n;

    unique_ptr<double_t[]> process_result_vector(new double_t[process_result_vector_size]);

    for (size_t i = 0; i < process_result_vector_size; i++)
    {
        process_result_vector[i] = 0.0;

        for (size_t j = 0; j < n; j++)
        {
            process_result_vector[i] += vector[j] * matrix_part_buffer[i * n + j];
        }
    }

    unique_ptr<int[]> result_val_per_process(new int[size]);

    unique_ptr<int[]> result_dispacepents(new int[size]);

    auto resultRowsLeft = m;

    for (size_t i = 0; i < size; i++)
    {
        auto rows_for_process = i + 1 == size ? resultRowsLeft : max_rows_for_process;

        result_val_per_process[i] = rows_for_process;

        result_dispacepents[i] = m - resultRowsLeft;

        resultRowsLeft -= rows_for_process;
    }

    MPI_Gatherv(process_result_vector.get(), process_result_vector_size, MPI_DOUBLE, result_vector.get(), result_val_per_process.get(), result_dispacepents.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << "Output matrix" << endl;

        LogVector(result_vector.get(), m);
    }

    MPI_Finalize();
}
