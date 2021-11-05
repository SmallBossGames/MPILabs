#include <iostream>
#include <mpi.h>
#include <string>
#include <sstream>

namespace mpi_labs::algorythms::cannon 
{
	using namespace std;

	int demo_program(int argc, char* argv[])
	{
		const int count_argyment = 4;

		MPI_Init(NULL, NULL);

		int rank;
		int size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		int mod = 0, err = 0;

		size_t m = 0, n = 0;

		if (rank == 0)
		{
			if (argc != count_argyment)
			{
				err = 1;
				cout << "\nWrong arguments!\nPrigram stopped" << endl;
			}
			else
			{
				mod = stoi(argv[1]);
				m = stoi(argv[2]);
				n = stoi(argv[3]);

				if (mod < 1 || mod > 3 || m < 1 || n < 1 || m < size || n < size
					|| (mod == 3 && (m != n || n % size != 0 || (sqrt(size) - (int)sqrt(size) > 0.001) || (sqrt(n) - (int)sqrt(n) > 0.001))))
				{
					err = 1;
					cout << "\nWrong arguments!\nPrigram stopped" << endl;
				}
				else
				{
					switch (mod)
					{
					case 1:
						cout << "\nProgram mod: Matrix-vector multiplication when dividing data by rows" << endl;
						break;
					case 2:
						cout << "\nProgram mod: Matrix-vector multiplication when dividing data by columns" << endl;
						break;
					case 3:
						cout << "\nProgram mod: Matrix-vector multiplication when dividing data by blocks" << endl;
						break;
					}
					cout << "Start matrix size:" << m << "x" << n << endl;
					cout << "Start vector size:" << n << endl;
				}
			}
		}

		MPI_Bcast(&err, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

		if (err == 1)
		{
			MPI_Finalize();
			return 0;
		}

		MPI_Bcast(&m, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&mod, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

		auto A = make_unique<double[]>(m * n);
		auto B = make_unique<double[]>(n);

		if (rank == 0) //инициализация матрицы и вектора
		{
			for (int i = 0; i < m * n; i++)
			{
				A[i] = rand() % 100;
				if (i < n)
				{
					B[i] = rand() % 100;
				}
			}
		}

		MPI_Bcast(B.get(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double start_timestamp = MPI_Wtime();

		if (mod == 1) //ленточное умножение по строкам
		{
			auto sendcounts = make_unique<int[]>(size); //число отправляемых элементов каждому процессу
			auto displs = make_unique<int[]>(size); //смещения в исходном массиве

			const auto step = m / size * n; //шаг смещения в исходном массиве матрицы
			for (int i = 0; i < size; i++)
			{
				displs[i] = step * i;
				sendcounts[i] = (i == size - 1) ? (step + (m % size * n)) : step;
			}

			auto buffer = make_unique<double[]> (sendcounts[rank]);

			MPI_Scatterv(A.get(), sendcounts.get(), displs.get(), MPI_DOUBLE, buffer.get(), sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

			for (int i = 0; i < sendcounts[rank];)
			{
				for (int j = 0; j < n; j++, i++)
				{
					buffer[i] *= B[j]; //операция умножения
				}
			}

			MPI_Gatherv(buffer.get(), sendcounts[rank], MPI_DOUBLE, A.get(), sendcounts.get(), displs.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		else if (mod == 2) //леточное умножение по столбцам
		{
			auto At = make_unique<double[]>(n * m);

			//транспонируем матрицу A
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					At[j + (i * m)] = A[j * n + i];
				}
			}

			auto sendcounts = make_unique<int[]>(size); //число отправляемых элементов каждому процессу
			auto displs = make_unique<int[]>(size); //смещения в исходном массиве

			const auto step = (n / size) * m; //шаг смещения в исходном массиве матрицы
			for (int i = 0; i < size; i++)
			{
				displs[i] = step * i;
				sendcounts[i] = (i == size - 1) ? (step + (n % size * m)) : step;
			}

			auto buffer = make_unique<double[]>(sendcounts[rank]);

			MPI_Scatterv(At.get(), sendcounts.get(), displs.get(), MPI_DOUBLE, buffer.get(), sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

			for (int i = 0; i < sendcounts[rank]; i++)
			{
				buffer[i] *= B[i / m + rank * (n / size)];// операция умножения
			}

			MPI_Gatherv(buffer.get(), sendcounts[rank], MPI_DOUBLE, At.get(), sendcounts.get(), displs.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//транспонируем матрицу A
			for (int i = 0; i < n; i++) 
			{
				for (int j = 0; j < m; j++)
				{
					A[j * n + i] = At[j + (i * m)];
				}
			}
		}
		else if (mod == 3) // блочное умножение
		{
			int h = sqrt(size);
			int step = m / h;
			
			auto buffer = make_unique<double[]>(step * step);

			if (rank == 0)
			{
				for (int k = 1; k < size; k++)
				{
					for (int j = 0; j < step; j++) 
					{
						for (int i = 0; i < step; i++) 
						{
							buffer[i + (j * step)] = A[i + (j * m) + (step * (k % h + m * (k / h)))];
						}
					}
					MPI_Send(buffer.get(), step* step, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
				}
				for (int j = 0; j < step; j++)
				{
					for (int i = 0; i < step; i++)
					{
						buffer[i + (j * step)] = A[i + (j * m)];
					}
				}
			}
			else
			{
				MPI_Status status;
				MPI_Recv(buffer.get(), step * step, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			}
			for (int i = 0; i < step * step;)
			{
				for (int j = 0; j < step; j++, i++)
				{
					buffer[i] *= B[j + step * (rank % h)]; //операция умножения
				}
			}
			if (rank == 0)
			{
				for (int j = 0; j < step; j++) 
				{
					for (int i = 0; i < step; i++)
					{
						A[i + (j * m)] = buffer[i + (j * step)];
					}
				}
					
				MPI_Status status;
				for (int k = 1; k < size; k++)
				{
					MPI_Recv(buffer.get(), step* step, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
					for (int j = 0; j < step; j++)
					{
						for (int i = 0; i < step; i++) {
							A[i + (j * m) + (step * (k % h + m * (k / h)))] = buffer[i + (j * step)];
						}
					}
				}
			}
			else
			{
				MPI_Send(buffer.get(), step* step, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			}
		}
		double end_timestamp = MPI_Wtime();
		if (rank == 0)
		{
			cout << "\nNumper of processes: " << size << endl;
			cout << "Result Matrix Size: " << m << '*' << n << endl;
			cout << "Time: " << (end_timestamp - start_timestamp) << endl;
			/*/ вывод результата:
			cout << "\nResult:" << endl;
			for (int j = 0; j < m; j++)
			{
			for (int i = 0; i < n; i++)
			cout << A[i+j*n] << " ";
			cout << endl;
			17
			}
			//*/
		}
		MPI_Finalize();
	}
}