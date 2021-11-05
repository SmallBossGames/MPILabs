#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

#include <Linear.h>
#include <Fox.h>
#include <Cannon.h>

void main(int argc, char* argv[]) 
{
	mpi_labs::algorythms::cannon::demo_program(argc, argv);
}
