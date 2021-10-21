# CMakeList.txt: проект CMake для MPILabs; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
cmake_minimum_required (VERSION 3.8)

# Добавьте источник в исполняемый файл этого проекта.
add_executable (MPILabs "MPILabs.cpp" "MPILabs.h" "Algorythm/Linear.h" "Algorythm/Linear.cpp" "Algorythm/Fox.h" "Algorythm/Fox.cpp")

# TODO: Добавьте тесты и целевые объекты, если это необходимо.
target_include_directories(MPILabs PRIVATE "Algorythm")

target_include_directories (MPILabs PUBLIC 
	"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
	"C:/Program Files (x86)/Microsoft SDKs/MPI/Include/x64"
)

target_link_directories(MPILabs PUBLIC
	"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64"
)

target_link_libraries(MPILabs PUBLIC 
	msmpi
	msmpifec
	msmpifmc
)