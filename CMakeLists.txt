cmake_minimum_required(VERSION 3.1.0)
project (YEEFDTD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
find_package(MPI REQUIRED)

set(CMAKE_BUILD_TYPE Debug)

include_directories(include ${MPI_INCLUDE_PATH})

set (PROJECT_LINK_LIBS "${PROJECT_LINK_LIBS} " )
link_directories(/usr/lib)

file(GLOB SOURCES "src/*.cpp")

#${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} hello

add_executable(main ${SOURCES})

target_link_libraries(main ${MPI_CXX_LIBRARIES})

