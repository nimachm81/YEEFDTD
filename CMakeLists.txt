cmake_minimum_required(VERSION 3.1.0)
project (YEEFDTD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
find_package(MPI REQUIRED)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

include_directories(include test ${MPI_INCLUDE_PATH})

set(PROJECT_LINK_LIBS ${PROJECT_LINK_LIBS} "libmyBlas.a" )
link_directories("lib" "/usr/lib")

file(GLOB SOURCES "src/*.cpp" "test/*.cpp")

#${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} hello

add_executable(main ${SOURCES})

target_link_libraries(main ${PROJECT_LINK_LIBS} ${MPI_CXX_LIBRARIES} -lpthread)

