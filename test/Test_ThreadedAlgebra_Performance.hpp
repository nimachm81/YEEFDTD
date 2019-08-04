

#include <cmath>
#include <iostream>
#include <any>
#include <typeinfo>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"
#include "ThreadedAlgebra.h"

void Test_threaded_algebra_performance() {

    std::cout << "Test: " << std::endl;

    srand (static_cast <unsigned> (time(0)));

    std::array<std::size_t, 3> shape{100,200,200};

    NumberArray3D<FPNumber> A(shape, 0.0);
    NumberArray3D<FPNumber> B(shape, 0.0);
    NumberArray3D<FPNumber> C(shape, 0.0);

    FPNumber*** x_A = A.GetArrayData();
    FPNumber*** x_B = B.GetArrayData();
    FPNumber*** x_C = C.GetArrayData();

    for(std::size_t i = 0; i < shape[0]; ++i) {
        for(std::size_t j = 0; j < shape[1]; ++j) {
            for(std::size_t k = 0; k < shape[2]; ++k) {
                x_A[i][j][k] = static_cast <FPNumber> (rand()) / static_cast <FPNumber> (RAND_MAX);
                x_B[i][j][k]  = static_cast <FPNumber> (rand()) / static_cast <FPNumber> (RAND_MAX);
                x_C[i][j][k]  = static_cast <FPNumber> (rand()) / static_cast <FPNumber> (RAND_MAX);
            }
        }
    }

    ThreadedAlgebra ta(4);

    FPNumber a = 0.05;
    int N = 5;

    auto t_start_0 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < N; ++i) {
        ta.Get_y_pe_ax(A, a, B);
    }
    auto t_end_0 = std::chrono::high_resolution_clock::now();
    auto duration_0 = std::chrono::duration_cast<std::chrono::microseconds>(t_end_0 - t_start_0);
    std::cout << "Threaded duration (ms) :    " << duration_0.count() << std::endl;


    auto t_start_1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < N; ++i) {
        A.Add_aX(a, B);
    }
    auto t_end_1 = std::chrono::high_resolution_clock::now();
    auto duration_1 = std::chrono::duration_cast<std::chrono::microseconds>(t_end_1 - t_start_1);
    std::cout << "NonThreaded duration (ms) : "  << duration_1.count() << std::endl;

    std::cout << "speed up: " << (float)duration_1.count()/duration_0.count() << std::endl;

}
