
#include <cmath>
#include <iostream>
#include <any>
#include <typeinfo>
#include <stdlib.h>
#include <time.h>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"
#include "ThreadedAlgebra.h"

void Test_threaded_algebra() {

    std::cout << "Test: " << std::endl;

    srand (static_cast <unsigned> (time(0)));

    std::array<std::size_t, 3> shape{50,40,30};

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

    std::array<std::size_t, 3> sliceShape{4, 5, 3};
    std::array<std::size_t, 3> ind_st_A1{1, 3, 4};
    std::array<std::size_t, 3> stride_A1{2, 3, 2};
    std::array<std::size_t, 3> ind_end_A1{ind_st_A1[0] + stride_A1[0]*sliceShape[0],
                                          ind_st_A1[1] + stride_A1[1]*sliceShape[1],
                                          ind_st_A1[2] + stride_A1[2]*sliceShape[2]
                                       };
    NumberArray3D<FPNumber> A1 = A.GetSlice(ind_st_A1, ind_end_A1, stride_A1);

    std::array<std::size_t, 3> ind_st_B1{1, 3, 4};
    std::array<std::size_t, 3> stride_B1{1, 2, 2};
    std::array<std::size_t, 3> ind_end_B1{ind_st_B1[0] + stride_B1[0]*sliceShape[0],
                                          ind_st_B1[1] + stride_B1[1]*sliceShape[1],
                                          ind_st_B1[2] + stride_B1[2]*sliceShape[2]
                                          };
    NumberArray3D<FPNumber> B1 = B.GetSlice(ind_st_B1, ind_end_B1, stride_B1);

    std::array<std::size_t, 3> ind_st_C1{2, 4, 3};
    std::array<std::size_t, 3> stride_C1{2, 3, 2};
    std::array<std::size_t, 3> ind_end_C1{ind_st_C1[0] + stride_C1[0]*sliceShape[0],
                                          ind_st_C1[1] + stride_C1[1]*sliceShape[1],
                                          ind_st_C1[2] + stride_C1[2]*sliceShape[2]
                                          };
    NumberArray3D<FPNumber> C1 = C.GetSlice(ind_st_C1, ind_end_C1, stride_C1);

    NumberArray3D<FPNumber> A2 = A1;
    FPNumber a = static_cast <FPNumber> (rand()) / static_cast <FPNumber> (RAND_MAX);

    ThreadedAlgebra ta(4);

    std::cout << "-----------------------------------------------------------" << std::endl;
    ta.Get_y_pe_ax(A1, a, B1);
    A1.Print();
    A2.Add_aX(a, B1);
    A2.Print();
    (A1 - A2).Print();


};

