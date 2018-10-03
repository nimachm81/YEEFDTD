 //


#include <iostream>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"


#include <any>
#include <typeinfo>
void test_multidim_array() {
    std::cout << "Test." << std::endl;

    std::array<std::size_t, 3> shape{2,3,4};
    NumberArray3D<RealNumber> A(shape, 2.0);
    std::cout << " A : " << &A << std::endl;
    A.Print();

    NumberArray3D<RealNumber> B(shape, 3.0);
    std::cout << " B : " << &B << std::endl;
    B.Print();

    A += B;
    A.Print();

    A /= B;
    A.Print();

    A *= 2.0;
    A.Print();

    A.GetSlice({0, 0, 1}, {2, 2, 2}).Print();

    NumberArray3D<RealNumber> C = (A + B);
    std::cout << " C : " << &C << std::endl;
    C.Print();

    A = B;
    std::cout << " A : " << &A << std::endl;
    A.Print();

    std::any a = std::make_any<NumberArray3D<RealNumber>>(shape, 4.01);
    std::cout << " a : " << a.type().name() << std::endl;

    std::cout << " a : " << typeid(NumberArray3D<RealNumber>).name() << std::endl;
    std::any_cast<NumberArray3D<RealNumber>>(a).Print();
}

#include <tuple>
void test_void_pointer() {
    std::tuple<int, double> p1(1, 2.0);
    std::cout << std::get<0>(p1) << "  " << std::get<1>(p1) << std::endl;

    void* vp_p1 = &p1;
    std::tuple<int, double>* p2 = static_cast<std::tuple<int, double>*>(vp_p1);
    std::cout << std::get<0>(*p2) << "  " << std::get<1>(*p2) << std::endl;
}

#include "Vector.hpp"
void test_vector3() {
    Vector3<double> v1(1.0, 2.0, 3.0);
    std::cout << "v1: " << v1 << std::endl;
    Vector3<double> v2(4.0, 5.0, 7.5);
    std::cout << "v2: " << v2 << std::endl;

    Vector3<double> v = v1 + v2;
    std::cout << "v1 + v2: " << v << std::endl;
}

#include <any>
#include <tuple>
#include <utility>
#include <vector>
auto* test_container() {
    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>> p0({1,2,3}, {4,5,6});
    std::cout << p0.first[0] << " " << p0.first[1] << " " << p0.first[2] << std::endl;

    std::array<int, 3> arr0{3, 4, 5};
    std::array<int, 3> arr1{7, 8, 9};
    std::pair<std::array<int, 3>, std::array<int, 3>> p1(arr0, arr1);
    std::cout << p1.first[0] << " " << p1.first[1] << " " << p1.first[2] << std::endl;

    std::vector<std::string> v{"this", "is", "a", "string"};
    std::cout << v[0] << v[1] << std::endl;

    std::tuple<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>, std::vector<std::string>> tp(p0, v);
    std::cout << std::get<0>(tp).first[0] << std::endl;

    double* p_d = new double[10];
    p_d[0] = 10.12;
    void* p_v = static_cast<void*>(p_d);
    std::pair<int, void*> p2(2, p_v);
    std::cout << static_cast<double*>(p2.second)[0] << "  " << static_cast<double*>(p2.second)[1] << std::endl;
    delete[] p_d;

    std::any a = p0;
    auto r_p0 = std::any_cast<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>>(a);
    std::cout << "any cast: " << r_p0.first[0] << " " << r_p0.first[1] << " " << r_p0.first[2] << std::endl;


    auto tp2 = new std::tuple<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,
            std::vector<std::string>>(p0, v);
    return tp2;
}

#include "YeeGrid.h"
void test_yeegrid() {
    std::size_t nz = 10;
    std::size_t indJ = 5;
    RealNumber dt = 0.1;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 0.0, 1.0};
    std::array<std::size_t, 3> nCells{1, 0, nz};
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianPointSource("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/, 0.0, 0.0, 0.0);
    void* E_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 1},      // ind_start_A
        {1, 1, nz},   // ind_end_A
        "E",            // arrayA_name
        0,              // arrayA_component
        {1.0, -1.0},    // bValues
        {"H", "H"},     // arrayC_names
        {1, 1},         // arraC_components
        {{0, 0, 1}, {0, 0, 0}}     // arrayC_indsStart
    );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "E",
        0,
        {1.0},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* H_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, 1, nz},
        "H",
        1,
        {1.0, -1.0},
        {"E", "E"},
        {0, 0},
        {{0, 0, 1}, {0, 0, 0}}
    );
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params);
    yee.SetIterationSequence({"J-update", "E-update", "E-J-update", "H-update"});
    yee.ApplyUpdateInstructions(10);
}

#include "GaussianGridArrayManipulator.h"
int main(int argc, char** argv) {
    test_yeegrid();
    //test_multidim_array();
}


