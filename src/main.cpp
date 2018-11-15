 //


#include <cmath>
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

#include <fstream>
#include "MultiDimArrayFileIO.hpp"
void test_file_write() {
    std::ofstream myfile;
    myfile.open("example.bin", std::ios::out | std::ios::binary);
    if(myfile.is_open()) {
        int* i = new int[3];
        i[0] = 1; i[1] = 13;
        std::array<std::size_t, 3> shape{3,4,5};
        NumberArray3D<RealNumber> A(shape, 2.01);
        NumberArray3D<RealNumber> B = A.GetSlice({1, 1, 1}, {2, 3, 4});
        double* d = new double[4];
        d[3] = 32.4;

        myfile.write((char*)i, sizeof(int)*3);
        B.WriteArrayDataToFile(&myfile);
        myfile.write((char*)d, sizeof(double)*4);

        myfile.close();
        delete[] i;
        delete[] d;
    }
}

void test_file_read() {
    std::streampos objsize;
    std::ifstream myfile;
    myfile.open("example.bin", std::ios::in | std::ios::binary);
    if(myfile.is_open()) {
        int* i = new int[3];
        double* d = new double[4];
        //myfile.seekg(0, std::ios::beg);
        std::array<std::size_t, 3> shape{3,4,5};
        NumberArray3D<RealNumber> A(shape, 0.0);
        NumberArray3D<RealNumber> B = A.GetSlice({1, 1, 1}, {2, 3, 4});

        myfile.read((char*)i, sizeof(int)*3);
        B.ReadArrayDataFromFile(&myfile);
        myfile.read((char*)d, sizeof(double)*4);
        myfile.close();
        std::cout << "i[1] expected : 13, got :" << i[1] << std::endl;
        std::cout << "d[3] expected : 32.4, got :" << d[3] << std::endl;
        A.Print();
        delete[] i;
        delete[] d;
    }
}

void test_file_read_backwards() {
    std::streampos objsize;
    std::ifstream myfile;
    myfile.open("example.bin", std::ios::in | std::ios::binary);
    if(myfile.is_open()) {
        int* i = new int[3];
        double* d = new double[4];
        myfile.seekg(-sizeof(double)*4, std::ios::end);
        myfile.read((char*)d, sizeof(double)*4);
        myfile.seekg(0, std::ios::beg);
        myfile.read((char*)i, sizeof(int)*3);
        myfile.close();
        std::cout << "i[1] expected : 13, got :" << i[1] << std::endl;
        std::cout << "d[3] expected : 32.4, got :" << d[3] << std::endl;
        delete[] i;
        delete[] d;
    }
}

#include "GaussianGridArrayManipulator.h"
#include "FDInstructionFactory.h"
#include "YeeGrid.h"
void test_yeegrid_1d() {
    std::size_t nz = 300;
    std::size_t indJ = nz/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 0.1, 10.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz*0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* E_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "E",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "E",
        0,
        {-1.0*dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );

    void* H_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params);
    yee.SetIterativeSequence({"J-update", "E-update", "E-J-update", "H-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplyIterativeInstructions(600);
}

void test_yeegrid_2d() {
    std::size_t nz = 200;
    std::size_t ny = 200;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{1, ny, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dy = (r1[1] - r0[1])/ny;
    RealNumber dt = dz/std::sqrt(2.0)*0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, indyJx, indzJx}, {1, 0, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    void* Ex_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ex_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, indyJx, indzJx},
        {1, indyJx + 1, indzJx + 1},
        "E",
        0,
        {-1.0*dt/(dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Hy_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny + 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hz_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny, nz + 1},
            "H",
            -1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("Ex-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hy_update_params);
    yee.AddUpdateInstruction("Ex-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hz_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_update_params);
    yee.AddUpdateInstruction("Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_update_params);
    yee.SetIterativeSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update", "E-J-update", "Hy-update", "Hz-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyIterativeInstructions(401);
}

void test_yeegrid_3d() {
    std::size_t nx = 100;
    std::size_t nz = 100;
    std::size_t ny = 100;
    std::size_t indxJx = nx/2;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{10.0, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{nx, ny, nz};
    RealNumber dx = (r1[0] - r0[0])/nx;
    RealNumber dy = (r1[1] - r0[1])/ny;
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz/std::sqrt(3.0)*0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {indxJx, indyJx, indzJx}, {1, 0, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            1.0 /*modulation frequency*/, M_PI/2.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* Ex_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ex_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hx_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hz_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hx_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hy_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {indxJx, indyJx, indzJx},
        {indxJx + 1, indyJx + 1, indzJx + 1},
        "E",
        0,
        {-1.0*dt/(dx*dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Hx_Ey_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hx_Ez_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hy_Ex_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hy_Ez_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hz_Ex_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -1.0*dt,
            "E");
    void* Hz_Ey_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("Ex-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hy_update_params);
    yee.AddUpdateInstruction("Ex-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hz_update_params);
    yee.AddUpdateInstruction("Ey-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hx_update_params);
    yee.AddUpdateInstruction("Ey-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hz_update_params);
    yee.AddUpdateInstruction("Ez-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hx_update_params);
    yee.AddUpdateInstruction("Ez-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hy_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("Hx-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ey_update_params);
    yee.AddUpdateInstruction("Hx-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ez_update_params);
    yee.AddUpdateInstruction("Hy-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ex_update_params);
    yee.AddUpdateInstruction("Hy-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ez_update_params);
    yee.AddUpdateInstruction("Hz-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ex_update_params);
    yee.AddUpdateInstruction("Hz-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ey_update_params);
    yee.SetIterativeSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update",
                                          "Ey-Hx-update", "Ey-Hz-update",
                                          "Ez-Hx-update", "Ez-Hy-update",
                                          "E-J-update",
                                          "Hx-Ey-update", "Hx-Ez-update",
                                          "Hy-Ex-update", "Hy-Ez-update",
                                          "Hz-Ex-update", "Hz-Ey-update"});
    yee.AddGridElementView("E-x", "E", 0, {indxJx, 0, 0}, {indxJx + 1, ny + 1, nz + 1});
    yee.AddGridElementView("H-y", "H", 1, {indxJx + int(nx/4), 0, 0}, {indxJx + int(nx/4) + 1, ny + 1, nz});
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyIterativeInstructions(150);
}

void test_yeegrid_dielectric_1d() {
    std::size_t nz = 1000;
    std::size_t indJ = nz/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 0.1, 10.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz*0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("EpsilonInv", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("EpsilonUpdater", "EpsilonInv", 0 /*x*/,
            {r0[0] - 0.1, r0[1] - 0.1, 7.0} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 9.0} /*box corner 1*/,
            {0.0, 0.0, 0.5} /*edge thickness*/,
            1.0/4.0 /*insideValue*/, 1.0 /*outsideValue*/);

    void* D_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "D",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* E_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 1},
        {1, 1, nz},
        "E",
        0,
        {1.0},          // bValues
        {"EpsilonInv"},
        {0},
        {{0, 0, 1}},
        {"D"},
        {0},
        {{0, 0, 1}}
    );
    void* D_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "D",
        0,
        {-1.0*dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* H_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    void* Epsilon_update_params = yee.ConstructParams_A_equal_func_r_t(
        "EpsilonUpdater"
    );
    yee.AddUpdateInstruction("Epsilon-update", FDInstructionCode::A_equal_func_r_t, Epsilon_update_params);
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("D-update", FDInstructionCode::A_plusequal_sum_b_C, D_update_params);
    yee.AddUpdateInstruction("D-J-update", FDInstructionCode::A_plusequal_sum_b_C, D_J_update_params);
    yee.AddUpdateInstruction("E-update", FDInstructionCode::A_equal_sum_bB_C, E_update_params);
    yee.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params);
    yee.SetSingleRunSequence({"Epsilon-update"});
    yee.SetIterativeSequence({"J-update", "D-update", "D-J-update", "E-update", "H-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("Eps", "EpsilonInv", 0);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("Eps", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplySingleRunInstructions();
    yee.ApplyIterativeInstructions(600);
}

void test_yeegrid_dielectric_pml_1d() {
    std::size_t nz = 600;
    std::size_t indJ = nz/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 0.1, 15.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz*0.95;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("sig-E", ElementType::EdgeE);
    yee.AddEntireGridElement("EpsilonInv", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddEntireGridElement("sig-H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("EpsilonUpdater", "EpsilonInv", 0 /*x*/,
            {r0[0] - 0.1, r0[1] - 0.1, 9.0} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 15.0} /*box corner 1*/,
            {0.0, 0.0, 0.0} /*edge thickness*/,
            1.0/4.0 /*insideValue*/, 1.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEUpdater", "sig-E", 2 /*z*/,
            {r0[0] - 0.1, r0[1] - 0.1, 4.01} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 11.01} /*box corner 1*/,
            {0.0, 0.0, 1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +1.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHUpdater", "sig-H", 2 /*z*/,
            {r0[0] - 0.1, r0[1] - 0.1, 4.01} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 11.01} /*box corner 1*/,
            {0.0, 0.0, 1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +1.0 /*outsideValue*/);

    void* D_curlH_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "D",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* D_sigz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 1},
        {1, 1, nz},
        "D",
        0,
        {-1.0*dt},          // bValues
        {"sig-E"},
        {2},
        {{0, 0, 1}},
        {"D"},
        {0},
        {{0, 0, 1}}
    );
    void* E_D_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 1},
        {1, 1, nz},
        "E",
        0,
        {1.0},          // bValues
        {"EpsilonInv"},
        {0},
        {{0, 0, 1}},
        {"D"},
        {0},
        {{0, 0, 1}}
    );
    void* D_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "D",
        0,
        {-1.0*dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* H_curlE_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* H_sigz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, 1, nz},
        "H",
        1,
        {-1.0*dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"H"},
        {1},
        {{0, 0, 0}}
    );
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    void* Epsilon_update_params = yee.ConstructParams_A_equal_func_r_t(
        "EpsilonUpdater"
    );
    void* SigE_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEUpdater"
    );
    void* SigH_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHUpdater"
    );
    yee.AddUpdateInstruction("Epsilon-update", FDInstructionCode::A_equal_func_r_t, Epsilon_update_params);
    yee.AddUpdateInstruction("SigE-update", FDInstructionCode::A_equal_func_r_t, SigE_update_params);
    yee.AddUpdateInstruction("SigH-update", FDInstructionCode::A_equal_func_r_t, SigH_update_params);
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("D-curlH-update", FDInstructionCode::A_plusequal_sum_b_C, D_curlH_update_params);
    yee.AddUpdateInstruction("D-sigz-update", FDInstructionCode::A_plusequal_sum_bB_C, D_sigz_update_params);
    yee.AddUpdateInstruction("D-J-update", FDInstructionCode::A_plusequal_sum_b_C, D_J_update_params);
    yee.AddUpdateInstruction("E-D-update", FDInstructionCode::A_equal_sum_bB_C, E_D_update_params);
    yee.AddUpdateInstruction("H-curlE-update", FDInstructionCode::A_plusequal_sum_b_C, H_curlE_update_params);
    yee.AddUpdateInstruction("H-sigz-update", FDInstructionCode::A_plusequal_sum_bB_C, H_sigz_update_params);
    yee.SetSingleRunSequence({"Epsilon-update", "SigE-update", "SigH-update"});
    yee.SetIterativeSequence({"J-update",
                              "D-sigz-update", "D-curlH-update", "D-J-update", "E-D-update",
                              "H-sigz-update", "H-curlE-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("Eps", "EpsilonInv", 0);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("Eps", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplySingleRunInstructions();
    yee.ApplyIterativeInstructions(600);
}

void test_yeegrid_2d_pml() {
    std::size_t nz = 200;
    std::size_t ny = 200;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{1, ny, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dy = (r1[1] - r0[1])/ny;
    RealNumber dt = dz/std::sqrt(2.0)*0.95;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("F", ElementType::EdgeE);
    yee.AddEntireGridElement("dF", ElementType::EdgeE);
    yee.AddEntireGridElement("sig-E", ElementType::EdgeE);

    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddEntireGridElement("B", ElementType::EdgeH);
    yee.AddEntireGridElement("G", ElementType::EdgeH);
    yee.AddEntireGridElement("dG", ElementType::EdgeH);
    yee.AddEntireGridElement("sig-H", ElementType::EdgeH);

    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, indyJx, indzJx}, {1, 0, 0});

    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEzUpdater", "sig-E", 2 /*z*/,
            {r0[0] - 0.1, r0[1] - 0.1, 2.01} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 8.01} /*box corner 1*/,
            {0.0, 0.0, 1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEyUpdater", "sig-E", 1 /*y*/,
            {r0[0] - 0.1, 2.01, r0[2] - 0.1} /*box corner 0*/, {r1[0] + 0.1, 8.01, r1[2] + 0.1} /*box corner 1*/,
            {0.0, 1.0, 0.0} /*edge thickness*/,
            0.0 /*insideValue*/, +2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHzUpdater", "sig-H", 2 /*z*/,
            {r0[0] - 0.1, r0[1] - 0.1, 2.01} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 8.01} /*box corner 1*/,
            {0.0, 0.0, 1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHyUpdater", "sig-H", 1 /*y*/,
            {r0[0] - 0.1, 2.01, r0[2] - 0.1} /*box corner 0*/, {r1[0] + 0.1, 8.01, r1[2] + 0.1} /*box corner 1*/,
            {0.0, 1.0, 0.0} /*edge thickness*/,
            0.0 /*insideValue*/, +2.0 /*outsideValue*/);
    void* dFx_sigzFx_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 1, 1},
        {1, ny, nz},
        "dF",
        0,
        {-1.0},          // bValues
        {"sig-E"},
        {2},
        {{0, 1, 1}},
        {"F"},
        {0},
        {{0, 1, 1}}
    );
    void* dFx_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "dF",              // name of A
                +1.0,             // b value
                "H"                  // name of C
                );
    void* dFx_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "dF",              // name of A
                +1.0,             // b value
                "H"                  // name of C
                );
    void* dFx_Jx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, indyJx, indzJx},
        {1, indyJx + 1, indzJx + 1},
        "dF",    //
        0,
        {-1.0/(dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Fx_dFx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "F",
        0,
        {+1.0*dt},          // bValues
        {"dF"},
        {0},
        {{0, 1, 1}}
    );
    void* Dx_dFx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "D",
        0,
        {+1.0*dt},          // bValues
        {"dF"},
        {0},
        {{0, 1, 1}}
    );
    void* Dx_sigyDx_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 1, 1},
        {1, ny, nz},
        "D",
        0,
        {-1.0*dt},          // bValues
        {"sig-E"},
        {1},            // y component
        {{0, 1, 1}},
        {"D"},
        {0},
        {{0, 1, 1}}
    );
    void* Ex_Dx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "E",
        0,
        {+1.0},          // bValues
        {"D"},
        {0},
        {{0, 1, 1}}
    );
    void* dGy_Ex_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny + 1, nz},
            "dG",
            -1.0,
            "E");
    void* dGz_Ex_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny, nz + 1},
            "dG",
            -1.0,
            "E");
    void* dGz_sigyGz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "dG",
        2,
        {-1.0},          // bValues
        {"sig-H"},
        {1},
        {{0, 0, 0}},
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    void* Bz_sigzGz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "B",
        2,
        {+1.0*dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    void* Gz_dGz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "G",
        2,
        {+1.0*dt},          // bValues
        {"dG"},
        {2},
        {{0, 0, 0}}
    );
    void* Gy_dGy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "G",
        1,
        {+1.0*dt},          // bValues
        {"dG"},
        {1},
        {{0, 0, 0}}
    );
    void* Bz_dGz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "B",
        2,
        {+1.0*dt},          // bValues
        {"dG"},
        {2},
        {{0, 0, 0}}
    );
    void* By_sigzBy_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz},    // TODO: ny+1 should be done separately
        "B",
        1,
        {-1.0*dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"B"},
        {1},
        {{0, 0, 0}}
    );
    void* By_dGy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "B",
        1,
        {+1.0*dt},          // bValues
        {"dG"},
        {1},
        {{0, 0, 0}}
    );
    void* Hy_By_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "H",
        1,
        {+1.0},          // bValues
        {"B"},
        {1},
        {{0, 0, 0}}
    );
    void* Hz_Bz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "H",
        2,
        {+1.0},          // bValues
        {"B"},
        {2},
        {{0, 0, 0}}
    );
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    void* SigEz_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEzUpdater"
    );
    void* SigEy_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEyUpdater"
    );
    void* SigHz_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHzUpdater"
    );
    void* SigHy_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHyUpdater"
    );
    void* Ex_Fx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "E",
        0,
        {+1.0},          // bValues
        {"F"},
        {0},
        {{0, 1, 1}}
    );
    void* Hy_Gy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "H",
        1,
        {+1.0},          // bValues
        {"G"},
        {1},
        {{0, 0, 0}}
    );
    void* Hz_Gz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "H",
        2,
        {+1.0},          // bValues
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("sigEz-update", FDInstructionCode::A_equal_func_r_t, SigEz_update_params);
    yee.AddUpdateInstruction("sigEy-update", FDInstructionCode::A_equal_func_r_t, SigEy_update_params);
    yee.AddUpdateInstruction("sigHz-update", FDInstructionCode::A_equal_func_r_t, SigHz_update_params);
    yee.AddUpdateInstruction("sigHy-update", FDInstructionCode::A_equal_func_r_t, SigHy_update_params);

    yee.AddUpdateInstruction("dFx-Jx-update", FDInstructionCode::A_plusequal_sum_b_C, dFx_Jx_update_params);
    yee.AddUpdateInstruction("dFx-Hy-update", FDInstructionCode::A_equal_sum_b_C, dFx_Hy_update_params);
    yee.AddUpdateInstruction("dFx-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, dFx_Hz_update_params);
    yee.AddUpdateInstruction("dFx-sigz-Fx-update", FDInstructionCode::A_plusequal_sum_bB_C, dFx_sigzFx_update_params);
    yee.AddUpdateInstruction("Fx-dFx-update", FDInstructionCode::A_plusequal_sum_b_C, Fx_dFx_update_params);
    yee.AddUpdateInstruction("Dx-dFx-update", FDInstructionCode::A_plusequal_sum_b_C, Dx_dFx_update_params);
    yee.AddUpdateInstruction("Dx-sigy-Dx-update", FDInstructionCode::A_plusequal_sum_bB_C, Dx_sigyDx_update_params);
    yee.AddUpdateInstruction("Ex-Dx-update", FDInstructionCode::A_equal_sum_b_C, Ex_Dx_update_params);

    yee.AddUpdateInstruction("dGy-Ex-update", FDInstructionCode::A_equal_sum_b_C, dGy_Ex_update_params);
    yee.AddUpdateInstruction("dGz-Ex-update", FDInstructionCode::A_equal_sum_b_C, dGz_Ex_update_params);
    yee.AddUpdateInstruction("dGz-sigy-Gz-update", FDInstructionCode::A_plusequal_sum_bB_C, dGz_sigyGz_update_params);
    yee.AddUpdateInstruction("Gz-dGz-update", FDInstructionCode::A_plusequal_sum_b_C, Gz_dGz_update_params);
    yee.AddUpdateInstruction("Gy-dGy-update", FDInstructionCode::A_plusequal_sum_b_C, Gy_dGy_update_params);
    yee.AddUpdateInstruction("Bz-dGz-update", FDInstructionCode::A_plusequal_sum_b_C, Bz_dGz_update_params);
    yee.AddUpdateInstruction("Bz-sigz-Gz-update", FDInstructionCode::A_plusequal_sum_bB_C, Bz_sigzGz_update_params);
    yee.AddUpdateInstruction("By-dGy-update", FDInstructionCode::A_plusequal_sum_b_C, By_dGy_update_params);
    yee.AddUpdateInstruction("By-sigz-By-update", FDInstructionCode::A_plusequal_sum_bB_C, By_sigzBy_update_params);
    yee.AddUpdateInstruction("Hy-By-update", FDInstructionCode::A_equal_sum_b_C, Hy_By_update_params);
    yee.AddUpdateInstruction("Hz-Bz-update", FDInstructionCode::A_equal_sum_b_C, Hz_Bz_update_params);

    yee.AddUpdateInstruction("Ex-Fx-update", FDInstructionCode::A_equal_sum_b_C, Ex_Fx_update_params);
    yee.AddUpdateInstruction("Hy-Gy-update", FDInstructionCode::A_equal_sum_b_C, Hy_Gy_update_params);
    yee.AddUpdateInstruction("Hz-Gz-update", FDInstructionCode::A_equal_sum_b_C, Hz_Gz_update_params);

    yee.SetSingleRunSequence({"sigEz-update", "sigEy-update", "sigHy-update", "sigHz-update"});
    yee.SetIterativeSequence({"J-update",
                              "dFx-Hy-update", "dFx-Jx-update", "dFx-Hz-update", "dFx-sigz-Fx-update",
                              "Fx-dFx-update",
                              "Dx-dFx-update", "Dx-sigy-Dx-update",
                              "Ex-Dx-update",
                              "dGy-Ex-update", "dGz-Ex-update", "dGz-sigy-Gz-update",
                              "Gz-dGz-update", "Gy-dGy-update",
                              "Bz-dGz-update", "Bz-sigz-Gz-update",
                              "By-dGy-update", "By-sigz-By-update",
                              "Hy-By-update", "Hz-Bz-update"
                              });
    /*yee.SetIterativeSequence({"J-update",
                              "dFx-Hy-update", "dFx-Jx-update", "dFx-Hz-update", "dFx-sigz-Fx-update",
                              "Fx-dFx-update",
                              "Ex-Fx-update",
                              "dGy-Ex-update", "dGz-Ex-update", "dGz-sigy-Gz-update",
                              "Gz-dGz-update", "Gy-dGy-update",
                              "Hy-Gy-update", "Hz-Gz-update"
                              });*/
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("sigEz", "sig-E", 2);
    yee.AddFullGridElementView("sigEy", "sig-E", 1);
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("sigEz", 10);
    yee.SetDataStoreRate("sigEy", 10);
    yee.ApplySingleRunInstructions();
    yee.ApplyIterativeInstructions(401);
}


void test_yeegrid_3d_pml() {
    std::size_t nx = 100;
    std::size_t nz = 100;
    std::size_t ny = 100;
    std::size_t indxJx = nx/2;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{10.0, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{nx, ny, nz};
    RealNumber dx = (r1[0] - r0[0])/nx;
    RealNumber dy = (r1[1] - r0[1])/ny;
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz/std::sqrt(3.0)*0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {indxJx, indyJx, indzJx}, {1, 0, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            1.0 /*modulation frequency*/, M_PI/2.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* Ex_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ex_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hx_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hz_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hx_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hy_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +1.0*dt,             // b value
                "H"                  // name of C
                );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {indxJx, indyJx, indzJx},
        {indxJx + 1, indyJx + 1, indzJx + 1},
        "E",
        0,
        {-1.0*dt/(dx*dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Hx_Ey_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hx_Ez_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hy_Ex_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hy_Ez_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -1.0*dt,
            "E");
    void* Hz_Ex_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -1.0*dt,
            "E");
    void* Hz_Ey_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("Ex-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hy_update_params);
    yee.AddUpdateInstruction("Ex-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hz_update_params);
    yee.AddUpdateInstruction("Ey-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hx_update_params);
    yee.AddUpdateInstruction("Ey-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hz_update_params);
    yee.AddUpdateInstruction("Ez-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hx_update_params);
    yee.AddUpdateInstruction("Ez-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hy_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("Hx-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ey_update_params);
    yee.AddUpdateInstruction("Hx-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ez_update_params);
    yee.AddUpdateInstruction("Hy-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ex_update_params);
    yee.AddUpdateInstruction("Hy-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ez_update_params);
    yee.AddUpdateInstruction("Hz-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ex_update_params);
    yee.AddUpdateInstruction("Hz-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ey_update_params);
    yee.SetIterativeSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update",
                                          "Ey-Hx-update", "Ey-Hz-update",
                                          "Ez-Hx-update", "Ez-Hy-update",
                                          "E-J-update",
                                          "Hx-Ey-update", "Hx-Ez-update",
                                          "Hy-Ex-update", "Hy-Ez-update",
                                          "Hz-Ex-update", "Hz-Ey-update"});
    yee.AddGridElementView("E-x", "E", 0, {indxJx, 0, 0}, {indxJx + 1, ny + 1, nz + 1});
    yee.AddGridElementView("H-y", "H", 1, {indxJx + int(nx/4), 0, 0}, {indxJx + int(nx/4) + 1, ny + 1, nz});
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyIterativeInstructions(150);
}

#include "YeeGridCollection.h"
void test_yeegrid_1d_collection() {
    YeeGridCollection yeeCollection;
    std::size_t yeeLeftInd = yeeCollection.AddGrid();
    std::size_t yeeRightInd = yeeCollection.AddGrid();

    std::size_t nz = 300/2;
    std::size_t indJ = nz/2;
    std::array<RealNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<RealNumber, 3> r1{0.1, 0.0, 5.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    RealNumber dz = (r1[2] - r0[2])/nz;
    RealNumber dt = dz*0.99;
    //------ left side
    YeeGrid3D& yee_l = yeeCollection.GetGrid(yeeLeftInd);
    YeeGrid3D& yee_r = yeeCollection.GetGrid(yeeRightInd);

    yee_l.SetCornerCoordinates(r0, r1);
    yee_l.SetNumOfCells(nCells);
    yee_l.SetTimeResolution(dt);
    yee_l.AddEntireGridElement("E", ElementType::EdgeE);
    yee_l.AddEntireGridElement("H", ElementType::EdgeH);
    yee_l.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee_l.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* E_update_params_l = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee_l,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "E",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* E_J_update_params_l = yee_l.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "E",
        0,
        {-1.0*dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );

    void* H_update_params_l = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee_l,
            {0, 0, 0},
            {1, 1, nz - 1},     // except the last H element
            "H",
            -1.0*dt,
            "E");
    void* H_update_params_nz_inside_l = yee_l.ConstructParams_A_plusequal_sum_b_C(
            {0, 0, nz - 1},  // indStart
            {1, 1, nz},  // indEnd
            "H",        // A name
            1,
            {+1.0*dt/dz},     // b values
            {"E"},           // C names
            {0},
            {{0, 0, nz - 1}}
            );
    void* H_update_params_nz_outside_l = yee_l.ConstructParams_A_plusequal_sum_b_C_neighbor(
            &yee_r,     // neighbor grid
            {0, 0, nz - 1},  // indStart
            {1, 1, nz},  // indEnd
            "H",        // A name
            1,
            {-1.0*dt/dz},     // b values
            {"E"},           // C names
            {0},
            {{0, 0, 0}}
            );
    void* J_update_params_l = yee_l.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee_l.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params_l);
    yee_l.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params_l);
    yee_l.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params_l);
    yee_l.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params_l);
    yee_l.AddUpdateInstruction("H-update-last-inside", FDInstructionCode::A_plusequal_sum_b_C,
                                H_update_params_nz_inside_l);
    yee_l.AddUpdateInstruction("H-update-last-outside", FDInstructionCode::A_plusequal_sum_b_C_neighbor,
                                H_update_params_nz_outside_l);

    //-------- right side
    r0[2] += 5.0;
    r1[2] += 5.0;
    yee_r.SetCornerCoordinates(r0, r1);
    yee_r.SetNumOfCells(nCells);
    yee_r.SetTimeResolution(dt);
    yee_r.AddEntireGridElement("E", ElementType::EdgeE);
    yee_r.AddEntireGridElement("H", ElementType::EdgeH);

    void* E_update_params_r = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee_l,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "E",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* E_update_params_0_inside_r = yee_r.ConstructParams_A_plusequal_sum_b_C(
            {0, 0, 0},  // indStart
            {1, 1, 1},  // indEnd
            "E",        // A name
            0,
            {-1.0*dt/dz},     // b values
            {"H"},           // C names
            {1},
            {{0, 0, 0}}
            );
    void* E_update_params_0_outside_r = yee_r.ConstructParams_A_plusequal_sum_b_C_neighbor(
            &yee_l,     // neighbor grid
            {0, 0, 0},  // indStart
            {1, 1, 1},  // indEnd
            "E",        // A name
            0,
            {+1.0*dt/dz},     // b values
            {"H"},           // C names
            {1},
            {{0, 0, nz - 1}}
            );

    void* H_update_params_r = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee_l,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -1.0*dt,
            "E");

    yee_r.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params_r);
    yee_r.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params_r);
    yee_r.AddUpdateInstruction("E-update-first-inside", FDInstructionCode::A_plusequal_sum_b_C,
                                E_update_params_0_inside_r);
    yee_r.AddUpdateInstruction("E-update-first-outside", FDInstructionCode::A_plusequal_sum_b_C_neighbor,
                                E_update_params_0_outside_r);



    yee_l.AddInstructionSequence("iterative-E-update", {"J-update", "E-update", "E-J-update"});
    yee_l.AddInstructionSequence("iterative-H-update", {"H-update", "H-update-last-inside", "H-update-last-outside"});

    yee_r.AddInstructionSequence("iterative-E-update", {"E-update", "E-update-first-inside", "E-update-first-outside"});
    yee_r.AddInstructionSequence("iterative-H-update", {"H-update"});


    yee_l.AddFullGridElementView("E-x-l", "E", 0);
    yee_l.AddFullGridElementView("H-y-l", "H", 1);
    yee_l.SetDataStoreRate("E-x-l", 1);
    yee_l.SetDataStoreRate("H-y-l", 1);
    yee_l.DeleteOlderViewFiles();

    yee_r.AddFullGridElementView("E-x-r", "E", 0);
    yee_r.AddFullGridElementView("H-y-r", "H", 1);
    yee_r.SetDataStoreRate("E-x-r", 1);
    yee_r.SetDataStoreRate("H-y-r", 1);
    yee_r.DeleteOlderViewFiles();

    yeeCollection.RunInstructionsPeriodically(0, 600, {"iterative-E-update", "iterative-H-update"});
}


#include "ParameterExtractor.h"
void test_read_json() {
    ParameterExtractor paramExtractor("instructions/MaxwellYee1D_processed.json");

    auto simulationType = paramExtractor.GetStringProperty("simulationType");

    ParameterExtractor dimensionsExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.dimensions"));
    auto r0 = dimensionsExtractor.Get3VecRealProperty("r0");
    auto r1 = dimensionsExtractor.Get3VecRealProperty("r1");
    auto nCells = dimensionsExtractor.Get3VecUintProperty("nCells");

    std::cout << "simulationType : " << simulationType << std::endl;
    std::cout << "r0 : " << r0[0] << " " << r0[1] << " " << r0[2] << std::endl;
    std::cout << "r1 : " << r1[0] << " " << r1[1] << " " << r1[2] << std::endl;
    std::cout << "nCells : " << nCells[0] << " " << nCells[1] << " " << nCells[2] << std::endl;

    ParameterExtractor entireGridArrayExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.entireGridArrays"));
    ParameterExtractor entireGridArrayExtractorFirst(entireGridArrayExtractor.GetSubTreeByIndex(1));
    std::string gridName = entireGridArrayExtractorFirst.GetStringProperty("name");
    std::cout << "GridArrayName : " << gridName << std::endl;

    ParameterExtractor updateSequenceExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.updateSequences"));
    for(std::size_t i = 0; i < updateSequenceExtractor.GetSize(); ++i) {
        ParameterExtractor updateSequenceExtractor_i(updateSequenceExtractor.GetSubTreeByIndex(i));
        auto update_sequence_name = updateSequenceExtractor_i.GetStringProperty("name");
        std::cout << "update sequence : " << update_sequence_name <<  std::endl;
        ParameterExtractor updateSequence_sequenceExtractor(updateSequenceExtractor_i.GetSubTreeRootNode("sequence"));
        for(std::size_t j = 0; j < updateSequence_sequenceExtractor.GetSize(); ++j) {
            ParameterExtractor updateSequence_sequenceExtractor_j(updateSequence_sequenceExtractor.GetSubTreeByIndex(j));
            auto update_name_j = updateSequence_sequenceExtractor_j.GetStringProperty("");
            std::cout << update_name_j <<  std::endl;
        }
    }
}

#include "ParameterExtractor.h"
#include "boost/lexical_cast.hpp"
#include "ParamFileTranslator.h"
void test_run_fdtd_1d_from_json() {
    RealNumber z0 = 0.0;
    RealNumber z1 = 10.0;
    std::size_t nz = 300;
    RealNumber dz = (z1 - z0)/nz;
    RealNumber stabilityFactor = 0.99;
    RealNumber dt = dz*stabilityFactor;
    RealNumber z_j = 5.0;
    std::size_t indJ = std::round(z_j/dz);
    std::size_t numOfTimeSamples = 600;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_z0_\"", boost::lexical_cast<std::string>(z0)},
            {"\"_z1_\"", boost::lexical_cast<std::string>(z1)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_dz_\"", boost::lexical_cast<std::string>(dz)},
            {"\"_dt_\"", boost::lexical_cast<std::string>(dt)},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(dt/dz)},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(-dt/dz)},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(z_j)},
            {"\"_indJ_\"", boost::lexical_cast<std::string>(indJ)},
            {"\"_indJ_p1_\"", boost::lexical_cast<std::string>(indJ + 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D.json",
                "instructions/MaxwellYee1D_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee1D_processed.json");
    fileTranslator.Translate();
}

void test_run_fdtd_2d_from_json() {
    RealNumber y0 = 0.0;
    RealNumber y1 = 10.0;
    RealNumber z0 = 0.0;
    RealNumber z1 = 10.0;
    std::size_t ny = 200;
    std::size_t nz = 200;
    RealNumber dy = (y1 - y0)/ny;
    RealNumber dz = (z1 - z0)/nz;
    RealNumber stabilityFactor = 0.99;
    RealNumber dt = 1.0/std::sqrt(1.0/(dy*dy) + 1.0/(dz*dz))*stabilityFactor;
    RealNumber y_j = (y0 + y1)/2;
    RealNumber z_j = (z0 + z1)/2;
    std::size_t indyJ = std::round(y_j/dy);
    std::size_t indzJ = std::round(z_j/dz);
    std::size_t numOfTimeSamples = 401;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_y0_\"", boost::lexical_cast<std::string>(y0)},
            {"\"_y1_\"", boost::lexical_cast<std::string>(y1)},
            {"\"_z0_\"", boost::lexical_cast<std::string>(z0)},
            {"\"_z1_\"", boost::lexical_cast<std::string>(z1)},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(dy)},
            {"\"_dz_\"", boost::lexical_cast<std::string>(dz)},
            {"\"_dt_\"", boost::lexical_cast<std::string>(dt)},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(dt/dy)},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(-dt/dy)},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(dt/dz)},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(-dt/dz)},
            {"\"_m_dt_dydz_\"", boost::lexical_cast<std::string>(-dt/(dy*dz))},
            {"\"_y_j_\"", boost::lexical_cast<std::string>(y_j)},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(z_j)},
            {"\"_indyJ_\"", boost::lexical_cast<std::string>(indyJ)},
            {"\"_indyJ_p1_\"", boost::lexical_cast<std::string>(indyJ + 1)},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D.json",
                "instructions/MaxwellYee2D_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee2D_processed.json");
    fileTranslator.Translate();
}

int main(int argc, char** argv) {
    //test_yeegrid_1d();
    // test_read_json();
    test_run_fdtd_2d_from_json();
}


