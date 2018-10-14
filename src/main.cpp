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
    std::array<RealNumber, 3> r1{0.1, 0.0, 10.0};
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
    yee.SetIterationSequence({"J-update", "E-update", "E-J-update", "H-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplyUpdateInstructions(600);
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
    yee.SetIterationSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update", "E-J-update", "Hy-update", "Hz-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyUpdateInstructions(401);
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
    yee.SetIterationSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update",
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
    yee.ApplyUpdateInstructions(150);
}

void test_yeegrid_dielectric_1d() {
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
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("EpsilonInv", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("EpsilonUpdater", "EpsilonInv", 0 /*x*/,
            {r0[0] - 0.1, r0[1] - 0.1, 8.0} /*box corner 0*/, {r1[0] + 0.1, r1[1] + 0.1, 9.0} /*box corner 1*/,
            1.0/4.0 /*insideValue*/, 1.0 /*outsideValue*/);

    void* D_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "D",        // A name
            1.0*dt,     // b value
            "H"         // C name
            );
    void* E_update_params = yee.ConstructParams_A_plusequal_sum_B_C(
        {0, 0, 1},
        {1, 1, nz},
        "E",
        0,
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
    yee.AddUpdateInstruction("E-update", FDInstructionCode::A_equal_sum_B_C, E_update_params);
    yee.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params);
    yee.SetIterationSequence({"Epsilon-update", "J-update", "D-update", "D-J-update", "E-update", "H-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("Eps", "EpsilonInv", 0);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("Eps", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplyUpdateInstructions(600);
}

int main(int argc, char** argv) {
    test_yeegrid_dielectric_1d();
}


