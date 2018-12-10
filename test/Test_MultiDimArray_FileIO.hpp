

#include <fstream>
#include "MultiDimArrayFileIO.hpp"

void test_file_write() {
    std::ofstream myfile;
    myfile.open("example.bin", std::ios::out | std::ios::binary);
    if(myfile.is_open()) {
        int* i = new int[3];
        i[0] = 1; i[1] = 13;
        std::array<std::size_t, 3> shape{3,4,5};
        NumberArray3D<FPNumber> A(shape, 2.01);
        NumberArray3D<FPNumber> B = A.GetSlice({1, 1, 1}, {2, 3, 4});
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
        NumberArray3D<FPNumber> A(shape, 0.0);
        NumberArray3D<FPNumber> B = A.GetSlice({1, 1, 1}, {2, 3, 4});

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



