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


int main(int argc, char** argv) {
    test_multidim_array();

}


