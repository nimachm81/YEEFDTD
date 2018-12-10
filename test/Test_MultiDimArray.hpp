

#include <cmath>
#include <iostream>
#include <any>
#include <typeinfo>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"


void test_multidim_array() {
    std::cout << "Test." << std::endl;

    std::array<std::size_t, 3> shape{2,3,4};
    NumberArray3D<FPNumber> A(shape, 2.0);
    std::cout << " A : " << &A << std::endl;
    A.Print();

    NumberArray3D<FPNumber> B(shape, 3.0);
    std::cout << " B : " << &B << std::endl;
    B.Print();

    A += B;
    A.Print();

    A /= B;
    A.Print();

    A *= 2.0;
    A.Print();

    A.GetSlice({0, 0, 1}, {2, 2, 2}).Print();

    NumberArray3D<FPNumber> C = (A + B);
    std::cout << " C : " << &C << std::endl;
    C.Print();

    A = B;
    std::cout << " A : " << &A << std::endl;
    A.Print();

    std::any a = std::make_any<NumberArray3D<FPNumber>>(shape, 4.01);
    std::cout << " a : " << a.type().name() << std::endl;

    std::cout << " a : " << typeid(NumberArray3D<FPNumber>).name() << std::endl;
    std::any_cast<NumberArray3D<FPNumber>>(a).Print();
}


