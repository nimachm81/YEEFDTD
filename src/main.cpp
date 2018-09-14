 // 

#include <iostream>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"

int main(int argc, char** argv) {
    std::cout << "Test." << std::endl;
    
    std::array<std::size_t, 3> shape{2,3,4};
    NumberArray3D<RealNumber> A(shape, 2.0);
    A.Print();
    
    NumberArray3D<RealNumber> B(shape, 3.0);
    B.Print();
    
    A += B;
    A.Print();
    
    A /= B;
    A.Print();

    A *= 2.0;
    A.Print();

    A.GetSlice({0, 0, 1}, {2, 2, 2}).Print();

}


