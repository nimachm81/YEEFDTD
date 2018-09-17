
#ifndef FDTD_MULTIDIMARRAY_H_
#define FDTD_MULTIDIMARRAY_H_


#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <array>       //std::array
#include <iostream>

#include "MultiDimArrayAllocator.hpp"
#include "MultiDimArrayPrinting.hpp"

template <typename T>
class NumberArray3D {
    private:
    std::array<std::size_t, 3> indStart;
    std::array<std::size_t, 3> shape;
    T*** arrayData = nullptr;
    bool isSlice = true;

    public:
    NumberArray3D(std::array<std::size_t, 3> shape, T initValue) : isSlice(false), shape(shape), indStart{0, 0, 0} {
        arrayData = Create3DNumberArray(shape, initValue);
    }

    // The copy constructor shares the Array pointer with the array being copied. To create a fresh copy use
    // the explicite copy methods.
    NumberArray3D(const NumberArray3D& numArray) : isSlice(true), shape(numArray.GetShape()),
            indStart(numArray.GetIndStart()) {
        arrayData = numArray.GetArrayData();
        std::cout << "inside the copy constructor" << std::endl;
    }

    NumberArray3D(const NumberArray3D&& numArray) :
            isSlice(numArray.IsSlice()),
            shape(std::move(numArray.GetShape())),
            indStart(std::move(numArray.GetIndStart())) {
        std::cout << "inside the move constructor" << std::endl;
        arrayData = std::move(numArray.GetArrayData());
    }

    NumberArray3D(T*** arraySlice, std::array<std::size_t, 3> shape, std::array<std::size_t, 3> indStart) :
            shape(shape), indStart(indStart), isSlice(true) {
        arrayData = arraySlice;
    }

    ~NumberArray3D() {
        if(!isSlice) {
            // when indStart is different than {0, 0, 0} the shape of the underlying allocated array should be
            // adjusted:
            std::array<std::size_t, 3> shapeOfArrayData{shape[0] + indStart[0],
                                                        shape[1] + indStart[1],
                                                        shape[2] + indStart[2] };
            Delete3DNumberArray(arrayData, shapeOfArrayData);
            std::cout << " Deleting : " << this << std::endl;
        }
    }

    T*** GetArrayData() const {
        return arrayData;
    }

    const std::array<std::size_t, 3>& GetIndStart() const {
        return indStart;
    }

    const std::array<std::size_t, 3>& GetShape() const {
        return shape;
    }

    bool IsSlice() const {
        return isSlice;
    }

    NumberArray3D GetSlice(std::array<std::size_t, 3> indStart_slice, std::array<std::size_t, 3> indEnd_slice) {
        std::array<std::size_t, 3> shape_slice;
        for(std::size_t i = 0; i < 3; ++i) {
            shape_slice[i] = indEnd_slice[i] - indStart_slice[i];
        }
        // indStart_slice_total: with respect to the origin of the allocated array
        std::array<std::size_t, 3> indStart_slice_total{indStart_slice[0] + indStart[0],
                                                        indStart_slice[1] + indStart[1],
                                                        indStart_slice[2] + indStart[2]};
        return NumberArray3D(arrayData, shape_slice, indStart_slice_total);
    }

    //------------------------------ Operators ---------------------------------------

    NumberArray3D& operator=(NumberArray3D& numArray) {
        if(this != &numArray) {
            this->~NumberArray3D();
        }
        isSlice = true;
        indStart = numArray.GetIndStart();
        shape = numArray.GetShape();
        arrayData = numArray.GetArrayData();
        return *this;
    }

    friend NumberArray3D operator+(const NumberArray3D& numArrA, const NumberArray3D& numArrB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        const std::array<std::size_t, 3>& b_shape = numArrB.GetShape();
        T*** b_arrayData = numArrB.GetArrayData();
        const std::array<std::size_t, 3>& b_indStart = numArrB.GetIndStart();
        std::size_t b_ind0 = b_indStart[0];
        std::size_t b_ind1 = b_indStart[1];
        std::size_t b_ind2 = b_indStart[2];

        assert( a_shape == b_shape );
        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] +
                            b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator+(const NumberArray3D& numArrA, const T numB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] + numB;
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    NumberArray3D& operator+=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] +=
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator+=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] += rhs;
                }
            }
        }
        return *this;
    }

    friend NumberArray3D operator*(const NumberArray3D& numArrA, const NumberArray3D& numArrB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        const std::array<std::size_t, 3>& b_shape = numArrB.GetShape();
        T*** b_arrayData = numArrB.GetArrayData();
        const std::array<std::size_t, 3>& b_indStart = numArrB.GetIndStart();
        std::size_t b_ind0 = b_indStart[0];
        std::size_t b_ind1 = b_indStart[1];
        std::size_t b_ind2 = b_indStart[2];

        assert( a_shape == b_shape );
        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] *
                            b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator*(const NumberArray3D& numArrA, const T numB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] * numB;
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    NumberArray3D& operator*=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *=
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator*=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= rhs;
                }
            }
        }
        return *this;
    }

    friend NumberArray3D operator/(const NumberArray3D& numArrA, const NumberArray3D& numArrB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        const std::array<std::size_t, 3>& b_shape = numArrB.GetShape();
        T*** b_arrayData = numArrB.GetArrayData();
        const std::array<std::size_t, 3>& b_indStart = numArrB.GetIndStart();
        std::size_t b_ind0 = b_indStart[0];
        std::size_t b_ind1 = b_indStart[1];
        std::size_t b_ind2 = b_indStart[2];

        assert( a_shape == b_shape );
        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] /
                            b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator/(const NumberArray3D& numArrA, const T numB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] / numB;
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator/(const T numB, const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    c_arrayData[i0][i1][i2] =
                            numB / a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2];
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    NumberArray3D& operator/=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /=
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator/=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];

        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /= rhs;
                }
            }
        }
        return *this;
    }


    //------------------------------------------------------------------------------------

    void Print() {
        Print3DNumberArray(shape, arrayData);
    }

};


#endif  // FDTD_MULTIDIMARRAY_H_
