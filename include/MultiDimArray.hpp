
#ifndef FDTD_MULTIDIMARRAY_H_
#define FDTD_MULTIDIMARRAY_H_


#include <cmath>
#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <complex>
#include <array>       //std::array
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <thread>

#include "MultiDimArrayAllocator.hpp"
#include "MultiDimArrayPrinting.hpp"
#include "MultiDimArrayFileIO.hpp"
#include "MultiDimArrayBufferIO.hpp"

#include "UtilityFunctions.hpp"

template <typename T>
class NumberArray3D {
    private:
    std::array<std::size_t, 3> indStart;
    std::array<std::size_t, 3> shape;
    std::array<std::size_t, 3> stride;
    T*** arrayData = nullptr;
    bool isSlice = true;
    bool isSingleStride = false;

    public:
    NumberArray3D() : isSlice(false), shape{0, 0, 0}, indStart{0, 0, 0}, stride{1, 1, 1}, isSingleStride(true) {}
    NumberArray3D(std::array<std::size_t, 3> shape, T initValue) : isSlice(false), shape(shape),
                                                                   indStart{0, 0, 0}, stride{1, 1, 1},
                                                                   isSingleStride(true) {
        arrayData = Create3DNumberArray(shape, initValue);
    }

    // The copy constructor shares the Array pointer with the array being copied. To create a fresh copy use
    // the explicite copy methods.
    NumberArray3D(const NumberArray3D& numArray) : isSlice(true), shape(numArray.GetShape()),
            indStart(numArray.GetIndStart()), stride(numArray.GetStride()), isSingleStride(numArray.IsSingleStride()) {
        arrayData = numArray.GetArrayData();
        //std::cout << "inside the copy constructor" << std::endl;
    }

    NumberArray3D(const NumberArray3D&& numArray) :
            isSlice(numArray.IsSlice()),
            shape(std::move(numArray.GetShape())),
            indStart(std::move(numArray.GetIndStart())),
            stride(std::move(numArray.GetStride())),
            isSingleStride(numArray.IsSingleStride()) {
        arrayData = numArray.GetArrayData();

        if( !numArray.IsSlice() ) {
            assert(false);  // it is not clear what happens if numArray is not a slice and it is destroyed
        }
    }

    NumberArray3D(T*** arrayDataOfAnother3DArray, std::array<std::size_t, 3> shape, std::array<std::size_t, 3> indStart,
                  std::array<std::size_t, 3> stride = {1, 1, 1}) :
            shape(shape), indStart(indStart), isSlice(true), stride(stride) {
        arrayData = arrayDataOfAnother3DArray;
        if(stride[0] == 1 && stride[1] == 1 && stride[2] == 1) {
            isSingleStride = true;
        } else {
            isSingleStride = false;
        }
    }

    ~NumberArray3D() {
        if(!isSlice) {
            // when indStart is different than {0, 0, 0} the shape of the underlying allocated array should be
            // adjusted:
            std::array<std::size_t, 3> actualShapeOfArrayData{shape[0]*stride[0] + indStart[0],
                                                              shape[1]*stride[1] + indStart[1],
                                                              shape[2]*stride[2] + indStart[2] };
            Delete3DNumberArray(arrayData, actualShapeOfArrayData);
        }
    }

    void ReInitialize(std::array<std::size_t, 3> shape, T initValue) {
        this->~NumberArray3D();
        this->isSlice = false;
        this->shape = shape;
        this->indStart = std::array<std::size_t, 3>{0, 0, 0};
        this->stride = std::array<std::size_t, 3>{1, 1, 1};
        this->arrayData = Create3DNumberArray(shape, initValue);
        this->isSingleStride = true;
    }

    T*** GetArrayData() const {
        return arrayData;
    }

    std::array<std::size_t, 3>& GetIndStart() {
        return indStart;
    }

    const std::array<std::size_t, 3>& GetIndStart() const {
        return indStart;
    }

    std::array<std::size_t, 3>& GetShape() {
        return shape;
    }

    const std::array<std::size_t, 3>& GetShape() const {
        return shape;
    }

    std::array<std::size_t, 3>& GetStride() {
        return stride;
    }

    const std::array<std::size_t, 3>& GetStride() const {
        return stride;
    }

    bool IsSlice() const {
        return isSlice;
    }

    bool IsSingleStride() const {
        return isSingleStride;
    }

    NumberArray3D GetSlice(std::array<std::size_t, 3> indStart_slice, std::array<std::size_t, 3> indEnd_slice,
                           std::array<std::size_t, 3> stride_slice = {1, 1, 1}) {
        // The slice array includes indStart_slice but it does not include indEnd_slice
        assert(indEnd_slice[0] <= shape[0] && indEnd_slice[1] <= shape[1] && indEnd_slice[2] <= shape[2]);
        std::array<std::size_t, 3> shape_slice;
        for(std::size_t i = 0; i < 3; ++i) {
            assert(stride_slice[i] > 0);
            shape_slice[i] = (indEnd_slice[i] - 1 - indStart_slice[i]) / stride_slice[i] + 1;
            assert(shape_slice[i] >= 0);
        }
        // indStart_slice_total: with respect to the origin of the allocated array
        std::array<std::size_t, 3> indStart_slice_total{indStart_slice[0]*stride[0] + indStart[0],
                                                        indStart_slice[1]*stride[1] + indStart[1],
                                                        indStart_slice[2]*stride[2] + indStart[2]};
        std::array<std::size_t, 3> stride_slice_total{stride[0]*stride_slice[0],
                                                      stride[1]*stride_slice[1],
                                                      stride[2]*stride_slice[2]};

        return NumberArray3D(arrayData, shape_slice, indStart_slice_total, stride_slice_total);
    }

    NumberArray3D& MakeThisASliceOf(const NumberArray3D& rhs) {
        if(this != &rhs) {
            this->~NumberArray3D();
            isSlice = true;
            shape = rhs.GetShape();
            indStart = rhs.GetIndStart();
            stride = rhs.GetStride();
            isSingleStride = rhs.IsSingleStride();

            arrayData = rhs.GetArrayData();
        }
        return *this;
    }

    NumberArray3D& Add_aX_threaded(T a, NumberArray3D& x, const std::size_t N_thread = 4) {
        assert( x.GetShape() == shape );

        std::size_t n_th_x, n_th_y;
        if(N_thread%2 == 0) {
            n_th_x = 2;
            n_th_y = N_thread / n_th_x;
        } else if(N_thread%3 == 0) {
            n_th_x = 3;
            n_th_y = N_thread / n_th_x;
        }

        if( shape[0] < n_th_x || shape[1] < n_th_y ) {
            return Add_aX(a, x);
        } else {
            std::size_t n0 = shape[0];
            std::size_t n1 = shape[1];
            std::size_t n2 = shape[2];

            std::size_t dn0 = n0/n_th_x;
            std::size_t dn1 = n1/n_th_y;

            std::size_t ind0_i, ind1_j;
            std::size_t n0_i, n1_j;

            //std::size_t ind2_st = 0;

            std::vector<std::thread> workers;
            std::vector<NumberArray3D> thisSlices;
            std::vector<NumberArray3D> xSlices;

            for(std::size_t i = 0; i < n_th_x; ++i) {
                ind0_i = i * dn0;
                n0_i = dn0;
                if(i == n_th_x - 1) {
                    n0_i = n0 - i * dn0;
                }
                for(std::size_t j = 0; j < n_th_y; ++j) {
                    ind1_j = j * dn1;
                    n1_j = dn1;
                    if(j == n_th_y - 1) {
                        n1_j = n1 - j * dn1;
                    }

                    //std::cout << ind0_i << " " << ind1_j << std::endl;

                    thisSlices.push_back(
                                GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                         std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                        )
                                        );
                    xSlices.push_back(
                                x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                }
            }

            for(std::size_t i = 0; i < thisSlices.size(); ++i) {
                workers.push_back(std::thread(&NumberArray3D<T>::Add_aX, &(thisSlices[i]),
                                  a, std::ref(xSlices[i])));
            }

            for(auto& worker : workers) {
                worker.join();
            }

            return *this;
        }
    }

    NumberArray3D& Equate_aX_threaded(T a, NumberArray3D& x, const std::size_t N_thread = 4) {
        assert( x.GetShape() == shape );

        std::size_t n_th_x, n_th_y;
        if(N_thread%2 == 0) {
            n_th_x = 2;
            n_th_y = N_thread / n_th_x;
        } else if(N_thread%3 == 0) {
            n_th_x = 3;
            n_th_y = N_thread / n_th_x;
        }

        if( shape[0] < n_th_x || shape[1] < n_th_y ) {
            return Equate_aX(a, x);
        } else {
            std::size_t n0 = shape[0];
            std::size_t n1 = shape[1];
            std::size_t n2 = shape[2];

            std::size_t dn0 = n0/n_th_x;
            std::size_t dn1 = n1/n_th_y;

            std::size_t ind0_i, ind1_j;
            std::size_t n0_i, n1_j;

            //std::size_t ind2_st = 0;

            std::vector<std::thread> workers;
            std::vector<NumberArray3D> thisSlices;
            std::vector<NumberArray3D> xSlices;

            for(std::size_t i = 0; i < n_th_x; ++i) {
                ind0_i = i * dn0;
                n0_i = dn0;
                if(i == n_th_x - 1) {
                    n0_i = n0 - i * dn0;
                }
                for(std::size_t j = 0; j < n_th_y; ++j) {
                    ind1_j = j * dn1;
                    n1_j = dn1;
                    if(j == n_th_y - 1) {
                        n1_j = n1 - j * dn1;
                    }

                    //std::cout << ind0_i << " " << ind1_j << std::endl;

                    thisSlices.push_back(
                                GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                         std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                        )
                                        );
                    xSlices.push_back(
                                x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                }
            }

            for(std::size_t i = 0; i < thisSlices.size(); ++i) {
                workers.push_back(std::thread(&NumberArray3D<T>::Equate_aX, &(thisSlices[i]),
                                  a, std::ref(xSlices[i])));
            }

            for(auto& worker : workers) {
                worker.join();
            }

            return *this;
        }
    }

    NumberArray3D& Add_aXY_threaded(T a, NumberArray3D& x, NumberArray3D& y, const std::size_t N_thread = 4) {
        assert( x.GetShape() == shape );

        std::size_t n_th_x, n_th_y;
        if(N_thread%2 == 0) {
            n_th_x = 2;
            n_th_y = N_thread / n_th_x;
        } else if(N_thread%3 == 0) {
            n_th_x = 3;
            n_th_y = N_thread / n_th_x;
        }

        if( shape[0] < n_th_x || shape[1] < n_th_y ) {
            return Add_aXY(a, x, y);
        } else {
            std::size_t n0 = shape[0];
            std::size_t n1 = shape[1];
            std::size_t n2 = shape[2];

            std::size_t dn0 = n0/n_th_x;
            std::size_t dn1 = n1/n_th_y;

            std::size_t ind0_i, ind1_j;
            std::size_t n0_i, n1_j;

            //std::size_t ind2_st = 0;

            std::vector<std::thread> workers;
            std::vector<NumberArray3D> thisSlices;
            std::vector<NumberArray3D> xSlices;
            std::vector<NumberArray3D> ySlices;

            for(std::size_t i = 0; i < n_th_x; ++i) {
                ind0_i = i * dn0;
                n0_i = dn0;
                if(i == n_th_x - 1) {
                    n0_i = n0 - i * dn0;
                }
                for(std::size_t j = 0; j < n_th_y; ++j) {
                    ind1_j = j * dn1;
                    n1_j = dn1;
                    if(j == n_th_y - 1) {
                        n1_j = n1 - j * dn1;
                    }

                    //std::cout << ind0_i << " " << ind1_j << std::endl;

                    thisSlices.push_back(
                                GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                         std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                        )
                                        );
                    xSlices.push_back(
                                x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                    ySlices.push_back(
                                y.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                }
            }

            for(std::size_t i = 0; i < thisSlices.size(); ++i) {
                workers.push_back(std::thread(&NumberArray3D<T>::Add_aXY, &(thisSlices[i]),
                                  a, std::ref(xSlices[i]), std::ref(ySlices[i])));
            }

            for(auto& worker : workers) {
                worker.join();
            }

            return *this;
        }
    }

    NumberArray3D& Equate_aXY_threaded(T a, NumberArray3D& x, NumberArray3D& y, const std::size_t N_thread = 4) {
        assert( x.GetShape() == shape );

        std::size_t n_th_x, n_th_y;
        if(N_thread%2 == 0) {
            n_th_x = 2;
            n_th_y = N_thread / n_th_x;
        } else if(N_thread%3 == 0) {
            n_th_x = 3;
            n_th_y = N_thread / n_th_x;
        }

        if( shape[0] < n_th_x || shape[1] < n_th_y ) {
            return Equate_aXY(a, x, y);
        } else {
            std::size_t n0 = shape[0];
            std::size_t n1 = shape[1];
            std::size_t n2 = shape[2];

            std::size_t dn0 = n0/n_th_x;
            std::size_t dn1 = n1/n_th_y;

            std::size_t ind0_i, ind1_j;
            std::size_t n0_i, n1_j;

            //std::size_t ind2_st = 0;

            std::vector<std::thread> workers;
            std::vector<NumberArray3D> thisSlices;
            std::vector<NumberArray3D> xSlices;
            std::vector<NumberArray3D> ySlices;

            for(std::size_t i = 0; i < n_th_x; ++i) {
                ind0_i = i * dn0;
                n0_i = dn0;
                if(i == n_th_x - 1) {
                    n0_i = n0 - i * dn0;
                }
                for(std::size_t j = 0; j < n_th_y; ++j) {
                    ind1_j = j * dn1;
                    n1_j = dn1;
                    if(j == n_th_y - 1) {
                        n1_j = n1 - j * dn1;
                    }

                    //std::cout << ind0_i << " " << ind1_j << std::endl;

                    thisSlices.push_back(
                                GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                         std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                        )
                                        );
                    xSlices.push_back(
                                x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                    ySlices.push_back(
                                y.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                           std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                           )
                                     );
                }
            }

            for(std::size_t i = 0; i < thisSlices.size(); ++i) {
                workers.push_back(std::thread(&NumberArray3D<T>::Equate_aXY, &(thisSlices[i]),
                                  a, std::ref(xSlices[i]), std::ref(ySlices[i])));
            }

            for(auto& worker : workers) {
                worker.join();
            }

            return *this;
        }
    }


    NumberArray3D& Add_aX(T a, const NumberArray3D& x) {
        assert( x.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** x_arrayData = x.GetArrayData();
        const std::array<std::size_t, 3>& x_indStart = x.GetIndStart();
        std::size_t x_ind0 = x_indStart[0];
        std::size_t x_ind1 = x_indStart[1];
        std::size_t x_ind2 = x_indStart[2];
        const std::array<std::size_t, 3>& x_stride = x.GetStride();
        std::size_t x_s0 = x_stride[0];
        std::size_t x_s1 = x_stride[1];
        std::size_t x_s2 = x_stride[2];

        T* arrayData_i0i1;
        T*  x_arrayData_i0i1;

        bool isUnitStride = IsSingleStride() && x.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0][ind1 + i1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0][x_ind1 + i1] + x_ind2;
#ifdef __USE_AVX__
                    y_pe_ax(n2, &a, x_arrayData_i0i1, arrayData_i0i1);
#else
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2] += a * x_arrayData_i0i1[i2];
                    }
#endif
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0*s0][ind1 + i1*s1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0*x_s0][x_ind1 + i1*x_s1] + x_ind2;
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2*s2] += a * x_arrayData_i0i1[i2*x_s2];
                    }
                }
            }
        }
        return *this;
    }


    NumberArray3D& Equate_aX(T a, const NumberArray3D& x) {
        assert( x.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** x_arrayData = x.GetArrayData();
        const std::array<std::size_t, 3>& x_indStart = x.GetIndStart();
        std::size_t x_ind0 = x_indStart[0];
        std::size_t x_ind1 = x_indStart[1];
        std::size_t x_ind2 = x_indStart[2];
        const std::array<std::size_t, 3>& x_stride = x.GetStride();
        std::size_t x_s0 = x_stride[0];
        std::size_t x_s1 = x_stride[1];
        std::size_t x_s2 = x_stride[2];

        T* arrayData_i0i1;
        T*  x_arrayData_i0i1;

        bool isUnitStride = IsSingleStride() && x.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0][ind1 + i1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0][x_ind1 + i1] + x_ind2;
#ifdef __USE_AVX__
                    y_e_ax(n2, &a, x_arrayData_i0i1, arrayData_i0i1);
#else
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2] = a * x_arrayData_i0i1[i2];
                    }
#endif
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0*s0][ind1 + i1*s1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0*x_s0][x_ind1 + i1*x_s1] + x_ind2;
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2*s2] = a * x_arrayData_i0i1[i2*x_s2];
                    }
                }
            }
        }
        return *this;
    }

    NumberArray3D& Add_aXY(T a, const NumberArray3D& x, const NumberArray3D& y) {
        assert( x.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** x_arrayData = x.GetArrayData();
        const std::array<std::size_t, 3>& x_indStart = x.GetIndStart();
        std::size_t x_ind0 = x_indStart[0];
        std::size_t x_ind1 = x_indStart[1];
        std::size_t x_ind2 = x_indStart[2];
        const std::array<std::size_t, 3>& x_stride = x.GetStride();
        std::size_t x_s0 = x_stride[0];
        std::size_t x_s1 = x_stride[1];
        std::size_t x_s2 = x_stride[2];

        T*** y_arrayData = y.GetArrayData();
        const std::array<std::size_t, 3>& y_indStart = y.GetIndStart();
        std::size_t y_ind0 = y_indStart[0];
        std::size_t y_ind1 = y_indStart[1];
        std::size_t y_ind2 = y_indStart[2];
        const std::array<std::size_t, 3>& y_stride = y.GetStride();
        std::size_t y_s0 = y_stride[0];
        std::size_t y_s1 = y_stride[1];
        std::size_t y_s2 = y_stride[2];

        T* arrayData_i0i1;
        T*  x_arrayData_i0i1;
        T*  y_arrayData_i0i1;

        bool isUnitStride = IsSingleStride() && x.IsSingleStride() && y.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0][ind1 + i1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0][x_ind1 + i1] + x_ind2;
                    y_arrayData_i0i1 = y_arrayData[y_ind0 + i0][y_ind1 + i1] + y_ind2;
#ifdef __USE_AVX__
                    z_pe_axy(n2, &a, x_arrayData_i0i1, y_arrayData_i0i1, arrayData_i0i1);
#else
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2] += a * x_arrayData_i0i1[i2] * y_arrayData_i0i1[i2];
                    }
#endif
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0*s0][ind1 + i1*s1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0*x_s0][x_ind1 + i1*x_s1] + x_ind2;
                    y_arrayData_i0i1 = y_arrayData[y_ind0 + i0*y_s0][y_ind1 + i1*y_s1] + y_ind2;
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2*s2] += a * x_arrayData_i0i1[i2*x_s2] * y_arrayData_i0i1[i2*y_s2];
                    }
                }
            }
        }
        return *this;
    }

    NumberArray3D& Equate_aXY(T a, const NumberArray3D& x, const NumberArray3D& y) {
        assert( x.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** x_arrayData = x.GetArrayData();
        const std::array<std::size_t, 3>& x_indStart = x.GetIndStart();
        std::size_t x_ind0 = x_indStart[0];
        std::size_t x_ind1 = x_indStart[1];
        std::size_t x_ind2 = x_indStart[2];
        const std::array<std::size_t, 3>& x_stride = x.GetStride();
        std::size_t x_s0 = x_stride[0];
        std::size_t x_s1 = x_stride[1];
        std::size_t x_s2 = x_stride[2];

        T*** y_arrayData = y.GetArrayData();
        const std::array<std::size_t, 3>& y_indStart = y.GetIndStart();
        std::size_t y_ind0 = y_indStart[0];
        std::size_t y_ind1 = y_indStart[1];
        std::size_t y_ind2 = y_indStart[2];
        const std::array<std::size_t, 3>& y_stride = y.GetStride();
        std::size_t y_s0 = y_stride[0];
        std::size_t y_s1 = y_stride[1];
        std::size_t y_s2 = y_stride[2];

        T* arrayData_i0i1;
        T*  x_arrayData_i0i1;
        T*  y_arrayData_i0i1;

        bool isUnitStride = IsSingleStride() && x.IsSingleStride() && y.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0][ind1 + i1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0][x_ind1 + i1] + x_ind2;
                    y_arrayData_i0i1 = y_arrayData[y_ind0 + i0][y_ind1 + i1] + y_ind2;
#ifdef __USE_AVX__
                    z_e_axy(n2, &a, x_arrayData_i0i1, y_arrayData_i0i1, arrayData_i0i1);
#else
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2] = a * x_arrayData_i0i1[i2] * y_arrayData_i0i1[i2];
                    }
#endif
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    arrayData_i0i1 = arrayData[ind0 + i0*s0][ind1 + i1*s1] + ind2;
                    x_arrayData_i0i1 = x_arrayData[x_ind0 + i0*x_s0][x_ind1 + i1*x_s1] + x_ind2;
                    y_arrayData_i0i1 = y_arrayData[y_ind0 + i0*y_s0][y_ind1 + i1*y_s1] + y_ind2;
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData_i0i1[i2*s2] = a * x_arrayData_i0i1[i2*x_s2] * y_arrayData_i0i1[i2*y_s2];
                    }
                }
            }
        }
        return *this;
    }

    //------------------------------ Operators ---------------------------------------
    //-------------------------- automatic substitution (start) ----------------------


    NumberArray3D& operator=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        const std::array<std::size_t, 3>& rhs_stride = rhs.GetStride();
        std::size_t rhs_s0 = rhs_stride[0];
        std::size_t rhs_s1 = rhs_stride[1];
        std::size_t rhs_s2 = rhs_stride[2];

        bool isUnitStride = IsSingleStride() && rhs.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] =
                                rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] =
                                rhs_arrayData[rhs_ind0 + i0*rhs_s0][rhs_ind1 + i1*rhs_s1][rhs_ind2 + i2*rhs_s2];
                    }
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        if(isSingleStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = rhs;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] = rhs;
                    }
                }
            }
        }
        return *this;
    }


    NumberArray3D& operator+=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        const std::array<std::size_t, 3>& rhs_stride = rhs.GetStride();
        std::size_t rhs_s0 = rhs_stride[0];
        std::size_t rhs_s1 = rhs_stride[1];
        std::size_t rhs_s2 = rhs_stride[2];

        bool isUnitStride = IsSingleStride() && rhs.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] +=
                                rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] +=
                                rhs_arrayData[rhs_ind0 + i0*rhs_s0][rhs_ind1 + i1*rhs_s1][rhs_ind2 + i2*rhs_s2];
                    }
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
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        if(isSingleStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] += rhs;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] += rhs;
                    }
                }
            }
        }
        return *this;
    }


    NumberArray3D& operator-=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        const std::array<std::size_t, 3>& rhs_stride = rhs.GetStride();
        std::size_t rhs_s0 = rhs_stride[0];
        std::size_t rhs_s1 = rhs_stride[1];
        std::size_t rhs_s2 = rhs_stride[2];

        bool isUnitStride = IsSingleStride() && rhs.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] -=
                                rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] -=
                                rhs_arrayData[rhs_ind0 + i0*rhs_s0][rhs_ind1 + i1*rhs_s1][rhs_ind2 + i2*rhs_s2];
                    }
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator-=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        if(isSingleStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] -= rhs;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] -= rhs;
                    }
                }
            }
        }
        return *this;
    }


    NumberArray3D& operator*=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        const std::array<std::size_t, 3>& rhs_stride = rhs.GetStride();
        std::size_t rhs_s0 = rhs_stride[0];
        std::size_t rhs_s1 = rhs_stride[1];
        std::size_t rhs_s2 = rhs_stride[2];

        bool isUnitStride = IsSingleStride() && rhs.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *=
                                rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] *=
                                rhs_arrayData[rhs_ind0 + i0*rhs_s0][rhs_ind1 + i1*rhs_s1][rhs_ind2 + i2*rhs_s2];
                    }
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
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        if(isSingleStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= rhs;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] *= rhs;
                    }
                }
            }
        }
        return *this;
    }


    NumberArray3D& operator/=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape );
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        const std::array<std::size_t, 3>& rhs_stride = rhs.GetStride();
        std::size_t rhs_s0 = rhs_stride[0];
        std::size_t rhs_s1 = rhs_stride[1];
        std::size_t rhs_s2 = rhs_stride[2];

        bool isUnitStride = IsSingleStride() && rhs.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /=
                                rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] /=
                                rhs_arrayData[rhs_ind0 + i0*rhs_s0][rhs_ind1 + i1*rhs_s1][rhs_ind2 + i2*rhs_s2];
                    }
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
        std::size_t s0 = stride[0];
        std::size_t s1 = stride[1];
        std::size_t s2 = stride[2];

        if(isSingleStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /= rhs;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        arrayData[ind0 + i0*s0][ind1 + i1*s1][ind2 + i2*s2] /= rhs;
                    }
                }
            }
        }
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        const std::array<std::size_t, 3>& b_stride = numArrB.GetStride();
        std::size_t b_s0 = b_stride[0];
        std::size_t b_s1 = b_stride[1];
        std::size_t b_s2 = b_stride[2];


        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride() && numArrB.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] +
                                b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] +
                                b_arrayData[b_ind0 + i0*b_s0][b_ind1 + i1*b_s1][b_ind2 + i2*b_s2];
                    }
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] + numB;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] + numB;
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator+(const T numB, const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB + a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB + a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2];
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }


    friend NumberArray3D operator-(const NumberArray3D& numArrA, const NumberArray3D& numArrB) {
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        const std::array<std::size_t, 3>& b_stride = numArrB.GetStride();
        std::size_t b_s0 = b_stride[0];
        std::size_t b_s1 = b_stride[1];
        std::size_t b_s2 = b_stride[2];


        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride() && numArrB.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] -
                                b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] -
                                b_arrayData[b_ind0 + i0*b_s0][b_ind1 + i1*b_s1][b_ind2 + i2*b_s2];
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator-(const NumberArray3D& numArrA, const T numB) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] - numB;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] - numB;
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator-(const T numB, const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB - a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB - a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2];
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        const std::array<std::size_t, 3>& b_stride = numArrB.GetStride();
        std::size_t b_s0 = b_stride[0];
        std::size_t b_s1 = b_stride[1];
        std::size_t b_s2 = b_stride[2];


        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride() && numArrB.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] *
                                b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] *
                                b_arrayData[b_ind0 + i0*b_s0][b_ind1 + i1*b_s1][b_ind2 + i2*b_s2];
                    }
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] * numB;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] * numB;
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    friend NumberArray3D operator*(const T numB, const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB * a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB * a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2];
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        const std::array<std::size_t, 3>& b_stride = numArrB.GetStride();
        std::size_t b_s0 = b_stride[0];
        std::size_t b_s1 = b_stride[1];
        std::size_t b_s2 = b_stride[2];


        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride() && numArrB.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] /
                                b_arrayData[b_ind0 + i0][b_ind1 + i1][b_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] /
                                b_arrayData[b_ind0 + i0*b_s0][b_ind1 + i1*b_s1][b_ind2 + i2*b_s2];
                    }
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2] / numB;
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2] / numB;
                    }
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

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        bool isUnitStride = numArrA.IsSingleStride();
        if(isUnitStride) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB / a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2];
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] =
                                numB / a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2];
                    }
                }
            }
        }
        return numArrC;     // TODO: define move constructors
    }

    //-------------------------- automatic substitution (end) ----------------------


    T& operator[](const std::array<std::size_t, 3> indx) {
        if(isSingleStride) {
            return arrayData[indStart[0] + indx[0]][indStart[1] + indx[1]][indStart[2] + indx[2]];
        } else{
            return arrayData[indStart[0] + indx[0]*stride[0]]
                            [indStart[1] + indx[1]*stride[1]]
                            [indStart[2] + indx[2]*stride[2]];
        }
    }

    const T& operator[](const std::array<std::size_t, 3> indx) const {
        if(isSingleStride) {
            return arrayData[indStart[0] + indx[0]][indStart[1] + indx[1]][indStart[2] + indx[2]];
        } else{
            return arrayData[indStart[0] + indx[0]*stride[0]]
                            [indStart[1] + indx[1]*stride[1]]
                            [indStart[2] + indx[2]*stride[2]];
        }
    }

    //------------------------- in-place math functions ----------------------------------

//    NumberArray3D& SetToNumber(const T num) {
//        std::size_t n0 = shape[0];
//        std::size_t n1 = shape[1];
//        std::size_t n2 = shape[2];
//        std::size_t ind0 = indStart[0];
//        std::size_t ind1 = indStart[1];
//        std::size_t ind2 = indStart[2];
//
//        for(std::size_t i0 = 0; i0 < n0; ++i0) {
//            for(std::size_t i1 = 0; i1 < n1; ++i1) {
//                for(std::size_t i2 = 0; i2 < n2; ++i2) {
//                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = num;
//                }
//            }
//        }
//        return *this;
//    }

    //------------------------- special functions -----------------------------------
    //-------------------------- automatic substitution (start) ----------------------


    static NumberArray3D exp(const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        if(numArrA.IsSingleStride()) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::exp(a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2]);
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::exp(a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2]);
                    }
                }
            }
        }
        return numArrC;
    }


    static NumberArray3D cos(const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        if(numArrA.IsSingleStride()) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::cos(a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2]);
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::cos(a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2]);
                    }
                }
            }
        }
        return numArrC;
    }


    static NumberArray3D sin(const NumberArray3D& numArrA) {
        const std::array<std::size_t, 3>& a_shape = numArrA.GetShape();
        T*** a_arrayData = numArrA.GetArrayData();
        const std::array<std::size_t, 3>& a_indStart = numArrA.GetIndStart();
        std::size_t a_ind0 = a_indStart[0];
        std::size_t a_ind1 = a_indStart[1];
        std::size_t a_ind2 = a_indStart[2];

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        const std::array<std::size_t, 3>& a_stride = numArrA.GetStride();
        std::size_t a_s0 = a_stride[0];
        std::size_t a_s1 = a_stride[1];
        std::size_t a_s2 = a_stride[2];

        NumberArray3D numArrC(a_shape, 0);  // TODO: use an uninitialized array
        T*** c_arrayData = numArrC.GetArrayData();

        if(numArrA.IsSingleStride()) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::sin(a_arrayData[a_ind0 + i0][a_ind1 + i1][a_ind2 + i2]);
                    }
                }
            }
        } else {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        c_arrayData[i0][i1][i2] = std::sin(a_arrayData[a_ind0 + i0*a_s0][a_ind1 + i1*a_s1][a_ind2 + i2*a_s2]);
                    }
                }
            }
        }
        return numArrC;
    }

    //-------------------------- automatic substitution (end) ----------------------

    //---------------------   meshgrid

    static NumberArray3D GetMeshGrid(const std::array<std::size_t, 3>& a_shape,
                                     const std::array<FPNumber, 3>& r_min, const std::array<FPNumber, 3>& r_max,
                                     int direction) {
        NumberArray3D numArrA(a_shape, 0);
        T*** a_arrayData = numArrA.GetArrayData();

        std::size_t n0 = a_shape[0];
        std::size_t n1 = a_shape[1];
        std::size_t n2 = a_shape[2];

        std::array<FPNumber, 3> d_max;
        for(int i = 0; i < 3; ++i) {
            d_max[i] = (r_max[i] - r_min[i])/(FPNumber)(a_shape[i] - 1);
        }

        std::array<std::size_t, 3> i;
        std::size_t& i0 = i[0];
        std::size_t& i1 = i[1];
        std::size_t& i2 = i[2];
        for(i0 = 0; i0 < n0; ++i0) {
            for(i1 = 0; i1 < n1; ++i1) {
                for(i2 = 0; i2 < n2; ++i2) {
                    a_arrayData[i0][i1][i2] = r_min[direction] + (FPNumber)(i[direction])*d_max[direction];
                }
            }
        }
        return numArrA;
    }

    //-----------------------------  bufferIO ----------------
    void WriteArrayDataToBuffer(FPNumber* buffer,
                                std::size_t bufferSize,     // maximum buffer size
                                std::size_t indStartBuffer  // write starting this index
                                ) {
        assert(arrayData != nullptr);
        std::size_t totalElements = shape[0]*shape[1]*shape[2];
        assert(indStartBuffer + totalElements <= bufferSize);
        FPNumber* bufferStart = &(buffer[indStartBuffer]);
        ndarray::buffer::Write3D(bufferStart, shape, indStart, arrayData);
    }

    void ReadArrayDatafromBuffer(FPNumber* buffer,
                                 std::size_t bufferSize,     // maximum buffer size
                                 std::size_t indStartBuffer  // read starting this index
                                 ) {
        assert(arrayData != nullptr);
        std::size_t totalElements = shape[0]*shape[1]*shape[2];
        assert(indStartBuffer + totalElements <= bufferSize);
        FPNumber* bufferStart = &(buffer[indStartBuffer]);
        ndarray::buffer::Read3D(bufferStart, shape, indStart, arrayData);
    }

    //------------------------------  ostream  -----------

    friend std::ostream& operator<<(std::ostream& out, const NumberArray3D& numArr)
    {
        Print3DNumberArrayToOstream(out, numArr.GetShape(), numArr.GetArrayData());
        return out;
    }

    //------------------------------------------------------------------------------------

    void Print() {
        Print3DNumberArray(shape, arrayData);
    }

    //-----------------------------------------------------------------------------------
    void WriteArrayDataToFile(std::ofstream* fileOut,
                              bool writeShape = false,
                              bool writeDataTypeSize = false) {
        assert(arrayData != nullptr);
        int dataTypeCode = UtilityFunctions::GetDatatypeNumericalCode<T>();

        return Write3DNumberArrayData(fileOut, shape, indStart, arrayData, dataTypeCode, writeShape, writeDataTypeSize);
    }

    void WriteArrayDataToMemory(char* buffer,
                              std::size_t& bufferInd,  //start writing to buffer from this index
                              bool writeShape = false,
                              bool writeDataTypeSize = false) {
        assert(arrayData != nullptr);
        int dataTypeCode = UtilityFunctions::GetDatatypeNumericalCode<T>();

        return Write3DNumberArrayDataToMemory(buffer, bufferInd,
                    shape, indStart, arrayData, dataTypeCode, writeShape, writeDataTypeSize);
    }

    // This function should take into account all the possible preambles that could be added during output writes
    std::size_t GetMaxDataSizeInBytes() {
        return  sizeof(int)             // datatype code
                + sizeof(std::size_t)   // datatype size
                + 3*sizeof(std::size_t)     // shape
                + shape[0]*shape[1]*shape[2]*sizeof(T);     // data array
    }

    void ReadArrayDataFromFile(std::ifstream* fileIn) {
        assert(arrayData != nullptr);
        return Read3DNumberArrayData(fileIn, shape, indStart, arrayData);
    }
};


#endif  // FDTD_MULTIDIMARRAY_H_
