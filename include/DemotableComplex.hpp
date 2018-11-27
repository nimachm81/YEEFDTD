

#ifndef FDTD_DEMOTABELCOMPLEX_H_
#define FDTD_DEMOTABELCOMPLEX_H_


template<typename T>
class DemotableComplex : public std::complex<T> {
    public:
    DemotableComplex() : std::complex<T>() {}
    DemotableComplex(T r) : std::complex<T>(r) {}
    DemotableComplex(T r, T i) : std::complex<T>(r, i) {}
    DemotableComplex(std::complex<T> c) : std::complex<T>(c) {}
    explicit operator T() const {return this->real();}
};


#endif  // FDTD_DEMOTABELCOMPLEX_H_
