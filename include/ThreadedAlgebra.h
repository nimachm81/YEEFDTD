#ifndef _THREADED_ALGEBRA_
#define _THREADED_ALGEBRA_

#include <vector>
#include <array>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "MultiDimArray.hpp"

enum class TAInstructionCode {
    y_pe_ax,
    y_e_ax,
    z_pe_axy,
    z_e_axy
};

class ThreadedAlgebra {
    public:
    ThreadedAlgebra(std::size_t n_threads = 2);
    ~ThreadedAlgebra();

    void MapThreadIds();
    void CalculateForever();

    void SetArgsFlagsAndSignalStart();
    void WaitForResults();

    void Get_y_pe_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);
    void Get_y_e_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);
    void Get_z_pe_axy(NumberArray3D<FPNumber>& z, FPNumber a, NumberArray3D<FPNumber>& x, NumberArray3D<FPNumber>& y);
    void Get_z_e_axy(NumberArray3D<FPNumber>& z, FPNumber a, NumberArray3D<FPNumber>& x, NumberArray3D<FPNumber>& y);

    bool Set_y_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);
    bool Set_y_e_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);
    bool Set_y_e_or_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);
    bool Set_z_pe_axy_args(NumberArray3D<FPNumber>& z, FPNumber a, NumberArray3D<FPNumber>& x, NumberArray3D<FPNumber>& y);
    bool Set_z_e_axy_args(NumberArray3D<FPNumber>& z, FPNumber a, NumberArray3D<FPNumber>& x, NumberArray3D<FPNumber>& y);
    bool Set_z_e_or_pe_axy_args(NumberArray3D<FPNumber>& z, FPNumber a, NumberArray3D<FPNumber>& x, NumberArray3D<FPNumber>& y);

    void Terminate();

    private:
    const std::size_t numThreads;
    std::atomic<bool> terminateAll = false;
    std::atomic<int> numOfThreadsProcessed = -1;
    std::vector<std::atomic<int>> startOperation;

    std::mutex sharedMutex;

    std::vector<std::thread> threads;
    std::map<std::thread::id, std::size_t> idMap;
    std::size_t lastId = 0;

    TAInstructionCode instruction;

    FPNumber numA;
    std::vector<NumberArray3D<FPNumber>> arrASlices;
    std::vector<NumberArray3D<FPNumber>> arrBSlices;
    std::vector<NumberArray3D<FPNumber>> arrCSlices;
};

#endif // _THREADED_ALGEBRA_
