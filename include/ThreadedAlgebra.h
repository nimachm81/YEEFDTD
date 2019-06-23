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

class ThreadedAlgebra {
    public:
    ThreadedAlgebra(std::size_t n_threads = 4);
    ~ThreadedAlgebra();

    void MapThreadIds();
    void CalculateForever();

    void Get_y_pe_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);

    void Set_y_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x);

    void Terminate();

    private:
    std::size_t numThreads;
    std::condition_variable endCondition;
    std::vector<bool> argsAreReady;
    std::vector<std::condition_variable> startConditions;
    std::vector<std::mutex> mutexs;
    std::mutex sharedMutex;
    std::atomic<bool> terminateAll = false;
    std::atomic<std::size_t> numOfThreadsProcessed = 0;

    std::vector<std::thread> threads;
    std::map<std::thread::id, std::size_t> idMap;
    std::size_t lastId = 0;

    FPNumber numA;
    std::vector<NumberArray3D<FPNumber>> arrASlices;
    std::vector<NumberArray3D<FPNumber>> arrBSlices;
};

#endif // _THREADED_ALGEBRA_
