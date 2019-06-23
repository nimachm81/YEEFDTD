
#include "ThreadedAlgebra.h"

ThreadedAlgebra::ThreadedAlgebra(std::size_t n_threads) : startConditions(n_threads),
                                                          mutexs(n_threads)
                                                          {
    numThreads = n_threads;

    for(std::size_t i = 0; i < numThreads; ++i) {
        threads.push_back(std::thread(&ThreadedAlgebra::CalculateForever, this));
        argsAreReady.push_back(false);
    }
}

ThreadedAlgebra::~ThreadedAlgebra() {
    std::cout << "inside the Dtor." << std::endl;
    Terminate();
}

void ThreadedAlgebra::MapThreadIds() {
    std::unique_lock<std::mutex> threadlock(sharedMutex);
    idMap[std::this_thread::get_id()] = lastId;
    lastId++;
    threadlock.unlock();
}

void ThreadedAlgebra::CalculateForever() {
    MapThreadIds();

    std::size_t threadIndex = idMap[std::this_thread::get_id()];

    std::unique_lock<std::mutex> threadlock(mutexs[threadIndex]);
    threadlock.unlock();

    while(true) {
        threadlock.lock();
        startConditions[threadIndex].wait(threadlock, [&,this]{return argsAreReady[threadIndex];});

        if(terminateAll) {
            threadlock.unlock();
            break;
        }

        NumberArray3D<FPNumber>& y = arrASlices[threadIndex];
        FPNumber a = numA;
        NumberArray3D<FPNumber>& x = arrBSlices[threadIndex];

        y.Add_aX(a, x);

        {
            std::lock_guard<std::mutex> lk(sharedMutex);
            argsAreReady[threadIndex] = false;
            numOfThreadsProcessed++;

            std::cout << "numOfThreadsProcessed : " << numOfThreadsProcessed << std::endl;
        }

        threadlock.unlock();
        endCondition.notify_one();

        assert(!terminateAll);
    }
}

void ThreadedAlgebra::Get_y_pe_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {

    Set_y_pe_ax_args(y, a, x);

    {
        std::lock_guard<std::mutex> lk(sharedMutex);
        for (std::size_t i = 0; i < numThreads; ++i) {
            argsAreReady[i] = true;
            startConditions[i].notify_one();
        }
        std::cout << "Args are set!" << std::endl;
    }


    {
        std::unique_lock<std::mutex> lk(sharedMutex);
        endCondition.wait(lk, [this]{return numOfThreadsProcessed >= numThreads;});

        for (std::size_t i = 0; i < numThreads; ++i) {
            assert(argsAreReady[i] == false);
        }
        if(numOfThreadsProcessed > numThreads) {
            assert(false);
        }
        numOfThreadsProcessed = 0;
        std::cout << "Results are ready." << std::endl;
    }
}

void ThreadedAlgebra::Set_y_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {
    arrASlices.clear();
    arrBSlices.clear();
    numA = a;

    std::array<std::size_t, 3> shape = y.GetShape();
    assert( x.GetShape() == shape );

    std::size_t n_th_x, n_th_y;
    if(numThreads%2 == 0) {
        n_th_x = 2;
        n_th_y = numThreads / n_th_x;
    } else if(numThreads%3 == 0) {
        n_th_x = 3;
        n_th_y = numThreads / n_th_x;
    }

    if( shape[0] < n_th_x || shape[1] < n_th_y ) {
        y.Add_aX(a, x);
        return;
    } else {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];

        std::size_t dn0 = n0/n_th_x;
        std::size_t dn1 = n1/n_th_y;

        std::size_t ind0_i, ind1_j;
        std::size_t n0_i, n1_j;


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

                arrASlices.push_back(
                          y.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                     std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                    )
                                    );
                arrBSlices.push_back(
                            x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                       std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                       )
                                 );
            }
        }
    }
}


void ThreadedAlgebra::Terminate() {
    {
        std::lock_guard<std::mutex> lk(sharedMutex);
        terminateAll = true;
        for (std::size_t i = 0; i < numThreads; ++i) {
            argsAreReady[i] = true;
            startConditions[i].notify_one();
        }
    }

    std::cout << "joining threads: " << std::endl;
    //std::this_thread::sleep_for(std::chrono::seconds(3));
    for(std::size_t i = 0; i < numThreads; ++i) {
        threads[i].join();
        std::cout << "thread " << i << " joined" << std::endl;
    }
    std::cout << "Done. " << std::endl;
}
