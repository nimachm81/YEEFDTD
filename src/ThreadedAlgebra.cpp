
#include "ThreadedAlgebra.h"

ThreadedAlgebra::ThreadedAlgebra(std::size_t n_threads) : numThreads(n_threads),
                                                          startOperation(n_threads)
                                                          {
    for(std::size_t i = 0; i < numThreads; ++i) {
        threads.push_back(std::thread(&ThreadedAlgebra::CalculateForever, this));

        startOperation[i].store(0, std::memory_order::memory_order_seq_cst);
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

    const std::size_t threadIndex = idMap[std::this_thread::get_id()];
    std::atomic<int>& startIt = startOperation[threadIndex];

    while(true) {
        do {

        } while(!startIt.load(std::memory_order::memory_order_seq_cst));

        if(terminateAll.load(std::memory_order::memory_order_seq_cst)) {
            break;
        }

        if(instruction == TAInstructionCode::y_pe_ax) {
            NumberArray3D<FPNumber>& y = arrASlices[threadIndex];
            FPNumber a = numA;
            NumberArray3D<FPNumber>& x = arrBSlices[threadIndex];

            y.Add_aX(a, x);
        } else if(instruction == TAInstructionCode::y_e_ax) {
            NumberArray3D<FPNumber>& y = arrASlices[threadIndex];
            FPNumber a = numA;
            NumberArray3D<FPNumber>& x = arrBSlices[threadIndex];

            y.Equate_aX(a, x);
        } else if(instruction == TAInstructionCode::z_pe_axy) {
            NumberArray3D<FPNumber>& z = arrASlices[threadIndex];
            FPNumber a = numA;
            NumberArray3D<FPNumber>& x = arrBSlices[threadIndex];
            NumberArray3D<FPNumber>& y = arrCSlices[threadIndex];

            z.Add_aXY(a, x, y);
        } else if(instruction == TAInstructionCode::z_e_axy) {
            NumberArray3D<FPNumber>& z = arrASlices[threadIndex];
            FPNumber a = numA;
            NumberArray3D<FPNumber>& x = arrBSlices[threadIndex];
            NumberArray3D<FPNumber>& y = arrCSlices[threadIndex];

            z.Equate_aXY(a, x, y);
        }

        startIt.store(0, std::memory_order::memory_order_seq_cst);
        numOfThreadsProcessed.fetch_add(1, std::memory_order::memory_order_seq_cst);
        //std::cout << numOfThreadsProcessed << std::endl;
    }
}

void ThreadedAlgebra::SetArgsFlagsAndSignalStart() {
    for (std::size_t i = 0; i < numThreads; ++i) {
        startOperation[i].store(1, std::memory_order::memory_order_seq_cst);
    }
}

void ThreadedAlgebra::WaitForResults() {
    int numProcessed;
    do {
        numProcessed = numOfThreadsProcessed.load(std::memory_order::memory_order_seq_cst);
    } while(numProcessed < numThreads);
    numOfThreadsProcessed.store(0, std::memory_order::memory_order_seq_cst);
}

void ThreadedAlgebra::Get_y_pe_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {

    bool arraySizeLargeEnough = Set_y_pe_ax_args(y, a, x);

    if(arraySizeLargeEnough) {
        SetArgsFlagsAndSignalStart();
        WaitForResults();
    } else {
        y.Add_aX(a, x);
        auto shape = x.GetShape();
        //std::cout << " did not parallelize .. " << shape[0] << " " << shape[1] << " " << shape[2] << std::endl;
    }
}

void ThreadedAlgebra::Get_y_e_ax(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {

    bool arraySizeLargeEnough = Set_y_e_ax_args(y, a, x);

    if(arraySizeLargeEnough) {
        SetArgsFlagsAndSignalStart();
        WaitForResults();
    } else {
        y.Equate_aX(a, x);
        //std::cout << " did not parallelize .. " << std::endl;
    }
}

void ThreadedAlgebra::Get_z_pe_axy(NumberArray3D<FPNumber>& z,
                                   FPNumber a,
                                   NumberArray3D<FPNumber>& x,
                                   NumberArray3D<FPNumber>& y) {

    bool arraySizeLargeEnough = Set_z_pe_axy_args(z, a, x, y);

    if(arraySizeLargeEnough) {
        SetArgsFlagsAndSignalStart();
        WaitForResults();
    } else {
        z.Add_aXY(a, x, y);
        //std::cout << " did not parallelize .. " << std::endl;
    }
}

void ThreadedAlgebra::Get_z_e_axy(NumberArray3D<FPNumber>& z,
                                  FPNumber a,
                                  NumberArray3D<FPNumber>& x,
                                  NumberArray3D<FPNumber>& y) {

    bool arraySizeLargeEnough = Set_z_e_axy_args(z, a, x, y);

    if(arraySizeLargeEnough) {
        SetArgsFlagsAndSignalStart();
        WaitForResults();
    } else {
        z.Equate_aXY(a, x, y);
        //std::cout << " did not parallelize .. " << std::endl;
    }
}

bool ThreadedAlgebra::Set_y_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {
    instruction = TAInstructionCode::y_pe_ax;
    return Set_y_e_or_pe_ax_args(y, a, x);
}

bool ThreadedAlgebra::Set_y_e_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {
    instruction = TAInstructionCode::y_e_ax;
    return Set_y_e_or_pe_ax_args(y, a, x);
}

bool ThreadedAlgebra::Set_y_e_or_pe_ax_args(NumberArray3D<FPNumber>& y, FPNumber a, NumberArray3D<FPNumber>& x) {
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
    } else {
        n_th_x = 1;
        n_th_y = numThreads;
    }

    if( shape[0] < n_th_x || shape[1] < n_th_y ) {
        return false;

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

        return true;
    }
}

bool ThreadedAlgebra::Set_z_pe_axy_args(NumberArray3D<FPNumber>& z,
                                        FPNumber a,
                                        NumberArray3D<FPNumber>& x,
                                        NumberArray3D<FPNumber>& y) {
    instruction = TAInstructionCode::z_pe_axy;
    return Set_z_e_or_pe_axy_args(z, a, x, y);
}

bool ThreadedAlgebra::Set_z_e_axy_args(NumberArray3D<FPNumber>& z,
                                       FPNumber a,
                                       NumberArray3D<FPNumber>& x,
                                       NumberArray3D<FPNumber>& y) {
    instruction = TAInstructionCode::z_e_axy;
    return Set_z_e_or_pe_axy_args(z, a, x, y);
}


bool ThreadedAlgebra::Set_z_e_or_pe_axy_args(NumberArray3D<FPNumber>& z,
                                             FPNumber a,
                                             NumberArray3D<FPNumber>& x,
                                             NumberArray3D<FPNumber>& y
                                             ) {
    arrASlices.clear();
    arrBSlices.clear();
    arrCSlices.clear();
    numA = a;

    std::array<std::size_t, 3> shape = z.GetShape();
    assert( x.GetShape() == shape && y.GetShape() == shape );

    std::size_t n_th_x, n_th_y;
    if(numThreads%2 == 0) {
        n_th_x = 2;
        n_th_y = numThreads / n_th_x;
    } else if(numThreads%3 == 0) {
        n_th_x = 3;
        n_th_y = numThreads / n_th_x;
    }

    if( shape[0] < n_th_x || shape[1] < n_th_y ) {
        return false;   // use single thread

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
                         z.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                    std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                    )
                                    );
                arrBSlices.push_back(
                         x.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                    std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                    )
                                    );
                arrCSlices.push_back(
                         y.GetSlice(std::array<std::size_t, 3>{ind0_i,        ind1_j,        0},
                                    std::array<std::size_t, 3>{ind0_i + n0_i, ind1_j + n1_j, n2}
                                    )
                                    );
            }
        }

        return true;    // use multiple thread
    }
}


void ThreadedAlgebra::Terminate() {
    {
        terminateAll.store(true, std::memory_order::memory_order_seq_cst);
        for (std::size_t i = 0; i < numThreads; ++i) {
            startOperation[i].store(1, std::memory_order::memory_order_seq_cst);
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
