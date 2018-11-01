
#ifndef FDTD_YEEGRIDCOLLECTIONPARALLEL_H_
#define FDTD_YEEGRIDCOLLECTIONPARALLEL_H_

#include <cstddef>
#include <string>
#include <vector>

#include "NumberTypes.h"
#include "YeeGrid.h"
#include "YeeGridCollection.h"

class YeeGridCollectionParallel : public YeeGridCollection {
    public:
    ~YeeGridCollectionParallel();

    void AddInBuffer(std::string name, std::size_t bufferSize);
    void AddOutBuffer(std::string name, std::size_t bufferSize);


    private:
    std::unordered_map<std::string, RealNumber*> inBuffers;
    std::unordered_map<std::string, RealNumber*> outBuffers;
};


#endif // FDTD_YEEGRIDCOLLECTIONPARALLEL_H_


