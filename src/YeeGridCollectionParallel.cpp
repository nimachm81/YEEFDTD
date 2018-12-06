
#include "YeeGridCollectionParallel.h"


YeeGridCollectionParallel::~YeeGridCollectionParallel() {
    for(auto it = inBuffers.begin(); it != inBuffers.end(); ++it) {
        if(it->second != nullptr) {
            delete[] it->second;
            it->second = nullptr;
        }
    }
    for(auto it = outBuffers.begin(); it != inBuffers.end(); ++it) {
        if(it->second != nullptr) {
            delete[] it->second;
            it->second = nullptr;
        }
    }
}

void YeeGridCollectionParallel::AddInBuffer(std::string name, std::size_t bufferSize) {
    FPNumber* buffer = new FPNumber[bufferSize];
    inBuffers[name] = buffer;
}

void YeeGridCollectionParallel::AddOutBuffer(std::string name, std::size_t bufferSize) {
    FPNumber* buffer = new FPNumber[bufferSize];
    outBuffers[name] = buffer;
}


