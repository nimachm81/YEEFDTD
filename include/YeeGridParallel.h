
#ifndef FDTD_YEEGRIDPARALLEL_H_
#define FDTD_YEEGRIDPARALLEL_H_

#include "YeeGrid.h"

class YeeGrid3DParallel : public YeeGrid3D {
    public:
    void SetInBuffers(std::unordered_map<std::string, RealNumber*>* buffer);
    void SetOutBuffers(std::unordered_map<std::string, RealNumber*>* buffer);

    private:
    std::unordered_map<std::string, RealNumber*>* inBuffers;
    std::unordered_map<std::string, RealNumber*>* outBuffers;
};

#endif // FDTD_YEEGRIDPARALLEL_H_
