
#include "YeeGridParallel.h"


void YeeGrid3DParallel::SetInBuffers(std::unordered_map<std::string, RealNumber*>* buffers) {
    inBuffers = buffers;
}

void YeeGrid3DParallel::SetOutBuffers(std::unordered_map<std::string, RealNumber*>* buffers) {
    outBuffers = buffers;
}


