
#ifndef FDTD_YEEGRIDPARALLEL_H_
#define FDTD_YEEGRIDPARALLEL_H_

#include "FDInstructionCode.h"
#include "YeeGrid.h"

class YeeGrid3DParallel : public YeeGrid3D {
    public:
    ~YeeGrid3DParallel();

    void SetInBuffers(std::unordered_map<std::string, RealNumber*>* buffer);
    void SetOutBuffers(std::unordered_map<std::string, RealNumber*>* buffer);

    void* ConstructParams_A_to_buffer(
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    RealNumber* buffer,
                                    std::size_t bufferSize,
                                    std::size_t ind_start_buffer
                                    ) const;

    void ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params);

    private:
    std::unordered_map<std::string, RealNumber*>* inBuffers = nullptr;
    std::unordered_map<std::string, RealNumber*>* outBuffers = nullptr;
};

#endif // FDTD_YEEGRIDPARALLEL_H_
