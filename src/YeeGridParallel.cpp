
#include "YeeGridParallel.h"


YeeGrid3DParallel::~YeeGrid3DParallel() {
    for(std::string updateName : iterativeSequence) {
        auto& instructCode_param_pair = updateInstructions[updateName];
        auto& instructionCode = instructCode_param_pair.first;
        void* params = instructCode_param_pair.second;
        if(instructionCode == FDInstructionCode::A_to_buffer ||
                  instructionCode == FDInstructionCode::A_from_buffer) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                        std::string,    // 1
                        int,    // 2
                        FPNumber*,   // 3
                        std::size_t,    // 4
                        std::size_t     // 5
                    >*
                >(params);
            delete params_tuple;
        }
    }
}

void YeeGrid3DParallel::SetInBuffers(std::unordered_map<std::string, FPNumber*>* buffers) {
    inBuffers = buffers;
}

void YeeGrid3DParallel::SetOutBuffers(std::unordered_map<std::string, FPNumber*>* buffers) {
    outBuffers = buffers;
}

void* YeeGrid3DParallel::ConstructParams_A_to_buffer(
                            std::array<std::size_t, 3> ind_start_A,
                            std::array<std::size_t, 3> ind_end_A,
                            std::string arrayA_name,
                            int arrayA_component,
                            FPNumber* buffer,
                            std::size_t bufferSize,
                            std::size_t ind_start_buffer
                            ) const {
    auto* params_tuple =
        new std::tuple<
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
            std::string,    // 1
            int,    // 2
            FPNumber*,   // 3
            std::size_t,    // 4
            std::size_t     // 5
        >(
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>(ind_start_A, ind_end_A),  // 0
            arrayA_name,        // 1
            arrayA_component,   // 2
            buffer,             // 3
            bufferSize,         // 4
            ind_start_buffer    // 5
        );
     return static_cast<void*>(params_tuple);
}


void YeeGrid3DParallel::ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params) {
    YeeGrid3D::ApplyUpdateInstruction(instructionCode, params);
    if(instructionCode == FDInstructionCode::A_to_buffer ||
            instructionCode == FDInstructionCode::A_from_buffer) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    FPNumber*,    // 3
                    std::size_t,    // 4
                    std::size_t     // 5
                >*
            >(params);
        std::array<std::size_t, 3>& ind_start_A = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end_A = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        FPNumber* buffer = std::get<3>(params_tuple);
        std::size_t bufferSize = std::get<4>(params_tuple);
        std::size_t ind_start_buffer = std::get<5>(params_tuple);

        NumberArray3D<FPNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<FPNumber> arrayASlice = arrayA.GetSlice(ind_start_A, ind_end_A);

        if(instructionCode == FDInstructionCode::A_to_buffer) {
            arrayASlice.WriteArrayDataToBuffer(buffer, bufferSize, ind_start_buffer);
        } else if(instructionCode == FDInstructionCode::A_from_buffer) {
            arrayASlice.ReadArrayDatafromBuffer(buffer, bufferSize, ind_start_buffer);
        }
    }
}
