
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t
#include <string>
#include <unordered_map>
#include <utility>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"
#include "FDInstructionCode.h"
#include "ElementType.h"
#include "GridArrayManipulator.h"


class YeeGrid3D {
    public:
    YeeGrid3D(std::array<std::size_t, 3>& nCells);
    ~YeeGrid3D();

    void SetCornerCoordinates(std::array<RealNumber, 3> r_0, std::array<RealNumber, 3> r_1);
    const std::array<RealNumber, 3>& GetCornerR0() const;
    const std::array<RealNumber, 3>& GetCornerR1() const;
    const std::array<std::size_t, 3>& GetNumberOfCells() const;

    // an element that spans over the entire Yee grid
    void AddEntireGridElement(const std::string name, ElementType elemType);
    // an element that spans over a fraction of the Yee grid
    void AddPartialGridElement(const std::string name, ElementType elemType
        ,std::array<std::size_t, 3> indOrigin     // the element start on this index of the background grid
        ,std::array<std::size_t, 3> numCells      // number of cells covered by the element
        );
    YeeGridData3D& GetGridElement(const std::string name);

    void AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params);
    void SetIterationSequence(std::vector<std::string> sequence);
    void ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params);
    void ApplyUpdateInstructions(std::size_t numIterations);

    void* ConstructParams_A_plusequal_sum_b_C(
                                    std::array<std::size_t, 3> ind_start,
                                    std::array<std::size_t, 3> ind_end,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<RealNumber> bValues,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<int, 3>> arrayC_indexShifts
                                    );

    void AddGaussianPointSource(int direction, std::array<std::size_t, 3> index, RealNumber amplitude,
            RealNumber t_center, RealNumber t_decay, RealNumber modulationFrequecy, RealNumber modulatioPhase);

    private:
    std::array<RealNumber, 3> r_0;  // coordinates of the lower left corner
    std::array<RealNumber, 3> r_1;  // coordinates of the upper right corner
    std::array<std::size_t, 3> nCells;      // number of Yee cells
    std::unordered_map<std::string, std::unique_ptr<YeeGridData3D>> gridElements;
    std::unordered_map<std::string, std::unique_ptr<GridArrayManipulator>> gridArrayManipulators;
    std::unordered_map<std::string, std::pair<FDInstructionCode, void*>> instructions;  // field update instructions
    std::vector<std::string> iterationSequence;     // sequence in which to apply the field update instructions

};


#endif  // FDTD_YEEGRID_H_


