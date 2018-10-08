
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
    YeeGrid3D() { };
    ~YeeGrid3D();

    void SetCornerCoordinates(std::array<RealNumber, 3> r_0, std::array<RealNumber, 3> r_1);
    void SetNumOfCells(std::array<std::size_t, 3>& nCells);
    void SetTimeResolution(const RealNumber dt);
    const std::array<RealNumber, 3>& GetCornerR0() const;
    const std::array<RealNumber, 3>& GetCornerR1() const;
    const std::array<std::size_t, 3>& GetNumberOfCells() const;

    // an element that spans over the entire Yee grid
    void AddEntireGridElement(const std::string name, ElementType elemType);
    // an element that spans over a fraction of the Yee grid
    void AddPartialGridElement(const std::string name, ElementType elemType
        ,std::array<std::size_t, 3> startCell     // the element start on this cell index of the background grid
        ,std::array<std::size_t, 3> numCells      // number of cells covered by the element
        );
    YeeGridData3D& GetGridElement(const std::string name);

    void AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params);
    void SetIterationSequence(std::vector<std::string> sequence);
    void ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params);
    void ApplyUpdateInstructions(std::size_t numIterations);

    void* ConstructParams_A_plusequal_sum_b_C(
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<RealNumber> bValues,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                    );
    void* ConstructParams_A_plusequal_sum_B_C(
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<std::string> arrayB_names,
                                    std::vector<int> arrayB_components,
                                    std::vector<std::array<std::size_t, 3>> arrayB_indsStart,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                    );
    void* ConstructParams_A_equal_func_r_t(
                                    std::string gridManipulator_name
                                    );

    void AddGaussianPointSource(const std::string name,
            const std::string gridDataName,     // name of the gridDataObject whose data is manipulated by this point source
            int direction, RealNumber amplitude,
            RealNumber t_center, RealNumber t_decay, RealNumber modulationFrequecy,
            RealNumber modulatioPhase, RealNumber timeOffsetFraction);

    void PrintAllGridData();

    void AddGridElementView(std::string gridElemViewName,   // name of the gridView
                            std::string gridElemName , int gridElemComponent,   // name of the gridElement and its x,y,z component
                            std::array<std::size_t, 3> indStart, std::array<std::size_t, 3> indEnd // slice start and end
                            );
    void AddFullGridElementView(std::string gridElemViewName,   // name of the gridView
                            std::string gridElemName , int gridElemComponent   // name of the gridElement and its x,y,z component
                            );
    void WriteGridDataToFile(std::string fileName, std::string gridElemViewName);
    void WriteAllGridElemViewsToFile();
    void DeleteOlderViewFiles();
    void SetDataStoreRate(std::size_t saveEveryNSammples);

    private:
    std::size_t timeIndex;
    RealNumber dt;                  // time resolution
    std::array<RealNumber, 3> dr;   // spatial resolution
    std::array<RealNumber, 3> r_0{0, 0, 0};  // coordinates of the lower left corner
    std::array<RealNumber, 3> r_1{0, 0, 0};  // coordinates of the upper right corner
    std::array<std::size_t, 3> nCells;      // number of Yee cells
    std::unordered_map<std::string, std::shared_ptr<YeeGridData3D>> gridElements;
    std::unordered_map<std::string, std::shared_ptr<GridArrayManipulator>> gridArrayManipulators;
    std::unordered_map<std::string, std::pair<FDInstructionCode, void*>> instructions;  // field update instructions
    std::vector<std::string> iterationSequence;     // sequence in which to apply the field update instructions
    std::unordered_map<std::string, NumberArray3D<RealNumber>> gridElementViews;  // a slice of gridElements for printing to output
    std::size_t saveDataEveryNTimeSamples = 1;
};


#endif  // FDTD_YEEGRID_H_


