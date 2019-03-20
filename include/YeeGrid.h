
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"
#include "FDInstructionCode.h"
#include "ElementType.h"
#include "GridArrayManipulator.h"
#include "GridElementView.h"
#include "DiscretePointsGAMDataUpdater.h"
#include "ParticleEmitter.h"
#include "Geometry.h"
#include "VectorField.h"

class YeeGrid3D {
    public:
    YeeGrid3D() { };
    //YeeGrid3D(YeeGrid3D& other) = default;
    YeeGrid3D(YeeGrid3D&& other) = default;
    YeeGrid3D& operator=(YeeGrid3D&& other) = default;
    ~YeeGrid3D();

    void SetCornerCoordinates(std::array<FPNumber, 3> r_0, std::array<FPNumber, 3> r_1);
    void SetNumOfCells(std::array<std::size_t, 3>& nCells);
    void SetTimeResolution(const FPNumber dt);
    void SetTimeIndex(const std::size_t ind);
    const std::array<FPNumber, 3>& GetCornerR0() const;
    const std::array<FPNumber, 3>& GetCornerR1() const;
    const std::array<std::size_t, 3>& GetNumberOfCells() const;
    FPNumber GetTimeResolution() const;
    FPNumber GetSpaceResolution(int i) const;
    const std::array<FPNumber, 3>& GetSpaceResolution() const;

    // a grid array that spans over the entire Yee grid
    void AddEntireGridElement(const std::string name, ElementType elemType);
    // a grid array that spans over a fraction of the Yee grid
    void AddPartialGridElement(const std::string name, ElementType elemType
        ,std::array<std::size_t, 3> startCell     // the element start on this cell index of the background grid
        ,std::array<std::size_t, 3> numCells      // number of cells covered by the element
        );
    YeeGridData3D& GetGridElement(const std::string name);

    void AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params);
    void SetIterativeSequence(std::vector<std::string> sequence);
    void SetSingleRunSequence(std::vector<std::string> sequence);
    void AddInstructionSequence(std::string name, std::vector<std::string> sequence);
    void ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params);
    void ApplyIterativeInstructions(std::size_t numIterations);
    void ApplyIterativeInstructionsOnce();
    void ApplySingleRunInstructions();
    void ApplyInstructionsOnce(std::string name);
    void ApplyInstructions(std::string name, std::size_t timeIndStart, std::size_t timeIndEnd, bool writeToFile = true);

    void* ConstructParams_A_plusequal_sum_b_C(
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<FPNumber> bValues,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                    ) const;
    void* ConstructParams_A_plusequal_sum_bB_C(
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<FPNumber> bValues,
                                    std::vector<std::string> arrayB_names,
                                    std::vector<int> arrayB_components,
                                    std::vector<std::array<std::size_t, 3>> arrayB_indsStart,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                    ) const;
    void* ConstructParams_A_equal_func_r_t(
                                    std::string gridManipulator_name
                                    ) const;

    void* ConstructParams_A_plusequal_sum_b_C_neighbor(
                                    YeeGrid3D* neighborGrid,
                                    std::array<std::size_t, 3> ind_start_A,
                                    std::array<std::size_t, 3> ind_end_A,
                                    std::string arrayA_name,
                                    int arrayA_component,
                                    std::vector<FPNumber> bValues,
                                    std::vector<std::string> arrayC_names,
                                    std::vector<int> arrayC_components,
                                    std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                    ) const;

    std::array<FPNumber, 3> GetCoordinatesOfFirstElementOfGridDataArray(const std::string& gridDataName, int direction);

    void AddGaussianGridArrayManipulator(const std::string name,
            const std::string gridDataName,     // name of the gridDataObject whose data is manipulated by this point source
            int direction, FPNumber amplitude,
            FPNumber t_center, FPNumber t_decay, FPNumber modulationFrequecy,
            FPNumber modulatioPhase, FPNumber timeOffsetFraction
            );
    void AddGaussianGridArrayManipulator(const std::string name,
            const std::string gridDataName,     // name of the gridDataObject whose data is manipulated by this point source
            int direction,
            std::array<std::size_t, 3> indStart,    // The manipulator only operates on a slice of the array starting from indStart
            std::array<std::size_t, 3> indEnd,      // and ending at indEnd-1
            FPNumber amplitude,
            FPNumber t_center, FPNumber t_decay, FPNumber modulationFrequecy,
            FPNumber modulatioPhase, FPNumber timeOffsetFraction
            );
    void AddSpatialCubeGridArrayManipulator(const std::string name,
            const std::string gridDataName,     // name of the gridDataObject whose data is manipulated by this point source
            int direction, // 0:x-component 1:y-component 2-z-component
            std::array<FPNumber, 3> boxCornerR0, std::array<FPNumber, 3> boxCornerR1, // corners of the cube
            std::array<FPNumber, 3> edgeThickness,    // thickness of the smooth edge
            FPNumber insideValue, FPNumber outsideValue // inside the cube set the array value to insideValue and...
            );
    void AddSpaceTimeCubeGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            std::array<FPNumber, 4> boxCornerR0, std::array<FPNumber, 4> boxCornerR1,   // space-time corners of the hypercube, index [3]->time
            std::array<FPNumber, 4> edgeThickness,      // space-time edge thickness. [0..2] -> space, [3] -> time
            FPNumber insideValue, FPNumber outsideValue,
            FPNumber timeOffsetFraction
            );
    void AddGaussianSpaceTimeGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            FPNumber amplitude,
            std::array<FPNumber, 4> st_center,
            std::array<FPNumber, 4> st_decay_rate,
            std::array<FPNumber, 4> st_modulationFrequecy,
            std::array<FPNumber, 4> st_modulatioPhase,
            FPNumber timeOffsetFraction
            );
    void AddPeriodicGaussianGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            FPNumber amplitude,
            std::array<FPNumber, 3> center,
            std::array<FPNumber, 3> decay_rate,
            std::array<FPNumber, 3> unitCellOrigin,
            std::array<std::array<FPNumber, 3>, 3> primitiveVectors
            );
    void AddSpherialShellGaussianGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            FPNumber amplitude,
            std::array<FPNumber, 3> centerPoint,
            FPNumber radius,
            FPNumber r_decay_rate,
            FPNumber r_modulationFrequecy,
            FPNumber r_modulatioPhase
            );
    void AddBivalueGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            const std::string geometryName,
            FPNumber valueInside,
            FPNumber valueOutside
            );
    void AddGaussianPlaneWaveGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            FPNumber amplitude,
            FPNumber t_center,
            FPNumber t_decay_rate,
            FPNumber t_modulationFrequecy,
            FPNumber t_modulatioPhase,
            FPNumber timeOffsetFraction
            );
    void AddRectPlaneWaveGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            FPNumber amplitude,
            FPNumber t_center,
            FPNumber t_rect_width,
            FPNumber t_edge_width,
            FPNumber t_modulationFrequecy,
            FPNumber t_modulatioPhase,
            FPNumber timeOffsetFraction
            );
    void AddDataTruncationGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            FPNumber minValue,
            FPNumber maxValue,
            bool truncateMin,
            bool truncateMax
            );

    void AddWedgeGeometry(const std::string name,
            const FPNumber wedgeAngle,
            const FPNumber tipRadius,
            const FPNumber apexToBaseDistance,
            const std::array<FPNumber, 3> apexPosition,
            const bool closeBase = true
            );

    void AddGaussianPlaneWaveVectorField(const std::string name,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            std::array<FPNumber, 3> amplitude,
            FPNumber t_center,
            FPNumber t_decayRate,
            FPNumber t_modulationFrequency,
            FPNumber t_modulationPhase
            );

    void AddRectPlaneWaveVectorField(const std::string name,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            std::array<FPNumber, 3> amplitude,
            FPNumber t_center,
            FPNumber t_rectWidth,
            FPNumber t_edgeWidth,
            FPNumber t_modulationFrequency,
            FPNumber t_modulationPhase
            );

    void AddManualChargedParticleEmitter(const std::string name,
            FPNumber particleCharge,
            FPNumber particleMasse,
            const std::vector<FPNumber>& emissionTimes,        // particles are emitted at these times
            const std::vector<FPNumber>& emissionNumbers,      // number of particles emitted at emissionTimes[i]
            const std::vector<std::array<FPNumber, 3>>& particlesInitialPositions,
            const std::vector<std::array<FPNumber, 3>>& particlesInitialVelocities
            );

    void AddChargedParticleEmitter(const std::string name,
            FPNumber particleCharge,
            FPNumber particleMasse,
            const std::string geometryName,
            int dimensions,
            FPNumber maxElemSize,
            const std::string eFieldName,
            const std::string analyticEFieldName,
            FPNumber unitLength,
            std::size_t numOfSubPoints = 0
            );

    void AddChargedParticlesTracer(const std::string name,
            const std::string eFieldName,
            const std::string bFieldName,
            const std::string analyticEFieldName,
            const std::string analyticBFieldName,
            const std::string srFieldName,
            const std::string particleEmitterName,
            const std::size_t numberOfReservedParticles,
            const std::size_t bunchSize,
            const std::string constrainingGeometryName,
            bool keepPointsInside
            );

    void AddDiscretePointsGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            const std::string dataUpdaterName,
            const std::string dataUpdaterDataName,      // name of the array inside the dataUpdater to associate with gridData
            int dataUpdaterDataDirection                // direction of the array inside dataUpdater
    );

    void PrintAllGridData();

    void AddGridElementView(std::string gridElemViewName,   // name of the gridView
                            std::string gridElemName , int gridElemComponent,   // name of the gridElement and its x,y,z component
                            std::array<std::size_t, 3> indStart, std::array<std::size_t, 3> indEnd // slice start and end
                            );
    void AddFullGridElementView(std::string gridElemViewName,   // name of the gridView
                            std::string gridElemName , int gridElemComponent   // name of the gridElement and its x,y,z component
                            );

    void WriteAllGridElemViewsToFile();
    void DeleteOlderViewFiles();
    void SetDataStoreRate(std::string gridElemViewName, std::size_t saveEveryNSammples);
    void CloseGridViewFiles();

    protected:
    std::size_t timeIndex = 0;
    FPNumber dt;                  // time resolution
    std::array<FPNumber, 3> dr;   // spatial resolution
    std::array<FPNumber, 3> r_0{0, 0, 0};  // coordinates of the lower left corner
    std::array<FPNumber, 3> r_1{0, 0, 0};  // coordinates of the upper right corner
    std::array<std::size_t, 3> nCells;      // number of Yee cells
    std::unordered_map<std::string, std::shared_ptr<YeeGridData3D>> gridElements;   // TODO: pointer is unnecessary here
    std::unordered_map<std::string, std::shared_ptr<GridArrayManipulator>> gridArrayManipulators;   // using pointer since GridArrayManipulator is virtual and cannot be instantiated
    std::unordered_map<std::string, std::pair<FDInstructionCode, void*>> updateInstructions;  // field update instructions
    std::vector<std::string> iterativeSequence;     // sequence in which to apply the field update instructions
    std::vector<std::string> singleRunSequence;     // sequence in which to apply the field update instructions
    std::unordered_map<std::string, std::vector<std::string>> instructionSequences;
    std::unordered_map<std::string, GridElementView> gridElementViews;  // a slice of gridElements for printing to output
    std::unordered_map<std::string, std::shared_ptr<DiscretePointsGAMDataUpdater>> gamDataUpdaters;  // data updaters for (some) grid array manipulators
    std::unordered_map<std::string, std::shared_ptr<Geometry>> geometries;
    std::unordered_map<std::string, std::shared_ptr<ParticleEmitter>> particleEmitters;
    std::unordered_map<std::string, std::shared_ptr<VectorField>> vectorFields;
};


#endif  // FDTD_YEEGRID_H_


