
#ifndef FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_
#define FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_

class DiscretePointsGAMDataUpdater {
    public:
    virtual ~DiscretePointsGAMDataUpdater() { };

    virtual void PointToDataPositions(std::vector<std::array<FPNumber, 3>>*& positions) = 0;

    virtual void PointToScalarData(std::vector<FPNumber>*& values,
                                   std::string dataName,
                                   int direction) = 0;

    virtual void PointToVectorData(std::vector<std::array<FPNumber, 3>>*& values,
                                 std::string dataNamen) = 0;

    virtual void UpdateGAMValues(const FPNumber t) = 0;
};

#endif // FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_
