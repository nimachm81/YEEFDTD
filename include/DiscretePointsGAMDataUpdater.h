
#ifndef FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_
#define FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_

class DiscretePointsGAMDataUpdater {
    public:
    virtual ~DiscretePointsGAMDataUpdater() { };

    virtual void AttachDataToGAMPositions(std::vector<std::array<FPNumber, 3>>*& positions) = 0;

    virtual void AttachScalarDataToGAMValues(std::vector<FPNumber>*& values,
                                             std::string dataName,
                                             int direction) = 0;

    virtual void AttachVectorDataToGAMValues(std::vector<std::array<FPNumber, 3>>*& values,
                                             std::string dataName,
                                             int direction) = 0;

    virtual void UpdateGAMValues(const FPNumber t) = 0;
};

#endif // FDTD_DISCRETEPOINTSGAMDATAUPDATER_H_
