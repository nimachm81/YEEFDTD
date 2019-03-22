
#include <cstddef>

#include "DiscretePointsVectorGridArrayManipulator.h"


void DiscretePointsVectorGridArrayManipulator::AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                                                        std::string dataName,
                                                        int direction
                                                        ) {
    dataUpdater = updater;
    updater->AttachDataToGAMPositions(positions);
    updater->AttachVectorDataToGAMValues(values, dataName, direction);

    assert(direction >= 0 && direction < 3);
    DiscretePointsVectorGridArrayManipulator::direction = direction;
}

FPNumber DiscretePointsVectorGridArrayManipulator::GetDataValueByInddex(std::size_t ind) {
    return (*values)[ind][direction];
}

