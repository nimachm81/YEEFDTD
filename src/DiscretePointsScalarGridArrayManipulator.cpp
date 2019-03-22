
#include <cstddef>

#include "DiscretePointsScalarGridArrayManipulator.h"


void DiscretePointsScalarGridArrayManipulator::AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                                                        std::string dataName,
                                                        int direction
                                                        ) {
    dataUpdater = updater;
    updater->AttachDataToGAMPositions(positions);
    updater->AttachScalarDataToGAMValues(values, dataName, direction);

    assert(direction >= 0 && direction < 3);
}

FPNumber DiscretePointsScalarGridArrayManipulator::GetDataValueByInddex(std::size_t ind) {
    return (*values)[ind];
}
