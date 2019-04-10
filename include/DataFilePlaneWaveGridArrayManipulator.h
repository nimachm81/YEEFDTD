
#ifndef FDTD_DATAFILELANEWAVEGRIDARRAYMANIPULATOR_H_
#define FDTD_DATAFILELANEWAVEGRIDARRAYMANIPULATOR_H_

#include "DataSampleFileReader.h"
#include "PlaneWaveGridArrayManipulator.h"

class DataFilePlaneWaveGridArrayManipulator : public DataSampleFileReader, public PlaneWaveGridArrayManipulator {
    public:
    virtual ~DataFilePlaneWaveGridArrayManipulator() { };

    virtual FPNumber GetNormalizedTemporalWaveform(FPNumber t);
};

#endif // DataFilePlaneWaveGridArrayManipulator



