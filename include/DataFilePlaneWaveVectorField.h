
#ifndef FDTD_DATAFILEPLANEWAVEVECTORFIELD_H_
#define FDTD_DATAFILEPLANEWAVEVECTORFIELD_H_

#include "DataSampleFileReader.h"
#include "PlaneWaveVectorField.h"

class DataFilePlaneWaveVectorField : public DataSampleFileReader, public PlaneWaveVectorField {
    public:
    virtual ~DataFilePlaneWaveVectorField() { };

    virtual FPNumber GetNormalizedWaveform(const FPNumber& t);
};

#endif // FDTD_DATAFILEPLANEWAVEVECTORFIELD_H_


