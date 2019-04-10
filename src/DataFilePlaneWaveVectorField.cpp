
#include <fstream>

#include "DataFilePlaneWaveVectorField.h"

FPNumber DataFilePlaneWaveVectorField::GetNormalizedWaveform(const FPNumber& t) {
    return DataSampleFileReader::GetWaveform(t);
}

