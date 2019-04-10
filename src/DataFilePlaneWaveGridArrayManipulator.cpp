
#include "DataFilePlaneWaveGridArrayManipulator.h"

FPNumber DataFilePlaneWaveGridArrayManipulator::GetNormalizedTemporalWaveform(FPNumber t) {
    return DataSampleFileReader::GetWaveform(t);
}


