
#ifndef FDTD_DATASAMPLEFILEREADER_H_
#define FDTD_DATASAMPLEFILEREADER_H_

#include <vector>

#include "NumberTypes.h"

class DataSampleFileReader {
    public:

    void SetTimeSamplesFromFile(std::string fileName);
    void SetFieldSamplesFromFile(std::string fileName);
    std::size_t GetLowerTimeIndex(FPNumber t);  // get ind such that timeSamples[ind] <= t < timeSamples[ind + 1]

    FPNumber GetWaveform(const FPNumber& t);

    protected:
    std::vector<FPNumber> timeSamples;
    std::vector<FPNumber> fieldSamples;
};

#endif // FDTD_DATASAMPLEFILEREADER_H_


