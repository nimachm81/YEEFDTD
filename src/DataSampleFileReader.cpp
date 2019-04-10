
#include <fstream>

#include "UtilityFunctions.hpp"
#include "DataSampleFileReader.h"

void DataSampleFileReader::SetTimeSamplesFromFile(std::string fileName) {
    UtilityFunctions::Read1DArrayFromFile<FPNumber>(fileName, timeSamples);
    for(std::size_t i = 0; i < timeSamples.size() - 1; ++i) {
        assert(timeSamples[i] < timeSamples[i + 1]);
        //std::cout << timeSamples[i] << " ";
    }
    //std::cout << std::endl;
}

void DataSampleFileReader::SetFieldSamplesFromFile(std::string fileName) {
    UtilityFunctions::Read1DArrayFromFile<FPNumber>(fileName, fieldSamples);
}

std::size_t DataSampleFileReader::GetLowerTimeIndex(FPNumber t) {
    for(std::size_t i = 0; i < timeSamples.size(); ++i) {
        if(timeSamples[i] > t) {
            return i - 1;
        }
    }
    assert(false);
}

FPNumber DataSampleFileReader::GetWaveform(const FPNumber& t) {
    std::size_t n_samples = timeSamples.size();
    FPNumber val = 0.0;
    if(t <= timeSamples[0]) {
        val = fieldSamples[0];
    } else if(t >= timeSamples[n_samples - 1]) {
        val =  fieldSamples[n_samples - 1];
    } else {
        std::size_t ind = GetLowerTimeIndex(t);
        FPNumber alpha = (t - timeSamples[ind]) / (timeSamples[ind + 1] - timeSamples[ind]);
        val = (1.0 - alpha)*fieldSamples[ind] + alpha*fieldSamples[ind + 1];
    }

    return val;
}

