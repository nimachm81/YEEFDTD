
#ifndef FDTD_GRIDELEMENTVIEW_H_
#define FDTD_GRIDELEMENTVIEW_H_

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"

class GridElementView {
    public:
    GridElementView();
    //~GridElementView();
    void SetSaveOnDiskFrequency(const std::size_t eachNSamples);
    std::size_t GetSaveOnDiskFrequency() const;
    void SetName(std::string name);
    void SetNumArray(const NumberArray3D<FPNumber>& numArray);
    void OpenFileToWrite();
    void CloseFile();
    void StoreData(std::size_t iterationIndex);
    void DeleteOlderFiles();

    private:
    std::size_t saveOnDiskFrequency = 1;
    NumberArray3D<FPNumber> numArray;   // it is initialized as a slice
    std::string viewName;
    std::string viewFolder;
    std::string viewFileExtension;

    std::ofstream file;
    bool fileIsOpen = false;
    std::size_t bufferSize;
    //std::unique_ptr<char> buffer = nullptr;
    char* buffer = nullptr;
    std::size_t bufferInd = 0;  // points to the next empty index
    std::size_t maxArraySizeInBytes = 0;
};



#endif // FDTD_GRIDELEMENTVIEW_H_
