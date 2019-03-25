
#ifndef FDTD_DATAVIEW_H_
#define FDTD_DATAVIEW_H_

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "NumberTypes.h"

class DataView {
    public:
    DataView();
    //~DataView();
    void SetSaveOnDiskFrequency(const std::size_t eachNSamples);
    std::size_t GetSaveOnDiskFrequency() const;
    void SetName(std::string name);
    void OpenFileToWrite();
    void CloseFile();
    void DeleteOlderFiles();

    void StoreData(std::size_t iterationIndex);

    virtual std::size_t GetMaxDataSizeInBytes() = 0;
    virtual void WriteDataToBuffer() = 0;
    virtual void WriteDataToFile() = 0;

    protected:
    std::size_t saveOnDiskFrequency = 1;
    std::string viewName;
    std::string viewFolder;
    std::string viewFileExtension;

    std::ofstream file;
    bool fileIsOpen = false;
    std::size_t bufferSize;

    //std::unique_ptr<char> buffer = nullptr;
    char* buffer = nullptr;
    std::size_t bufferInd = 0;  // points to the next empty index
};

#endif // FDTD_DATAVIEW_H_

