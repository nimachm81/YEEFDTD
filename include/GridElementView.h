
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
    void SetSaveOnDiskFrequency(const std::size_t eachNSamples);
    std::size_t GetSaveOnDiskFrequency() const;
    void SetName(std::string name);
    void SetNumArray(const NumberArray3D<RealNumber>& numArray);
    void OpenFileToWrite();
    void CloseFile();
    void StoreData(std::size_t iterationIndex);
    void DeleteOlderFiles();

    private:
    std::size_t saveOnDiskFrequency;
    std::unique_ptr<NumberArray3D<RealNumber>> numArray;
    std::string viewName;
    std::string viewFolder;
    std::string viewFileExtension;

    std::ofstream file;
    bool fileIsOpen = false;
};



#endif // FDTD_GRIDELEMENTVIEW_H_
