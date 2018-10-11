
#include <cmath>

#include "GridElementView.h"

GridElementView::GridElementView() {
    viewFolder = "data";
    viewFileExtension = ".data";
    saveOnDiskFrequency = 1;
}

void GridElementView::SetSaveOnDiskFrequency(const std::size_t eachNSamples) {
    saveOnDiskFrequency = eachNSamples;
}

std::size_t GridElementView::GetSaveOnDiskFrequency() const {
    return saveOnDiskFrequency;
}

void GridElementView::SetName(std::string name) {
    viewName = name;
}

void GridElementView::SetNumArray(const NumberArray3D<RealNumber>& numArray) {
    GridElementView::numArray = std::make_unique<NumberArray3D<RealNumber>>(numArray);
}

void GridElementView::OpenFileToWrite() {
    assert(!viewName.empty());
    std::string fileName = viewFolder + "/" + viewName + viewFileExtension;
    std::cout << "output file: " << fileName << std::endl;
    file.open(fileName.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    assert(file.is_open());
    fileIsOpen = true;
}

void GridElementView::CloseFile() {
    if(fileIsOpen) {
        file.close();
        fileIsOpen = false;
    }
}

void GridElementView::StoreData(std::size_t iterationIndex) {
    if(!fileIsOpen) {
        OpenFileToWrite();
    }
    if(std::remainder(iterationIndex, saveOnDiskFrequency) == 0) {
        numArray->WriteArrayDataToFile(&file, true, true);
    }
}

void GridElementView::DeleteOlderFiles() {
    assert(!viewName.empty());
    std::string fileName = viewFolder + "/" + viewName + viewFileExtension;
    std::cout << "Older data file: " << fileName << " will be overwritten." << std::endl;
    std::ifstream ifile(fileName.c_str());
    if(ifile) {
        int file_deleted = std::remove(fileName.c_str());
        assert(file_deleted == 0);
        std::cout << fileName << " deleted!" << std::endl;
    }
}

