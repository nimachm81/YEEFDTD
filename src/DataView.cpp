
#include <cmath>
#include <cassert>

#include "DataView.h"

DataView::DataView() {
    fileIsOpen = false;
    viewFolder = "data";
    viewFileExtension = ".data";
    saveOnDiskFrequency = 1;

    bufferSize = 1024*1024*200;

    //buffer = std::make_unique<char>(bufferSize);
    //buffer = new char[bufferSize];
}

//DataView::~DataView() {
//    CloseFile();
//    if(buffer != nullptr) {
//        delete[] buffer;
//        buffer = nullptr;
//    }
//}

void DataView::SetSaveOnDiskFrequency(const std::size_t eachNSamples) {
    saveOnDiskFrequency = eachNSamples;
}

std::size_t DataView::GetSaveOnDiskFrequency() const {
    return saveOnDiskFrequency;
}

void DataView::SetName(std::string name) {
    viewName = name;
}

void DataView::OpenFileToWrite() {
    if(!fileIsOpen) {
        assert(!viewName.empty());
        std::string fileName = viewFolder + "/" + viewName + viewFileExtension;
        std::cout << "output file: " << fileName << std::endl;
        file.open(fileName.c_str(), std::ios::out | std::ios::app | std::ios::binary);
        assert(file.is_open());
        fileIsOpen = true;

        if(buffer == nullptr) {
            buffer = new char[bufferSize];
        } else {
            assert(false);
        }
    } else {
        assert(false);
    }

    // TODO: control file synchronization and buffering
    //file.rdbuf()->pubsetbuf(buffer, bufferSize);
    //file.unsetf(std::ios_base::unitbuf);
}

void DataView::CloseFile() {
    if(fileIsOpen) {
        if(bufferInd > 0) {
            file.write(buffer, bufferInd);
            bufferInd = 0;
        }

        file.close();
        fileIsOpen = false;

        if(buffer != nullptr) {
            delete[] buffer;
        } else {
            assert(false);
        }
    }
}


void DataView::DeleteOlderFiles() {
    assert(!viewName.empty());
    std::string fileName = viewFolder + "/" + viewName + viewFileExtension;
    //std::cout << "Older data file: " << fileName << " will be overwritten." << std::endl;
    std::ifstream ifile(fileName.c_str());
    if(ifile) {
        int file_deleted = std::remove(fileName.c_str());
        assert(file_deleted == 0);
        std::cout << fileName << " already exists. It was Overwritten!" << std::endl;
    }
}


void DataView::StoreData(std::size_t iterationIndex) {
    if(!fileIsOpen) {
        OpenFileToWrite();
    }
    if(std::remainder(iterationIndex, saveOnDiskFrequency) == 0) {
        if(bufferInd + GetMaxDataSizeInBytes() < bufferSize) {

            WriteDataToBuffer();

        } else {
            if(bufferInd > 0) {
                file.write(buffer, bufferInd);
                bufferInd = 0;
            }

            WriteDataToFile();
        }
        if(verbosityLevel >= 1 && std::remainder(iterationIndex, 100*saveOnDiskFrequency) == 0) {
            std::cout << "Writing data at iterationIndex: " << iterationIndex << std::endl;
        }
    }
}
