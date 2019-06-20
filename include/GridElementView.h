
#ifndef FDTD_GRIDELEMENTVIEW_H_
#define FDTD_GRIDELEMENTVIEW_H_

#include "NumberTypes.h"
#include "MultiDimArray.hpp"
#include "DataView.h"

class GridElementView : public DataView {
    public:
    GridElementView(std::size_t bufferSize = 1024*1024*20);
    //~GridElementView();
    void SetNumArray(const NumberArray3D<FPNumber>& numArray);

    virtual std::size_t GetMaxDataSizeInBytes();
    virtual void WriteDataToBuffer();
    virtual void WriteDataToFile();

    private:
    NumberArray3D<FPNumber> numArray;   // it is initialized as a slice
    std::size_t maxArraySizeInBytes = 0;
};

#endif // FDTD_GRIDELEMENTVIEW_H_
