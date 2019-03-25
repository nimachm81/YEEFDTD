
#ifndef FDTD_DISCRETEVECTORDATAVIEW_H_
#define FDTD_DISCRETEVECTORDATAVIEW_H_

#include <string>
#include <vector>
#include <array>

#include "NumberTypes.h"
#include "DataView.h"

class DiscreteVectorDataView : public DataView {
    public:
    DiscreteVectorDataView();

    void SetVectorProperty(std::vector<std::array<FPNumber, 3>>* vecProp);

    virtual std::size_t GetMaxDataSizeInBytes();
    virtual void WriteDataToBuffer();
    virtual void WriteDataToFile();

    private:
    std::vector<std::array<FPNumber, 3>>* vectorProperty;

};

#endif // FDTD_DISCRETEVECTORDATAVIEW_H_


