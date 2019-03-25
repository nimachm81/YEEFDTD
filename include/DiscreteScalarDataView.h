
#ifndef FDTD_DISCRETESCALARDATAVIEW_H_
#define FDTD_DISCRETESCALARDATAVIEW_H_

#include <string>
#include <vector>
#include <array>

#include "NumberTypes.h"
#include "DataView.h"

class DiscreteScalarDataView : public DataView {
    public:
    DiscreteScalarDataView();

    void SetScalarProperty(std::vector<FPNumber>* scalarProp);

    virtual std::size_t GetMaxDataSizeInBytes();
    virtual void WriteDataToBuffer();
    virtual void WriteDataToFile();

    private:
    std::vector<FPNumber>* scalarProperty;

};

#endif // FDTD_DISCRETESCALARDATAVIEW_H_

