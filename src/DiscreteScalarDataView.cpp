
#include <cassert>

#include "UtilityFunctions.hpp"
#include "DiscreteScalarDataView.h"

DiscreteScalarDataView::DiscreteScalarDataView() {
    viewFileExtension = ".dsdata";
}

void DiscreteScalarDataView::SetScalarProperty(std::vector<FPNumber>* scalarProp) {
    scalarProperty = scalarProp;
}

std::size_t DiscreteScalarDataView::GetMaxDataSizeInBytes() {
    return scalarProperty->size() * sizeof(FPNumber) +
           sizeof(int) +                                    // datatypeCode
           sizeof(std::size_t) +                            // dataSize
           sizeof(std::size_t);                             // numOfPoints
}

void DiscreteScalarDataView::WriteDataToBuffer() {
    int dataypeCode = UtilityFunctions::GetDatatypeNumericalCode<FPNumber>();
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&dataypeCode, sizeof(int));

    std::size_t dataSize = sizeof(FPNumber);
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&dataSize, sizeof(std::size_t));

    std::size_t numPoints = scalarProperty->size();
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&numPoints, sizeof(std::size_t));

    std::vector<FPNumber>& scalarPropertyRef = *scalarProperty;
    for(std::size_t i = 0; i < numPoints; ++i) {
        UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)(&scalarPropertyRef[i]), sizeof(FPNumber));
    }
}

void DiscreteScalarDataView::WriteDataToFile() {
    int dataypeCode = UtilityFunctions::GetDatatypeNumericalCode<FPNumber>();
    file.write((char*)&dataypeCode, sizeof(int));

    std::size_t dataSize = sizeof(FPNumber);
    file.write((char*)&dataSize, sizeof(std::size_t));

    std::size_t numPoints = scalarProperty->size();
    file.write((char*)&numPoints, sizeof(std::size_t));

    std::vector<FPNumber>& scalarPropertyRef = *scalarProperty;
    for(std::size_t i = 0; i < numPoints; ++i) {
        file.write((char*)(&scalarPropertyRef[i]), sizeof(FPNumber));
    }
}
