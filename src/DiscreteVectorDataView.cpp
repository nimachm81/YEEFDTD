
#include <cassert>

#include "UtilityFunctions.hpp"
#include "DiscreteVectorDataView.h"

DiscreteVectorDataView::DiscreteVectorDataView() {
    viewFileExtension = ".dvdata";
}

void DiscreteVectorDataView::SetVectorProperty(std::vector<std::array<FPNumber, 3>>* vecProp) {
    vectorProperty = vecProp;
}

std::size_t DiscreteVectorDataView::GetMaxDataSizeInBytes() {
    return vectorProperty->size() * sizeof(FPNumber) * 3 +
           sizeof(int) +                                    // datatypeCode
           sizeof(std::size_t) +                            // dataSize
           sizeof(std::size_t);                             // numOfPoints
}

void DiscreteVectorDataView::WriteDataToBuffer() {
    int dataypeCode = UtilityFunctions::GetDatatypeNumericalCode<FPNumber>();
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&dataypeCode, sizeof(int));

    std::size_t dataSize = sizeof(FPNumber);
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&dataSize, sizeof(std::size_t));

    std::size_t numPoints = vectorProperty->size();
    UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)&numPoints, sizeof(std::size_t));

    std::vector<std::array<FPNumber, 3>>& vectorPropertyRef = *vectorProperty;
    for(std::size_t i = 0; i < numPoints; ++i) {
        for(int j = 0; j < 3; ++j) {
            UtilityFunctions::WriteToBuffer(buffer, bufferInd, (char*)(&vectorPropertyRef[i][j]), sizeof(FPNumber));
        }
    }
}

void DiscreteVectorDataView::WriteDataToFile() {
   int dataypeCode = UtilityFunctions::GetDatatypeNumericalCode<FPNumber>();
    file.write((char*)&dataypeCode, sizeof(int));

    std::size_t dataSize = sizeof(FPNumber);
    file.write((char*)&dataSize, sizeof(std::size_t));

    std::size_t numPoints = vectorProperty->size();
    file.write((char*)&numPoints, sizeof(std::size_t));

    std::vector<std::array<FPNumber, 3>>& vectorPropertyRef = *vectorProperty;
    for(std::size_t i = 0; i < numPoints; ++i) {
        for(int j = 0; j < 3; ++j) {
            file.write((char*)(&vectorPropertyRef[i][j]), sizeof(FPNumber));
        }
    }
}


