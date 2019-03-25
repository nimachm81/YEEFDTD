
#include <cmath>

#include "GridElementView.h"

GridElementView::GridElementView() {

    maxArraySizeInBytes = NumberArray3D<FPNumber>().GetMaxDataSizeInBytes();
}

void GridElementView::SetNumArray(const NumberArray3D<FPNumber>& numArray) {
    GridElementView::numArray.MakeThisASliceOf(numArray);
    maxArraySizeInBytes = GridElementView::numArray.GetMaxDataSizeInBytes();
}

std::size_t GridElementView::GetMaxDataSizeInBytes() {
    maxArraySizeInBytes;
}

void GridElementView::WriteDataToBuffer() {
    numArray.WriteArrayDataToMemory(buffer,
                                    bufferInd,
                                    true,    // writeShape
                                    true     // writeDataTypeSize
                                    );
}

void GridElementView::WriteDataToFile() {
    numArray.WriteArrayDataToFile(&file,
                                   true,    // writeShape
                                   true     // writeDataTypeSize
                                   );
}

