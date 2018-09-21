
#include <cassert>
#include <cmath>

#include "UniformTimeVaryingVector.h"

void UniformSingleDirectionGaussianTimeVaryingVector::SetDirection(const int direction) {
    assert(direction>=0 && direction<3);
    this->direction = direction;
}

void UniformSingleDirectionGaussianTimeVaryingVector::SetAmplitude(const RealNumber amplitude) {
    amplitude = amplitude;
}

void UniformSingleDirectionGaussianTimeVaryingVector::SetCenterTime(const RealNumber t_center) {
    t_center = t_center;
}

void UniformSingleDirectionGaussianTimeVaryingVector::SetDecayTime(const RealNumber t_decay) {
    t_decay = t_decay;
}


virtual void UniformSingleDirectionGaussianTimeVaryingVector::UpdateArray(YeeGridData3D& gridData, RealNumber t) {
    NumberArray3D<RealNumber>& arrayA = gridData.GetNumArray(direction);
    RealNumber gaussianValue = amplitude * std::exp(-(t - t_center)*(t - t_center) / (t_decay*t_decay));
    arrayA.SetToNumber(gaussianValue);
}

