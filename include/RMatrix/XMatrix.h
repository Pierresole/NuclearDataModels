#ifndef XMATRIX_H
#define XMATRIX_H

#include "LevelMatrix.h"
#include "SpinGroup.h"
#include <Eigen/Dense>

class XMatrix {
public:
    XMatrix(const LevelMatrix& A, const SpinGroup& spinGroup, double E);

    const Eigen::MatrixXcd& getMatrix() const { return X_; }

private:
    Eigen::MatrixXcd X_;
};

#endif // XMATRIX_H
