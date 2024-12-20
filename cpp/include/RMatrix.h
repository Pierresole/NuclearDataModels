#ifndef RMATRIX_H
#define RMATRIX_H

#include "SpinGroup.h"
#include <Eigen/Dense>

class RMatrix {
public:
    RMatrix(const SpinGroup& spinGroup, double E);

    const Eigen::MatrixXd& getMatrix() const { return R_; }

private:
    Eigen::MatrixXd R_;
};

#endif // RMATRIX_H
