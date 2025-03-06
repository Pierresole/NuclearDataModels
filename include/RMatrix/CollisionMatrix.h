#ifndef COLLISION_MATRIX_H
#define COLLISION_MATRIX_H

#include "XMatrix.h"
#include "SpinGroup.h"
#include <Eigen/Dense>

class CollisionMatrix {
public:
    CollisionMatrix(const XMatrix& X, const SpinGroup& spinGroup);

    const Eigen::MatrixXcd& getMatrix() const { return U_; }

private:
    Eigen::MatrixXcd U_;
};

#endif // COLLISION_MATRIX_H