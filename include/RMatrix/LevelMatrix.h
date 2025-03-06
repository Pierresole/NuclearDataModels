#ifndef LEVEL_MATRIX_H
#define LEVEL_MATRIX_H

#include <Eigen/Dense>

class LevelMatrix {
public:
    LevelMatrix(const Eigen::MatrixXcd& R, const Eigen::MatrixXcd& L);

    const Eigen::MatrixXcd& getMatrix() const { return A_; }

private:
    Eigen::MatrixXcd A_;
};

#endif // LEVEL_MATRIX_H
