#include "LevelMatrix.h"

LevelMatrix::LevelMatrix(const Eigen::MatrixXcd& R, const Eigen::MatrixXcd& L) {
    size_t size = R.rows();
    Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(size, size);
    
    // Solve for A using LU decomposition instead of direct inversion
    // Eigen::PartialPivLU<Eigen::MatrixXcd> lu(I - R * L);
    // A_ = lu.solve(I);
}
