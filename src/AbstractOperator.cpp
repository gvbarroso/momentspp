// ==== AbstractOperator.cpp ====
#include "AbstractOperator.hpp"

#include <fstream>

void AbstractOperator::printDeltaLDMat(const std::string& fileName)
{
    // 1) Open file
    std::ofstream matFile(fileName);
    if (!matFile.is_open())
        throw bpp::Exception("AbstractOperator::failed to open file: " + fileName);

    // 2) Ensure we have a transition matrix
    if (!transition_)
        throw bpp::Exception("AbstractOperator::attempted to print un-initialized transition matrix!");

    // 3) Extract an Eigen‐variant view of the transition
    auto eigenVar = transition_->toEigenMatrixVariant();

    // 4) Visit whichever underlying Eigen type it is
    std::visit(overloaded{
        [&](const auto& mat) {
            const int rows = mat.rows();
            const int cols = mat.cols();
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    matFile << mat.coeff(i, j);
                    if (j + 1 < cols) matFile << ",";
                }
                matFile << "\n";
            }
        }
    }, eigenVar);

    matFile.close();
}
