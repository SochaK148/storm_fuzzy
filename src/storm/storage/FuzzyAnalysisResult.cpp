#include <iostream>
#include <vector>

#include "storm/storage/TriangularFuzzyNumber.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/FuzzyAnalysisResult.h"

namespace storm {
namespace storage {

    storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> FuzzyAnalysisResult::getMatrix() const {
        return this->matrix;
    }

    Interval storm::storage::FuzzyAnalysisResult::getAlphaCut(storm::storage::TriangularFuzzyNumber const& t, double alpha){
        double i = alpha * (t.getPeak() - t.getLeftBound()) + t.getLeftBound();
        double j = alpha * (t.getRightBound() - t.getPeak()) + t.getPeak();
        return Interval(i, j);
    }

    bool storm::storage::FuzzyAnalysisResult::isFeasible(std::vector<std::vector<double>> const& m, double alpha){
        for (int row = 0; row < this->getMatrix().getRowCount(); ++row) {
            std::cout << "Row: " << row << std::endl;
            // std::cout << this->getMatrix().getRow(row) << std::endl;
            for (storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber>::const_iterator it = this->getMatrix().begin(row), ite = this->getMatrix().end(row); it != ite; ++it) {
                int column = it->getColumn();
                storm::storage::TriangularFuzzyNumber value = it->getValue();
                std::cout << "Column: " << column << ", Value: " << value << std::endl;
                Interval alphaCut = this->getAlphaCut(value, alpha);
                double crispValue = m[row][column];
                // TODO check if crispValue is within alphaCut
            }
        }
        return true; // placeholder
    }
    
}  // namespace storage
}  // namespace storm