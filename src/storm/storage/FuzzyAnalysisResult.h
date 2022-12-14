#ifndef STORM_STORAGE_FUZZYANALYSISRESULT_H_
#define STORM_STORAGE_FUZZYANALYSISRESULT_H_

#include "storm/storage/TriangularFuzzyNumber.h"
#include "storm/storage/SparseMatrix.h"

#include <iostream>
#include <vector>

namespace storm {
namespace storage {

    class FuzzyAnalysisResult{
        private:
            storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> matrix;
            
        public:
            FuzzyAnalysisResult(storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> m) : matrix(m) {}
            storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> getMatrix() const;
            Interval getAlphaCut(storm::storage::TriangularFuzzyNumber const& t, double alpha);
            bool isFeasible(std::vector<std::vector<double>> const& m, double alpha);
    };
    
}  // namespace storage
}  // namespace storm

#endif /* STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_ */