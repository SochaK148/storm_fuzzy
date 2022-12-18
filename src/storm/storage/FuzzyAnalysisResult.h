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
            std::vector<std::vector<storm::storage::TriangularFuzzyNumber>> matrixV;
            std::vector<std::vector<Interval>> matrixI;
        public:
            class member{
                private:
                    std::vector<std::vector<double>> matrix;
                    int n;
                    std::pair<int, int> idx;
                    double fitness;
                public:
                    member(std::vector<std::vector<double>> matrix, int n, std::pair<int, int> idx) : matrix(matrix), n(n), idx(idx) {}
                    std::vector<std::vector<double>> getMatrix();
                    double getFitness();
                    void updateFitness();
                    std::vector<std::vector<double>> matrixMul(std::vector<std::vector<double>> matrix, int n);
                    bool operator< (const member &other) const;
            };
            FuzzyAnalysisResult(storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> m) : matrix(m) {}
            storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> getMatrix() const;
            std::vector<std::vector<Interval>> getMatrixI() const;
            std::vector<std::vector<double>> randomCrisp();
            void calcMatrixI(double alpha);
            Interval getAlphaCut(storm::storage::TriangularFuzzyNumber const& t, double alpha);
            double restrictedMatrixMul(std::vector<std::vector<Interval>> intvlP, int steps, std::pair<int, int> idx, bool isMin);
            std::vector<member> initializePopulation(int populationSize, int steps, std::pair<int, int> idx);
            std::vector<FuzzyAnalysisResult::member> selectPopulation(std::vector<FuzzyAnalysisResult::member> population, int selectionSample, bool isMin);
            std::vector<FuzzyAnalysisResult::member> crossPopulation(std::vector<FuzzyAnalysisResult::member> population);
            std::vector<FuzzyAnalysisResult::member> mutatePopulation(std::vector<FuzzyAnalysisResult::member> population, int mutationRate);

            bool isFeasible(std::vector<std::vector<double>> const& m, double alpha);
    };
    
}  // namespace storage
}  // namespace storm

#endif /* STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_ */