#ifndef STORM_STORAGE_FUZZYANALYSISRESULT_H_
#define STORM_STORAGE_FUZZYANALYSISRESULT_H_

#include "storm/storage/FlexibleSparseMatrix.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/TriangularFuzzyNumber.h"

#include <iostream>
#include <vector>

namespace storm {
namespace storage {

class FuzzyAnalysisResult {
   private:
    FlexibleSparseMatrix<TriangularFuzzyNumber> matrix;

   public:
    class member {
       private:
        std::vector<std::vector<double>> matrix;
        int n;
        std::pair<int, int> idx;
        double fitness;
        bool reachability;

       public:
        member(std::vector<std::vector<double>> matrix, int n, std::pair<int, int> idx, bool reachability = true) : matrix(matrix), n(n), idx(idx), reachability(reachability) {}
        std::vector<std::vector<double>> getMatrix();
        void setMatrix(std::vector<std::vector<double>> matrix);
        int getN();
        std::pair<int, int> getIdx();
        double getFitness();
        void updateFitness();
        bool getReachability();
        bool operator<(const member& other) const;
        bool operator==(member& other) const;
    };
    FuzzyAnalysisResult(SparseMatrix<TriangularFuzzyNumber> m) : matrix(FlexibleSparseMatrix<TriangularFuzzyNumber>(m)) {}
    FlexibleSparseMatrix<TriangularFuzzyNumber>& getMatrix();
    std::pair<std::vector<std::vector<Interval>>, std::vector<std::pair<int, int>>> getIntervalMatrix(double alpha);
    std::vector<std::vector<double>> getCrispMatrix();
    bool isRegular();
    bool isAbsorbing();
    std::vector<std::vector<double>> randomCrisp(std::vector<std::vector<Interval>> intervalMatrix);
    Interval getAlphaCut(TriangularFuzzyNumber const& t, double alpha);
    std::vector<std::pair<double, double>> fuzzyGeneticAlgorithm(std::pair<int, int> idx, int steps, int acc, int populationSize = 1000, int threshold = 100, double selectionSample = 0.2,
                          double mutationRate = 0.01, bool timeBased = false, bool reachability = true);
    double restrictedMatrixMul(std::vector<std::vector<Interval>> intervalMatrix, std::vector<std::pair<int, int>> intervalIndices, int steps,
                               std::pair<int, int> idx, bool isMin, int populationSize = 1000, int generations = 100, double selectionSample = 0.2,
                               double mutationRate = 0.01, bool reachability = true);
    double timeBasedMatrixMul(std::vector<std::vector<Interval>> intervalMatrix, std::vector<std::pair<int, int>> intervalIndices, int steps,
                              std::pair<int, int> idx, bool isMin, int populationSize, int milliseconds, double selectionSample, double mutationRate, bool reachability = true);
    std::vector<member> initializePopulation(std::vector<std::vector<Interval>> intervalMatrix, int populationSize, int steps, std::pair<int, int> idx, bool reachability = true);
    std::vector<FuzzyAnalysisResult::member> selectPopulation(std::vector<FuzzyAnalysisResult::member> population, double selectionSample, bool isMin);
    std::vector<FuzzyAnalysisResult::member> crossPopulation(std::vector<FuzzyAnalysisResult::member> population, int populationSize,
                                                             std::vector<int> intervalRowIndices);
    std::vector<FuzzyAnalysisResult::member> mutatePopulation(std::vector<FuzzyAnalysisResult::member> population, double mutationRate,
                                                              std::vector<std::pair<int, int>> intervalIndices, std::vector<int> intervalRowIndices,
                                                              std::vector<std::vector<Interval>> intervalMatrix);
    std::vector<FuzzyAnalysisResult::member> crossParents(FuzzyAnalysisResult::member q1, FuzzyAnalysisResult::member q2, std::vector<int> intervalRowIndices);

    // bool isFeasible(std::vector<std::vector<double>> const& m, double alpha);
};
storm::Interval getFeasibleMutationRange(double mutationValue1, double mutationValue2, storm::Interval mutationRange1, storm::Interval mutationRange2);
double randomDouble(double min, double max);
int randomElement(std::vector<FuzzyAnalysisResult::member> population);
int randomIntFromList(std::vector<int> rows);
std::vector<std::vector<double>> matrixMul(std::vector<std::vector<double>> matrix, int n);
std::vector<int> randomAmountOfRandomRows(std::vector<int> rows);
FuzzyAnalysisResult::member mutateMember(FuzzyAnalysisResult::member member, std::vector<std::pair<int, int>> intervalIndices,
                                         std::vector<int> intervalRowIndices, std::vector<std::vector<Interval>> intervalMatrix);
std::vector<double> stationaryDistribution(std::vector<std::vector<double>> P);
}  // namespace storage
}  // namespace storm

#endif /* STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_ */