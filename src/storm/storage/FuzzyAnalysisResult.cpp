#include <iostream>
#include <vector>

#include "storm/storage/TriangularFuzzyNumber.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/FuzzyAnalysisResult.h"

namespace storm {
namespace storage {

    // Member

    std::vector<std::vector<double>> FuzzyAnalysisResult::member::getMatrix(){
        return this->matrix;
    }
    double FuzzyAnalysisResult::member::getFitness(){
        return this->fitness;
    }
    void FuzzyAnalysisResult::member::updateFitness(){
        this->fitness = matrixMul(matrix, n)[idx.first][idx.second];
    }
    std::vector<std::vector<double>> FuzzyAnalysisResult::member::matrixMul(std::vector<std::vector<double>> matrix, int n){
        // TODO
        return matrix;
    }
    bool FuzzyAnalysisResult::member::operator< (const FuzzyAnalysisResult::member &other) const {
        return fitness < other.fitness;
    }

    // Analysis result

    storm::storage::SparseMatrix<storm::storage::TriangularFuzzyNumber> FuzzyAnalysisResult::getMatrix() const {
        return this->matrix;
    }

    std::vector<std::vector<Interval>> FuzzyAnalysisResult::getMatrixI() const {
        return this->matrixI;
    }

    void FuzzyAnalysisResult::calcMatrixI(double alpha) {
        std::vector<std::vector<Interval>> mI;
        this->matrixI = mI;
    }

    Interval FuzzyAnalysisResult::getAlphaCut(storm::storage::TriangularFuzzyNumber const& t, double alpha){
        double i = alpha * (t.getPeak() - t.getLeftBound()) + t.getLeftBound();
        double j = alpha * (t.getRightBound() - t.getPeak()) + t.getPeak();
        return Interval(i, j);
    }

    double FuzzyAnalysisResult::restrictedMatrixMul(std::vector<std::vector<Interval>> intvlP, int steps, std::pair<int, int> idx, bool isMin){
        int populationSize = 1000;
        int generations = 10000;
        int selectionSample = 0.2;
        int mutationRate = 0.01;
        std::vector<member> population = this->initializePopulation(populationSize, steps, idx);
        for(int i = 0; i < generations; i++){
            population = this->mutatePopulation(this->crossPopulation(this->selectPopulation(population, selectionSample, isMin)), mutationRate);
        }
        return population[0].getFitness();
    }

    std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::initializePopulation(int populationSize, int steps, std::pair<int, int> idx){
        std::vector<member> population;
        for(int i = 0; i < populationSize; i++){
            std::vector<std::vector<double>> r = this->randomCrisp();
            member m(r, steps, idx);
            population.push_back(m);
        }
        return population;
    }

    std::vector<std::vector<double>> FuzzyAnalysisResult::randomCrisp(){
        std::vector<std::vector<Interval>> mi = this->getMatrixI();
        int size = mi.size();
        std::vector<std::vector<double>> r(size, std::vector<double>(size));
        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                // r[i][j] = fRand(mi[i][j].lower(), mi[i][j].upper());
            }
        }
        return r;
    }

    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }

    std::vector<FuzzyAnalysisResult::member> storm::storage::FuzzyAnalysisResult::selectPopulation(std::vector<FuzzyAnalysisResult::member> population, int selectionSample, bool isMin){
        std::sort(population.begin(), population.end());
        int n = population.size() * selectionSample;
        if(isMin){
            std::vector<FuzzyAnalysisResult::member> selected(population.begin(), population.begin() + n);
        } else {
            std::vector<FuzzyAnalysisResult::member> selected(population.end() - n, population.end());
        }
        // return selected;
        return population;
    }

    std::vector<FuzzyAnalysisResult::member> storm::storage::FuzzyAnalysisResult::crossPopulation(std::vector<FuzzyAnalysisResult::member> population){
        // TODO
        return population;
    }

    std::vector<FuzzyAnalysisResult::member> storm::storage::FuzzyAnalysisResult::mutatePopulation(std::vector<FuzzyAnalysisResult::member> population, int mutationRate){
        // TODO
        return population;
    }
    
    bool storm::storage::FuzzyAnalysisResult::isFeasible(std::vector<std::vector<double>> const& m, double alpha){
        for (SparseMatrix<TriangularFuzzyNumber>::index_type row = 0; row < this->getMatrix().getRowCount(); ++row) {
            std::cout << "Row: " << row << std::endl;
            std::cout << (this->getMatrix()).begin(row)->getColumn() << std::endl;
            std::cout << (this->getMatrix()).begin(row)->getValue() << std::endl;
            for (SparseMatrix<TriangularFuzzyNumber>::const_iterator it = (this->getMatrix()).begin(row); it != (this->getMatrix()).end(row); ++it) {
                SparseMatrix<TriangularFuzzyNumber>::index_type column = it->getColumn();
                SparseMatrix<TriangularFuzzyNumber>::value_type value = it->getValue();
                std::cout << "Column: " << column << ", Value: " << value << std::endl;
                // Interval alphaCut = this->getAlphaCut(value, alpha);
                // double crispValue = m[row][column];
                // TODO check if crispValue is within alphaCut
            }
        }
        return true; // placeholder
    }
    
}  // namespace storage
}  // namespace storm