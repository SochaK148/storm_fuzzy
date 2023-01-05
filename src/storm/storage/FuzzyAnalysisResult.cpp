#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "storm/storage/FuzzyAnalysisResult.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/TriangularFuzzyNumber.h"

namespace storm {
namespace storage {

// Member

std::vector<std::vector<double>> FuzzyAnalysisResult::member::getMatrix() {
    return this->matrix;
}
void FuzzyAnalysisResult::member::setMatrix(std::vector<std::vector<double>> matrix) {
    this->matrix = matrix;
}
int FuzzyAnalysisResult::member::getN() {
    return this->n;
}
std::pair<int, int> FuzzyAnalysisResult::member::getIdx() {
    return this->idx;
}
double& FuzzyAnalysisResult::member::getFitness() {
    return this->fitness;
}
void FuzzyAnalysisResult::member::updateFitness() {
    this->fitness = matrixMul(matrix, n)[idx.first][idx.second];
}
std::vector<std::vector<double>> FuzzyAnalysisResult::member::matrixMul(std::vector<std::vector<double>> matrix, int p) {
    const int n = matrix.size();
    const int m = matrix[0].size();
    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0));
    std::vector<std::vector<double>> temp = matrix;

    for (auto l = 0; l < p - 1; ++l) {
        for (auto j = 0; j < m; ++j) {
            for (auto k = 0; k < m; ++k) {
                for (auto i = 0; i < n; ++i) {
                    result[i][j] += temp[i][k] * matrix[k][j];
                }
            }
        }
        temp = result;
        for (auto& i : result) std::fill(i.begin(), i.end(), 0);
    }

    return temp;
}
bool FuzzyAnalysisResult::member::operator<(const FuzzyAnalysisResult::member& other) const {
    return fitness < other.fitness;
}
bool FuzzyAnalysisResult::member::operator==(FuzzyAnalysisResult::member& other) const {
    return matrix == other.getMatrix();
}

// Analysis result

FlexibleSparseMatrix<TriangularFuzzyNumber>& FuzzyAnalysisResult::getMatrix() {
    return matrix;
}

std::pair<std::vector<std::vector<Interval>>, std::vector<std::pair<int, int>>> FuzzyAnalysisResult::getIntervalMatrix(double alpha) {
    std::vector<std::vector<Interval>> intervalMatrix(this->getMatrix().getRowCount(), std::vector<Interval>(this->getMatrix().getColumnCount(), Interval(0)));
    std::vector<std::pair<int, int>> intervalIndices;

    for (FlexibleSparseMatrix<TriangularFuzzyNumber>::index_type row = 0; row < this->getMatrix().getRowCount(); ++row) {
        for (auto const element : this->getMatrix().getRow(row)) {
            int column = element.getColumn();
            TriangularFuzzyNumber value = element.getValue();
            if (value.getLeftBound() == 0 && value.getRightBound() == 0) {  // crisp number
                intervalMatrix[row][column] = Interval(value.getPeak());
            } else {
                intervalMatrix[row][column] = this->getAlphaCut(value, alpha);
                intervalIndices.push_back(std::pair(row, column));
            }
        }
    }

    return std::pair(intervalMatrix, intervalIndices);
}

Interval FuzzyAnalysisResult::getAlphaCut(TriangularFuzzyNumber const& t, double alpha) {
    double i = alpha * (t.getPeak() - t.getLeftBound()) + t.getLeftBound();
    double j = (1 - alpha) * (t.getRightBound() - t.getPeak()) + t.getPeak();
    return Interval(i, j);
}

void FuzzyAnalysisResult::reachableInSteps(std::pair<int, int> idx, int steps, int acc) {
    for (int i = 0; i < acc; i++) {
        double alpha = i / acc;
        std::cout << "alpha=" << alpha << " processing... " << std::endl;
        std::pair<std::vector<std::vector<Interval>>, std::vector<std::pair<int, int>>> matrixAndIndexes = this->getIntervalMatrix(alpha);
        std::vector<std::vector<Interval>> intervalMatrix = matrixAndIndexes.first;
        std::vector<std::pair<int, int>> intervalIndices = matrixAndIndexes.second;
        double left = this->restrictedMatrixMul(intervalMatrix, intervalIndices, steps, idx, true);
        std::cout << "alpha=" << alpha << ", left: " << left << std::endl;
        double right = this->restrictedMatrixMul(intervalMatrix, intervalIndices, steps, idx, false);
        std::cout << "alpha=" << alpha << ", right: " << right << std::endl;
    }
}

double FuzzyAnalysisResult::restrictedMatrixMul(std::vector<std::vector<Interval>> intervalMatrix, std::vector<std::pair<int, int>> intervalIndices, int steps,
                                                std::pair<int, int> idx, bool isMin) {
    int populationSize = 1000;
    int generations = 100;
    double selectionSample = 0.2;
    double mutationRate = 0.01;

    std::vector<int> intervalRowIndices;
    std::transform(intervalIndices.begin(), intervalIndices.end(), std::back_inserter(intervalRowIndices),
                   [](const std::pair<int, int>& p) { return p.first; });
    std::vector<FuzzyAnalysisResult::member> population = this->initializePopulation(intervalMatrix, populationSize, steps, idx);
    for (int i = 0; i < generations; i++) {
        std::cout << "Generation: " << i << std::endl;
        population =
            this->mutatePopulation(this->crossPopulation(this->selectPopulation(population, selectionSample, isMin), populationSize, intervalRowIndices),
                                   mutationRate, intervalIndices, intervalRowIndices, intervalMatrix);
    }
    return population[0].getFitness();
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::initializePopulation(std::vector<std::vector<Interval>> intervalMatrix, int populationSize,
                                                                                   int steps, std::pair<int, int> idx) {
    std::vector<FuzzyAnalysisResult::member> population;
    for (int i = 0; i < populationSize; i++) {
        std::vector<std::vector<double>> r = this->randomCrisp(intervalMatrix);
        member m(r, steps, idx);
        population.push_back(m);
    }
    return population;
}

double randomDouble(double min, double max) {
    std::uniform_real_distribution<double> unif(min, max);
    std::default_random_engine re;
    return unif(re);
}

std::vector<std::vector<double>> FuzzyAnalysisResult::randomCrisp(std::vector<std::vector<Interval>> intervalMatrix) {
    int rows = intervalMatrix.size();
    int columns = intervalMatrix[0].size();
    std::vector<std::vector<double>> r(rows, std::vector<double>(columns));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (intervalMatrix[i][j].isPointInterval()) {
                r[i][j] = intervalMatrix[i][j].lower();
            } else {
                r[i][j] = randomDouble(intervalMatrix[i][j].lower(), intervalMatrix[i][j].upper());
            }
        }
    }
    return r;
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::selectPopulation(std::vector<FuzzyAnalysisResult::member> population, double selectionSample,
                                                                               bool isMin) {
    for (auto m : population) {
        m.updateFitness();
    }
    std::sort(population.begin(), population.end());
    int n = population.size() * selectionSample;
    std::vector<FuzzyAnalysisResult::member> selected;
    if (isMin) {
        selected.insert(selected.begin(), population.begin(), population.begin() + n);
    } else {
        selected.insert(selected.begin(), population.end() - n, population.end());
    }
    std::cout << "Best 3: " << selected[0].getFitness() << ", " << selected[1].getFitness() << selected[2].getFitness() << std::endl;
    return selected;
}

int randomElement(std::vector<FuzzyAnalysisResult::member> population) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, population.size() - 1);
    return distr(gen);
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::crossPopulation(std::vector<FuzzyAnalysisResult::member> population, int populationSize,
                                                                              std::vector<int> intervalRowIndices) {
    std::vector<FuzzyAnalysisResult::member> crossed;
    while (crossed.size() < populationSize) {
        int q1 = randomElement(population);
        int q2 = q1;
        while (q2 == q1) {
            q2 = randomElement(population);
        }
        std::vector<FuzzyAnalysisResult::member> children = crossParents(population[q1], population[q2], intervalRowIndices);
        crossed.insert(crossed.end(), children.begin(), children.end());
    }
    return crossed;
}

int randomIntFromList(std::vector<int> rows) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, rows.size() - 1);
    // std::uniform_int_distribution<> unif(0,rows.size());
    // std::default_random_engine re;
    int result = distr(gen);
    return rows[result];
}

std::vector<int> randomAmountOfRandomRows(std::vector<int> rows) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, rows.size() - 1);
    int amount = distr(gen);
    std::vector<int> randomRows;
    for (int i = 0; i < amount; i++) {
        int row = randomIntFromList(rows);
        randomRows.push_back(row);
        auto it = std::find(rows.begin(), rows.end(), row);
        rows.erase(it);
    }
    return randomRows;
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::crossParents(FuzzyAnalysisResult::member q1, FuzzyAnalysisResult::member q2,
                                                                           std::vector<int> intervalRowIndices) {
    std::vector<int> randomRows = randomAmountOfRandomRows(intervalRowIndices);
    FuzzyAnalysisResult::member c1 = member(q1.getMatrix(), q1.getN(), q1.getIdx());
    FuzzyAnalysisResult::member c2 = member(q2.getMatrix(), q2.getN(), q2.getIdx());
    for (auto const randomRow : randomRows) {
        std::vector<std::vector<double>> c1New = c1.getMatrix();
        c1New[randomRow] = q2.getMatrix()[randomRow];
        c1.setMatrix(c1New);

        std::vector<std::vector<double>> c2New = c2.getMatrix();
        c2New[randomRow] = q1.getMatrix()[randomRow];
        c2.setMatrix(c2New);
    }
    return {c1, c2};
}

Interval getFeasibleMutationRange(int mutationValue1, int mutationValue2, Interval mutationRange1, Interval mutationRange2) {
    int left = std::min(mutationValue1 - mutationRange1.lower(), mutationRange2.upper() - mutationValue2);
    int right = std::min(mutationRange1.upper() - mutationValue1, mutationValue2 - mutationRange2.lower());
    return Interval(-left, right);
}

FuzzyAnalysisResult::member mutateMember(FuzzyAnalysisResult::member member, std::vector<std::pair<int, int>> intervalIndices,
                                         std::vector<int> intervalRowIndices, std::vector<std::vector<Interval>> intervalMatrix) {
    int row = randomIntFromList(intervalRowIndices);
    std::vector<int> elementsInRandomRow;
    for (auto const intervalIndex : intervalIndices) {
        if (intervalIndex.first == row) {
            elementsInRandomRow.push_back(intervalIndex.second);
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, elementsInRandomRow.size() - 1);
    // std::uniform_int_distribution<> unif(0,elementsInRandomRow.size());
    // std::default_random_engine re;
    int toMutate1 = elementsInRandomRow[distr(gen)];
    int toMutate2 = toMutate1;
    while (toMutate1 == toMutate2) {
        toMutate2 = elementsInRandomRow[distr(gen)];
    }

    std::vector<std::vector<double>> mutableMatrix = member.getMatrix();
    int mutationValue1 = mutableMatrix[row][toMutate1];
    int mutationValue2 = mutableMatrix[row][toMutate2];
    Interval mutationRange1 = intervalMatrix[row][toMutate1];
    Interval mutationRange2 = intervalMatrix[row][toMutate2];

    Interval feasibleMutationRange = getFeasibleMutationRange(mutationValue1, mutationValue2, mutationRange1, mutationRange2);
    double mutationOffset = randomDouble(feasibleMutationRange.lower(), feasibleMutationRange.upper());

    mutableMatrix[row][toMutate1] += mutationOffset;
    mutableMatrix[row][toMutate2] -= mutationOffset;
    member.setMatrix(mutableMatrix);

    return member;
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::mutatePopulation(std::vector<FuzzyAnalysisResult::member> population, double mutationRate,
                                                                               std::vector<std::pair<int, int>> intervalIndices,
                                                                               std::vector<int> intervalRowIndices,
                                                                               std::vector<std::vector<Interval>> intervalMatrix) {
    int amount = mutationRate * population.size();
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> indices(population.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), gen);

    for (int i = 0; i < amount; i++) {
        int index = indices[i];
        population[index] = mutateMember(population[index], intervalIndices, intervalRowIndices, intervalMatrix);
    }

    return population;
}

}  // namespace storage
}  // namespace storm