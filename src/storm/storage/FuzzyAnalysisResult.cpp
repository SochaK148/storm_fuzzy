#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <thread>

#include "storm/storage/FuzzyAnalysisResult.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/TriangularFuzzyNumber.h"

namespace storm {
namespace storage {

void print2DVector(const std::vector<std::vector<double>> &vec)
{
    for (const auto &row : vec)
    {
        for (const auto &element : row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

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
bool FuzzyAnalysisResult::member::getReachability() {
    return this->reachability;
}
double FuzzyAnalysisResult::member::getFitness() {
    return this->fitness;
}
void FuzzyAnalysisResult::member::updateFitness() {
    if(this->reachability)
        this->fitness = matrixMul(matrix, n)[idx.first][idx.second];
    else
        this->fitness = stationaryDistribution(matrix)[idx.second];
}

std::vector<double> stationaryDistribution(std::vector<std::vector<double>> P)
{
    int n = P.size();
    std::vector<double> w(n);
    for (int i = 0; i < n; i++) {
        w[i] = 1.0;
    }

    std::vector<double> prev_w(n);
    double prev_diff = 0;
    double diff = 1;
    while (diff > 1e-8) {
        prev_w = w;
        for (int i = 0; i < n; i++) {
            w[i] = 0;
            for (int j = 0; j < n; j++) {
                w[i] += prev_w[j] * P[j][i];
            }
        }

        prev_diff = diff;
        diff = 0;
        for (int i = 0; i < n; i++) {
            diff += abs(w[i] - prev_w[i]);
        }
    }

    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += w[i];
    }
    for (int i = 0; i < n; i++) {
        w[i] /= sum;
    }

    return w;
}

std::vector<std::vector<double>> matrixMul(std::vector<std::vector<double>> matrix, int p) {
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

bool FuzzyAnalysisResult::isRegular()
{
    std::vector<std::vector<double>> P = this->getCrispMatrix();
    int n = P.size();

    std::vector<int> reachable(n);
    reachable[0] = 1;
    for (int k = 0; k < n; k++)
    {
        std::vector<int> new_reachable(reachable);
        for (int i = 0; i < n; i++)
        {
            if (reachable[i] == 1)
            {
                for (int j = 0; j < n; j++)
                {
                    if (P[i][j] > 0)
                    {
                        new_reachable[j] = 1;
                    }
                }
            }
        }
        reachable = new_reachable;
    }
    bool irreducible = true;
    for (int i = 0; i < n; i++)
    {
        if (reachable[i] == 0)
        {
            irreducible = false;
            break;
        }
    }
    if (!irreducible)
        return false;

    std::vector<int> d(n);
    for (int i = 0; i < n; i++)
    {
        d[i] = 1;
        for (int k = 1;; k++)
        {
            if (P[i][i] > 0)
            {
                d[i] = k;
                break;
            }
            std::vector<double> p(n);
            for (int j = 0; j < n; j++)
            {
                for (int l = 0; l < n; l++)
                {
                    p[j] += P[l][j] * P[i][l];
                }
            }
            P[i] = p;
        }
    }
    bool aperiodic = true;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (std::__gcd(d[i], d[j]) != 1)
            {
                aperiodic = false;
                break;
            }
        }
        if (!aperiodic)
            break;
    }
    if (!aperiodic)
        return false;

    return true;
}

bool FuzzyAnalysisResult::isAbsorbing()
{
    std::vector<std::vector<double>> P = this->getCrispMatrix();

    int n = P.size();
    int k = 0;
    std::vector<int> absorbing_states;
    for (int i = 0; i < n; i++)
    {
        if (P[i][i] == 1)
        {
            absorbing_states.push_back(i);
            k++;
        }
    }
    if (k == 0)
    {
        return false;
    }

    for (int i = 0; i < n; i++)
    {
        bool reachable = false;
        for (int j = 0; j < absorbing_states.size(); j++)
        {
            if (i == absorbing_states[j])
            {
                reachable = true;
                break;
            }
        }
        if (!reachable)
        {
            std::vector<double> p(n);
            for (int k = 0; k < n; k++)
            {
                p[k] = P[i][k];
            }
            for (int k = 0; k < n; k++)
            {
                for (int l = 0; l < n; l++)
                {
                    p[k] += P[i][l] * P[l][k];
                }
            }
            for (int j = 0; j < absorbing_states.size(); j++)
            {
                if (p[absorbing_states[j]] > 0)
                {
                    reachable = true;
                    break;
                }
            }
            if (!reachable)
            {
                return false;
            }
        }
    }

    return true;
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

std::vector<std::vector<double>> FuzzyAnalysisResult::getCrispMatrix() {
    std::vector<std::vector<double>> crispMatrix(this->getMatrix().getRowCount(), std::vector<double>(this->getMatrix().getColumnCount()));
    
    for (FlexibleSparseMatrix<TriangularFuzzyNumber>::index_type row = 0; row < this->getMatrix().getRowCount(); ++row) {
        for (auto const element : this->getMatrix().getRow(row)) {
            int column = element.getColumn();
            TriangularFuzzyNumber value = element.getValue();
            crispMatrix[row][column] = value.getPeak();
        }
    }

    return crispMatrix;
}

Interval FuzzyAnalysisResult::getAlphaCut(TriangularFuzzyNumber const& t, double alpha) {
    double i = alpha * (t.getPeak() - t.getLeftBound()) + t.getLeftBound();
    double j = (1 - alpha) * (t.getRightBound() - t.getPeak()) + t.getPeak();
    return Interval(i, j);
}

std::vector<std::pair<double, double>> FuzzyAnalysisResult::fuzzyGeneticAlgorithm(std::pair<int, int> idx, int steps, int acc, int populationSize, int threshold,
                                                                             double selectionSample, double mutationRate, bool timeBased, bool reachability) {
    std::vector<std::pair<double, double>> result;
    for (int i = 0; i < acc; i++) {
        double alpha = (double)i / (double)acc;
        double left;
        double right;
        std::cout << "alpha=" << alpha << " processing... " << std::endl;
        std::pair<std::vector<std::vector<Interval>>, std::vector<std::pair<int, int>>> matrixAndIndexes = this->getIntervalMatrix(alpha);
        std::vector<std::vector<Interval>> intervalMatrix = matrixAndIndexes.first;
        std::vector<std::pair<int, int>> intervalIndices = matrixAndIndexes.second;
        if (timeBased) {
            left = this->timeBasedMatrixMul(intervalMatrix, intervalIndices, steps, idx, true, populationSize, threshold, selectionSample, mutationRate);
            std::cout << "alpha=" << alpha << ", left: " << left << std::endl;
            right =  this->timeBasedMatrixMul(intervalMatrix, intervalIndices, steps, idx, false, populationSize, threshold, selectionSample, mutationRate);
            std::cout << "alpha=" << alpha << ", right: " << right << std::endl;

        } else {
            left = this->restrictedMatrixMul(intervalMatrix, intervalIndices, steps, idx, true, populationSize, threshold, selectionSample, mutationRate);
            std::cout << "alpha=" << alpha << ", left: " << left << std::endl;
            right = this->restrictedMatrixMul(intervalMatrix, intervalIndices, steps, idx, false, populationSize, threshold, selectionSample, mutationRate);
            std::cout << "alpha=" << alpha << ", right: " << right << std::endl;
        }
        result.push_back({left, alpha});
        result.push_back({right, alpha});
    }
    result.push_back({matrixMul(this->getCrispMatrix(), steps)[idx.first][idx.second], 1.0});
    std::sort(result.begin(), result.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
        return a.first < b.first;
    });
    return result;
}

double FuzzyAnalysisResult::restrictedMatrixMul(std::vector<std::vector<Interval>> intervalMatrix, std::vector<std::pair<int, int>> intervalIndices, int steps,
                                                std::pair<int, int> idx, bool isMin, int populationSize, int generations, double selectionSample,
                                                double mutationRate, bool reachability) {
    std::vector<int> intervalRowIndices;
    std::transform(intervalIndices.begin(), intervalIndices.end(), std::back_inserter(intervalRowIndices),
                   [](const std::pair<int, int>& p) { return p.first; });
    std::sort(intervalRowIndices.begin(), intervalRowIndices.end());
    auto last = std::unique(intervalRowIndices.begin(), intervalRowIndices.end());
    intervalRowIndices.erase(last, intervalRowIndices.end());
    std::vector<FuzzyAnalysisResult::member> population = this->initializePopulation(intervalMatrix, populationSize, steps, idx, reachability);

    for (int i = 0; i < generations; i++) {
        std::cout << "Generation: " << i << std::endl;
        population =
            this->mutatePopulation(this->crossPopulation(this->selectPopulation(population, selectionSample, isMin), populationSize, intervalRowIndices),
                                   mutationRate, intervalIndices, intervalRowIndices, intervalMatrix);
    }

    population = this->selectPopulation(population, selectionSample, isMin);
    if (isMin)
        return population[0].getFitness();
    return population[population.size() - 1].getFitness();
}

double FuzzyAnalysisResult::timeBasedMatrixMul(std::vector<std::vector<Interval>> intervalMatrix, std::vector<std::pair<int, int>> intervalIndices, int steps,
                                               std::pair<int, int> idx, bool isMin, int populationSize, int milliseconds, double selectionSample,
                                               double mutationRate, bool reachability, int mutationRows) {

    std::ofstream log_file;
    std::string file_name = "time-results-5.txt";
    log_file.open(file_name, std::ofstream::out | std::ofstream::app);
    log_file << std::fixed << std::setprecision(10);

    std::vector<int> intervalRowIndices;
    std::transform(intervalIndices.begin(), intervalIndices.end(), std::back_inserter(intervalRowIndices),
                   [](const std::pair<int, int>& p) { return p.first; });
    std::sort(intervalRowIndices.begin(), intervalRowIndices.end());
    auto last = std::unique(intervalRowIndices.begin(), intervalRowIndices.end());
    intervalRowIndices.erase(last, intervalRowIndices.end());
    std::vector<FuzzyAnalysisResult::member> population = this->initializePopulation(intervalMatrix, populationSize, steps, idx, reachability);
    auto start = std::chrono::high_resolution_clock::now();
    auto end = start + std::chrono::milliseconds(milliseconds);
    int i = 0;
    double current_result;
    auto now = std::chrono::high_resolution_clock::now();
    while (now < end) {
        population = this->selectPopulation(population, selectionSample, isMin);
        current_result = isMin ? population[0].getFitness() : population[population.size() - 1].getFitness();
        population = this->mutatePopulation(this->crossPopulation(population, populationSize, intervalRowIndices), mutationRate, intervalIndices,
                                            intervalRowIndices, intervalMatrix, mutationRows);
        now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
        std::cout << "Generation: " << i++ << ", result: " << current_result << ", time: " << duration.count() << std::endl;
        log_file << "Generation: " << i++ << ", result: " << current_result << ", time: " << duration.count() << std::endl;
    }

    population = this->selectPopulation(population, selectionSample, isMin);
    if (isMin){
        //print2DVector(population[0].getMatrix());
        std::this_thread::sleep_for(std::chrono::milliseconds(5000));
        return population[0].getFitness();
    }
    //print2DVector(population[population.size() - 1].getMatrix());
    std::this_thread::sleep_for(std::chrono::milliseconds(5000));
    return population[population.size() - 1].getFitness();
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::initializePopulation(std::vector<std::vector<Interval>> intervalMatrix, int populationSize,
                                                                                   int steps, std::pair<int, int> idx, bool reachability) {
    std::vector<FuzzyAnalysisResult::member> population;
    for (int i = 0; i < populationSize; i++) {
        std::vector<std::vector<double>> r = this->randomCrisp(intervalMatrix);
        member m(r, steps, idx, reachability);
        population.push_back(m);
    }
    return population;
}

double randomDouble(double min, double max) {
    std::uniform_real_distribution<double> unif(min, max);
    std::random_device rd;
    std::default_random_engine re(rd());
    return unif(re);
}

std::vector<std::vector<double>> FuzzyAnalysisResult::randomCrisp(std::vector<std::vector<Interval>> intervalMatrix) {
    int rows = intervalMatrix.size();
    int columns = intervalMatrix[0].size();
    std::vector<std::vector<double>> crispMatrix;
    for (int i = 0; i < rows; i++) {
        std::vector<double> ranges_lengths;
        std::vector<double> random_values;
        for (auto r : intervalMatrix[i]) {
            double random_value;
            ranges_lengths.push_back(r.upper() - r.lower());
            random_value = randomDouble(r.lower(), r.upper());
            random_values.push_back(random_value);
        }
        double difference = 1.0 - std::accumulate(random_values.begin(), random_values.end(), 0.0);
        if (difference != 0.0) {
            std::vector<double> margins;
            std::vector<double> distributed_corrections;
            if (difference > 0.0) {
                for (int j = 0; j < random_values.size(); j++) {
                    margins.push_back(intervalMatrix[i][j].upper() - random_values[j]);
                }
            } else {
                for (int j = 0; j < random_values.size(); j++) {
                    margins.push_back(random_values[j] - intervalMatrix[i][j].lower());
                }
            }
            double margin_sum = std::accumulate(margins.begin(), margins.end(), 0.0);
            for (int j = 0; j < margins.size(); j++) {
                distributed_corrections.push_back(difference * (margins[j] / margin_sum));
                random_values[j] += distributed_corrections[j];
            }
        }
        crispMatrix.push_back(random_values);
    }
    return crispMatrix;
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::selectPopulation(std::vector<FuzzyAnalysisResult::member> population, double selectionSample,
                                                                               bool isMin) {
    std::cout << std::fixed;
    for (auto& m : population) {
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

    return selected;
}

int randomElement(std::vector<FuzzyAnalysisResult::member> population) {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_int_distribution<> distr(0, population.size() - 1);
    return distr(re);
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
    std::default_random_engine re(rd());
    std::uniform_int_distribution<> distr(0, rows.size() - 1);
    int result = distr(re);
    return rows[result];
}

std::vector<int> randomAmountOfRandomRows(std::vector<int> rows) {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_int_distribution<> distr(0, rows.size() - 1);
    int amount = distr(re);
    std::vector<int> randomRows;
    for (int i = 0; i < amount; i++) {
        int row = randomIntFromList(rows);
        randomRows.push_back(row);
        auto it = std::find(rows.begin(), rows.end(), row);
        rows.erase(it);
    }
    return randomRows;
}

std::vector<int> fixedAmountOfRandomRows(std::vector<int> rows, int amount) {
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
    FuzzyAnalysisResult::member c1 = member(q1.getMatrix(), q1.getN(), q1.getIdx(), q1.getReachability());
    FuzzyAnalysisResult::member c2 = member(q2.getMatrix(), q2.getN(), q2.getIdx(), q2.getReachability());
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

Interval getFeasibleMutationRange(double mutationValue1, double mutationValue2, Interval mutationRange1, Interval mutationRange2) {
    double left = std::min(mutationValue1 - mutationRange1.lower(), mutationRange2.upper() - mutationValue2);
    double right = std::min(mutationRange1.upper() - mutationValue1, mutationValue2 - mutationRange2.lower());
    return Interval(-left, right);
}

FuzzyAnalysisResult::member mutateMember(FuzzyAnalysisResult::member member, std::vector<std::pair<int, int>> intervalIndices,
                                         std::vector<int> intervalRowIndices, std::vector<std::vector<Interval>> intervalMatrix, int mutationRows) {
    
    std::vector<int> rowsToMutate = fixedAmountOfRandomRows(intervalRowIndices, mutationRows);

    for(int i = 0; i < rowsToMutate.size(); i++)
    {
        int row = rowsToMutate[i];
        std::vector<int> elementsInRandomRow;
        for (auto const intervalIndex : intervalIndices) {
            if (intervalIndex.first == row) {
                elementsInRandomRow.push_back(intervalIndex.second);
            }
        }

        std::random_device rd;
        std::default_random_engine re(rd());
        std::uniform_int_distribution<> distr(0, elementsInRandomRow.size() - 1);
        int toMutate1 = elementsInRandomRow[distr(re)];
        int toMutate2 = toMutate1;
        while (toMutate1 == toMutate2) {
            toMutate2 = elementsInRandomRow[distr(re)];
        }

        std::vector<std::vector<double>> mutableMatrix = member.getMatrix();
        double mutationValue1 = mutableMatrix[row][toMutate1];
        double mutationValue2 = mutableMatrix[row][toMutate2];
        Interval mutationRange1 = intervalMatrix[row][toMutate1];
        Interval mutationRange2 = intervalMatrix[row][toMutate2];

        Interval feasibleMutationRange = getFeasibleMutationRange(mutationValue1, mutationValue2, mutationRange1, mutationRange2);
        double mutationOffset = randomDouble(feasibleMutationRange.lower(), feasibleMutationRange.upper());

        mutableMatrix[row][toMutate1] += mutationOffset;
        mutableMatrix[row][toMutate2] -= mutationOffset;
        member.setMatrix(mutableMatrix);
    }
    return member;
}

std::vector<FuzzyAnalysisResult::member> FuzzyAnalysisResult::mutatePopulation(std::vector<FuzzyAnalysisResult::member> population, double mutationRate,
                                                                               std::vector<std::pair<int, int>> intervalIndices,
                                                                               std::vector<int> intervalRowIndices,
                                                                               std::vector<std::vector<Interval>> intervalMatrix, int mutationRows) {
    int amount = mutationRate * population.size();
    std::random_device rd;
    std::default_random_engine re(rd());

    std::vector<int> indices(population.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), re);

    for (int i = 0; i < amount; i++) {
        int index = indices[i];
        population[index] = mutateMember(population[index], intervalIndices, intervalRowIndices, intervalMatrix, mutationRows);
    }

    return population;
}

}  // namespace storage
}  // namespace storm