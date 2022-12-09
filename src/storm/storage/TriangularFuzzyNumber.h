#ifndef STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_
#define STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_

#include <iostream>
#include <string>

namespace storm {
namespace storage {
class TriangularFuzzyNumber{
    public:
        explicit TriangularFuzzyNumber(double v) : leftBound(v), peak(v), rightBound(v) {}

        TriangularFuzzyNumber() = default;

        TriangularFuzzyNumber(double lb, double p, double rb) : leftBound(lb), peak(p), rightBound(rb) {}

        double getLeftBound() const {
            return leftBound;
        }

        double getPeak() const {
            return peak;
        }

        double getRightBound() const {
            return rightBound;
        }

    private:
        double leftBound;
        double peak;
        double rightBound;

};

std::ostream& operator<<(std::ostream& os, TriangularFuzzyNumber const& i);
TriangularFuzzyNumber operator + (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator += (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i);
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator * (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator *= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator / (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator -= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator /= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator > (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator < (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator == (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator != (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
//TriangularFuzzyNumber parseTriangularFuzzyNumber(std::string const& stringRepr);
std::size_t hash_value(TriangularFuzzyNumber const& i);
}  // namespace storage
}  // namespace storm


#endif /* STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_ */