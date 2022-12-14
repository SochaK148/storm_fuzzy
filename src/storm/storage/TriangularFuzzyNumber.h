#ifndef STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_
#define STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_

#include <iostream>
#include <string>

namespace storm {
namespace storage {
class TriangularFuzzyNumber{
    private:
        double leftBound;
        double peak;
        double rightBound;
    public:
        TriangularFuzzyNumber(double v) : leftBound(0), peak(v), rightBound(0) {} // crisp number

        TriangularFuzzyNumber() = default;

        TriangularFuzzyNumber(double lb, double p, double rb) : leftBound(lb), peak(p), rightBound(rb) {}

        double getLeftBound() const;

        double getPeak() const;

        double getRightBound() const;

        void setLeftBound(double const& lb);

        void setPeak(double const& p);

        void setRightBound(double const& rb);

        TriangularFuzzyNumber& operator += (const TriangularFuzzyNumber& j);
        TriangularFuzzyNumber& operator *= (const TriangularFuzzyNumber& j);
        TriangularFuzzyNumber& operator -= (const TriangularFuzzyNumber& j);
        TriangularFuzzyNumber& operator /= (const TriangularFuzzyNumber& j);
};

std::ostream& operator<<(std::ostream& os, TriangularFuzzyNumber const& i);
TriangularFuzzyNumber operator + (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i);
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator * (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
TriangularFuzzyNumber operator / (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator > (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator < (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator == (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
bool operator != (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j);
//TriangularFuzzyNumber parseTriangularFuzzyNumber(std::string const& stringRepr);
std::size_t hash_value(TriangularFuzzyNumber const& i);
}  // namespace storage
}  // namespace storm


#endif /* STORM_STORAGE_TRIANGULARFUZZYNUMBER_H_ */