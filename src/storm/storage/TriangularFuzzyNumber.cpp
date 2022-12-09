#include "TriangularFuzzyNumber.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <boost/functional/hash.hpp>
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/macros.h"

namespace storm {
namespace storage {
/*
TriangularFuzzyNumber parseTriangularFuzzyNumber(std::string const& stringRepr) {
    return TriangularFuzzyNumber(0, 0, 0);
}
*/
TriangularFuzzyNumber operator + (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "+ not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator += (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "+= not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Unary - not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator - (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Binary - not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator * (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "* not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator *= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "*= not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator / (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "/ not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator -= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "-= not implemented for fuzzy numbers" );
}
TriangularFuzzyNumber operator /= (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "/= not implemented for fuzzy numbers" );
}
bool operator == (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    return i.getLeftBound() == j.getLeftBound() && i.getPeak() == j.getPeak() && i.getRightBound() == j.getRightBound();
}
bool operator != (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    return i.getLeftBound() != j.getLeftBound() || i.getPeak() != j.getPeak() || i.getRightBound() != j.getRightBound();
}
bool operator > (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    return i.getPeak() > j.getPeak();
}
bool operator < (const TriangularFuzzyNumber& i, const TriangularFuzzyNumber& j) {
    return i.getPeak() < j.getPeak();
}

std::size_t hash_value(TriangularFuzzyNumber const& i){
    boost::hash<double> hasher;
    return hasher(i.getPeak());
}
std::ostream& operator<<(std::ostream& os, TriangularFuzzyNumber const& i){
    os << '(' << i.getLeftBound() << '/' << i.getPeak() << '/' << i.getRightBound() << ')';
    return os;
}
}  // namespace storage
}  // namespace storm