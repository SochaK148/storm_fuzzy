#include "src/storage/expressions/SimpleValuation.h"

#include <boost/functional/hash.hpp>

namespace storm {
    namespace expressions {
        SimpleValuation::SimpleValuation(std::size_t booleanVariableCount, std::size_t integerVariableCount, std::size_t doubleVariableCount) : identifierToIndexMap(new std::unordered_map<std::string, uint_fast64_t>), booleanValues(booleanVariableCount), integerValues(integerVariableCount), doubleValues(doubleVariableCount) {
            // Intentionally left empty.
        }
        
        SimpleValuation::SimpleValuation(std::shared_ptr<std::unordered_map<std::string, uint_fast64_t>> identifierToIndexMap, std::vector<bool> booleanValues, std::vector<int_fast64_t> integerValues, std::vector<double> doubleValues) : identifierToIndexMap(identifierToIndexMap), booleanValues(booleanValues), integerValues(integerValues), doubleValues(doubleValues) {
            // Intentionally left empty.
        }
        
        bool SimpleValuation::operator==(SimpleValuation const& other) const {
            return this->identifierToIndexMap.get() == other.identifierToIndexMap.get() && this->booleanValues == other.booleanValues && this->integerValues == other.integerValues && this->doubleValues == other.doubleValues;
        }

        void SimpleValuation::setIdentifierIndex(std::string const& name, uint_fast64_t index) {
            (*this->identifierToIndexMap)[name] = index;
        }
        
        void SimpleValuation::setBooleanValue(std::string const& name, bool value) {
            this->booleanValues[this->identifierToIndexMap->at(name)] = value;
        }
        
        void SimpleValuation::setIntegerValue(std::string const& name, int_fast64_t value) {
            this->integerValues[this->identifierToIndexMap->at(name)] = value;
        }
        
        void SimpleValuation::setDoubleValue(std::string const& name, double value) {
            this->doubleValues[this->identifierToIndexMap->at(name)] = value;
        }
        
        bool SimpleValuation::getBooleanValue(std::string const& name) const {
            auto const& nameIndexPair = this->identifierToIndexMap->find(name);
            return this->booleanValues[nameIndexPair->second];
        }
        
        int_fast64_t SimpleValuation::getIntegerValue(std::string const& name) const {
            auto const& nameIndexPair = this->identifierToIndexMap->find(name);
            return this->integerValues[nameIndexPair->second];
        }
        
        double SimpleValuation::getDoubleValue(std::string const& name) const {
            auto const& nameIndexPair = this->identifierToIndexMap->find(name);
            return this->doubleValues[nameIndexPair->second];
        }
        
        std::ostream& operator<<(std::ostream& stream, SimpleValuation const& valuation) {
            stream << "valuation { bool[";
            for (uint_fast64_t i = 0; i < valuation.booleanValues.size() - 1; ++i) {
                stream << valuation.booleanValues[i] << ", ";
            }
            stream << valuation.booleanValues.back() << "] ints[";
            for (uint_fast64_t i = 0; i < valuation.integerValues.size() - 1; ++i) {
                stream << valuation.integerValues[i] << ", ";
            }
            stream << valuation.integerValues.back() << "] double[";
            for (uint_fast64_t i = 0; i < valuation.doubleValues.size() - 1; ++i) {
                stream << valuation.doubleValues[i] << ", ";
            }
            stream << valuation.doubleValues.back() << "] }";
            
            return stream;
        }
        
        std::size_t SimpleValuationPointerHash::operator()(SimpleValuation* valuation) const {
            size_t seed = 0;
            for (auto const& value : valuation->booleanValues) {
                boost::hash_combine<bool>(seed, value);
            }
            for (auto const& value : valuation->integerValues) {
                boost::hash_combine<int_fast64_t>(seed, value);
            }
            for (auto const& value : valuation->doubleValues) {
                boost::hash_combine<double>(seed, value);
            }
            return seed;
        }
        
        bool SimpleValuationPointerCompare::operator()(SimpleValuation* valuation1, SimpleValuation* valuation2) const {
            return *valuation1 == *valuation2;
        }
        
        bool SimpleValuationPointerLess::operator()(SimpleValuation* valuation1, SimpleValuation* valuation2) const {
            // Compare boolean variables.
            bool less = valuation1->booleanValues < valuation2->booleanValues;
            if (less) {
                return true;
            }
            less = valuation1->integerValues < valuation2->integerValues;
            if (less) {
                return true;
            }
            less = valuation1->doubleValues < valuation2->doubleValues;
            if (less) {
                return true;
            } else {
                return false;
            }
        }
    }
}