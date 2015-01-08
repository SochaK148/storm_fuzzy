#include "src/storage/prism/Constant.h"
#include "src/utility/macros.h"
#include "src/exceptions/IllegalFunctionCallException.h"

namespace storm {
    namespace prism {
        Constant::Constant(storm::expressions::Variable const& variable, storm::expressions::Expression const& expression, std::string const& filename, uint_fast64_t lineNumber) : LocatedInformation(filename, lineNumber), variable(variable), defined(true), expression(expression) {
            // Intentionally left empty.
        }
        
        Constant::Constant(storm::expressions::Variable const& variable, std::string const& filename, uint_fast64_t lineNumber) : LocatedInformation(filename, lineNumber), variable(variable), defined(false), expression() {
            // Intentionally left empty.
        }
        
        std::string const& Constant::getName() const {
            return this->variable.getName();
        }
        
        storm::expressions::Type const& Constant::getType() const {
            return this->getExpressionVariable().getType();
        }
        
        storm::expressions::Variable const& Constant::getExpressionVariable() const {
            return this->variable;
        }
        
        bool Constant::isDefined() const {
            return this->defined;
        }
        
        storm::expressions::Expression const& Constant::getExpression() const {
            STORM_LOG_THROW(this->isDefined(), storm::exceptions::IllegalFunctionCallException, "Unable to retrieve defining expression for undefined constant.");
            return this->expression;
        }
        
        Constant Constant::substitute(std::map<storm::expressions::Variable, storm::expressions::Expression> const& substitution) const {
            return Constant(variable, this->getExpression().substitute(substitution), this->getFilename(), this->getLineNumber());
        }
        
        std::ostream& operator<<(std::ostream& stream, Constant const& constant) {
            stream << "const " << constant.getExpressionVariable().getType() << " ";
            stream << constant.getName();
            if (constant.isDefined()) {
                stream << " = " << constant.getExpression();
            }
            stream << ";";
            return stream;
        }
    } // namespace prism
} // namespace storm