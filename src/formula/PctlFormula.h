/*
 * Pctlformula.h
 *
 *  Created on: 19.10.2012
 *      Author: Thomas Heinemann
 */

#ifndef STORM_FORMULA_PCTLFORMULA_H_
#define STORM_FORMULA_PCTLFORMULA_H_

#include <string>

namespace storm {

namespace formula {


//abstract
/*!
 * @brief
 * Abstract base class for PCTL formulas in general.
 *
 * @attention This class is abstract.
 * @note Formula classes do not have copy constructors. The parameters of the constructors are usually the subtrees, so
 * 	   the syntax conflicts with copy constructors for unary operators. To produce an identical object, the classes
 * 	   PctlPathFormula and PctlStateFormula offer the method clone().
 */
template <class T>
class PctlFormula {

public:
	/*!
	 * virtual destructor
	 */
	virtual ~PctlFormula() { }

	/*!
	 * @note This function is not implemented in this class.
	 * @returns a string representation of the formula
	 */
	virtual std::string toString() const = 0;
};

} //namespace formula

} //namespace storm

#endif /* STORM_FORMULA_PCTLFORMULA_H_ */
