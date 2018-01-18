#include "storm/settings/modules/TopologicalEquationSolverSettings.h"


#include "storm/settings/modules/CoreSettings.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Argument.h"
#include "storm/solver/SolverSelectionOptions.h"

#include "storm/storage/dd/DdType.h"

#include "storm/utility/macros.h"
#include "storm/exceptions/IllegalArgumentValueException.h"
#include "storm/exceptions/InvalidOptionException.h"

namespace storm {
    namespace settings {
        namespace modules {
            
            const std::string TopologicalEquationSolverSettings::moduleName = "topological";
            const std::string TopologicalEquationSolverSettings::underlyingEquationSolverOptionName = "eqsolver";
            const std::string TopologicalEquationSolverSettings::underlyingMinMaxMethodOptionName = "minmax";
            
            TopologicalEquationSolverSettings::TopologicalEquationSolverSettings() : ModuleSettings(moduleName) {
                std::vector<std::string> linearEquationSolver = {"gmm++", "native", "eigen", "elimination"};
                this->addOption(storm::settings::OptionBuilder(moduleName, underlyingEquationSolverOptionName, true, "Sets which solver is considered for solving the underlying equation systems.")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the used solver.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(linearEquationSolver)).setDefaultValueString("gmm++").build()).build());
                std::vector<std::string> minMaxSolvingTechniques = {"vi", "value-iteration", "pi", "policy-iteration", "linear-programming", "lp", "ratsearch", "qvi", "quick-value-iteration"};
                this->addOption(storm::settings::OptionBuilder(moduleName, underlyingMinMaxMethodOptionName, true, "Sets which minmax method is considered for solving the underlying minmax equation systems.")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the used min max method.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(minMaxSolvingTechniques)).setDefaultValueString("value-iteration").build()).build());
            }

            bool TopologicalEquationSolverSettings::isUnderlyingEquationSolverTypeSet() const {
                return this->getOption(underlyingEquationSolverOptionName).getHasOptionBeenSet();
            }
            
            bool TopologicalEquationSolverSettings::isUnderlyingEquationSolverTypeSetFromDefaultValue() const {
                return !this->getOption(underlyingEquationSolverOptionName).getHasOptionBeenSet() || this->getOption(underlyingEquationSolverOptionName).getArgumentByName("name").wasSetFromDefaultValue();
            }

            storm::solver::EquationSolverType  TopologicalEquationSolverSettings::getUnderlyingEquationSolverType() const {
                std::string equationSolverName = this->getOption(underlyingEquationSolverOptionName).getArgumentByName("name").getValueAsString();
                if (equationSolverName == "gmm++") {
                    return storm::solver::EquationSolverType::Gmmxx;
                } else if (equationSolverName == "native") {
                    return storm::solver::EquationSolverType::Native;
                } else if (equationSolverName == "eigen") {
                    return storm::solver::EquationSolverType::Eigen;
                } else if (equationSolverName == "elimination") {
                    return storm::solver::EquationSolverType::Elimination;
                }
                STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Unknown underlying equation solver '" << equationSolverName << "'.");
            }
            
            bool TopologicalEquationSolverSettings::isUnderlyingMinMaxMethodSet() const {
                return this->getOption(underlyingMinMaxMethodOptionName).getHasOptionBeenSet();
            }
            
            bool TopologicalEquationSolverSettings::isUnderlyingMinMaxMethodSetFromDefaultValue() const {
                return !this->getOption(underlyingMinMaxMethodOptionName).getHasOptionBeenSet() || this->getOption(underlyingMinMaxMethodOptionName).getArgumentByName("name").wasSetFromDefaultValue();
            }

            storm::solver::MinMaxMethod  TopologicalEquationSolverSettings::getUnderlyingMinMaxMethod() const {
                std::string minMaxEquationSolvingTechnique = this->getOption(underlyingMinMaxMethodOptionName).getArgumentByName("name").getValueAsString();
                 if (minMaxEquationSolvingTechnique == "value-iteration" || minMaxEquationSolvingTechnique == "vi") {
                    return storm::solver::MinMaxMethod::ValueIteration;
                } else if (minMaxEquationSolvingTechnique == "policy-iteration" || minMaxEquationSolvingTechnique == "pi") {
                    return storm::solver::MinMaxMethod::PolicyIteration;
                } else if (minMaxEquationSolvingTechnique == "linear-programming" || minMaxEquationSolvingTechnique == "lp") {
                    return storm::solver::MinMaxMethod::LinearProgramming;
                } else if (minMaxEquationSolvingTechnique == "ratsearch") {
                    return storm::solver::MinMaxMethod::RationalSearch;
                } else if (minMaxEquationSolvingTechnique == "quick-value-iteration" || minMaxEquationSolvingTechnique == "qvi") {
                    return storm::solver::MinMaxMethod::QuickValueIteration;
                }
                STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Unknown underlying equation solver '" << minMaxEquationSolvingTechnique << "'.");
            }
            
            bool TopologicalEquationSolverSettings::check() const {
                if (this->isUnderlyingEquationSolverTypeSet() && getUnderlyingEquationSolverType() == storm::solver::EquationSolverType::Topological) {
                    STORM_LOG_WARN("Underlying solver type of the topological solver can not be the topological solver.");
                    return false;
                }
                if (this->isUnderlyingMinMaxMethodSet() && getUnderlyingMinMaxMethod() == storm::solver::MinMaxMethod::Topological) {
                    STORM_LOG_WARN("Underlying minmax method of the topological solver can not be topological.");
                    return false;
                }
                return true;
            }

        } // namespace modules
    } // namespace settings
} // namespace storm


