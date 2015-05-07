#include "gtest/gtest.h"
#include "storm-config.h"

#include "src/logic/Formulas.h"
#include "src/utility/solver.h"
#include "src/modelchecker/prctl/SparseDtmcPrctlModelChecker.h"
#include "src/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "src/settings/SettingsManager.h"
#include "src/settings/SettingMemento.h"
#include "src/parser/AutoParser.h"

TEST(SparseDtmcPrctlModelCheckerTest, Die) {
    std::shared_ptr<storm::models::sparse::Model<double>> abstractModel = storm::parser::AutoParser::parseModel(STORM_CPP_BASE_PATH "/examples/dtmc/die/die.tra", STORM_CPP_BASE_PATH "/examples/dtmc/die/die.lab", "", STORM_CPP_BASE_PATH "/examples/dtmc/die/die.coin_flips.trans.rew");
    
    ASSERT_EQ(abstractModel->getType(), storm::models::ModelType::Dtmc);
    
    std::shared_ptr<storm::models::sparse::Dtmc<double>> dtmc = abstractModel->as<storm::models::sparse::Dtmc<double>>();
    
    ASSERT_EQ(dtmc->getNumberOfStates(), 13ull);
    ASSERT_EQ(dtmc->getNumberOfTransitions(), 20ull);
    
    storm::modelchecker::SparseDtmcPrctlModelChecker<double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::NativeLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("one");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult1 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult1[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("two");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult2 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult2[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("three");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult3 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult3[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    auto done = std::make_shared<storm::logic::AtomicLabelFormula>("done");
    auto reachabilityRewardFormula = std::make_shared<storm::logic::ReachabilityRewardFormula>(done);
    
    result = checker.check(*reachabilityRewardFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult4 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(3.6666650772094727, quantitativeResult4[0], storm::settings::nativeEquationSolverSettings().getPrecision());
}

TEST(SparseDtmcPrctlModelCheckerTest, Crowds) {
    std::shared_ptr<storm::models::sparse::Model<double>> abstractModel = storm::parser::AutoParser::parseModel(STORM_CPP_BASE_PATH "/examples/dtmc/crowds/crowds5_5.tra", STORM_CPP_BASE_PATH "/examples/dtmc/crowds/crowds5_5.lab", "", "");
    
    ASSERT_EQ(abstractModel->getType(), storm::models::ModelType::Dtmc);
    
    std::shared_ptr<storm::models::sparse::Dtmc<double>> dtmc = abstractModel->as<storm::models::sparse::Dtmc<double>>();
    
    ASSERT_EQ(8607ull, dtmc->getNumberOfStates());
    ASSERT_EQ(15113ull, dtmc->getNumberOfTransitions());
    
    storm::modelchecker::SparseDtmcPrctlModelChecker<double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::NativeLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observe0Greater1");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult1 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(0.33288205191646525, quantitativeResult1[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observeIGreater1");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult2 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(0.15222066094730619, quantitativeResult2[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observeOnlyTrueSender");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult3 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(0.32153900158185761, quantitativeResult3[0], storm::settings::nativeEquationSolverSettings().getPrecision());
}

TEST(SparseDtmcPrctlModelCheckerTest, SynchronousLeader) {
    std::shared_ptr<storm::models::sparse::Model<double>> abstractModel = storm::parser::AutoParser::parseModel(STORM_CPP_BASE_PATH "/examples/dtmc/synchronous_leader/leader4_8.tra", STORM_CPP_BASE_PATH "/examples/dtmc/synchronous_leader/leader4_8.lab", "", STORM_CPP_BASE_PATH "/examples/dtmc/synchronous_leader/leader4_8.pick.trans.rew");
    
    ASSERT_EQ(abstractModel->getType(), storm::models::ModelType::Dtmc);
    std::shared_ptr<storm::models::sparse::Dtmc<double>> dtmc = abstractModel->as<storm::models::sparse::Dtmc<double>>();
    
    ASSERT_EQ(12400ull, dtmc->getNumberOfStates());
    ASSERT_EQ(16495ull, dtmc->getNumberOfTransitions());
    
    storm::modelchecker::SparseDtmcPrctlModelChecker<double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::NativeLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult1 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(1.0, quantitativeResult1[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto trueFormula = std::make_shared<storm::logic::BooleanLiteralFormula>(true);
    auto boundedUntilFormula = std::make_shared<storm::logic::BoundedUntilFormula>(trueFormula, labelFormula, 20);
    
    result = checker.check(*boundedUntilFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult2 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(0.9999965911265462636, quantitativeResult2[0], storm::settings::nativeEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto reachabilityRewardFormula = std::make_shared<storm::logic::ReachabilityRewardFormula>(labelFormula);
    
    result = checker.check(*reachabilityRewardFormula);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeResult3 = result->asExplicitQuantitativeCheckResult<double>();
    
    EXPECT_NEAR(1.0448979589010925, quantitativeResult3[0], storm::settings::nativeEquationSolverSettings().getPrecision());
}