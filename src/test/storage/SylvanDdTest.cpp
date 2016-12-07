#include "gtest/gtest.h"
#include "storm-config.h"

#include "storm/adapters/CarlAdapter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/storage/dd/DdManager.h"
#include "storm/storage/dd/Add.h"
#include "storm/storage/dd/Odd.h"
#include "storm/storage/dd/DdMetaVariable.h"
#include "storm/settings/SettingsManager.h"

#include "storm/storage/SparseMatrix.h"

#include <memory>
#include <iostream>

TEST(SylvanDd, Constants) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    storm::dd::Add<storm::dd::DdType::Sylvan, double> zero;
    ASSERT_NO_THROW(zero = manager->template getAddZero<double>());
    
    EXPECT_EQ(0ul, zero.getNonZeroCount());
    EXPECT_EQ(1ul, zero.getLeafCount());
    EXPECT_EQ(1ul, zero.getNodeCount());
    EXPECT_EQ(0, zero.getMin());
    EXPECT_EQ(0, zero.getMax());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> one;
    ASSERT_NO_THROW(one = manager->template getAddOne<double>());
    
    EXPECT_EQ(0ul, one.getNonZeroCount());
    EXPECT_EQ(1ul, one.getLeafCount());
    EXPECT_EQ(1ul, one.getNodeCount());
    EXPECT_EQ(1, one.getMin());
    EXPECT_EQ(1, one.getMax());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> two;
    ASSERT_NO_THROW(two = manager->template getConstant<double>(2));
    
    EXPECT_EQ(0ul, two.getNonZeroCount());
    EXPECT_EQ(1ul, two.getLeafCount());
    EXPECT_EQ(1ul, two.getNodeCount());
    EXPECT_EQ(2, two.getMin());
    EXPECT_EQ(2, two.getMax());
}

TEST(SylvanDd, BddConstants) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    storm::dd::Bdd<storm::dd::DdType::Sylvan> zero;
    ASSERT_NO_THROW(zero = manager->getBddZero());
    
    EXPECT_EQ(0ul, zero.getNonZeroCount());
    EXPECT_EQ(1ul, zero.getLeafCount());
    EXPECT_EQ(1ul, zero.getNodeCount());
    
    storm::dd::Bdd<storm::dd::DdType::Sylvan> one;
    ASSERT_NO_THROW(one = manager->getBddOne());
    
    EXPECT_EQ(0ul, one.getNonZeroCount());
    EXPECT_EQ(1ul, one.getLeafCount());
    EXPECT_EQ(1ul, one.getNodeCount());
}

TEST(SylvanDd, BddExistAbstractRepresentative) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> zero;
    ASSERT_NO_THROW(zero = manager->getBddZero());
	storm::dd::Bdd<storm::dd::DdType::Sylvan> one;
    ASSERT_NO_THROW(one = manager->getBddOne());
	
	std::pair<storm::expressions::Variable, storm::expressions::Variable> x;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> y;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> z;
	ASSERT_NO_THROW(x = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(y = manager->addMetaVariable("y", 0, 1));
	ASSERT_NO_THROW(z = manager->addMetaVariable("z", 0, 1));
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(x.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(x.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(y.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(y.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ0 = manager->getEncoding(z.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ1 = manager->getEncoding(z.first, 1);
	
	// Abstract from FALSE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_false_x = zero.existsAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_false_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_false_x.getLeafCount());
    EXPECT_EQ(1ul, representative_false_x.getNodeCount());
	EXPECT_TRUE(representative_false_x == zero);
	
	// Abstract from TRUE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_x = one.existsAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_true_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_x.getLeafCount());
    EXPECT_EQ(2ul, representative_true_x.getNodeCount());
	EXPECT_TRUE(representative_true_x == bddX0);
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_xyz = one.existsAbstractRepresentative({x.first, y.first, z.first});
	EXPECT_EQ(0ul, representative_true_xyz.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_xyz.getLeafCount());
    EXPECT_EQ(4ul, representative_true_xyz.getNodeCount());
	EXPECT_TRUE(representative_true_xyz == ((bddX0 && bddY0) && bddZ0));

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1Y0Z0 = (bddX1 && bddY0) && bddZ0;
	EXPECT_EQ(1ul, bddX1Y0Z0.getNonZeroCount());
	EXPECT_EQ(1ul, bddX1Y0Z0.getLeafCount());
	EXPECT_EQ(4ul, bddX1Y0Z0.getNodeCount());

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_x = bddX1Y0Z0.existsAbstractRepresentative({x.first});
	EXPECT_EQ(1ul, representative_x.getNonZeroCount());
	EXPECT_EQ(1ul, representative_x.getLeafCount());
	EXPECT_EQ(4ul, representative_x.getNodeCount());
	EXPECT_TRUE(bddX1Y0Z0 == representative_x);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_y = bddX1Y0Z0.existsAbstractRepresentative({y.first});
	EXPECT_EQ(1ul, representative_y.getNonZeroCount());
	EXPECT_EQ(1ul, representative_y.getLeafCount());
	EXPECT_EQ(4ul, representative_y.getNodeCount());
	EXPECT_TRUE(bddX1Y0Z0 == representative_y);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_z = bddX1Y0Z0.existsAbstractRepresentative({z.first});
	EXPECT_EQ(1ul, representative_z.getNonZeroCount());
	EXPECT_EQ(1ul, representative_z.getLeafCount());
	EXPECT_EQ(4ul, representative_z.getNodeCount());
	EXPECT_TRUE(bddX1Y0Z0 == representative_z);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_xyz = bddX1Y0Z0.existsAbstractRepresentative({x.first, y.first, z.first});
	EXPECT_EQ(1ul, representative_xyz.getNonZeroCount());
	EXPECT_EQ(1ul, representative_xyz.getLeafCount());
	EXPECT_EQ(4ul, representative_xyz.getNodeCount());
	EXPECT_TRUE(bddX1Y0Z0 == representative_xyz);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0Y0Z0 = (bddX0 && bddY0) && bddZ0;
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1Y1Z1 = (bddX1 && bddY1) && bddZ1;

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddAllTrueOrAllFalse = bddX0Y0Z0 || bddX1Y1Z1;

	representative_x = bddAllTrueOrAllFalse.existsAbstractRepresentative({x.first});
	EXPECT_EQ(2ul, representative_x.getNonZeroCount());
	EXPECT_EQ(1ul, representative_x.getLeafCount());
	EXPECT_EQ(5ul, representative_x.getNodeCount());
	EXPECT_TRUE(bddAllTrueOrAllFalse == representative_x);

	representative_y = bddAllTrueOrAllFalse.existsAbstractRepresentative({y.first});
	EXPECT_EQ(2ul, representative_y.getNonZeroCount());
	EXPECT_EQ(1ul, representative_y.getLeafCount());
	EXPECT_EQ(5ul, representative_y.getNodeCount());
	EXPECT_TRUE(bddAllTrueOrAllFalse == representative_y);

	representative_z = bddAllTrueOrAllFalse.existsAbstractRepresentative({z.first});
	EXPECT_EQ(2ul, representative_z.getNonZeroCount());
	EXPECT_EQ(1ul, representative_z.getLeafCount());
	EXPECT_EQ(5ul, representative_z.getNodeCount());
	EXPECT_TRUE(bddAllTrueOrAllFalse == representative_z);

	representative_xyz = bddAllTrueOrAllFalse.existsAbstractRepresentative({x.first, y.first, z.first});
	EXPECT_EQ(1ul, representative_xyz.getNonZeroCount());
	EXPECT_EQ(1ul, representative_xyz.getLeafCount());
	EXPECT_EQ(4ul, representative_xyz.getNodeCount());
	EXPECT_TRUE(bddX0Y0Z0 == representative_xyz);
}

TEST(SylvanDd, AddMinAbstractRepresentative) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZero;
    ASSERT_NO_THROW(bddZero = manager->getBddZero());
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddOne;
    ASSERT_NO_THROW(bddOne = manager->getBddOne());

	storm::dd::Add<storm::dd::DdType::Sylvan, double> addZero;
    ASSERT_NO_THROW(addZero = manager->template getAddZero<double>());
	storm::dd::Add<storm::dd::DdType::Sylvan, double> addOne;
    ASSERT_NO_THROW(addOne = manager->template getAddOne<double>());

	std::pair<storm::expressions::Variable, storm::expressions::Variable> x;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> y;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> z;
	ASSERT_NO_THROW(x = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(y = manager->addMetaVariable("y", 0, 1));
	ASSERT_NO_THROW(z = manager->addMetaVariable("z", 0, 1));

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(x.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(x.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(y.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(y.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ0 = manager->getEncoding(z.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ1 = manager->getEncoding(z.first, 1);

	storm::dd::Add<storm::dd::DdType::Sylvan, double> complexAdd = 
          ((bddX1 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.4))
        + ((bddX1 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.7))
        + ((bddX1 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.3))
        + ((bddX1 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.3))
        + ((bddX0 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.9))
        + ((bddX0 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.5))
        + ((bddX0 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(1.0))
        + ((bddX0 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.0));

	// Abstract from FALSE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_false_x = addZero.minAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_false_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_false_x.getLeafCount());
    EXPECT_EQ(2ul, representative_false_x.getNodeCount());
	EXPECT_TRUE(representative_false_x == bddX0);

	// Abstract from TRUE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_x = addOne.minAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_true_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_x.getLeafCount());
    EXPECT_EQ(2ul, representative_true_x.getNodeCount());
	EXPECT_TRUE(representative_true_x == bddX0);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_xyz = addOne.minAbstractRepresentative({x.first, y.first, z.first});
	EXPECT_EQ(0ul, representative_true_xyz.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_xyz.getLeafCount());
    EXPECT_EQ(4ul, representative_true_xyz.getNodeCount());
	EXPECT_TRUE(representative_true_xyz == ((bddX0 && bddY0) && bddZ0));
	
	// Abstract x
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_x = complexAdd.minAbstractRepresentative({x.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_x = (
	     ((bddX0 && (bddY0 && bddZ0)))
	  || ((bddX1 && (bddY0 && bddZ1)))
	  || ((bddX0 && (bddY1 && bddZ0)))
	  || ((bddX1 && (bddY1 && bddZ1)))
	);
	EXPECT_EQ(4ul, representative_complex_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_x.getLeafCount());
    EXPECT_EQ(3ul, representative_complex_x.getNodeCount());
	EXPECT_TRUE(representative_complex_x == comparison_complex_x);

	// Abstract y
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_y = complexAdd.minAbstractRepresentative({y.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_y = (
	     ((bddX0 && (bddY0 && bddZ0)))
	  || ((bddX0 && (bddY1 && bddZ1)))
	  || ((bddX1 && (bddY0 && bddZ0)))
	  || ((bddX1 && (bddY0 && bddZ1)))
	);
	EXPECT_EQ(4ul, representative_complex_y.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_y.getLeafCount());
    EXPECT_EQ(5ul, representative_complex_y.getNodeCount());
	EXPECT_TRUE(representative_complex_y == comparison_complex_y);
	
	// Abstract z
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_z = complexAdd.minAbstractRepresentative({z.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_z = (
	     ((bddX0 && (bddY0 && bddZ0)))
	  || ((bddX0 && (bddY1 && bddZ0)))
	  || ((bddX1 && (bddY0 && bddZ0)))
	  || ((bddX1 && (bddY1 && bddZ1)))
	);
	EXPECT_EQ(4ul, representative_complex_z.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_z.getLeafCount());
    EXPECT_EQ(4ul, representative_complex_z.getNodeCount());
	EXPECT_TRUE(representative_complex_z == comparison_complex_z);
	
	// Abstract x, y, z
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_xyz = complexAdd.minAbstractRepresentative({x.first, y.first, z.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_xyz = (bddX0 && (bddY0 && bddZ0));
	EXPECT_EQ(1ul, representative_complex_xyz.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_xyz.getLeafCount());
    EXPECT_EQ(4ul, representative_complex_xyz.getNodeCount());
	EXPECT_TRUE(representative_complex_xyz == comparison_complex_xyz);
}

TEST(SylvanDd, AddMaxAbstractRepresentative) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZero;
    ASSERT_NO_THROW(bddZero = manager->getBddZero());
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddOne;
    ASSERT_NO_THROW(bddOne = manager->getBddOne());

	storm::dd::Add<storm::dd::DdType::Sylvan, double> addZero;
    ASSERT_NO_THROW(addZero = manager->template getAddZero<double>());
	storm::dd::Add<storm::dd::DdType::Sylvan, double> addOne;
    ASSERT_NO_THROW(addOne = manager->template getAddOne<double>());

	std::pair<storm::expressions::Variable, storm::expressions::Variable> x;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> y;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> z;
	ASSERT_NO_THROW(x = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(y = manager->addMetaVariable("y", 0, 1));
	ASSERT_NO_THROW(z = manager->addMetaVariable("z", 0, 1));

	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(x.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(x.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(y.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(y.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ0 = manager->getEncoding(z.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ1 = manager->getEncoding(z.first, 1);

	storm::dd::Add<storm::dd::DdType::Sylvan, double> complexAdd = 
          ((bddX1 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.4))
        + ((bddX1 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.7))
        + ((bddX1 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.3))
        + ((bddX1 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.3))
        + ((bddX0 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.9))
        + ((bddX0 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.5))
        + ((bddX0 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(1.0))
        + ((bddX0 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.0));

	// Abstract from FALSE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_false_x = addZero.maxAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_false_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_false_x.getLeafCount());
    EXPECT_EQ(2ul, representative_false_x.getNodeCount());
	EXPECT_TRUE(representative_false_x == bddX0);

	// Abstract from TRUE
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_x = addOne.maxAbstractRepresentative({x.first});
	EXPECT_EQ(0ul, representative_true_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_x.getLeafCount());
    EXPECT_EQ(2ul, representative_true_x.getNodeCount());
	EXPECT_TRUE(representative_true_x == bddX0);

	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_true_xyz = addOne.maxAbstractRepresentative({x.first, y.first, z.first});
	EXPECT_EQ(0ul, representative_true_xyz.getNonZeroCount());
    EXPECT_EQ(1ul, representative_true_xyz.getLeafCount());
    EXPECT_EQ(4ul, representative_true_xyz.getNodeCount());
	EXPECT_TRUE(representative_true_xyz == ((bddX0 && bddY0) && bddZ0));

	// Abstract x
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_x = complexAdd.maxAbstractRepresentative({x.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_x = (
	     ((bddX1 && (bddY0 && bddZ0)))
	  || ((bddX0 && (bddY0 && bddZ1)))
	  || ((bddX1 && (bddY1 && bddZ0)))
	  || ((bddX0 && (bddY1 && bddZ1)))
	);
	EXPECT_EQ(4ul, representative_complex_x.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_x.getLeafCount());
    EXPECT_EQ(3ul, representative_complex_x.getNodeCount());
	EXPECT_TRUE(representative_complex_x == comparison_complex_x);

	// Abstract y
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_y = complexAdd.maxAbstractRepresentative({y.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_y = (
	     ((bddX0 && (bddY1 && bddZ0)))
	  || ((bddX0 && (bddY0 && bddZ1)))
	  || ((bddX1 && (bddY1 && bddZ0)))
	  || ((bddX1 && (bddY1 && bddZ1)))
	);
	EXPECT_EQ(4ul, representative_complex_y.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_y.getLeafCount());
    EXPECT_EQ(5ul, representative_complex_y.getNodeCount());
	EXPECT_TRUE(representative_complex_y == comparison_complex_y);

	// Abstract z
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_z = complexAdd.maxAbstractRepresentative({z.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_z = (
	     ((bddX0 && (bddY0 && bddZ1)))
	  || ((bddX0 && (bddY1 && bddZ1)))
	  || ((bddX1 && (bddY0 && bddZ0)))
	  || ((bddX1 && (bddY1 && bddZ0)))
	);
	EXPECT_EQ(4ul, representative_complex_z.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_z.getLeafCount());
    EXPECT_EQ(3ul, representative_complex_z.getNodeCount());
	EXPECT_TRUE(representative_complex_z == comparison_complex_z);

	// Abstract x, y, z
	storm::dd::Bdd<storm::dd::DdType::Sylvan> representative_complex_xyz = complexAdd.maxAbstractRepresentative({x.first, y.first, z.first});
	storm::dd::Bdd<storm::dd::DdType::Sylvan> comparison_complex_xyz = (bddX0 && (bddY0 && bddZ1));
	EXPECT_EQ(1ul, representative_complex_xyz.getNonZeroCount());
    EXPECT_EQ(1ul, representative_complex_xyz.getLeafCount());
    EXPECT_EQ(4ul, representative_complex_xyz.getNodeCount());
	EXPECT_TRUE(representative_complex_xyz == comparison_complex_xyz);
}

TEST(SylvanDd, AddGetMetaVariableTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    ASSERT_NO_THROW(manager->addMetaVariable("x", 1, 9));
    EXPECT_EQ(2ul, manager->getNumberOfMetaVariables());
    
    ASSERT_THROW(manager->addMetaVariable("x", 0, 3), storm::exceptions::InvalidArgumentException);
    
    ASSERT_NO_THROW(manager->addMetaVariable("y", 0, 3));
    EXPECT_EQ(4ul, manager->getNumberOfMetaVariables());
    
    EXPECT_TRUE(manager->hasMetaVariable("x'"));
    EXPECT_TRUE(manager->hasMetaVariable("y'"));
    
    std::set<std::string> metaVariableSet = {"x", "x'", "y", "y'"};
    EXPECT_EQ(metaVariableSet, manager->getAllMetaVariableNames());
}

TEST(SylvanDd, EncodingTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Bdd<storm::dd::DdType::Sylvan> encoding;
    ASSERT_THROW(encoding = manager->getEncoding(x.first, 0), storm::exceptions::InvalidArgumentException);
    ASSERT_THROW(encoding = manager->getEncoding(x.first, 10), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(encoding = manager->getEncoding(x.first, 4));
    EXPECT_EQ(1ul, encoding.getNonZeroCount());
    
    // As a BDD, this DD has one only leaf, because there does not exist a 0-leaf, and (consequently) one node less
    // than the MTBDD.
    EXPECT_EQ(5ul, encoding.getNodeCount());
    EXPECT_EQ(1ul, encoding.getLeafCount());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> add;
    ASSERT_NO_THROW(add = encoding.template toAdd<double>());

    // As an MTBDD, the 0-leaf is there, so the count is actually 2 and the node count is 6.
    EXPECT_EQ(6ul, add.getNodeCount());
    EXPECT_EQ(2ul, add.getLeafCount());
}

#ifdef STORM_HAVE_CARL
TEST(SylvanDd, RationalFunctionLeaveReplacementNonVariable) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> zero;
    ASSERT_NO_THROW(zero = manager->template getAddZero<storm::RationalFunction>());
	
	std::map<storm::RationalFunctionVariable, std::pair<storm::expressions::Variable, std::pair<storm::RationalNumber, storm::RationalNumber>>> replacementMap;
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> zeroReplacementResult = zero.replaceLeaves(replacementMap);
	
	EXPECT_EQ(0ul, zeroReplacementResult.getNonZeroCount());
    EXPECT_EQ(1ul, zeroReplacementResult.getLeafCount());
    EXPECT_EQ(1ul, zeroReplacementResult.getNodeCount());
	EXPECT_TRUE(zeroReplacementResult == zero);
}

TEST(SylvanDd, RationalFunctionLeaveReplacementSimpleVariable) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> function;
	
	std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
	
	carl::StringParser parser;
	parser.setVariables({"x"});
	
	storm::RawPolynomial polyX = parser.parseMultivariatePolynomial<storm::RationalNumber>("x");	
	storm::RationalFunction variableX = storm::RationalFunction(storm::Polynomial(polyX, cache));
	
	ASSERT_NO_THROW(function = manager->template getConstant<storm::RationalFunction>(variableX));
	
	std::pair<storm::expressions::Variable, storm::expressions::Variable> xExpr;
	ASSERT_NO_THROW(xExpr = manager->addMetaVariable("x", 0, 1));
	
	std::map<storm::RationalFunctionVariable, std::pair<storm::expressions::Variable, std::pair<storm::RationalNumber, storm::RationalNumber>>> replacementMap;
	storm::RationalNumber rnOneThird = storm::RationalNumber(1) / storm::RationalNumber(3);
	storm::RationalNumber rnTwoThird = storm::RationalNumber(2) / storm::RationalNumber(3);
	replacementMap.insert(std::make_pair(parser.variables().find("x")->second, std::make_pair(xExpr.first, std::make_pair(rnOneThird, rnTwoThird))));
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> replacedAddSimpleX = function.replaceLeaves(replacementMap);
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(xExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(xExpr.first, 1);

	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> complexAdd = 
          (bddX0.template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnTwoThird)))
        + (bddX1.template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnOneThird)));
	
	EXPECT_EQ(2ul, replacedAddSimpleX.getNonZeroCount());
    EXPECT_EQ(2ul, replacedAddSimpleX.getLeafCount());
    EXPECT_EQ(3ul, replacedAddSimpleX.getNodeCount());
	EXPECT_TRUE(replacedAddSimpleX == complexAdd);
	
	storm::dd::Add<storm::dd::DdType::Sylvan, double> abstractedAddMax = replacedAddSimpleX.toDouble().maxAbstract({xExpr.first});
	storm::dd::Add<storm::dd::DdType::Sylvan, double> abstractedAddMin = replacedAddSimpleX.toDouble().minAbstract({xExpr.first});
	EXPECT_TRUE(abstractedAddMax == manager->template getConstant<double>(0.66666666666666666666));
	EXPECT_TRUE(abstractedAddMin == manager->template getConstant<double>(0.33333333333333333333));
}

TEST(SylvanDd, RationalFunctionLeaveReplacementTwoVariables) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> function;
	std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
	carl::StringParser parser;
	parser.setVariables({"x", "y"});
	
	storm::RationalFunction variableX(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("x"), cache));
	storm::RationalFunction variableY(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("y"), cache));
	ASSERT_NO_THROW(function = manager->template getConstant<storm::RationalFunction>(variableX * variableY));
	
	std::pair<storm::expressions::Variable, storm::expressions::Variable> xExpr;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> yExpr;
	ASSERT_NO_THROW(xExpr = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(yExpr = manager->addMetaVariable("y", 0, 1));
	
	std::map<storm::RationalFunctionVariable, std::pair<storm::expressions::Variable, std::pair<storm::RationalNumber, storm::RationalNumber>>> replacementMap;
	storm::RationalNumber rnOneThird = storm::RationalNumber(1) / storm::RationalNumber(3);
	storm::RationalNumber rnTwoThird = storm::RationalNumber(2) / storm::RationalNumber(3);
	storm::RationalNumber rnOne = storm::RationalNumber(1);
	storm::RationalNumber rnTen = storm::RationalNumber(10);
	replacementMap.insert(std::make_pair(parser.variables().find("x")->second, std::make_pair(xExpr.first, std::make_pair(rnOneThird, rnTwoThird))));
	replacementMap.insert(std::make_pair(parser.variables().find("y")->second, std::make_pair(yExpr.first, std::make_pair(rnOne, rnTen))));
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> replacedAdd = function.replaceLeaves(replacementMap);
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(xExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(xExpr.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(yExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(yExpr.first, 1);

	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> complexAdd = 
          ((bddX0 && bddY0).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnTwoThird * rnTen)))
        + ((bddX0 && bddY1).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnTwoThird)))
		+ ((bddX1 && bddY0).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnOneThird * rnTen)))
		+ ((bddX1 && bddY1).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(rnOneThird)));
	
	EXPECT_EQ(4ul, replacedAdd.getNonZeroCount());
    EXPECT_EQ(4ul, replacedAdd.getLeafCount());
    EXPECT_EQ(7ul, replacedAdd.getNodeCount());
	EXPECT_TRUE(replacedAdd == complexAdd);
}

TEST(SylvanDd, RationalFunctionCarlSubstituteTest) {
	std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
	carl::StringParser parser;
	parser.setVariables({"x", "y", "z"});
	storm::RationalFunction zHalfed = storm::RationalFunction(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("z"), cache), storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("2"), cache));
	storm::RationalFunction rationalFunction = storm::RationalFunction(storm::Polynomial(parser.parseMultivariatePolynomial<storm::RationalNumber>("2*x+x*y"), cache), storm::Polynomial(parser.parseMultivariatePolynomial<storm::RationalNumber>("1"), cache)) + zHalfed;
	
	std::map<storm::RationalFunctionVariable, storm::RationalNumber> replacement = {{parser.variables().find("x")->second, storm::RationalNumber(2)}};
	storm::RationalFunction subX = rationalFunction.substitute(replacement);
	
	storm::RationalFunction cmp = storm::RationalFunction(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("4+2*y"), cache), storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("1"), cache)) + zHalfed;
	EXPECT_EQ(subX, cmp);
	
	storm::RawPolynomial poly(parser.template parseMultivariatePolynomial<storm::RationalNumber>("2*x+x*y") + parser.template parseMultivariatePolynomial<storm::RationalNumber>("z") / storm::RationalNumber(2));
	storm::RawPolynomial polySub = poly.substitute(replacement);
	EXPECT_EQ(polySub, parser.template parseMultivariatePolynomial<storm::RationalNumber>("4+2*y") + parser.template parseMultivariatePolynomial<storm::RationalNumber>("z") / storm::RationalNumber(2));
}

TEST(SylvanDd, RationalFunctionLeaveReplacementComplexFunction) {
	std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> function;
	std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
	carl::StringParser parser;
	parser.setVariables({"x", "y", "z"});
	storm::RationalFunction zDivTwoY = storm::RationalFunction(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("z"), cache), storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("2*y"), cache));

	storm::RationalFunction rationalFunction = storm::RationalFunction(storm::Polynomial(parser.parseMultivariatePolynomial<storm::RationalNumber>("2*x+x*y"), cache), storm::Polynomial(parser.parseMultivariatePolynomial<storm::RationalNumber>("1"), cache)) + zDivTwoY;	
	ASSERT_NO_THROW(function = manager->template getConstant<storm::RationalFunction>(rationalFunction));
	
	std::pair<storm::expressions::Variable, storm::expressions::Variable> xExpr;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> yExpr;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> zExpr;
	ASSERT_NO_THROW(xExpr = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(yExpr = manager->addMetaVariable("y", 0, 1));
	ASSERT_NO_THROW(zExpr = manager->addMetaVariable("z", 0, 1));
	
	std::map<storm::RationalFunctionVariable, std::pair<storm::expressions::Variable, std::pair<storm::RationalNumber, storm::RationalNumber>>> replacementMap;
	storm::RationalNumber rnTwo(2);
	storm::RationalNumber rnThree(3);
	storm::RationalNumber rnFive(5);
	storm::RationalNumber rnSeven(7);
	storm::RationalNumber rnEleven(11);
	storm::RationalNumber rnThirteen(13);
	replacementMap.insert(std::make_pair(parser.variables().find("x")->second, std::make_pair(xExpr.first, std::make_pair(rnTwo, rnSeven))));
	replacementMap.insert(std::make_pair(parser.variables().find("y")->second, std::make_pair(yExpr.first, std::make_pair(rnThree, rnEleven))));
	replacementMap.insert(std::make_pair(parser.variables().find("z")->second, std::make_pair(zExpr.first, std::make_pair(rnFive, rnThirteen))));
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> replacedAdd = function.replaceLeaves(replacementMap);
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(xExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(xExpr.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(yExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(yExpr.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ0 = manager->getEncoding(zExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ1 = manager->getEncoding(zExpr.first, 1);

	auto f = [&](bool x, bool y, bool z) {
		storm::RationalNumber result(2);
		if (!x) {
			result *= rnSeven;
		} else {
			result *= rnTwo;
		}
		
		storm::RationalNumber partTwo(1);
		if (!x) {
			partTwo *= rnSeven;
		} else {
			partTwo *= rnTwo;
		}
		if (!y) {
			partTwo *= rnEleven;
		} else {
			partTwo *= rnThree;
		}
		
		storm::RationalNumber partThree(1);
		if (!z) {
			partThree *= rnThirteen;
		} else {
			partThree *= rnFive;
		}
		if (!y) {
			partThree /= storm::RationalNumber(2) * rnEleven;
		} else {
			partThree /= storm::RationalNumber(2) * rnThree;
		}
		
		return result + partTwo + partThree;
	};
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> complexAdd = 
          ((bddX0 && (bddY0 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(false, false, false))))
		+ ((bddX0 && (bddY0 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(false, false, true))))
		+ ((bddX0 && (bddY1 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(false, true, false))))
		+ ((bddX0 && (bddY1 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(false, true, true))))
		+ ((bddX1 && (bddY0 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(true, false, false))))
		+ ((bddX1 && (bddY0 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(true, false, true))))
		+ ((bddX1 && (bddY1 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(true, true, false))))
		+ ((bddX1 && (bddY1 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(f(true, true, true))));
	
	EXPECT_EQ(8ul, replacedAdd.getNonZeroCount());
    EXPECT_EQ(8ul, replacedAdd.getLeafCount());
    EXPECT_EQ(15ul, replacedAdd.getNodeCount());
	EXPECT_TRUE(replacedAdd == complexAdd);
}

TEST(SylvanDd, RationalFunctionConstants) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> zero;
    ASSERT_NO_THROW(zero = manager->template getAddZero<storm::RationalFunction>());
    
    EXPECT_EQ(0ul, zero.getNonZeroCount());
    EXPECT_EQ(1ul, zero.getLeafCount());
    EXPECT_EQ(1ul, zero.getNodeCount());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> one;
    ASSERT_NO_THROW(one = manager->template getAddOne<storm::RationalFunction>());
    
    EXPECT_EQ(0ul, one.getNonZeroCount());
    EXPECT_EQ(1ul, one.getLeafCount());
    EXPECT_EQ(1ul, one.getNodeCount());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> two;
	storm::RationalFunction constantTwo(2);
	
    ASSERT_NO_THROW(two = manager->template getConstant<storm::RationalFunction>(constantTwo));
    
    EXPECT_EQ(0ul, two.getNonZeroCount());
    EXPECT_EQ(1ul, two.getLeafCount());
    EXPECT_EQ(1ul, two.getNodeCount());
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> function;
	std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
	carl::StringParser parser;
	parser.setVariables({"x", "y", "z"});
	
	storm::RationalFunction partA = storm::RationalFunction(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("2*x+x*y"), cache));
	storm::RationalFunction partB = storm::RationalFunction(storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("z"), cache), storm::Polynomial(parser.template parseMultivariatePolynomial<storm::RationalNumber>("2*y"), cache));

	storm::RationalFunction rationalFunction = storm::RationalFunction(partA + partB);
	
    ASSERT_NO_THROW(function = manager->template getConstant<storm::RationalFunction>(rationalFunction));
    
    EXPECT_EQ(0ul, function.getNonZeroCount());
    EXPECT_EQ(1ul, function.getLeafCount());
    EXPECT_EQ(1ul, function.getNodeCount());
}

TEST(SylvanDd, RationalFunctionToDouble) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    
	std::pair<storm::expressions::Variable, storm::expressions::Variable> xExpr;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> yExpr;
	std::pair<storm::expressions::Variable, storm::expressions::Variable> zExpr;
	ASSERT_NO_THROW(xExpr = manager->addMetaVariable("x", 0, 1));
	ASSERT_NO_THROW(yExpr = manager->addMetaVariable("y", 0, 1));
	ASSERT_NO_THROW(zExpr = manager->addMetaVariable("z", 0, 1));
	
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX0 = manager->getEncoding(xExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddX1 = manager->getEncoding(xExpr.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY0 = manager->getEncoding(yExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddY1 = manager->getEncoding(yExpr.first, 1);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ0 = manager->getEncoding(zExpr.first, 0);
	storm::dd::Bdd<storm::dd::DdType::Sylvan> bddZ1 = manager->getEncoding(zExpr.first, 1);
	
	storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> complexAdd = 
          ((bddX0 && (bddY0 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(-1))))
		+ ((bddX0 && (bddY0 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(0))))
		+ ((bddX0 && (bddY1 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(1) / storm::RationalNumber(2))))
		+ ((bddX0 && (bddY1 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(1) / storm::RationalNumber(3))))
		+ ((bddX1 && (bddY0 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(100000))))
		+ ((bddX1 && (bddY0 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(3))))
		+ ((bddX1 && (bddY1 && bddZ0)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(4))))
		+ ((bddX1 && (bddY1 && bddZ1)).template toAdd<storm::RationalFunction>() * manager->template getConstant<storm::RationalFunction>(storm::RationalFunction(storm::RationalNumber(0))));
	EXPECT_EQ(6ul, complexAdd.getNonZeroCount());
    EXPECT_EQ(7ul, complexAdd.getLeafCount());
    EXPECT_EQ(14ul, complexAdd.getNodeCount());
	
	storm::dd::Add<storm::dd::DdType::Sylvan, double> doubleAdd = complexAdd.toDouble();
    
    EXPECT_EQ(6ul, doubleAdd.getNonZeroCount());
    EXPECT_EQ(7ul, doubleAdd.getLeafCount());
    EXPECT_EQ(14ul, doubleAdd.getNodeCount());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> comparisonAdd = 
		  ((bddX0 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(-1.0))
		+ ((bddX0 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.0))
		+ ((bddX0 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(0.5))
		+ ((bddX0 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.33333333333333333333))
		+ ((bddX1 && (bddY0 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(100000.0))
		+ ((bddX1 && (bddY0 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(3.0))
		+ ((bddX1 && (bddY1 && bddZ0)).template toAdd<double>() * manager->template getConstant<double>(4.0))
		+ ((bddX1 && (bddY1 && bddZ1)).template toAdd<double>() * manager->template getConstant<double>(0.0));

    EXPECT_TRUE(comparisonAdd == doubleAdd);
}

TEST(SylvanDd, RationalFunctionEncodingTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Bdd<storm::dd::DdType::Sylvan> encoding;
    ASSERT_THROW(encoding = manager->getEncoding(x.first, 0), storm::exceptions::InvalidArgumentException);
    ASSERT_THROW(encoding = manager->getEncoding(x.first, 10), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(encoding = manager->getEncoding(x.first, 4));
    EXPECT_EQ(1ul, encoding.getNonZeroCount());
    
    // As a BDD, this DD has one only leaf, because there does not exist a 0-leaf, and (consequently) one node less
    // than the MTBDD.
    EXPECT_EQ(5ul, encoding.getNodeCount());
    EXPECT_EQ(1ul, encoding.getLeafCount());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> add;
    ASSERT_NO_THROW(add = encoding.template toAdd<storm::RationalFunction>());

    // As an MTBDD, the 0-leaf is there, so the count is actually 2 and the node count is 6.
    EXPECT_EQ(6ul, add.getNodeCount());
    EXPECT_EQ(2ul, add.getLeafCount());
}

TEST(SylvanDd, RationalFunctionIdentityTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, storm::RationalFunction> identity;
    ASSERT_NO_THROW(identity = manager->getIdentity<storm::RationalFunction>(x.first));

    EXPECT_EQ(9ul, identity.getNonZeroCount());
    EXPECT_EQ(10ul, identity.getLeafCount());
    EXPECT_EQ(21ul, identity.getNodeCount());
}
#endif

TEST(SylvanDd, RangeTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x;
    ASSERT_NO_THROW(x = manager->addMetaVariable("x", 1, 9));

    storm::dd::Bdd<storm::dd::DdType::Sylvan> range;
    ASSERT_NO_THROW(range = manager->getRange(x.first));
 
    EXPECT_EQ(9ul, range.getNonZeroCount());
    EXPECT_EQ(1ul, range.getLeafCount());
    EXPECT_EQ(5ul, range.getNodeCount());
}

TEST(SylvanDd, DoubleIdentityTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);

    storm::dd::Add<storm::dd::DdType::Sylvan, double> identity;
    ASSERT_NO_THROW(identity = manager->getIdentity<double>(x.first));

    EXPECT_EQ(9ul, identity.getNonZeroCount());
    EXPECT_EQ(10ul, identity.getLeafCount());
    EXPECT_EQ(21ul, identity.getNodeCount());
}

TEST(SylvanDd, UintIdentityTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);

    storm::dd::Add<storm::dd::DdType::Sylvan, uint_fast64_t> identity;
    ASSERT_NO_THROW(identity = manager->getIdentity<uint_fast64_t>(x.first));

    EXPECT_EQ(9ul, identity.getNonZeroCount());
    EXPECT_EQ(10ul, identity.getLeafCount());
    EXPECT_EQ(21ul, identity.getNodeCount());
}

TEST(SylvanDd, OperatorTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    EXPECT_TRUE(manager->template getAddZero<double>() == manager->template getAddZero<double>());
    EXPECT_FALSE(manager->template getAddZero<double>() == manager->template getAddOne<double>());
    
    EXPECT_FALSE(manager->template getAddZero<double>() != manager->template getAddZero<double>());
    EXPECT_TRUE(manager->template getAddZero<double>() != manager->template getAddOne<double>());
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd1 = manager->template getAddOne<double>();
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd2 = manager->template getAddOne<double>();
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd3 = dd1 + dd2;
    storm::dd::Bdd<storm::dd::DdType::Sylvan> bdd;
    EXPECT_TRUE(dd3 == manager->template getConstant<double>(2));
    
    dd3 += manager->template getAddZero<double>();
    EXPECT_TRUE(dd3 == manager->template getConstant<double>(2));
    
    dd3 = dd1 * manager->template getConstant<double>(3);
    EXPECT_TRUE(dd3 == manager->template getConstant<double>(3));
    
    dd3 *= manager->template getConstant<double>(2);
    EXPECT_TRUE(dd3 == manager->template getConstant<double>(6));
    
    dd3 = dd1 - dd2;
    EXPECT_TRUE(dd3.isZero());
    
    dd3 -= manager->template getConstant<double>(-2);
    EXPECT_TRUE(dd3 == manager->template getConstant<double>(2));
    
    dd3 /= manager->template getConstant<double>(2);
    EXPECT_TRUE(dd3.isOne());
    
    bdd = !dd3.toBdd();
    EXPECT_TRUE(bdd.isZero());
    
    bdd = !bdd;
    EXPECT_TRUE(bdd.isOne());
    
    bdd = dd1.toBdd() || dd2.toBdd();
    EXPECT_TRUE(bdd.isOne());
    
    dd1 = manager->template getIdentity<double>(x.first);
    dd2 = manager->template getConstant<double>(5);
    
    bdd = dd1.equals(dd2);
    EXPECT_EQ(1ul, bdd.getNonZeroCount());
    
    storm::dd::Bdd<storm::dd::DdType::Sylvan> bdd2 = dd1.notEquals(dd2);
    EXPECT_TRUE(bdd2 == !bdd);
    
    bdd = dd1.less(dd2);
    EXPECT_EQ(11ul, bdd.getNonZeroCount());
    
    bdd = dd1.lessOrEqual(dd2);
    EXPECT_EQ(12ul, bdd.getNonZeroCount());
    
    bdd = dd1.greater(dd2);
    EXPECT_EQ(4ul, bdd.getNonZeroCount());
    
    bdd = dd1.greaterOrEqual(dd2);
    EXPECT_EQ(5ul, bdd.getNonZeroCount());
    
    dd3 = manager->getEncoding(x.first, 2).ite(dd2, dd1);
    bdd = dd3.less(dd2);
    EXPECT_EQ(10ul, bdd.getNonZeroCount());

    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd4 = dd3.minimum(dd1);
    dd4 *= manager->getEncoding(x.first, 2).template toAdd<double>();
    dd4 = dd4.sumAbstract({x.first});
    EXPECT_EQ(2, dd4.getValue());

    dd4 = dd3.maximum(dd1);
    dd4 *= manager->getEncoding(x.first, 2).template toAdd<double>();
    dd4 = dd4.sumAbstract({x.first});
    EXPECT_EQ(5, dd4.getValue());
    
    dd1 = manager->template getConstant<double>(0.01);
    dd2 = manager->template getConstant<double>(0.01 + 1e-6);
    EXPECT_TRUE(dd1.equalModuloPrecision(dd2, 1e-6, false));
    EXPECT_FALSE(dd1.equalModuloPrecision(dd2, 1e-6));
}

TEST(SylvanDd, AbstractionTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd1;
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd2;
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd3;
    storm::dd::Bdd<storm::dd::DdType::Sylvan> bdd;
    
    dd1 = manager->template getIdentity<double>(x.first);
    dd2 = manager->template getConstant<double>(5);
    bdd = dd1.equals(dd2);
    EXPECT_EQ(1ul, bdd.getNonZeroCount());
    ASSERT_THROW(bdd = bdd.existsAbstract({x.second}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(bdd = bdd.existsAbstract({x.first}));
    EXPECT_EQ(0ul, bdd.getNonZeroCount());
    EXPECT_EQ(1, bdd.template toAdd<double>().getMax());
    
    dd3 = dd1.equals(dd2).template toAdd<double>();
    dd3 *= manager->template getConstant<double>(3);
    EXPECT_EQ(1ul, dd3.getNonZeroCount());
    ASSERT_THROW(bdd = dd3.toBdd().existsAbstract({x.second}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(bdd = dd3.toBdd().existsAbstract({x.first}));
    EXPECT_TRUE(bdd.isOne());
    
    dd3 = dd1.equals(dd2).template toAdd<double>();
    dd3 *= manager->template getConstant<double>(3);
    ASSERT_THROW(dd3 = dd3.sumAbstract({x.second}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(dd3 = dd3.sumAbstract({x.first}));
    EXPECT_EQ(0ul, dd3.getNonZeroCount());
    EXPECT_EQ(3, dd3.getMax());
    
    dd3 = dd1.equals(dd2).template toAdd<double>();
    dd3 *= manager->template getConstant<double>(3);
    ASSERT_THROW(dd3 = dd3.minAbstract({x.second}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(dd3 = dd3.minAbstract({x.first}));
    EXPECT_EQ(0ul, dd3.getNonZeroCount());
    EXPECT_EQ(0, dd3.getMax());
    
    dd3 = dd1.equals(dd2).template toAdd<double>();
    dd3 *= manager->template getConstant<double>(3);
    ASSERT_THROW(dd3 = dd3.maxAbstract({x.second}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(dd3 = dd3.maxAbstract({x.first}));
    EXPECT_EQ(0ul, dd3.getNonZeroCount());
    EXPECT_EQ(3, dd3.getMax());
}

TEST(SylvanDd, SwapTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    std::pair<storm::expressions::Variable, storm::expressions::Variable> z = manager->addMetaVariable("z", 2, 8);
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd1;
    
    dd1 = manager->template getIdentity<double>(x.first);
    ASSERT_THROW(dd1 = dd1.swapVariables({std::make_pair(x.first, z.first)}), storm::exceptions::InvalidArgumentException);
    ASSERT_NO_THROW(dd1 = dd1.swapVariables({std::make_pair(x.first, x.second)}));
    EXPECT_TRUE(dd1 == manager->template getIdentity<double>(x.second));
}

TEST(SylvanDd, MultiplyMatrixTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd1 = manager->template getIdentity<double>(x.first).equals(manager->template getIdentity<double>(x.second)).template toAdd<double>();
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd2 = manager->getRange(x.second).template toAdd<double>();
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd3;
    dd1 *= manager->template getConstant<double>(2);
    
    ASSERT_NO_THROW(dd3 = dd1.multiplyMatrix(dd2, {x.second}));
    ASSERT_NO_THROW(dd3 = dd3.swapVariables({std::make_pair(x.first, x.second)}));
    EXPECT_TRUE(dd3 == dd2 * manager->template getConstant<double>(2));
}

TEST(SylvanDd, GetSetValueTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd1 = manager->template getAddOne<double>();
    ASSERT_NO_THROW(dd1.setValue(x.first, 4, 2));
    EXPECT_EQ(2ul, dd1.getLeafCount());
    
    std::map<storm::expressions::Variable, int_fast64_t> metaVariableToValueMap;
    metaVariableToValueMap.emplace(x.first, 1);
    EXPECT_EQ(1, dd1.getValue(metaVariableToValueMap));
    
    metaVariableToValueMap.clear();
    metaVariableToValueMap.emplace(x.first, 4);
    EXPECT_EQ(2, dd1.getValue(metaVariableToValueMap));
}

TEST(SylvanDd, AddIteratorTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    std::pair<storm::expressions::Variable, storm::expressions::Variable> y = manager->addMetaVariable("y", 0, 3);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd;
    ASSERT_NO_THROW(dd = manager->getRange(x.first).template toAdd<double>());
    
    storm::dd::AddIterator<storm::dd::DdType::Sylvan, double> it, ite;
    ASSERT_NO_THROW(it = dd.begin());
    ASSERT_NO_THROW(ite = dd.end());
    std::pair<storm::expressions::SimpleValuation, double> valuationValuePair;
    uint_fast64_t numberOfValuations = 0;
    while (it != ite) {
        ASSERT_NO_THROW(valuationValuePair = *it);
        ASSERT_NO_THROW(++it);
        ++numberOfValuations;
    }
    EXPECT_EQ(9ul, numberOfValuations);
    
    dd = manager->getRange(x.first).template toAdd<double>();
    dd = dd.notZero().ite(manager->template getAddOne<double>(), manager->template getAddOne<double>());
    ASSERT_NO_THROW(it = dd.begin());
    ASSERT_NO_THROW(ite = dd.end());
    numberOfValuations = 0;
    while (it != ite) {
        ASSERT_NO_THROW(valuationValuePair = *it);
        ASSERT_NO_THROW(++it);
        ++numberOfValuations;
    }
    EXPECT_EQ(16ul, numberOfValuations);
    
    ASSERT_NO_THROW(it = dd.begin(false));
    ASSERT_NO_THROW(ite = dd.end());
    numberOfValuations = 0;
    while (it != ite) {
        ASSERT_NO_THROW(valuationValuePair = *it);
        ASSERT_NO_THROW(++it);
        ++numberOfValuations;
    }
    EXPECT_EQ(1ul, numberOfValuations);
}

TEST(SylvanDd, AddOddTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> a = manager->addMetaVariable("a");
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd = manager->template getIdentity<double>(x.first);
    storm::dd::Odd odd;
    ASSERT_NO_THROW(odd = dd.createOdd());
    EXPECT_EQ(9ul, odd.getTotalOffset());
    EXPECT_EQ(12ul, odd.getNodeCount());
    
    std::vector<double> ddAsVector;
    ASSERT_NO_THROW(ddAsVector = dd.toVector());
    EXPECT_EQ(9ul, ddAsVector.size());
    for (uint_fast64_t i = 0; i < ddAsVector.size(); ++i) {
        EXPECT_TRUE(i+1 == ddAsVector[i]);
    }
    
    // Create a non-trivial matrix.
    dd = manager->template getIdentity<double>(x.first).equals(manager->template getIdentity<double>(x.second)).template toAdd<double>() * manager->getRange(x.first).template toAdd<double>();
    dd += manager->getEncoding(x.first, 1).template toAdd<double>() * manager->getRange(x.second).template toAdd<double>() + manager->getEncoding(x.second, 1).template toAdd<double>() * manager->getRange(x.first).template toAdd<double>();
    
    // Create the ODDs.
    storm::dd::Odd rowOdd;
    ASSERT_NO_THROW(rowOdd = manager->getRange(x.first).template toAdd<double>().createOdd());
    storm::dd::Odd columnOdd;
    ASSERT_NO_THROW(columnOdd = manager->getRange(x.second).template toAdd<double>().createOdd());
    
    // Try to translate the matrix.
    storm::storage::SparseMatrix<double> matrix;
    ASSERT_NO_THROW(matrix = dd.toMatrix({x.first}, {x.second}, rowOdd, columnOdd));
    
    EXPECT_EQ(9ul, matrix.getRowCount());
    EXPECT_EQ(9ul, matrix.getColumnCount());
    EXPECT_EQ(25ul, matrix.getNonzeroEntryCount());
    
    dd = manager->getRange(x.first).template toAdd<double>() * manager->getRange(x.second).template toAdd<double>() * manager->getEncoding(a.first, 0).ite(dd, dd + manager->template getConstant<double>(1));
    ASSERT_NO_THROW(matrix = dd.toMatrix({a.first}, rowOdd, columnOdd));
    EXPECT_EQ(18ul, matrix.getRowCount());
    EXPECT_EQ(9ul, matrix.getRowGroupCount());
    EXPECT_EQ(9ul, matrix.getColumnCount());
    EXPECT_EQ(106ul, matrix.getNonzeroEntryCount());
}

TEST(SylvanDd, BddOddTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> manager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> a = manager->addMetaVariable("a");
    std::pair<storm::expressions::Variable, storm::expressions::Variable> x = manager->addMetaVariable("x", 1, 9);
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> dd = manager->template getIdentity<double>(x.first);
    storm::dd::Bdd<storm::dd::DdType::Sylvan> bdd = dd.notZero();
    storm::dd::Odd odd;
    ASSERT_NO_THROW(odd = bdd.createOdd());
    EXPECT_EQ(9ul, odd.getTotalOffset());
    EXPECT_EQ(5ul, odd.getNodeCount());

    std::vector<double> ddAsVector;
    ASSERT_NO_THROW(ddAsVector = dd.toVector());
    EXPECT_EQ(9ul, ddAsVector.size());
    for (uint_fast64_t i = 0; i < ddAsVector.size(); ++i) {
        EXPECT_EQ(i+1, ddAsVector[i]);
    }
    
    storm::dd::Add<storm::dd::DdType::Sylvan, double> vectorAdd = storm::dd::Add<storm::dd::DdType::Sylvan, double>::fromVector(*manager, ddAsVector, odd, {x.first});
    
    // Create a non-trivial matrix.
    dd = manager->template getIdentity<double>(x.first).equals(manager->template getIdentity<double>(x.second)).template toAdd<double>() * manager->getRange(x.first).template toAdd<double>();
    dd += manager->getEncoding(x.first, 1).template toAdd<double>() * manager->getRange(x.second).template toAdd<double>() + manager->getEncoding(x.second, 1).template toAdd<double>() * manager->getRange(x.first).template toAdd<double>();
    
    // Create the ODDs.
    storm::dd::Odd rowOdd;
    ASSERT_NO_THROW(rowOdd = manager->getRange(x.first).createOdd());
    storm::dd::Odd columnOdd;
    ASSERT_NO_THROW(columnOdd = manager->getRange(x.second).createOdd());
    
    // Try to translate the matrix.
    storm::storage::SparseMatrix<double> matrix;
    ASSERT_NO_THROW(matrix = dd.toMatrix({x.first}, {x.second}, rowOdd, columnOdd));
    
    EXPECT_EQ(9ul, matrix.getRowCount());
    EXPECT_EQ(9ul, matrix.getColumnCount());
    EXPECT_EQ(25ul, matrix.getNonzeroEntryCount());
    
    dd = manager->getRange(x.first).template toAdd<double>() * manager->getRange(x.second).template toAdd<double>() * manager->getEncoding(a.first, 0).ite(dd, dd + manager->template getConstant<double>(1));
    ASSERT_NO_THROW(matrix = dd.toMatrix({a.first}, rowOdd, columnOdd));
    EXPECT_EQ(18ul, matrix.getRowCount());
    EXPECT_EQ(9ul, matrix.getRowGroupCount());
    EXPECT_EQ(9ul, matrix.getColumnCount());
    EXPECT_EQ(106ul, matrix.getNonzeroEntryCount());
}

TEST(SylvanDd, BddToExpressionTest) {
    std::shared_ptr<storm::dd::DdManager<storm::dd::DdType::Sylvan>> ddManager(new storm::dd::DdManager<storm::dd::DdType::Sylvan>());
    std::pair<storm::expressions::Variable, storm::expressions::Variable> a = ddManager->addMetaVariable("a");
    std::pair<storm::expressions::Variable, storm::expressions::Variable> b = ddManager->addMetaVariable("b");
    
    storm::dd::Bdd<storm::dd::DdType::Sylvan> bdd = ddManager->getBddOne();
    bdd &= ddManager->getEncoding(a.first, 1);
    bdd |= ddManager->getEncoding(b.first, 0);
    
    std::shared_ptr<storm::expressions::ExpressionManager> manager = std::make_shared<storm::expressions::ExpressionManager>();
    storm::expressions::Variable c = manager->declareBooleanVariable("c");
    storm::expressions::Variable d = manager->declareBooleanVariable("d");
    
    auto result = bdd.toExpression(*manager);
}
