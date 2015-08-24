#include "gtest/gtest.h"
#include "storm-config.h"
#include "src/storage/BitVector.h"
#include "src/utility/vector.h"

TEST(VectorTest, sum_if) {
    std::vector<double> a = {1.0, 2.0, 4.0, 8.0, 16.0};
    storm::storage::BitVector f1(5, {2,4});
    storm::storage::BitVector f2(5, {3,4});
    
    ASSERT_EQ(20.0, storm::utility::vector::sum_if(a, f1));
    ASSERT_EQ(24.0, storm::utility::vector::sum_if(a, f2));
    
}
