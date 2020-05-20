#pragma once

// Include these utility headers so we can access utility function from Eigen.
#include "storm/utility/constants.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

// Finally include the parts of Eigen we need.
// Make sure to include our patched version of Eigen (and not a pre-installed one e.g. located at /usr/include)
#include <resources/3rdparty/StormEigen/Eigen/Dense>
#include <resources/3rdparty/StormEigen/Eigen/Sparse>
#include <resources/3rdparty/StormEigen/unsupported/Eigen/IterativeSolvers>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
