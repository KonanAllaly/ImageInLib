#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef GAUSSIAN_DISTRIBUTION
#define GAUSSIAN_DISTRIBUTION

#include "common_functions.h"

// Random number with Gaussian distribution (Box-Muller transform)
dataType generateRandNormal(dataType mean, dataType stddev);

#endif // !GAUSSIAN_DISTRIBUTION

#ifdef __cplusplus
}
#endif