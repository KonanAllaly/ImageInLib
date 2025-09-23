
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "gaussian_distribution.h"

#define M_PI 3.14159265358979323846

// Random number with Gaussian distribution (Box-Muller transform)
dataType generateRandNormal(dataType mean, dataType stddev) {
    dataType u1 = (dataType)((rand() + 1.0) / (RAND_MAX + 1.0));
    dataType u2 = (dataType)((rand() + 1.0) / (RAND_MAX + 1.0));
    return mean + stddev * sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}
