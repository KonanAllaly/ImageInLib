#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include "common_functions.h"
#pragma once

	typedef struct {
		dataType min, max;
	} extremeValues;

	extremeValues findExtremeValues(dataType* imageDataPtr, const size_t height, const size_t width, const BoundingBox2D box);

	bool computeImageHistogram(dataType* imageDataPtr, const size_t height, const size_t width, const size_t bins, size_t* histogram);

	bool histogramEqualization(dataType* imageDataPtr, const size_t height, const size_t width);

	bool CLACHE(dataType* imageDataPtr, const size_t height, const size_t width, const size_t nbDivision);

#ifdef __cplusplus
}
#endif