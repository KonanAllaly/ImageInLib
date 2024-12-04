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

	dataType getImageMeanValue(dataType* imageDataPtr, const size_t length, const size_t width);

	dataType getImageMeanValue3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height);

	bool computeImageHistogram(dataType* imageDataPtr, const size_t height, const size_t width, const size_t bins, size_t* histogram);

	size_t findHistogramPeaks(dataType* imageDataPtr, const size_t height, const size_t width, const size_t bins, const char* store_file);

	size_t histogramPeaks3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t bins, const char* store_file);

	bool histogramEqualization(dataType* imageDataPtr, const size_t height, const size_t width);

	bool CLACHE(dataType* imageDataPtr, const size_t height, const size_t width, const size_t nbDivision);

#ifdef __cplusplus
}
#endif