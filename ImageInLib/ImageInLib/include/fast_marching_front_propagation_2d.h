#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef FAST_MARCHING_FRONT_PROPAGATION_2D
#define FAST_MARCHING_FRONT_PROPAGATION_2D

#include <stdbool.h>
#include "common_functions.h"

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching;

	typedef struct {
		dataType K; //edge detection coef
		dataType thres;//edge detector threshold
		dataType eps; //path smothing parameter
		double radius;
	} Potential_Parameters;

	typedef enum
	{
		PARTIAL_FRONT_PROPAGATION = 1,
		FRONT_PROPAGATION = 2,
		DOUBLE_FRONT_PROPAGATION = 3,
		KEY_POINT_DETECTION = 4,
		DISTANCE_MAP = 5
	} PropagationType;

	void fastMarchingFrontPropagation(void* inputImageData, void* actionPtr, void* potential, void* endPoints, const PropagationType pType);

	bool partialFrontPropagation2D(Image_Data2D inputImage, dataType* action, dataType* potential, Point2D* endPoint);

#endif // !FAST_MARCHING_FRONT_PROPAGATION_2D

#ifdef __cplusplus
}
#endif
