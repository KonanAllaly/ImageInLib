#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef FAST_MARCHING_FRONT_PROPAGATION
#define FAST_MARCHING_FRONT_PROPAGATION
#include "common_function.h"

	typedef enum
	{
		PARTIAL_FRONT_PROPAGATION = 1,
		FRONT_PROPAGATION = 2,
		DOUBLE_FRONT_PROPAGATION = 3,
		KEY_POINT_DETECTION = 4,
		DISTANCE_MAP = 5
	} PropagationType;

	void fastMarchingFrontPropagation(void* inputImageData, void* actionPtr, void* potentialPtr, void* endPoints, const PropagationType pType);

	//bool partialFrontPropagation2D(Image_Data2D inputImage, dataType* action, dataType* potential, Point2D* endPoint);
	//void computePotential2D(dataType* imageDataPtr, dataType* potentialPtr, const size_t length, const size_t width, Point2D* endPoints, const dataType epsilon);

#endif // !FAST_MARCHING_FRONT_PROPAGATION

#ifdef __cplusplus
}
#endif
