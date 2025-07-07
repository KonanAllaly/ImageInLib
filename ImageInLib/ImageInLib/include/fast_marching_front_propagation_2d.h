#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef FAST_MARCHING_FRONT_PROPAGATION_2D
#define FAST_MARCHING_FRONT_PROPAGATION_2D

#include "common_functions.h"

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

	dataType upwindFiniteDifference2dX(dataType* action, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y);

	dataType upwindFiniteDifference2dY(dataType* action, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y);

	dataType solve2dQuadratic(dataType dx, dataType dy, dataType p, PixelSpacing h);

	void fastMarchingFrontPropagation2D(void* inputImageData, void* actionPtr, void* potential, void* endPoints, const PropagationType pType);

	bool partialFrontPropagation2D(Image_Data2D inputImage, dataType* action, dataType* potential, Point2D* endPoint);

#endif // !FAST_MARCHING_FRONT_PROPAGATION_2D

#ifdef __cplusplus
}
#endif
