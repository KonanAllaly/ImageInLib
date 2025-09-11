#ifdef __cplusplus
extern "C" {
#endif

#include "../src/filter_params.h"
#pragma once

	typedef struct
	{
		dataType** east; // East coefficient pointer
		dataType** west; // West coefficient pointer
		dataType** north; // North coefficient pointer
		dataType** south; // South coefficient pointer
		dataType** top; // Top coefficient pointer
		dataType** bottom; // Bottom coefficient pointer
	} Pointers_Neighbours;

	dataType getMinInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t x, const size_t y, const size_t z);

	dataType getMaxInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t x, const size_t y, const size_t z);

	bool normOfGradientReducedDiamondCells(Image_Data inputImageData, Pointers_Neighbours vGrad);

#ifdef __cplusplus
}
#endif
