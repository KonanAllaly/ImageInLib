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

	bool normOfGradientReducedDiamondCells(Image_Data inputImageData, Pointers_Neighbours vGrad);

#ifdef __cplusplus
}
#endif
