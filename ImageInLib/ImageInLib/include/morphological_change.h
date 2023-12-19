#ifdef __cplusplus
extern "C" {
#endif
#pragma once
#include "common_functions.h"
#include <stdbool.h>
#include "../src/filter_params.h"

	//Erosion removes floating pixels and thin lines so that only substantive objects remain
	bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool erosion3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	//Dilation makes objects more visible and fills in small holes in objects.
	bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	bool dilatation3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background);

	/// <summary>
	/// Perform 2D image erosion. Turn foreground pixels with background pixel to background pixels
	/// </summary>
	/// <param name="imageDataPtr">Binary input image</param>
	/// <param name="length">image Length</param>
	/// <param name="width">image Width</param>
	/// <param name="object">forground pixel</param>
	/// <param name="background">background pixel</param>
	/// <returns></returns>
	bool erosion2D(dataType* imageDataPtr, const size_t length, const size_t width, dataType object, dataType background);

	bool dilatation2D(dataType* imageDataPtr, const size_t length, const size_t width, dataType object, dataType background);
#ifdef __cplusplus
}
#endif
