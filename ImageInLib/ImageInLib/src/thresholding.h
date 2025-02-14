#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"

	//3D function
	/* image3DPtr is pointer to the 3D image, xDim, yDim and zDim are dimensions of x, y, and z respectively.
	* inputbgvalue is the value the voxel value (image3DPtr[k][i]) will be set to if the voxel value is
	less or equal to the threshold value
	* inputfgvalue is the value the voxel value (image3DPtr[k][i]) will be set to if the voxel value is
	greater than the threshold value
	*/
	bool thresholding3dFunctionUC(unsigned char ** image3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, const unsigned char thresholdvalue, const unsigned char inputbgvalue,
		const unsigned char inputfgvalue);

	bool thresholding3dFunctionD(dataType ** image3DPtr, const size_t xDim, const size_t yDim,
		const size_t zDim, const dataType thresholdvalue, const dataType inputbgvalue,
		const dataType inputfgvalue);

	/*
	* The function compute the histogram of a 3D given image
	* histoPtr contains the frequence of appearence of each bin members
	* image3DPtr : input image data
	* xDim, yDim, zDim : image dimensions
	* binCount : number of bins (to divide the intensity range)
	*/
	bool computeHistogram(dataType** image3DPtr, size_t* histoPtr, const size_t xDim, const size_t yDim, const size_t zDim, const size_t binCount);

	bool thresholding3dFunctionN(dataType** image3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, dataType backGround, dataType forGround);

	bool thresholdingOTSU(dataType ** image3DPtr, const size_t length, const size_t width, const size_t height, dataType background, dataType forground);

	bool iterativeThreshold(dataType** image3DPtr, const size_t length, const size_t width, const size_t height, dataType initial_threshold, dataType background, dataType foreground);

	//2D function
	bool thresholding2DFunction(dataType* image2DPtr, const size_t xDim, const size_t yDim, dataType thres_min, dataType thres_max);

	bool thresholdingOTSU2D(dataType* image2DPtr, const size_t length, const size_t width, dataType background, dataType foreground);

	dataType localOTSU2D(dataType* image2DPtr, const size_t length, const size_t width, Point2D seed, double radius, dataType background, dataType foreground);

#ifdef __cplusplus
}
#endif
