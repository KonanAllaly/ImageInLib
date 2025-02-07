#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/filter_params.h"
#include "../src/heat_equation.h"

	typedef struct {
		size_t x, y;
		int label;
	}point2dLabelling;

	typedef struct {
		size_t x, y, z;
		//bool isThePointProcessed;
		int label;
	} point3dLabelling;

	//for 2D image using 2D arrays
	bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, int object);

	//For 3D Images
	bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, const size_t height, dataType object);

	/// <summary>
	/// count foreground pixel in binary image
	/// </summary>
	/// <param name="binaryImagePtr">input binary image</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="height">image height</param>
	/// <param name="object">foreground pixel value</param>
	/// <returns>number of foreground pixels</returns>
	int countNumberOfRegionsCells(dataType** binaryImagePtr, const size_t length, const size_t width, const size_t height, dataType object);

	//=======================

	bool labeling(dataType* imageDataPtr, int* labelArray, bool* statusArray, const size_t length, const size_t width, dataType object);

	bool regionGrowing(dataType** maskThreshold, dataType** segmentedImage, const size_t length, const size_t width, const size_t height, Point3D pSeed);

	bool lungsSegmentation(Image_Data imageData, dataType** segmentedImage, Point3D pSeed, dataType tMin, dataType tMax);

#ifdef __cplusplus
}
#endif
