#ifdef __cplusplus
extern "C" {
#endif

	/*
	* Author: Markjoe Olunna UBA
	* Purpose: ImageInLife project - 4D Image Segmentation Methods
	* Language:  C
	*/
#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "vtk_params.h"

	bool load3dDataArrayD(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr);

	bool load3dDataArrayRAW(dataType ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, LoadDataType dType);

	bool load3dDataArrayVTK(unsigned char ** imageDataPtr, const size_t imageLength, const size_t imageWidth,
		const size_t imageHeight, unsigned char * pathPtr, VTK_Header_Lines * lines);

	//==================================
	//Load 2D .pgm (ascii) image
	bool load2dPGM(dataType* imageDataPtr, const size_t xDim, const size_t yDim, const char* pathPtr);

	//==================================

	bool load2dArrayRAW(dataType* imageDataPtr, const size_t length, const size_t width, const char* pathPtr, LoadDataType dType);

	//==================================

	bool loadListof3dPoints(Image_Data image, Curve3D* pCurve, const char* filePath);

	bool changeToRealworldCord(Image_Data image, Curve3D* pCurve);

#ifdef __cplusplus
}
#endif