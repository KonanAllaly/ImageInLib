#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/data_storage.h"
#include "../src/heat_equation.h"
#include "../src/filter_params.h"
#include "../src/segmentation3D_subsurf.h"

	typedef struct
	{
		dataType* East;
		dataType* West;
		dataType* North;
		dataType* South;
	} neighPtrs;

	void initialize2dArrayWithZero(dataType* arrayPtr, const size_t height, const size_t width);

	dataType getMinInNeighborhood(dataType* imageDataPtr, const size_t height, const size_t width, const size_t i, const size_t j);

	dataType getMaxInNeighborhood(dataType* imageDataPtr, const size_t height, const size_t width, const size_t i, const size_t j);

	dataType l2norm(dataType* arrayPtr1, dataType* arrayPtr2, const size_t height, const size_t width, dataType h);

	bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width);

	bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, Point2D* center, dataType v, dataType R);

	bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width);

	bool computeNormOfGradientDiamondCells(dataType* imageDataPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h);

	bool subsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

	bool gsubsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

	bool gsubsurf_iioe(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

#ifdef __cplusplus
}
#endif