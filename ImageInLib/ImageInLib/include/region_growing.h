#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef REGION_GROWING
#define REGION_GROWING

#include "labeling.h"

	bool regionGrowing2D(Image_Data2D inputImageData, dataType* segment, Point2D seed, dataType radius, dataType foreGroundValue);

	bool regionGrowing3D(Image_Data inputImageData, dataType** segment, Point3D seed, dataType radius, dataType foreGroundValue);

	void regionGrowingSegmentation(void* pInputImageData, void* segmentationResult, void* seed, const dataType radius, const dataType foreGroundValue, const pDimension dim);

#endif // !REGION_GROWING

#ifdef __cplusplus
}
#endif