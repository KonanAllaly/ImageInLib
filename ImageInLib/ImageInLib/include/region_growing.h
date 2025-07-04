#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef REGION_GROWING
#define REGION_GROWING

#include "labeling.h"

	typedef enum
	{
		RG_2D = 1,
		RG_3D,
	} dimRG;

	bool regionGrowing2D(Image_Data2D inputImageData, dataType* segment, dataType foreGroundValue);

	bool regionGrowing3D(Image_Data inputImageData, dataType** segment, dataType foreGroundValue);

	void regionGrowingSegmentation(void* pInputImageData, void* segmentationResult, void* seed, dataType foreGroundValue, const dimRG dimension);

#endif // !REGION_GROWING

#ifdef __cplusplus
}
#endif