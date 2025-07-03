#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef LABELING
#define LABELING

#include <stdbool.h>
#include "common_functions.h"

	typedef enum
	{
		LABEL_2D_ARRAY = 1,
		LABEL_3D_ARRAY,
	} labelingApproach;

	LinkedPoint* pushPointToList(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double x, const double y);

	void popFirstEltList(LinkedCurve* linked_curve);

	bool labeling2D(Image_Data2D inputImageData, dataType* segment, dataType foreGroundValue);
	
	bool labeling3D(Image_Data inputImageData, dataType** segment, dataType foreGroundValue);
	
	void connectedComponentsLabeling(void* pInputImageData, void* segmentationResult, dataType foreGroundValue, const labelingApproach dimension);

#endif // !LABELING

#ifdef __cplusplus
}
#endif
