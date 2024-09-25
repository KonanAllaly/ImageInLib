#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"
#include "segmentation.h"

	bool lagrangeanExplicit3DCurveSegmentation(Image_Data inputImage3D,
		const Lagrangean3DSegmentationParameters* pSegmentationParams, unsigned char* pOutputPathPtr, Curve3D* resultSegmentation);

#ifdef __cplusplus
}
#endif
