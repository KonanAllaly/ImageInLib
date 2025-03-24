#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "common_functions.h"
#include "heat_equation.h"
#include "filter_params.h"


	//3D function
	bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters, unsigned char* outputPathPtr);

	bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** edgeGradientPtr, Gradient_Pointers VPtrs,
		Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters);

	bool generalizedGaussSeidelCoefficients(Image_Data segmentationData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs,
		Gradient_Pointers VPtrs, Segmentation_Parameters segParameters);

	bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Image_Data segmentationData,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs);

#ifdef __cplusplus
}
#endif