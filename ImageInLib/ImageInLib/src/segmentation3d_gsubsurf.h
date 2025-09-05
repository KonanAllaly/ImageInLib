#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "segmentation3D_subsurf.h"

	dataType getMinInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k);

	dataType getMaxInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k);

	bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters, unsigned char* outputPathPtr);

	bool generalizedGaussSeidelCoefficients(Image_Data segmentationData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs,
		Gradient_Pointers VPtrs, Segmentation_Parameters segParameters);

	bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Image_Data segmentationData,
		Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs);

	bool computeNormOfGradientDiamondCell3D(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Coefficient_Pointers nGrad);

	bool generalizedSubsurf_iioe(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms);

#ifdef __cplusplus
}
#endif