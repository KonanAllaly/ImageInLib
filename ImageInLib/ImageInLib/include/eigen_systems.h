#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdlib.h>
#include "common_functions.h"

bool getHessianMatrix3D(dataType** imageDataPtr, dataType** hessianMatrix, const size_t length, const size_t width, const size_t height, const size_t ind_x, const size_t ind_y, const size_t ind_z, const VoxelSpacing fVolume);

bool eigenSystems(dataType** a, size_t n, dataType d[], dataType** v, size_t* nrot);

bool sortEigenValues(dataType* eigen_values, const size_t dim);

bool generateGaussianKernel(dataType** kernel, const size_t kernelSize, const dataType sigma);

bool gaussianFiltering(dataType** imageDataPtr, dataType** filteredImageDataPtr, const size_t length, const size_t width, const size_t height, const dataType sigma);

bool multiscaleFiltering(Image_Data inputImageData, dataType** vesselnessImage, dataType sigma);

//============================================

void reflection2DB(dataType* toReflectImage, size_t imageLength, size_t imageWidth, size_t p);

bool generateGaussianKernel2D(dataType* kernel, const size_t kernelSize, const dataType sigma);

bool gaussianSmoothing2D(dataType* imageDataPtr, dataType* smoothImageDataPtr, const size_t length, const size_t width, const dataType sigma);

#ifdef __cplusplus
}

#endif
