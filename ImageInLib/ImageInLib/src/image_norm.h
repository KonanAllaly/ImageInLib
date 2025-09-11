#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h"

	double norm3dDataArrayD(dataType ** dataArray3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h);

	bool l2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, unsigned char * pathPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType h,
		size_t step);

	dataType timespacel2norm3dDataArrayD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, dataType h);

	dataType l2normD(dataType ** dataArray3DPtr1, dataType ** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, double h);

	dataType l2normRectangularGrid(dataType** dataArray3DPtr1, dataType** dataArray3DPtr2, const size_t xDim, const size_t yDim, const size_t zDim, VoxelSpacing h);

#ifdef __cplusplus
}
#endif
