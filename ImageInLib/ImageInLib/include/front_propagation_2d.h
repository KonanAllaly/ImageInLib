#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef FRONT_PROPAGATION_2D
#define FRONT_PROPAGATION_2D

#include <stdbool.h>
#include "common_function.h"

	bool partialFrontPropagation2D(Image_Data2D inputImageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* seed);

#endif // !FRONT_PROPAGATION_2D

#ifdef __cplusplus
}
#endif