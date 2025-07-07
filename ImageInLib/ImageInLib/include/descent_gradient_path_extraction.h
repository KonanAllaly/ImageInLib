#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef DESCENT_GRADIENT_PATH_EXTRACTION
#define DESCENT_GRADIENT_PATH_EXTRACTION

#include "common_functions.h"

	typedef struct {
		dataType tau; // used in the descent gradient
		size_t max_iteration; // maximal iteration to stop descent gradient
		dataType tolerance; //minimal distance to stop
	} Path_Parameters;

	bool shortestPath2d(Image_Data2D actionMapStr, Point2D* seedPoints, Path_Parameters parameters, unsigned char* pathPtr);


#endif // !DESCENT_GRADIENT_PATH_EXTRACTION

#ifdef __cplusplus
}
#endif