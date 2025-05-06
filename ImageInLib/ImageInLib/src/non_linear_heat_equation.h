#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

// INCLUDEs
#include <stdbool.h>
#include <stddef.h>
#include "common_functions.h" // ImageData structs, Reflection functions
#include "filter_params.h"
// MACROs

// STRUCTs
// FUNCTION PROTOTYPES
/*explicitParameters
* Function To Perform Heat Explicit Scheme
*/
	bool nonLinearHeatExplicitScheme(Image_Data inputImageData, Filter_Parameters explicitParameters);
	/*
	* Function To Perform Heat Gauss-Seidel Method Implicit Scheme
	*/
	bool nonLinearHeatImplicitScheme(Image_Data inputImageData, Filter_Parameters implicitParameters);

	bool geodesicMeanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters);

	bool meanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters);

	//=================

	bool geodesicMeanCurvatureRectangularTimeStep(Image_Data inputImageData, Filtering_Parameters filterParameters);

	bool geodesicMeanCurvature2D(Image_Data2D inputImage, Filter_Parameters filtering_parameters);

#endif // !HEAT_QUATION_H

#ifdef __cplusplus
}
#endif
