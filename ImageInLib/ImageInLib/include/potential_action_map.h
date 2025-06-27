#pragma once

#include "common_functions.h"
#include "front_propagation.h"

typedef struct {  
    dataType tau; // used in the descent gradient  
    size_t max_iteration; // maximal iteration to stop descent gradient  
    dataType tolerance; // minimal distance to stop  
} Path_Parameters;

bool computePotential2D(Image_Data2D imageData, dataType* potentialFuncPtr, Point2D* seedPoints, Potential_Parameters parameters);