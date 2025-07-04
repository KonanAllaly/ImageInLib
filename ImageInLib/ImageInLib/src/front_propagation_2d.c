#include "front_propagation.h"

bool partialFrontPropagation2D(Image_Data2D inputImageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* endPoints)
{
	if(inputImageData.imageDataPtr == NULL || actionMapPtr == NULL 
		|| potentialPtr == NULL || endPoints == NULL)
	{
		return false;
	}

	return true;
}