#include "fast_marching_front_propagation_2d.h"

void fastMarchingFrontPropagation(void* inputImageData, void* actionPtr, void* potentialPtr, void* endPoints, const PropagationType pType)
{
	switch (pType)
	{
	case PARTIAL_FRONT_PROPAGATION:
	{
		Image_Data2D inputImage = *(Image_Data2D*)inputImageData;
		dataType* action = (dataType*)actionPtr;
		dataType* potential = (dataType*)potentialPtr;
		Point2D* endPoint = (Point2D*)endPoints;
		partialFrontPropagation2D(inputImage, action, potential, endPoint);
		break;
	}
	default:
		break;
	}
}

bool partialFrontPropagation2D(Image_Data2D inputImage, dataType* action, dataType* potential, Point2D* endPoint)
{
	if(inputImage.imageDataPtr == NULL || action == NULL || potential == NULL || endPoint == NULL)
	{
		return false;
	}

	return true;
}