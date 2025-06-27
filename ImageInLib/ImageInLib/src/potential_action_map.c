#include "potential_action_map.h"

/*
bool computePotential2D(Image_Data2D imageData, dataType* potentialPtr, Point2D* endPoints, Potential_Parameters parameters)
{
	
	const size_t length = imageData.height;
	const size_t width = imageData.width;
	const size_t dim2D = length * width;

	if (imageData.imageDataPtr == NULL || potentialPtr == NULL || endPoints == NULL)
		return false;

	dataType value_end_pt1 = imageData.imageDataPtr[x_new((size_t)endPoints[0].x, (size_t)endPoints[0].y, length)];
	dataType value_end_pt2 = imageData.imageDataPtr[x_new((size_t)endPoints[1].x, (size_t)endPoints[1].y, length)];
	dataType seedVal = (dataType)(value_end_pt1 + value_end_pt2) / 2.0);

	//Compute absolute difference from the seed value
	for (size_t i = 0; i < dim2D; i++) 
	{
		potentialPtr[i] = abs(seedVal - imageData.imageDataPtr[i]);
	}

	//find the max potential
	dataType max_potential = 0;
	for (size_t i = 0; i < dim2D; i++)
	{
		if (max_potential < potentialPtr[i]) {
			max_potential = potentialPtr[i];
		}
	}

	//Normalization
	for (size_t i = 0; i < dim2D; i++)
	{
		potentialPtr[i] = parameters.eps + potentialPtr[i] / max_potential;
	}

	return true;
}
*/