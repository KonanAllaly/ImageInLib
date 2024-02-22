
#include <stdio.h>
#include <stdlib.h>
#include "thresholding.h"
#include "common_functions.h"

bool thresholding2dFunctionD(dataType * image2DPtr, const size_t xDim, const size_t yDim,
	const dataType thresholdvalue, const dataType inputbgvalue, const dataType inputfgvalue);
bool thresholding2dFunctionUC(unsigned char * image2DPtr, const size_t xDim, const size_t yDim,
	const unsigned char thresholdvalue, const unsigned char inputbgvalue, const unsigned char inputfgvalue);


bool thresholding3dFunctionUC(unsigned char ** image3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const unsigned char thresholdvalue, const unsigned char inputbgvalue,
	const unsigned char inputfgvalue)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		thresholding2dFunctionUC(image3DPtr[k], xDim, yDim, thresholdvalue, inputbgvalue, inputfgvalue);
	}

	return true;
}

bool thresholding3dFunctionD(dataType ** image3DPtr, const size_t xDim, const size_t yDim,
	const size_t zDim, const dataType thresholdvalue, const dataType inputbgvalue,
	const dataType inputfgvalue)
{
	size_t k;//loop counter for z dimension

			 //checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		thresholding2dFunctionD(image3DPtr[k], xDim, yDim, thresholdvalue, inputbgvalue, inputfgvalue);
	}

	return true;
}

//function for thresholding of 2D array with some constant value.
bool thresholding2dFunctionUC(unsigned char * image2DPtr, const size_t xDim, const size_t yDim,
	const unsigned char thresholdvalue, const unsigned char inputbgvalue, const unsigned char inputfgvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	//checks for successful allocation of memory
	if (image2DPtr == NULL)
		return false;

	// Thresholding of 2D slice with thresholdvalue
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] <= thresholdvalue)
			image2DPtr[i] = inputbgvalue;
		else
			image2DPtr[i] = inputfgvalue;
	}


	return true;
}
//function for thresholding of 2D array with some constant value.
bool thresholding2dFunctionD(dataType * image2DPtr, const size_t xDim, const size_t yDim,
	const dataType thresholdvalue, const dataType inputbgvalue, const dataType inputfgvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	//checks for successful allocation of memory
	if (image2DPtr == NULL)
		return false;

	// Thresholding of 2D slice with thresholdvalue
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] <= thresholdvalue)
			image2DPtr[i] = inputfgvalue;
		else
			image2DPtr[i] = inputbgvalue;
	}


	return true;
}

bool thresholding3dFunctionN(dataType** image3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType thres_min, dataType thres_max, dataType backGround, dataType foreGround) {
	size_t i, k;

	if (image3DPtr == NULL)
		return false;

	if (thres_min != thres_max) {
		for (k = 0; k < zDim; k++) {
			for (i = 0; i < xDim * yDim; i++) {
				if (image3DPtr[k][i] >= thres_min && image3DPtr[k][i] <= thres_max) {
					image3DPtr[k][i] = foreGround;
				}
				else {
					image3DPtr[k][i] = backGround;
				}
			}
		}
	}
	else {
		for (k = 0; k < zDim; k++) {
			for (i = 0; i < xDim * yDim; i++) {
				if (image3DPtr[k][i] <= thres_min) {
					image3DPtr[k][i] = foreGround;
				}
				else {
					image3DPtr[k][i] = backGround;
				}
			}
		}
	}
	
}

bool thresholdingOTSU(dataType** image3DPtr, const size_t length, const size_t width, const size_t height, dataType background, dataType foreground) {

	size_t i, k, dim2D = length * width;

	//Find min and max data
	dataType min_data = 1000000, max_data = -10000000;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] < min_data) {
				min_data = image3DPtr[k][i];
			}
			if (image3DPtr[k][i] > max_data) {
				max_data = image3DPtr[k][i];
			}
		}
	}

	//Total number of classes
	size_t totalClass = (size_t)(max_data - min_data + 1);

	//Declare arrays
	int* histogram = (int*)malloc(totalClass * sizeof(int));
	dataType* Proba = (dataType*)malloc(totalClass * sizeof(dataType));
	dataType* interClassVariance = (dataType*)malloc(totalClass * sizeof(dataType));
	if (histogram == NULL || Proba == NULL || interClassVariance == NULL) {
		return false;
	}

	//initialize the array
	for (i = 0; i < totalClass; i++) {
		histogram[i] = 0;
		Proba[i] = 0;
		interClassVariance[i] = 0;
	}
	
	//compute histogram
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			histogram[(size_t)image3DPtr[k][i]]++;
		}
	}

	//Number of cells
	size_t numberOfCells = 0;
	for (i = 0; i < totalClass; i++) {
		numberOfCells = numberOfCells + histogram[i];
	}

	dataType sumProba = 0;
	for (i = 0; i < totalClass; i++) {
		Proba[i] = (dataType)histogram[i] / (dataType)numberOfCells;
	}

	size_t T;
	dataType sigma = 0.0, sum_weight = 0.0;
	dataType weightClass_b, meanClass_b, varClass_b;
	dataType weightClass_f, meanClass_f, varClass_f;

	for (T = 0; T < totalClass; T++) {
		
		weightClass_b = 0;
		meanClass_b = 0;
		varClass_b = 0;
		weightClass_f = 0;
		meanClass_f = 0;
		varClass_f = 0;

		//compute class weight
		for (i = 0; i <= T; i++) {
			weightClass_b = weightClass_b + Proba[i];
		}

		for (i = T + 1; i < totalClass; i++) {
			weightClass_f = weightClass_f + Proba[i];
		}
		sum_weight = weightClass_b + weightClass_f;

		//compute class mean
		for (i = 0; i <= T; i++) {
			meanClass_b = meanClass_b + i * Proba[i];
		}
		meanClass_b = meanClass_b / weightClass_b;

		for (i = T + 1; i < totalClass; i++) {
			meanClass_f = meanClass_f + i * Proba[i];
		}
		meanClass_f = meanClass_f / weightClass_f;

		//compute class variance
		for (i = 0; i <= T; i++) {
			varClass_b = varClass_b + (i - meanClass_b) * (i - meanClass_b) * Proba[i];
		}
		varClass_b = varClass_b / weightClass_b;

		for (i = T + 1; i < totalClass; i++) {
			varClass_f = varClass_f + (i - meanClass_f) * (i - meanClass_f) * Proba[i];
		}
		varClass_f = varClass_f / weightClass_f;

		//compute inter-class variance
		sigma = weightClass_b * varClass_b + weightClass_f * varClass_f;
		interClassVariance[T] = sigma;
	}

	dataType minVariance = 1e10, optimalThresholdValue = 0.0;
	for (T = 0; T < totalClass; T++) {
		if (interClassVariance[T] != 0 && interClassVariance[T] < minVariance) {
			minVariance = interClassVariance[T];
			optimalThresholdValue = (dataType)T;
		}
	}

	printf("optimal threshold value = %f \n", optimalThresholdValue);

	//Threshold
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] < optimalThresholdValue) {
				image3DPtr[k][i] = background;
			}
			else {
				image3DPtr[k][i] = foreground;
			}
		}
	}

	free(histogram);
	free(Proba);
	free(interClassVariance);

	return true;
}

bool iterativeThreshold(dataType** image3DPtr, const size_t length, const size_t width, const size_t height, dataType background, dataType foreground) {

	size_t i, k, dim2D = length * width;

	//Find min and max data
	dataType min_data = 1000000, max_data = -10000000;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] < min_data) {
				min_data = image3DPtr[k][i];
			}
			if (image3DPtr[k][i] > max_data) {
				max_data = image3DPtr[k][i];
			}
		}
	}

	dataType optimalThresholdValue = 500; //(min_data + max_data) / 2;

	int cpt1 = 0, cpt2 = 0;
	dataType s1 = 0.0, s2 = 0.0, mean1 = 0.0, mean2 = 0.0;
	do {
		s1 = 0.0; s1 = 0.0;
		cpt2 = 0; cpt1 = 0;
		for (k = 0; k < height; k++) {
			for (i = 0; i < dim2D; i++) {
				if (image3DPtr[k][i] < optimalThresholdValue) {
					cpt1++;
					s1 = s1 + image3DPtr[k][i];
				}
				else {
					cpt2++;
					s2 = s2 + image3DPtr[k][i];
				}
			}
		}
		mean1 = s1 / cpt1;
		mean2 = s2 / cpt2;
		optimalThresholdValue = (mean1 + mean2) / 2;
	} while (mean1 == mean2);

	printf("optimal threshold value = %f \n", optimalThresholdValue);

	//Threshold
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] < optimalThresholdValue) {
				image3DPtr[k][i] = foreground;
			}
			else {
				image3DPtr[k][i] = background;
			}
		}
	}

	return true;
}

//2D function
bool thresholding2DFunction(dataType* image2DPtr, const size_t xDim, const size_t yDim, dataType thres_min, dataType thres_max) {
	if (image2DPtr == NULL)
		return false;

	size_t i, dim2D = xDim * yDim;
	if (thres_min == thres_max) {
		dataType thres = thres_min;
		for (i = 0; i < dim2D; i++) {
			if (image2DPtr[i] < thres) {
				image2DPtr[i] = 1.0;
			}
			else {
				image2DPtr[i] = 0.0;
			}
		}
	}
	else {
		for (i = 0; i < dim2D; i++) {
			if (image2DPtr[i] >= thres_min && image2DPtr[i] <= thres_max) {
				image2DPtr[i] = 1.0;
			}
			else {
				image2DPtr[i] = 0.0;
			}
		}
	}
	return true;
}