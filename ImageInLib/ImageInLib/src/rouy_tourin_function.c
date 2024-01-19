/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_initialization.h"
#include "distance_function.h"
#include "common_functions.h"

dataType min0(dataType x, dataType y) {
	if (y - x > 0) 
		return pow(x - y, 2);
	else 
		return 0;
}

bool rouyTourinFunction_3D(dataType ** distance3DPtr, dataType ** image3DPtr, dataType tolerance,
	const size_t xDim, const size_t yDim, const size_t zDim, dataType tau, const dataType h)
{
	if (distance3DPtr == NULL || image3DPtr == NULL)
		return false;

	size_t k, i, i_n, iter;//i and k are loop counters. i_n is also a loop counter given by 

	size_t rowDim_ext = xDim + 2;
	size_t columnDim_ext = yDim + 2;
	size_t sliceDim_ext = zDim + 2;

	size_t origRowDim = xDim;
	size_t origColumnDim = yDim;

	size_t sliceBound = (rowDim_ext - 1)* columnDim_ext;
	size_t k_o, i_o, j_o, j;
	const size_t noIteration = (size_t)(sqrt((double)(xDim * xDim + yDim * yDim + zDim * zDim)) / tau);

	const size_t dim2D_ext = rowDim_ext * columnDim_ext;
	dataType tauu = tau / h, mass = 10.0;
	dataType current_dist;

	dataType ** zPtrtemp = (dataType **)malloc(sizeof(dataType*) * sliceDim_ext);

	//checks if the memory was allocated
	if (zPtrtemp == NULL)
		return false;

	for (i = 0; i < sliceDim_ext; i++)
	{
		zPtrtemp[i] = (dataType *)malloc(sizeof(dataType) * dim2D_ext);

		//checks if the memory was allocated
		if (zPtrtemp[i] == NULL)
			return false;
	}

	// filling the 3D array with number=0
	initialize3dArrayD(distance3DPtr, xDim, yDim, zDim, 0);

	iter = 0;

	// implementation of Rouy Tourin scheme in 3D
	while (mass > tolerance && iter < noIteration)
	{
		mass = 0;
		iter++;

		copyDataToExtendedArea(distance3DPtr, zPtrtemp, zDim, origRowDim, origColumnDim);
		//reflection of zPtrtemp 
		reflection3D(zPtrtemp, sliceDim_ext, rowDim_ext, columnDim_ext);

		//computation of 3D distance
		k_o = 0;
		for (k = 1; k <= zDim; k++, k_o++)
		{
			i_o = 0;
			for (i = 1; i <= xDim; i++, i_o++)//row loop
			{
				j_o = 0;
				for (j = 1; j <= yDim; j++, j_o++)// column loop i_n = i(row) + j(column) * rowDim
				{
					i_n = i + j * rowDim_ext;
					if (image3DPtr[k_o][(i_o + j_o * xDim)])
					{
						current_dist = zPtrtemp[k][i_n];

						distance3DPtr[k_o][(i_o + j_o * xDim)] = zPtrtemp[k][i_n] + tau - tauu *
							((dataType)sqrt(
								max(pow(min(zPtrtemp[k][i_n - 1] - current_dist, 0), 2),
									pow(min(zPtrtemp[k][i_n + 1] - current_dist, 0), 2)) +
								max(pow(min(zPtrtemp[k][i_n - rowDim_ext] - current_dist, 0), 2),
									pow(min(zPtrtemp[k][i_n + rowDim_ext] - current_dist, 0), 2)) +
								max(pow(min(zPtrtemp[k - 1][i_n] - current_dist, 0), 2),
									pow(min(zPtrtemp[k + 1][i_n] - current_dist, 0), 2))));

						mass += (dataType)pow(distance3DPtr[k_o][(i_o + j_o * xDim)] - zPtrtemp[k][i_n], 2);
					}
				}
			}
		}

		mass = (dataType)sqrt(mass);
	}

	for (i = 0; i < sliceDim_ext; i++)
		free(zPtrtemp[i]);

	free(zPtrtemp);
	return true;
}

bool rouyTourin2D(dataType* imageDataPtr, dataType* distancePtr, const size_t length, const size_t width, dataType tolerance, dataType tau, dataType h) {

	if (imageDataPtr == NULL || distancePtr == NULL)
		return false;

	dataType mass = 10.0;
	size_t i, j, xd, count = 0, max_iter = 50000;
	size_t i_ext, j_ext, xd_ext;

	size_t length_ext = length + 2, width_ext = width + 2;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);

	for (i = 0; i < length_ext * width_ext; i++){
		extendedArray[i] = 0.0;
	}
	
	//copyDataTo2dExtendedArea(imageDataPtr, extendedArray, length, width);
	//reflection2D(extendedArray, length_ext, width_ext);

	while (mass > tolerance && count < max_iter) {
		
		count++;
		mass = 0.0;
		copyDataTo2dExtendedArea(distancePtr, extendedArray, length, width);
		reflection2D(extendedArray, length_ext, width_ext);

		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				
				xd = x_new(i, j, length);
				xd_ext = x_new(i_ext, j_ext, length_ext);

				if (imageDataPtr[xd] == 1.0) {
					distancePtr[xd] = extendedArray[xd_ext] + tau - (tau / h) * sqrt( max(min(extendedArray[x_new(i_ext, j_ext - 1, length_ext)] - extendedArray[x_new(i_ext, j_ext, length_ext)], 0), min(extendedArray[x_new(i_ext, j_ext + 1, length_ext)] - extendedArray[x_new(i_ext, j_ext, length_ext)], 0.0) ) +
						                                                   max(min(extendedArray[x_new(i_ext - 1, j_ext, length_ext)] - extendedArray[x_new(i_ext, j_ext, length_ext)], 0.0), min(extendedArray[x_new(i_ext + 1, j_ext, length_ext)] - extendedArray[x_new(i_ext, j_ext, length_ext)], 0.0)) );
				}
				else {
					distancePtr[xd] = 0.0;
				}
				mass += pow(extendedArray[xd_ext] - distancePtr[xd], 2);
			}
		}
		mass = sqrt(mass);

	}

	free(extendedArray);
	return true;
}

bool bruteForceDistanceMap(dataType* imageDataPtr, dataType* distancePtr, const size_t length, const size_t width) {
	if (imageDataPtr == NULL || distancePtr == NULL)
		return false;

	size_t i, j, xd, in, jn, xdn;
	size_t min_distance;
	Point2D p1, p2;
	double dist = 0.0;

	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {

			xd = x_new(i, j, length);
			p1.x = i; 
			p1.y = j;

			min_distance = 100000000.0;

			for (in = 0; in < length; in++) {
				for (jn = 0; jn < width; jn++) {

					xdn = x_new(in, jn, length);
					p2.x = in; 
					p2.y = jn;

					if (imageDataPtr[xdn] == 1.0) {
						dist = getPoint2DDistance(p1, p2);
						if (min_distance > dist) {
							min_distance = dist;
						}
					}

				}
			}

			distancePtr[xd] = min_distance;
		}
	}

	return true;
}

