#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_initialization.h"
#include "distance_function.h"
#include "common_functions.h"

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
	const size_t noIteration = 1000;//(size_t)(sqrt((double)(xDim * xDim + yDim * yDim + zDim * zDim)) / tau);

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

bool RouyTourinDistanceMapRectangularGrid(Image_Data ctImageData, dataType** distancePtr, dataType tolerance, size_t max_iteration, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distancePtr == NULL)
		return false;

	size_t height = ctImageData.height;
	size_t length = ctImageData.length;
	size_t width = ctImageData.width;
	size_t i, j, k, x;

	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = length + 2;
	size_t i_ext, j_ext, k_ext, x_ext;

	dataType** previousSolution = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) {
		previousSolution[k] = (dataType*)malloc(sizeof(dataType)* length_ext * width_ext);
		if (previousSolution[k] == NULL) {
			return false;
		}
	}
	if (previousSolution == NULL) {
		return false;
	}

	double mass = 1.0;
	dataType hx = ctImageData.spacing.sx;
	dataType hy = ctImageData.spacing.sy;
	dataType hz = ctImageData.spacing.sz;
	dataType value = 0.0;

	dataType hx_2 = 1.0 / (hx * hx);
	dataType hy_2 = 1.0 / (hy * hy);
	dataType hz_2 = 1.0 / (hz * hz);

	dataType tau = hx * hy * hz / (2.0 * sqrt(hx * hx + hy * hy + hz * hz));

	size_t count_iteration = 0;

	while (mass > tolerance && count_iteration < max_iteration) {
		copyDataToExtendedArea(distancePtr, previousSolution, height, length, width);
		reflection3D(previousSolution, height_ext, length_ext, width_ext);
		count_iteration++;
		mass = 0.0;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					if (ctImageData.imageDataPtr[k][x_new(i, j, length)] == foregroundValue) {
						value = previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)];
						distancePtr[k][x_new(i, j, length)] = (dataType)(value + tau - tau * sqrt(
							  hx_2 * max(pow(min(previousSolution[k_ext][x_new(i_ext - 1, j_ext, length_ext)] - value, 0), 2),
								         pow(min(previousSolution[k_ext][x_new(i_ext + 1, j_ext, length_ext)] - value, 0), 2))
							+ hy_2 * max(pow(min(previousSolution[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - value, 0), 2), 
								         pow(min(previousSolution[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - value, 0), 2))
							+ hz_2 * max(pow(min(previousSolution[k_ext - 1][x_new(i_ext, j_ext, length_ext)] - value, 0), 2),
								         pow(min(previousSolution[k_ext + 1][x_new(i_ext, j_ext, length_ext)] - value, 0), 2))));
						//Compute the mass
						mass += pow(previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)] - distancePtr[k][x_new(i, j, length)], 2);
					}
				}
			}
		}
		mass = sqrt(mass);
	}

	for (k = 0; k < height_ext; k++) {
		free(previousSolution[k]);
	}
	free(previousSolution);

	return true;
}