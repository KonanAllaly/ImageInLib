#include <stdlib.h>
#include<stdio.h>
#include <stdbool.h> // Boolean function bool
#include <math.h> // Maths functions i.e. pow, sin, cos
#include "heat_equation.h"
#include "data_initialization.h"

// Local Function Prototype

// Functions for 3D images

void heatExplicitScheme(Image_Data toExplicitImage, const Filter_Parameters explicitParameters)
{
	size_t k, i, j;
	dataType hh = explicitParameters.h * explicitParameters.h;
	dataType tau = explicitParameters.timeStepSize;

	// Perform Reflection of the tempPtr
	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	// Less the borders because in the loops we add back the border p
	size_t height = toExplicitImage.height;
	size_t length = toExplicitImage.length;
	size_t width = toExplicitImage.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;

	// Create temp Image Data holder for Previous time step data - with extended boundary because of boundary condition

	dataType** tempPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++)
	{
		tempPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		if (tempPtr[k] == NULL)
			return;
	}
	if (tempPtr == NULL)
		return;

	initialize3dArrayD(tempPtr, length_ext, width_ext, height_ext, 0.0);

	size_t k_ext, j_ext, i_ext;
	//size_t sliceBound = (length_ext - 1) * width_ext;
	//size_t i_d = 0;

	copyDataToExtendedArea(toExplicitImage.imageDataPtr, tempPtr, height, length, width);
	reflection3D(tempPtr, height_ext, length_ext, width_ext);

	size_t x, x_ext, t;

	const dataType coeff = tau / hh;

	// The Explicit Scheme Evaluation
	for (t = 0; t < explicitParameters.timeStepsNum; t++) {
		
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Explicit formula
					toExplicitImage.imageDataPtr[k][x] = (dataType)((1.0 - 6.0 * coeff) * tempPtr[k_ext][x_ext]
						+ coeff * (tempPtr[k_ext][x_new(i_ext - 1, j, length_ext)]
							+ tempPtr[k_ext][x_new(i_ext + 1, j, length_ext)]
							+ tempPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ tempPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
							+ tempPtr[k_ext + 1][x_ext]
							+ tempPtr[k_ext - 1][x_ext]));
				}
			}
		}

		initialize3dArrayD(tempPtr, length_ext, width_ext, height_ext, 0);
		copyDataToExtendedArea(toExplicitImage.imageDataPtr, tempPtr, height, length, width);
		reflection3D(tempPtr, height_ext, length_ext, width_ext);

	}

	// Freeing Memory after use
	for (k = 0; k < height_ext; k++)
	{
		free(tempPtr[k]);
	}
	free(tempPtr);
}

void heatImplicitScheme(Image_Data toImplicitImage, const Filter_Parameters implicitParameters)
{
	if (toImplicitImage.imageDataPtr == NULL)
		return;

	size_t k, i, j, x, steps = implicitParameters.maxNumberOfSolverIteration, p = implicitParameters.p;
	dataType hh = implicitParameters.h * implicitParameters.h;
	dataType coeff = implicitParameters.timeStepSize / hh;
	size_t k_ext, j_ext, i_ext, x_ext;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error = 0.0, sor = 0.0;
	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	// Less the borders because in the loops we add back the border p
	size_t height = toImplicitImage.height, length = toImplicitImage.length, width = toImplicitImage.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	// Create temp Image Data holder for Previous time step data
	dataType** tempPtr = (dataType**)malloc(sizeof(dataType*) * (height_ext));
	dataType** currentPtr = (dataType**)malloc(sizeof(dataType*) * (height_ext));
	for (i = 0; i < height_ext; i++)
	{
		tempPtr[i] = malloc(sizeof(dataType) * length_ext * width_ext);
		currentPtr[i] = malloc(sizeof(dataType) * length_ext * width_ext);
		if (tempPtr[i] == NULL || currentPtr[i] == NULL)
			return;
	}
	if (tempPtr == NULL || currentPtr == NULL)
		return;

	//Initialize the arrays to avoid unwanted values
	initialize3dArrayD(tempPtr, length_ext, width_ext, height_ext, 0.0);
	initialize3dArrayD(currentPtr, length_ext, width_ext, height_ext, 0.0);
	
	//// Preparation
	////We use the manual copy because the current function doesn't work for non squared images
	//copyDataToExtendedArea(toImplicitImage.imageDataPtr, tempPtr, height, length, width);
	//copyDataToExtendedArea(toImplicitImage.imageDataPtr, currentPtr, height, length, width);
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				tempPtr[k_ext][x_ext] = toImplicitImage.imageDataPtr[k][x];
				currentPtr[k_ext][x_ext] = toImplicitImage.imageDataPtr[k][x];
			}
		}
	}
	
	//Perform Reflection of the tempPtr
	reflection3D(tempPtr, height_ext, length_ext, width_ext);
	reflection3D(currentPtr, height_ext, length_ext, width_ext);
	
	// The Gauss-Seidel Implicit Scheme
	size_t z = 0;
	for(size_t t = 0; t < implicitParameters.timeStepsNum; t++){
		
		z = 0; // Steps counter
		
		do
		{
			z = z + 1;
			// Gauss-Seidel Method
			for (k = 0, k_ext = 1; k < height; k++, k_ext++){
				for (i = 0, i_ext = 1; i < length; i++, i_ext++){
					for (j = 0, j_ext = 1; j < width; j++, j_ext++){
						// 2D to 1D representation for i, j
						x_ext = x_new(i_ext, j_ext, length_ext);
						// Begin Gauss-Seidel Formula Evaluation
						sor = (dataType)((tempPtr[k_ext][x_ext] + coeff * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]
							+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)]
							+ currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
							+ currentPtr[k_ext + 1][x_ext] + currentPtr[k_ext - 1][x_ext])) / (1 + 6.0 * coeff));
						// Gauss-Seidel
						currentPtr[k_ext][x_ext] = currentPtr[k_ext][x_ext] + implicitParameters.omega_c * (sor - currentPtr[k_ext][x_ext]);
					}
				}
			}

			// Error Evaluation
			error = 0.0; // Initialize
			for (k = 0, k_ext = 1; k < height; k++, k_ext++){
				for (i = 0, i_ext = 1; i < length; i++, i_ext++){
					for (j = 0, j_ext = 1; j < width; j++, j_ext++){
						// 2D to 1D representation for i, j
						x_ext = x_new(i_ext, j_ext, length_ext);
						// Begin Error Calculation
						error += (dataType)(pow(currentPtr[k_ext][x_ext] * (1 + 6.0 * coeff)
							- coeff * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]
								+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] + currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
								+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] + currentPtr[k_ext + 1][x_ext]
								+ currentPtr[k_ext - 1][x_ext]) - tempPtr[k_ext][x_ext], 2));
					}
				}
			}

		} while (error > implicitParameters.tolerance && z < steps);

		// Copy current to tempPtr before next time step
		copyDataToAnotherArray(currentPtr, tempPtr, height_ext, length_ext, width_ext);
	}
	printf("The number of iterations is %zd\n", z);
	printf("Error is %e\n", error);

	// Copy back to original after filtering
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				toImplicitImage.imageDataPtr[k][x_new(i, j, length)] = currentPtr[k_ext][x_new(i_ext, j_ext, length_ext)];
			}
		}
	}
	
	// Freeing Memory after use
	for (i = 0; i < height_ext; i++)
	{
		free(tempPtr[i]);
		free(currentPtr[i]);
	}
	free(tempPtr);
	free(currentPtr);
}

// Functions for 2D images

void heat2dExplicitScheme(Image_Data2D imageData, const Filter_Parameters explicitParameters)
{
	if (imageData.imageDataPtr == NULL)
		return;

	size_t i, j, n, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = explicitParameters.timeStepSize, hh = explicitParameters.h;
	dataType coef_tau = tau / hh;

	dataType* temporaryPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (temporaryPtr == NULL)
		return;

	initialize2dArrayD(temporaryPtr, height_ext, width_ext, 0.0);

	copyDataTo2dExtendedArea(imageData.imageDataPtr, temporaryPtr, height, width);
	reflection2D(temporaryPtr, height_ext, width_ext);

	for (n = 0; n < explicitParameters.timeStepsNum; n++) {

		for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				imageData.imageDataPtr[x_new(i, j, height)] = (1 - 4 * coef_tau) * temporaryPtr[x_new(i_ext, j_ext, height_ext)] + coef_tau * (temporaryPtr[x_new(i_ext - 1, j_ext, height_ext)] +
					temporaryPtr[x_new(i_ext + 1, j_ext, height_ext)] + temporaryPtr[x_new(i_ext, j_ext - 1, height_ext)] + temporaryPtr[x_new(i_ext, j_ext + 1, height_ext)]);
			}
		}

		initialize2dArrayD(temporaryPtr, height_ext, width_ext, 0.0);
		copyDataTo2dExtendedArea(imageData.imageDataPtr, temporaryPtr, height, width);
		reflection2D(temporaryPtr, height_ext, width_ext);
	}

	free(temporaryPtr);
}

void heatImplicit2dScheme(Image_Data2D imageData, const Filter_Parameters implicitParameters)
{
	size_t i, j, n, i_ext, j_ext, currentIndx;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = implicitParameters.timeStepSize, hh = implicitParameters.h * implicitParameters.h;
	dataType tol = implicitParameters.tolerance, omega = implicitParameters.omega_c;
	dataType coeff = tau / hh;
	size_t maxIteration = implicitParameters.maxNumberOfSolverIteration;

	dataType* previous_solution = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* gauss_seidel_solution = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	initialize2dArrayD(previous_solution, height_ext, width_ext, 0.0);
	initialize2dArrayD(gauss_seidel_solution, height_ext, width_ext, 0.0);

	copyDataTo2dExtendedArea(imageData.imageDataPtr, previous_solution, height, width);
	copyDataTo2dExtendedArea(imageData.imageDataPtr, gauss_seidel_solution, height, width);
	reflection2D(previous_solution, height_ext, width_ext);
	reflection2D(gauss_seidel_solution, height_ext, width_ext);

	size_t cpt = 0;
	dataType error = 0.0, gauss_seidel_coef = 0.0;
	for (n = 0; n < implicitParameters.timeStepsNum; n++) {

		cpt = 0;
		do {
			cpt = cpt + 1;

			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					currentIndx = x_new(i_ext, j_ext, height_ext);
					gauss_seidel_coef = (dataType)((previous_solution[x_new(i_ext, j_ext, height_ext)] + coeff * (gauss_seidel_solution[x_new(i_ext - 1, j_ext, height_ext)] +
						gauss_seidel_solution[x_new(i_ext + 1, j_ext, height_ext)] + gauss_seidel_solution[x_new(i_ext, j_ext - 1, height_ext)] + gauss_seidel_solution[x_new(i_ext, j_ext + 1, height_ext)])) /
						(1 + 4 * coeff));
					gauss_seidel_solution[currentIndx] = gauss_seidel_solution[currentIndx] + omega * (gauss_seidel_coef - gauss_seidel_solution[currentIndx]);
				}
			}

			error = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					currentIndx = x_new(i_ext, j_ext, height_ext);
					error += (dataType)(pow(gauss_seidel_solution[currentIndx] * (1 + 4.0 * coeff) - coeff * (gauss_seidel_solution[x_new(i_ext - 1, j_ext, height_ext)] + gauss_seidel_solution[x_new(i_ext + 1, j_ext, height_ext)] +
						gauss_seidel_solution[x_new(i_ext, j_ext - 1, height_ext)] + gauss_seidel_solution[x_new(i_ext, j_ext + 1, height_ext)]) - previous_solution[currentIndx], 2));
				}
			}

		} while (cpt < maxIteration && error > tol);

		copyDataToAnother2dArray(gauss_seidel_solution, previous_solution, height_ext, width_ext);
	}

	//Copy back
	copyDataTo2dReducedArea(imageData.imageDataPtr, gauss_seidel_solution, height, width);
	
	free(previous_solution);
	free(gauss_seidel_solution);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////

void heatImplicitRectangularScheme(Image_Data toImplicitImage, const Filter_Parameters implicitParameters)
{
	if (toImplicitImage.imageDataPtr == NULL)
		return;

	size_t k, i, j, x, steps = implicitParameters.maxNumberOfSolverIteration;
	dataType hx = toImplicitImage.spacing.sx;
	dataType hy = toImplicitImage.spacing.sy;
	dataType hz = toImplicitImage.spacing.sz;
	dataType tau = implicitParameters.timeStepSize;
	size_t k_ext, j_ext, i_ext, x_ext;

	dataType hx_2 = hx * hx;
	dataType hy_2 = hy * hy;
	dataType hz_2 = hz * hz;
	dataType coeff = 2.0 * tau * (1.0 / hx_2 + 1.0 / hy_2 + 1.0 / hz_2);

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error = 0.0, sor = 0.0;
	size_t height = toImplicitImage.height, length = toImplicitImage.length, width = toImplicitImage.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	// Create temp Image Data holder for Previous time step data
	dataType** tempPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** currentPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (i = 0; i < height_ext; i++)
	{
		tempPtr[i] = malloc(sizeof(dataType) * length_ext * width_ext);
		currentPtr[i] = malloc(sizeof(dataType) * length_ext * width_ext);
		if (tempPtr[i] == NULL || currentPtr[i] == NULL)
			return;
	}
	if (tempPtr == NULL || currentPtr == NULL)
		return;

	//Initialize the arrays to avoid unwanted values
	initialize3dArrayD(tempPtr, length_ext, width_ext, height_ext, 0.0);
	initialize3dArrayD(currentPtr, length_ext, width_ext, height_ext, 0.0);

	//// Preparation
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				tempPtr[k_ext][x_ext] = toImplicitImage.imageDataPtr[k][x];
				currentPtr[k_ext][x_ext] = toImplicitImage.imageDataPtr[k][x];
			}
		}
	}

	//Perform Reflection of the tempPtr
	reflection3D(tempPtr, height_ext, length_ext, width_ext);
	reflection3D(currentPtr, height_ext, length_ext, width_ext);

	// The Gauss-Seidel Implicit Scheme
	size_t z = 0;
	for (size_t t = 0; t < implicitParameters.timeStepsNum; t++) {

		z = 0; // Steps counter

		do
		{
			z = z + 1;
			// Gauss-Seidel Method
			for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
				for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
					for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
						x_ext = x_new(i_ext, j_ext, length_ext);
						// Begin Gauss-Seidel Formula Evaluation
						sor = (dataType)((tempPtr[k_ext][x_ext] + tau * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)] / hx_2
							+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] / hx_2
							+ currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] / hy_2
							+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] / hy_2
							+ currentPtr[k_ext + 1][x_ext] / hz_2 + currentPtr[k_ext - 1][x_ext] / hz_2)) / (1 + coeff));
						// Gauss-Seidel
						currentPtr[k_ext][x_ext] = currentPtr[k_ext][x_ext] + implicitParameters.omega_c * (sor - currentPtr[k_ext][x_ext]);
					}
				}
			}

			// Error Evaluation
			error = 0.0; // Initialize
			for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
				for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
					for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
						// 2D to 1D representation for i, j
						x_ext = x_new(i_ext, j_ext, length_ext);
						// Begin Error Calculation
						error += (dataType)pow(currentPtr[k_ext][x_ext] * (1 + coeff)
							- tau * (currentPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)] / hx_2
								+ currentPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] / hx_2 + currentPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] / hy_2
								+ currentPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] / hy_2 + currentPtr[k_ext + 1][x_ext] / hz_2
								+ currentPtr[k_ext - 1][x_ext] / hz_2) - tempPtr[k_ext][x_ext], 2);
					}
				}
			}

		} while (error > implicitParameters.tolerance && z < steps);

		// Copy current to tempPtr before next time step
		copyDataToAnotherArray(currentPtr, tempPtr, height_ext, length_ext, width_ext);
	}
	printf("The number of iterations is %zd\n", z);
	printf("Error is %e\n", error);

	// Copy back to original after filtering
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				toImplicitImage.imageDataPtr[k][x_new(i, j, length)] = currentPtr[k_ext][x_new(i_ext, j_ext, length_ext)];
			}
		}
	}

	// Freeing Memory after use
	for (i = 0; i < height_ext; i++)
	{
		free(tempPtr[i]);
		free(currentPtr[i]);
	}
	free(tempPtr);
	free(currentPtr);
}