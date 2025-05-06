#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "setting_boundary_values.h"
#include "common_functions.h"
#include "filter_params.h"
#include "eigen_systems.h"

// Local Function Prototype

bool geodesicMeanCurvatureTimeStep(Image_Data inputImageData, Filter_Parameters filterParameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL)
		return false;

	size_t k, i, j;
	dataType hh = filterParameters.h * filterParameters.h;
	dataType tau = filterParameters.timeStepSize;
	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error, gauss_seidel;

	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t k_ext, j_ext, i_ext;
	dataType ux, uy, uz, orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	size_t x; //x = x_new(i, j, length);
	size_t x_ext; //x_ext = x_new(i_ext, j_ext, length_ext);
	size_t z; // Steps counter

	const dataType coef_tauh = tau / hh;
	dataType  u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType  orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType voxel_coef, average_face_coef;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;

	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType** prevSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	// Create temporary Image Data holder for Current time step data - with extended boundary because of boundary condition
	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	/* Create tempporary Image Data holder for calculation of diffusion coefficients on presmoothed image
	- with extended boundary because of boundary condition*/
	dataType** presmoothed_coefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	//checks if the memory was allocated
	if (prevSolPtr == NULL || gauss_seidelPtr == NULL || presmoothed_coefPtr == NULL)
		return false;

	for (k = 0; k < height_ext; k++)
	{
		presmoothed_coefPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		gauss_seidelPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		prevSolPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		//checks if the memory was allocated
		if (presmoothed_coefPtr[k] == NULL || gauss_seidelPtr[k] == NULL || prevSolPtr[k] == NULL)
			return false;
	}

	/* Create tempporary Image Data holder for diffusion coefficients
	- with extended boundary because of boundary condition*/
	dataType** presmoot_e_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_w_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_n_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_s_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_t_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_b_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (presmoot_e_coefPtr == NULL || presmoot_w_coefPtr == NULL || presmoot_n_coefPtr == NULL ||
		presmoot_s_coefPtr == NULL || presmoot_t_coefPtr == NULL || presmoot_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		presmoot_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (presmoot_e_coefPtr[k] == NULL || presmoot_w_coefPtr[k] == NULL || presmoot_n_coefPtr[k] == NULL ||
			presmoot_s_coefPtr[k] == NULL || presmoot_t_coefPtr[k] == NULL || presmoot_b_coefPtr[k] == NULL)
			return false;
	}

	dataType** orig_e_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_w_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_n_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_s_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_t_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_b_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (orig_e_coefPtr == NULL || orig_w_coefPtr == NULL || orig_n_coefPtr == NULL || orig_s_coefPtr == NULL ||
		orig_t_coefPtr == NULL || orig_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		orig_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (orig_e_coefPtr[k] == NULL || orig_w_coefPtr[k] == NULL || orig_n_coefPtr[k] == NULL || orig_s_coefPtr[k] == NULL
			|| orig_t_coefPtr[k] == NULL || orig_b_coefPtr[k] == NULL)
			return false;
	}

	dataType** coefPtr_e = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_w = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_n = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_s = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_t = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_b = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (coefPtr_e == NULL || coefPtr_w == NULL || coefPtr_n == NULL || coefPtr_s == NULL || coefPtr_t == NULL ||
		coefPtr_b == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		coefPtr_e[k] = malloc(sizeof(float) * length * width);
		coefPtr_w[k] = malloc(sizeof(float) * length * width);
		coefPtr_n[k] = malloc(sizeof(float) * length * width);
		coefPtr_s[k] = malloc(sizeof(float) * length * width);
		coefPtr_t[k] = malloc(sizeof(float) * length * width);
		coefPtr_b[k] = malloc(sizeof(float) * length * width);
		//checks if the memory was allocated
		if (coefPtr_e[k] == NULL || coefPtr_w[k] == NULL || coefPtr_n[k] == NULL || coefPtr_s[k] == NULL ||
			coefPtr_t[k] == NULL || coefPtr_b[k] == NULL)
			return false;
	}

	presmoothingData.imageDataPtr = presmoothed_coefPtr;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.imageDataPtr, prevSolPtr, height, length, width);
	//copy data to extended area which will be used in each Gauss Seidel iteration
	copyDataToExtendedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);
	//copy data to extended area which will be used for calculation of diffusion coefficients
	copyDataToExtendedArea(inputImageData.imageDataPtr, presmoothed_coefPtr, height, length, width);

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(presmoothed_coefPtr, height_ext, length_ext, width_ext);
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);
	reflection3D(gauss_seidelPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	//heatExplicitScheme(presmoothingData, filterParameters);
	heatImplicitScheme(presmoothingData, filterParameters);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = presmoothed_coefPtr[k_ext][x_ext];
				uN = presmoothed_coefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = presmoothed_coefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = presmoothed_coefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				uW = presmoothed_coefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				uNW = presmoothed_coefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = presmoothed_coefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = presmoothed_coefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = presmoothed_coefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = presmoothed_coefPtr[kminus1][x_ext];
				TuN = presmoothed_coefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = presmoothed_coefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = presmoothed_coefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				TuW = presmoothed_coefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				TuNW = presmoothed_coefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = presmoothed_coefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = presmoothed_coefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = presmoothed_coefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = presmoothed_coefPtr[kplus1][x_ext];
				BuN = presmoothed_coefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = presmoothed_coefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = presmoothed_coefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				BuW = presmoothed_coefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				BuNW = presmoothed_coefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = presmoothed_coefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = presmoothed_coefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = presmoothed_coefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//values of voxels in the extended data container for the original image
				orig_u = prevSolPtr[k_ext][x_ext];
				orig_uN = prevSolPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				orig_uS = prevSolPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				orig_uE = prevSolPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_uW = prevSolPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_uNW = prevSolPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				orig_uNE = prevSolPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				orig_uSE = prevSolPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				orig_uSW = prevSolPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				orig_Tu = prevSolPtr[kminus1][x_ext];
				orig_TuN = prevSolPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				orig_TuS = prevSolPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				orig_TuE = prevSolPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_TuW = prevSolPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_TuNW = prevSolPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				orig_TuNE = prevSolPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				orig_TuSE = prevSolPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				orig_TuSW = prevSolPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				orig_Bu = prevSolPtr[kplus1][x_ext];
				orig_BuN = prevSolPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				orig_BuS = prevSolPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				orig_BuE = prevSolPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_BuW = prevSolPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_BuNW = prevSolPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				orig_BuNE = prevSolPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				orig_BuSE = prevSolPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				orig_BuSW = prevSolPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / filterParameters.h;
				uy = (dataType)(((uN + uNE) - (uS + uSE)) / (4.0 * filterParameters.h));
				uz = (dataType)(((Tu + TuE) - (Bu + BuE)) / (4.0 * filterParameters.h));
				presmoot_e_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				// Calculation of coefficients in west direction
				ux = (uW - u) / filterParameters.h;
				uy = (dataType)(((uNW + uN) - (uSW + uS)) / (4.0 * filterParameters.h));
				uz = (dataType)(((TuW + Tu) - (BuW + Bu)) / (4.0 * filterParameters.h));
				presmoot_w_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				// Calculation of coefficients in north direction
				ux = (dataType)(((uNE + uE) - (uNW + uW)) / (4.0 * filterParameters.h));
				uy = (uN - u) / filterParameters.h;
				uz = (dataType)(((TuN + Tu) - (BuN + Bu)) / (4.0 * filterParameters.h));
				presmoot_n_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				// Calculation of coefficients in south direction
				ux = (dataType)(((uE + uSE) - (uW + uSW))
					/ (4.0 * filterParameters.h));
				uy = (uS - u) / filterParameters.h;
				uz = (dataType)(((TuS + Tu) - (BuS + Bu))
					/ (4.0 * filterParameters.h));
				presmoot_s_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				// Calculation of coefficients in top direction
				ux = (dataType)(((TuE + uE) - (TuW + uW))
					/ (4.0 * filterParameters.h));
				uy = (dataType)(((TuN + uN) - (TuS + uS))
					/ (4.0 * filterParameters.h));
				uz = (Tu - u) / filterParameters.h;
				presmoot_t_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				// Calculation of coefficients in bottom direction
				ux = (dataType)(((BuW + uW) - (BuE + uE))
					/ (4.0 * filterParameters.h));
				uy = (dataType)(((BuN + uN) - (BuS + uS))
					/ (4.0 * filterParameters.h));
				uz = (Bu - u) / filterParameters.h;
				presmoot_b_coefPtr[k][x] = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), filterParameters.edge_detector_coefficient);

				//calculation of coefficients in the original image data
				// Calculation of coefficients in east direction
				orig_ux = (orig_uE - orig_u) / filterParameters.h;
				orig_uy = (dataType)(((orig_uN + orig_uNE) - (orig_uS + orig_uSE))
					/ (4.0 * filterParameters.h));
				orig_uz = (dataType)(((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE))
					/ (4.0 * filterParameters.h));
				orig_e_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / filterParameters.h;
				orig_uy = (dataType)(((orig_uNW + orig_uN) - (orig_uSW + orig_uS))
					/ (4.0 * filterParameters.h));
				orig_uz = (dataType)(((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu))
					/ (4.0 * filterParameters.h));
				orig_w_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = (dataType)(((orig_uNE + orig_uE) - (orig_uNW + orig_uW))
					/ (4.0 * filterParameters.h));
				orig_uy = (orig_uN - orig_u) / filterParameters.h;
				orig_uz = (dataType)(((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu))
					/ (4.0 * filterParameters.h));
				orig_n_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = (dataType)(((orig_uE + orig_uSE) - (orig_uW + orig_uSW))
					/ (4.0 * filterParameters.h));
				orig_uy = (orig_uS - orig_u) / filterParameters.h;
				orig_uz = (dataType)(((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu))
					/ (4.0 * filterParameters.h));
				orig_s_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = (dataType)(((orig_TuE + orig_uE) - (orig_TuW + orig_uW))
					/ (4.0 * filterParameters.h));
				orig_uy = (dataType)(((orig_TuN + orig_uN) - (orig_TuS + orig_uS))
					/ (4.0 * filterParameters.h));
				orig_uz = (orig_Tu - orig_u) / filterParameters.h;
				orig_t_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = (dataType)(((orig_BuW + orig_uW) - (orig_BuE + orig_uE))
					/ (4.0 * filterParameters.h));
				orig_uy = (dataType)(((orig_BuN + orig_uN) - (orig_BuS + orig_uS))
					/ (4.0 * filterParameters.h));
				orig_uz = (orig_Bu - orig_u) / filterParameters.h;
				orig_b_coefPtr[k][x] = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + filterParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e_coefPtr[k][x] + orig_w_coefPtr[k][x] + orig_n_coefPtr[k][x] + orig_s_coefPtr[k][x]
					+ orig_t_coefPtr[k][x] + orig_b_coefPtr[k][x]) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + filterParameters.eps2);

				/* evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				image at each voxel face and reciprocal of norm of gradient of image at each voxel face*/
				coefPtr_e[k][x] = (dataType)(voxel_coef * presmoot_e_coefPtr[k][x] * (1.0 / orig_e_coefPtr[k][x]));//east coefficient
				coefPtr_w[k][x] = (dataType)(voxel_coef * presmoot_w_coefPtr[k][x] * (1.0 / orig_w_coefPtr[k][x]));//west coefficient
				coefPtr_n[k][x] = (dataType)(voxel_coef * presmoot_n_coefPtr[k][x] * (1.0 / orig_n_coefPtr[k][x]));//north coefficient
				coefPtr_s[k][x] = (dataType)(voxel_coef * presmoot_s_coefPtr[k][x] * (1.0 / orig_s_coefPtr[k][x]));//south coefficient
				coefPtr_t[k][x] = (dataType)(voxel_coef * presmoot_t_coefPtr[k][x] * (1.0 / orig_t_coefPtr[k][x]));//top coefficient
				coefPtr_b[k][x] = (dataType)(voxel_coef * presmoot_b_coefPtr[k][x] * (1.0 / orig_b_coefPtr[k][x]));//bottom coefficient
			}
		}
	}

	// The Implicit Scheme Evaluation
	z = 0;
	do
	{
		z = z + 1;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Begin Gauss-Seidel Formula Evaluation
					gauss_seidel = (prevSolPtr[k_ext][x_ext] + coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coef_tauh * (coefPtr_e[k][x] + coefPtr_w[k][x] + coefPtr_n[k][x]
							+ coefPtr_s[k][x] + coefPtr_t[k][x] + coefPtr_b[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						filterParameters.omega_c*(gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
				}
			}
		}

		// Error Evaluation
		error = 0.0; // Initialize
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					error += (dataType)pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (coefPtr_e[k][x]
						+ coefPtr_w[k][x] + coefPtr_n[k][x] + coefPtr_s[k][x]
						+ coefPtr_t[k][x] + coefPtr_b[k][x]))
						- coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
							+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
							+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSolPtr[k_ext][x_ext], 2);
				}
			}
		}
	} while (error > filterParameters.tolerance && z < filterParameters.maxNumberOfSolverIteration);
	printf("The number of iterations is %zd\n", z);
	printf("Error is %e\n", error);
	//printf("Step is %zd\n", filterParameters.timeStepsNum);

	//Copy the current time step to original data holder after timeStepsNum
	copyDataToReducedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);

	// Freeing Memory after use
	for (k = 0; k < height_ext; k++)
	{
		free(presmoothed_coefPtr[k]);
		free(gauss_seidelPtr[k]);
		free(prevSolPtr[k]);
	}
	free(presmoothed_coefPtr);
	free(gauss_seidelPtr);
	free(prevSolPtr);

	// free _coefPtr pointers
	for (k = 0; k < height; k++)
	{
		free(presmoot_e_coefPtr[k]);
		free(presmoot_w_coefPtr[k]);
		free(presmoot_n_coefPtr[k]);
		free(presmoot_s_coefPtr[k]);
		free(presmoot_t_coefPtr[k]);
		free(presmoot_b_coefPtr[k]);
	}
	free(presmoot_e_coefPtr);
	free(presmoot_w_coefPtr);
	free(presmoot_n_coefPtr);
	free(presmoot_s_coefPtr);
	free(presmoot_t_coefPtr);
	free(presmoot_b_coefPtr);

	// free orig_ _coefPtr pointers
	for (k = 0; k < height; k++)
	{
		free(orig_e_coefPtr[k]);
		free(orig_w_coefPtr[k]);
		free(orig_n_coefPtr[k]);
		free(orig_s_coefPtr[k]);
		free(orig_t_coefPtr[k]);
		free(orig_b_coefPtr[k]);
	}
	free(orig_e_coefPtr);
	free(orig_w_coefPtr);
	free(orig_n_coefPtr);
	free(orig_s_coefPtr);
	free(orig_t_coefPtr);
	free(orig_b_coefPtr);

	// free coefPtr_ pointers
	for (k = 0; k < height; k++)
	{
		free(coefPtr_e[k]);
		free(coefPtr_w[k]);
		free(coefPtr_n[k]);
		free(coefPtr_s[k]);
		free(coefPtr_t[k]);
		free(coefPtr_b[k]);
	}
	free(coefPtr_e);
	free(coefPtr_w);
	free(coefPtr_n);
	free(coefPtr_s);
	free(coefPtr_t);
	free(coefPtr_b);

	return true;
}

bool geodesicMeanCurvatureRectangularTimeStep(Image_Data inputImageData, Filtering_Parameters filterParameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL)
		return false;

	size_t k, i, j;
	dataType hhh = filterParameters.hx * filterParameters.hy * filterParameters.hz;
	dataType tau = filterParameters.timeStepSize;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error, gauss_seidel;

	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t k_ext, j_ext, i_ext;

	dataType ux, uy, uz, orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	size_t x; //x = x_new(i, j, length);
	size_t x_ext; //x_ext = x_new(i_ext, j_ext, length_ext);
	size_t count_step; // Steps counter

	const dataType coef_tauh = tau / hhh;
	dataType  u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType  orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType voxel_coef, average_face_coef;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;

	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType** prevSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	// Create temporary Image Data holder for Current time step data - with extended boundary because of boundary condition
	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	//Create tempporary Image Data holder for calculation of diffusion coefficients on presmoothed image
	//- with extended boundary because of boundary condition

	dataType** presmoothed_coefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	//checks if the memory was allocated
	if (prevSolPtr == NULL || gauss_seidelPtr == NULL || presmoothed_coefPtr == NULL)
		return false;

	for (k = 0; k < height_ext; k++)
	{
		presmoothed_coefPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		gauss_seidelPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		prevSolPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		//checks if the memory was allocated
		if (presmoothed_coefPtr[k] == NULL || gauss_seidelPtr[k] == NULL || prevSolPtr[k] == NULL)
			return false;
	}

	//Create tempporary Image Data holder for diffusion coefficients
	//- with extended boundary because of boundary condition

	dataType** presmoot_e_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_w_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_n_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_s_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_t_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** presmoot_b_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (presmoot_e_coefPtr == NULL || presmoot_w_coefPtr == NULL || presmoot_n_coefPtr == NULL ||
		presmoot_s_coefPtr == NULL || presmoot_t_coefPtr == NULL || presmoot_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		presmoot_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		presmoot_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (presmoot_e_coefPtr[k] == NULL || presmoot_w_coefPtr[k] == NULL || presmoot_n_coefPtr[k] == NULL ||
			presmoot_s_coefPtr[k] == NULL || presmoot_t_coefPtr[k] == NULL || presmoot_b_coefPtr[k] == NULL)
			return false;
	}

	dataType** orig_e_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_w_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_n_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_s_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_t_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** orig_b_coefPtr = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (orig_e_coefPtr == NULL || orig_w_coefPtr == NULL || orig_n_coefPtr == NULL || orig_s_coefPtr == NULL ||
		orig_t_coefPtr == NULL || orig_b_coefPtr == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		orig_e_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_w_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_n_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_s_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_t_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		orig_b_coefPtr[k] = malloc(sizeof(dataType) * length * width);
		//checks if the memory was allocated
		if (orig_e_coefPtr[k] == NULL || orig_w_coefPtr[k] == NULL || orig_n_coefPtr[k] == NULL || orig_s_coefPtr[k] == NULL
			|| orig_t_coefPtr[k] == NULL || orig_b_coefPtr[k] == NULL)
			return false;
	}

	dataType** coefPtr_e = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_w = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_n = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_s = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_t = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** coefPtr_b = (dataType**)malloc(sizeof(dataType*) * height);
	//checks if the memory was allocated
	if (coefPtr_e == NULL || coefPtr_w == NULL || coefPtr_n == NULL || coefPtr_s == NULL || coefPtr_t == NULL ||
		coefPtr_b == NULL)
		return false;

	for (k = 0; k < height; k++)
	{
		coefPtr_e[k] = malloc(sizeof(float) * length * width);
		coefPtr_w[k] = malloc(sizeof(float) * length * width);
		coefPtr_n[k] = malloc(sizeof(float) * length * width);
		coefPtr_s[k] = malloc(sizeof(float) * length * width);
		coefPtr_t[k] = malloc(sizeof(float) * length * width);
		coefPtr_b[k] = malloc(sizeof(float) * length * width);
		//checks if the memory was allocated
		if (coefPtr_e[k] == NULL || coefPtr_w[k] == NULL || coefPtr_n[k] == NULL || coefPtr_s[k] == NULL ||
			coefPtr_t[k] == NULL || coefPtr_b[k] == NULL)
			return false;
	}

	Image_Data presmoothingData;
	presmoothingData.height = height_ext;
	presmoothingData.length = length_ext;
	presmoothingData.width = width_ext;
	presmoothingData.imageDataPtr = presmoothed_coefPtr;

	//copy data to extended area which will be used in each time step
	copyDataToExtendedArea(inputImageData.imageDataPtr, prevSolPtr, height, length, width);
	//copy data to extended area which will be used in each Gauss Seidel iteration
	copyDataToExtendedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);
	//copy data to extended area which will be used for calculation of diffusion coefficients
	copyDataToExtendedArea(inputImageData.imageDataPtr, presmoothed_coefPtr, height, length, width);

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(presmoothed_coefPtr, height_ext, length_ext, width_ext);
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);
	reflection3D(gauss_seidelPtr, height_ext, length_ext, width_ext);

	//perfom presmoothing
	heatExplicitRectangularScheme(presmoothingData, filterParameters);

	dataType hx = filterParameters.hx;
	dataType hy = filterParameters.hy;
	dataType hz = filterParameters.hz;
	dataType K = filterParameters.edge_detector_coefficient;
	dataType norm_gradient;

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = presmoothed_coefPtr[k_ext][x_ext];
				uN = presmoothed_coefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = presmoothed_coefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = presmoothed_coefPtr[k_ext][x_new(iplus1, j_ext, length_ext)];
				uW = presmoothed_coefPtr[k_ext][x_new(iminus1, j_ext, length_ext)];
				uNW = presmoothed_coefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = presmoothed_coefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = presmoothed_coefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = presmoothed_coefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				
				Tu = presmoothed_coefPtr[kminus1][x_ext];
				TuN = presmoothed_coefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = presmoothed_coefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = presmoothed_coefPtr[kminus1][x_new(iplus1, j_ext, length_ext)];
				TuW = presmoothed_coefPtr[kminus1][x_new(iminus1, j_ext, length_ext)];
				TuNW = presmoothed_coefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = presmoothed_coefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = presmoothed_coefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = presmoothed_coefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				
				Bu = presmoothed_coefPtr[kplus1][x_ext];
				BuN = presmoothed_coefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = presmoothed_coefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = presmoothed_coefPtr[kplus1][x_new(iplus1, j_ext, length_ext)];
				BuW = presmoothed_coefPtr[kplus1][x_new(iminus1, j_ext, length_ext)];
				BuNW = presmoothed_coefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = presmoothed_coefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = presmoothed_coefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = presmoothed_coefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//values of voxels in the extended data container for the original image
				orig_u = prevSolPtr[k_ext][x_ext];
				orig_uN = prevSolPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				orig_uS = prevSolPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				orig_uE = prevSolPtr[k_ext][x_new(iplus1, j_ext, length_ext)];
				orig_uW = prevSolPtr[k_ext][x_new(iminus1, j_ext, length_ext)];
				orig_uNW = prevSolPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				orig_uNE = prevSolPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				orig_uSE = prevSolPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				orig_uSW = prevSolPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				
				orig_Tu = prevSolPtr[kminus1][x_ext];
				orig_TuN = prevSolPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				orig_TuS = prevSolPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				orig_TuE = prevSolPtr[kminus1][x_new(iplus1, j_ext, length_ext)];
				orig_TuW = prevSolPtr[kminus1][x_new(iminus1, j_ext, length_ext)];
				orig_TuNW = prevSolPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				orig_TuNE = prevSolPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				orig_TuSE = prevSolPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				orig_TuSW = prevSolPtr[kminus1][x_new(iminus1, jplus1, length_ext)];

				orig_Bu = prevSolPtr[kplus1][x_ext];
				orig_BuN = prevSolPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				orig_BuS = prevSolPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				orig_BuE = prevSolPtr[kplus1][x_new(iplus1, j_ext, length_ext)];
				orig_BuW = prevSolPtr[kplus1][x_new(iminus1, j_ext, length_ext)];
				orig_BuNW = prevSolPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				orig_BuNE = prevSolPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				orig_BuSE = prevSolPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				orig_BuSW = prevSolPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / hx;
				uy = (dataType)(((uN + uNE) - (uS + uSE)) / (4.0 * hy));
				uz = (dataType)(((Tu + TuE) - (Bu + BuE)) / (4.0 * hz));
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_e_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				// Calculation of coefficients in west direction
				ux = (uW - u) / hx;
				uy = (dataType)(((uNW + uN) - (uSW + uS)) / (4.0 * hy));
				uz = (dataType)(((TuW + Tu) - (BuW + Bu)) / (4.0 * hz));
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_w_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				// Calculation of coefficients in north direction
				ux = (dataType)(((uNE + uE) - (uNW + uW)) / (4.0 * hx));
				uy = (uN - u) / hy;
				uz = (dataType)(((TuN + Tu) - (BuN + Bu)) / (4.0 * hz));
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_n_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				// Calculation of coefficients in south direction
				ux = (dataType)(((uE + uSE) - (uW + uSW)) / (4.0 * hx));
				uy = (uS - u) / hy;
				uz = (dataType)(((TuS + Tu) - (BuS + Bu)) / (4.0 * hz));
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_s_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				// Calculation of coefficients in top direction
				ux = (dataType)(((TuE + uE) - (TuW + uW)) / (4.0 * hx));
				uy = (dataType)(((TuN + uN) - (TuS + uS)) / (4.0 * hy));
				uz = (Tu - u) / hz;
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_t_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				// Calculation of coefficients in bottom direction
				ux = (dataType)(((BuW + uW) - (BuE + uE)) / (4.0 * hx));
				uy = (dataType)(((BuN + uN) - (BuS + uS)) / (4.0 * hy));
				uz = (Bu - u) / hz;
				norm_gradient = ux * ux + uy * uy + uz * uz;
				presmoot_b_coefPtr[k][x] = gradientFunction(norm_gradient, K);

				//calculation of coefficients in the original image data
				
				// Calculation of coefficients in east direction
				orig_ux = (orig_uE - orig_u) / hx;
				orig_uy = (dataType)(((orig_uN + orig_uNE) - (orig_uS + orig_uSE)) / (4.0 * hy));
				orig_uz = (dataType)(((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE)) / (4.0 * hz));
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_e_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / hx;
				orig_uy = (dataType)(((orig_uNW + orig_uN) - (orig_uSW + orig_uS)) / (4.0 * hy));
				orig_uz = (dataType)(((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu)) / (4.0 * hz));
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_w_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = (dataType)(((orig_uNE + orig_uE) - (orig_uNW + orig_uW)) / (4.0 * hx));
				orig_uy = (orig_uN - orig_u) / hy;
				orig_uz = (dataType)(((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu)) / (4.0 * hz));
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_n_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = (dataType)(((orig_uE + orig_uSE) - (orig_uW + orig_uSW)) / (4.0 * hx));
				orig_uy = (orig_uS - orig_u) / hy;
				orig_uz = (dataType)(((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu)) / (4.0 * hz));
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_s_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = (dataType)(((orig_TuE + orig_uE) - (orig_TuW + orig_uW)) / (4.0 * hx));
				orig_uy = (dataType)(((orig_TuN + orig_uN) - (orig_TuS + orig_uS)) / (4.0 * hy));
				orig_uz = (orig_Tu - orig_u) / hz;
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_t_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = (dataType)(((orig_BuW + orig_uW) - (orig_BuE + orig_uE)) / (4.0 * hx));
				orig_uy = (dataType)(((orig_BuN + orig_uN) - (orig_BuS + orig_uS)) / (4.0 * hy));
				orig_uz = (orig_Bu - orig_u) / hz;
				norm_gradient = orig_ux * orig_ux + orig_uy * orig_uy + orig_uz * orig_uz;
				orig_b_coefPtr[k][x] = (dataType)sqrt(norm_gradient + filterParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e_coefPtr[k][x] + orig_w_coefPtr[k][x] + orig_n_coefPtr[k][x] + orig_s_coefPtr[k][x]
					+ orig_t_coefPtr[k][x] + orig_b_coefPtr[k][x]) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + filterParameters.eps2);

				//evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				//image at each voxel face and reciprocal of norm of gradient of image at each voxel face

				coefPtr_e[k][x] = (dataType)(voxel_coef * presmoot_e_coefPtr[k][x] * (hy * hz / hx) * (1.0 / orig_e_coefPtr[k][x]));//east coefficient
				coefPtr_w[k][x] = (dataType)(voxel_coef * presmoot_w_coefPtr[k][x] * (hy * hz / hx) * (1.0 / orig_w_coefPtr[k][x]));//west coefficient
				coefPtr_n[k][x] = (dataType)(voxel_coef * presmoot_n_coefPtr[k][x] * (hx * hz / hy) * (1.0 / orig_n_coefPtr[k][x]));//north coefficient
				coefPtr_s[k][x] = (dataType)(voxel_coef * presmoot_s_coefPtr[k][x] * (hx * hz / hy) * (1.0 / orig_s_coefPtr[k][x]));//south coefficient
				coefPtr_t[k][x] = (dataType)(voxel_coef * presmoot_t_coefPtr[k][x] * (hx * hy / hz) * (1.0 / orig_t_coefPtr[k][x]));//top coefficient
				coefPtr_b[k][x] = (dataType)(voxel_coef * presmoot_b_coefPtr[k][x] * (hx * hy / hz) * (1.0 / orig_b_coefPtr[k][x]));//bottom coefficient
			}
		}
	}
	
	// The Implicit Scheme Evaluation
	count_step = 0;
	do
	{
		count_step++;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Begin Gauss-Seidel Formula Evaluation
					gauss_seidel = (prevSolPtr[k_ext][x_ext] + coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coef_tauh * (coefPtr_e[k][x] + coefPtr_w[k][x] + coefPtr_n[k][x]
							+ coefPtr_s[k][x] + coefPtr_t[k][x] + coefPtr_b[k][x]));

					// SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] +
						filterParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);
				}
			}
		}

		// Error Evaluation
		error = 0.0; // Initialize

		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					error += (dataType)pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (coefPtr_e[k][x]
						+ coefPtr_w[k][x] + coefPtr_n[k][x] + coefPtr_s[k][x] + coefPtr_t[k][x] + coefPtr_b[k][x]))
						- coef_tauh * ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)])
							+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)])
							+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSolPtr[k_ext][x_ext], 2);
				}
			}
		}

	} while (error > filterParameters.tolerance && count_step < filterParameters.maxNumberOfSolverIteration);
	
	printf("The number of iterations is %zd\n", count_step);
	printf("Error is %e\n", error);
	//printf("Step is %zd\n", filterParameters.timeStepsNum);

	//Copy the current time step to original data holder after timeStepsNum
	copyDataToReducedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);

	// Freeing Memory after use
	for (k = 0; k < height_ext; k++)
	{
		free(presmoothed_coefPtr[k]);
		free(gauss_seidelPtr[k]);
		free(prevSolPtr[k]);
	}
	free(presmoothed_coefPtr);
	free(gauss_seidelPtr);
	free(prevSolPtr);

	// free _coefPtr pointers
	for (k = 0; k < height; k++)
	{
		free(presmoot_e_coefPtr[k]);
		free(presmoot_w_coefPtr[k]);
		free(presmoot_n_coefPtr[k]);
		free(presmoot_s_coefPtr[k]);
		free(presmoot_t_coefPtr[k]);
		free(presmoot_b_coefPtr[k]);
	}
	free(presmoot_e_coefPtr);
	free(presmoot_w_coefPtr);
	free(presmoot_n_coefPtr);
	free(presmoot_s_coefPtr);
	free(presmoot_t_coefPtr);
	free(presmoot_b_coefPtr);

	// free orig_ _coefPtr pointers
	for (k = 0; k < height; k++)
	{
		free(orig_e_coefPtr[k]);
		free(orig_w_coefPtr[k]);
		free(orig_n_coefPtr[k]);
		free(orig_s_coefPtr[k]);
		free(orig_t_coefPtr[k]);
		free(orig_b_coefPtr[k]);
	}
	free(orig_e_coefPtr);
	free(orig_w_coefPtr);
	free(orig_n_coefPtr);
	free(orig_s_coefPtr);
	free(orig_t_coefPtr);
	free(orig_b_coefPtr);

	// free coefPtr_ pointers
	for (k = 0; k < height; k++)
	{
		free(coefPtr_e[k]);
		free(coefPtr_w[k]);
		free(coefPtr_n[k]);
		free(coefPtr_s[k]);
		free(coefPtr_t[k]);
		free(coefPtr_b[k]);
	}
	free(coefPtr_e);
	free(coefPtr_w);
	free(coefPtr_n);
	free(coefPtr_s);
	free(coefPtr_t);
	free(coefPtr_b);

	return true;
}

bool geodesicMeanCurvature2D(Image_Data2D inputImage, Filter_Parameters filtering_parameters)
{

	if (inputImage.imageDataPtr == NULL) {
		return false;
	}

	const size_t length = inputImage.height, width = inputImage.width;
	PixelSpacing spacing = inputImage.spacing;
	dataType eps2 = filtering_parameters.eps2;

	size_t i, j, i_ext, j_ext;
	const size_t length_ext = length + 2, width_ext = width + 2;
	size_t dim2D = length * width, dim2D_ext = length_ext * width_ext;
	dataType tau = filtering_parameters.timeStepSize, h = filtering_parameters.h;
	dataType tol = filtering_parameters.tolerance, omega = filtering_parameters.omega_c;
	dataType coef_edge_detector = filtering_parameters.edge_detector_coefficient, eps = filtering_parameters.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;
	size_t maxIter = filtering_parameters.maxNumberOfSolverIteration;

	dataType* smothedImage = (dataType*)malloc(dim2D * sizeof(dataType));
	Point2D iOrigin = { 0.0, 0.0 };
	Image_Data2D imageData = { length, width, smothedImage, iOrigin, spacing };

	copyDataToAnother2dArray(inputImage.imageDataPtr, smothedImage, length, width);

	//heat2dExplicitScheme(imageData, filtering_parameters);
	//heatImplicit2dScheme(imageData, filtering_parameters);
	gaussianSmoothing2D(inputImage.imageDataPtr, smothedImage, length, width, 1.0);

	dataType* previous_Solution = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	dataType* gauss_Seidel_Sol = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	dataType* extended_image_data = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	if (previous_Solution == NULL || gauss_Seidel_Sol == NULL) {
		return false;
	}

	copyDataTo2dExtendedArea(smothedImage, extended_image_data, length, width);
	reflection2D(extended_image_data, length_ext, width_ext);

	dataType* uNorth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* uSouth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* uEast = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* uWest = (dataType*)malloc(dim2D * sizeof(dataType));
	if (uNorth == NULL || uSouth == NULL || uEast == NULL || uWest == NULL)
		return false;

	dataType* gNorth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* gSouth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* gEast = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* gWest = (dataType*)malloc(dim2D * sizeof(dataType));
	if (gNorth == NULL || gSouth == NULL || gEast == NULL || gWest == NULL)
		return false;

	dataType* coefNorth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* coefSouth = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* coefEast = (dataType*)malloc(dim2D * sizeof(dataType));
	dataType* coefWest = (dataType*)malloc(dim2D * sizeof(dataType));
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	dataType uP, uN, uNW, uNE, uS, uSW, uSE, uW, uE;
	dataType ux, uy, current_value, u_average, avg_norm_gardient;

	for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

			size_t iplus = i_ext + 1;
			size_t iminus = i_ext - 1;
			size_t jplus = j_ext + 1;
			size_t jminus = j_ext - 1;

			size_t currentIndx = x_new(i, j, length);
			uP = extended_image_data[x_new(i_ext, j_ext, length_ext)];
			uE = extended_image_data[x_new(iplus, j_ext, length_ext)];
			uW = extended_image_data[x_new(iminus, j_ext, length_ext)];
			uN = extended_image_data[x_new(i_ext, jminus, length_ext)];
			uS = extended_image_data[x_new(i_ext, jplus, length_ext)];
			uNE = extended_image_data[x_new(iplus, jminus, length_ext)];
			uNW = extended_image_data[x_new(iminus, jminus, length_ext)];
			uSE = extended_image_data[x_new(iplus, jplus, length_ext)];
			uSW = extended_image_data[x_new(iminus, jplus, length_ext)];

			//East
			ux = (uE - uP) / h;
			uy = (uNE + uN - uS - uSE) / (4.0 * h);
			current_value = ux * ux + uy * uy;
			uEast[currentIndx] = sqrt(current_value + eps2);
			gEast[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//West
			ux = (uP - uW) / h;
			uy = (uNW + uN - uSW - uS) / (4.0 * h);
			current_value = ux * ux + uy * uy;
			uWest[currentIndx] = sqrt(current_value + eps2);
			gEast[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * h);
			uy = (uN - uP) / h;
			current_value = ux * ux + uy * uy;
			uNorth[currentIndx] = (current_value + eps2);
			gNorth[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * h);
			uy = (uP - uS) / h;
			current_value = ux * ux + uy * uy;
			uSouth[currentIndx] = (current_value + eps2);
			gSouth[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			u_average = (dataType)((uEast[currentIndx] + uWest[currentIndx] + uNorth[currentIndx] + uSouth[currentIndx]) / 4.0);
			avg_norm_gardient = sqrt(u_average * u_average + eps);

			coefEast[currentIndx] = (dataType)(coef_tau * avg_norm_gardient * gEast[currentIndx] * (1.0 / uEast[currentIndx]));
			coefWest[currentIndx] = (dataType)(coef_tau * avg_norm_gardient * gWest[currentIndx] * (1.0 / uWest[currentIndx]));
			coefNorth[currentIndx] = (dataType)(coef_tau * avg_norm_gardient * gNorth[currentIndx] * (1.0 / uNorth[currentIndx]));
			coefSouth[currentIndx] = (dataType)(coef_tau * avg_norm_gardient * gSouth[currentIndx] * (1.0 / uSouth[currentIndx]));

		}
	}

	copyDataToAnother2dArray(extended_image_data, gauss_Seidel_Sol, length_ext, width_ext);
	copyDataToAnother2dArray(extended_image_data, previous_Solution, length_ext, width_ext);

	//copyDataTo2dExtendedArea(inputImage.imageDataPtr, gauss_Seidel_Sol, length, width);
	//reflection2D(gauss_Seidel_Sol, length_ext, width_ext);
	//copyDataTo2dExtendedArea(inputImage.imageDataPtr, previous_Solution, length, width);
	//reflection2D(previous_Solution, length_ext, width_ext);

	size_t cpt = 0;
	dataType error = 0.0;
	do {
		cpt++;

		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

				size_t currentIndx = x_new(i, j, length);
				size_t currentIndx_ext = x_new(i_ext, j_ext, length_ext);

				gauss_seidel_coef = (dataType)((previous_Solution[currentIndx_ext] + coefEast[currentIndx] * gauss_Seidel_Sol[x_new(i_ext + 1, j_ext, length_ext)] + coefNorth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext - 1, length_ext)]
					+ coefWest[currentIndx] * gauss_Seidel_Sol[x_new(i_ext - 1, j_ext, length_ext)] + coefSouth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext + 1, length_ext)])
					/ (1 + coefEast[currentIndx] + coefNorth[currentIndx] + coefWest[currentIndx] + coefSouth[currentIndx]));

				gauss_Seidel_Sol[currentIndx_ext] = gauss_Seidel_Sol[currentIndx_ext] + omega * (gauss_seidel_coef - gauss_Seidel_Sol[currentIndx_ext]);
			}
		}

		error = 0.0;
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

				size_t currentIndx = x_new(i, j, length);
				size_t currentIndx_ext = x_new(i_ext, j_ext, length_ext);

				error += (dataType)(pow((1 + coefEast[currentIndx] + coefNorth[currentIndx] + coefWest[currentIndx] + coefSouth[currentIndx]) * gauss_Seidel_Sol[currentIndx_ext]
					- (coefEast[currentIndx] * gauss_Seidel_Sol[x_new(i_ext + 1, j_ext, length_ext)] + coefNorth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext - 1, length_ext)]
						+ coefWest[currentIndx] * gauss_Seidel_Sol[x_new(i_ext - 1, j_ext, length_ext)] + coefSouth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext + 1, length_ext)]) - previous_Solution[currentIndx_ext], 2) * h * h);
			}
		}

	} while (cpt < maxIter && error > tol);

	//Copy back
	copyDataTo2dReducedArea(inputImage.imageDataPtr, gauss_Seidel_Sol, length, width);

	free(uNorth);
	free(uSouth);
	free(uEast);
	free(uWest);

	free(gNorth);
	free(gSouth);
	free(gEast);
	free(gWest);

	free(coefNorth);
	free(coefSouth);
	free(coefEast);
	free(coefWest);

	free(previous_Solution);
	free(gauss_Seidel_Sol);
	free(extended_image_data);
}