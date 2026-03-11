#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool

#include "file.h"
#include "data_storage.h"

#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "setting_boundary_values.h"
#include "common_functions.h"
#include "filter_params.h"
#include "conmon_filtering.h"


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
	dataType tau = filtering_parameters.timeStepSize;
	dataType hx = spacing.sx, hy = spacing.sy;
	dataType hx2 = hx * hx, hy2 = hy * hy;
	dataType tol = filtering_parameters.tolerance, omega = filtering_parameters.omega_c;
	dataType coef_edge_detector = filtering_parameters.edge_detector_coefficient;
	dataType gauss_seidel_coef = 0.0;
	size_t maxIter = filtering_parameters.maxNumberOfSolverIteration;

	dataType* smothedImage = (dataType*)malloc(dim2D * sizeof(dataType));
	Image_Data2D imageData = { length, width, smothedImage, inputImage.origin, spacing };

	copyDataToAnother2dArray(inputImage.imageDataPtr, smothedImage, length, width);

	heatImplicit2dScheme(imageData, filtering_parameters);
	//gaussianSmoothing2D(inputImage.imageDataPtr, smothedImage, length, width, 1.0);

	Storage_Flags flags = { false, false };
	const char storing_path[] = "C:/Users/Konan Allaly/Documents/Tests/output/smoothed.raw";
	store2dRawData(imageData.imageDataPtr, length, width, storing_path, flags);

	dataType* previous_Solution = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	dataType* gauss_Seidel_Sol = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	dataType* extended_image_data = (dataType*)malloc(dim2D_ext * sizeof(dataType));
	if (previous_Solution == NULL || gauss_Seidel_Sol == NULL) {
		return false;
	}

	copyDataTo2dExtendedArea(smothedImage, extended_image_data, length, width);
	reflection2D(extended_image_data, length_ext, width_ext);

	copyDataToAnother2dArray(extended_image_data, gauss_Seidel_Sol, length_ext, width_ext);
	copyDataToAnother2dArray(extended_image_data, previous_Solution, length_ext, width_ext);

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
			ux = (uE - uP) / hx;
			uy = (uNE + uN - uS - uSE) / (4.0 * hy);
			current_value = ux * ux + uy * uy;
			uEast[currentIndx] = sqrt(current_value + eps2);
			gEast[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//West
			ux = (uP - uW) / hx;
			uy = (uNW + uN - uSW - uS) / (4.0 * hy);
			current_value = ux * ux + uy * uy;
			uWest[currentIndx] = sqrt(current_value + eps2);
			gWest[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * hx);
			uy = (uN - uP) / hy;
			current_value = ux * ux + uy * uy;
			uNorth[currentIndx] = sqrt(current_value + eps2);
			gNorth[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * hx);
			uy = (uP - uS) / hy;
			current_value = ux * ux + uy * uy;
			uSouth[currentIndx] = sqrt(current_value + eps2);
			gSouth[currentIndx] = gradientFunction(current_value, coef_edge_detector);

			u_average = (dataType)((uEast[currentIndx] + uWest[currentIndx] + uNorth[currentIndx] + uSouth[currentIndx]) / 4.0);
			avg_norm_gardient = sqrt(u_average * u_average + eps2);

			coefEast[currentIndx] = (dataType)(tau * avg_norm_gardient * gEast[currentIndx] * (1.0 / (hx2 * uEast[currentIndx])));
			coefWest[currentIndx] = (dataType)(tau * avg_norm_gardient * gWest[currentIndx] * (1.0 / (hx2 * uWest[currentIndx])));
			coefNorth[currentIndx] = (dataType)(tau * avg_norm_gardient * gNorth[currentIndx] * (1.0 / (hy2 * uNorth[currentIndx])));
			coefSouth[currentIndx] = (dataType)(tau * avg_norm_gardient * gSouth[currentIndx] * (1.0 / (hy2 * uSouth[currentIndx])));

		}
	}

	const char edge_detector_path[] = "C:/Users/Konan Allaly/Documents/Tests/output/edge_East.raw";
	store2dRawData(gEast, length, width, edge_detector_path, flags);

	size_t cpt = 0;
	dataType error = 0.0;
	
	do {
		cpt++;

		for (i = 0, i_ext = 1; i < length; i++, i_ext++) 
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
			{

				size_t currentIndx = x_new(i, j, length);
				size_t currentIndx_ext = x_new(i_ext, j_ext, length_ext);

				gauss_seidel_coef = (dataType)((previous_Solution[currentIndx_ext] + 
					coefEast[currentIndx] * gauss_Seidel_Sol[x_new(i_ext + 1, j_ext, length_ext)] + 
					coefNorth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext - 1, length_ext)]+ 
					coefWest[currentIndx] * gauss_Seidel_Sol[x_new(i_ext - 1, j_ext, length_ext)] + 
					coefSouth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext + 1, length_ext)])
					/ (1 + coefEast[currentIndx] + coefNorth[currentIndx] + coefWest[currentIndx] + coefSouth[currentIndx]));

				gauss_Seidel_Sol[currentIndx_ext] = gauss_Seidel_Sol[currentIndx_ext] + omega * (gauss_seidel_coef - gauss_Seidel_Sol[currentIndx_ext]);
			}
		}

		error = 0.0;
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) 
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
			{

				size_t currentIndx = x_new(i, j, length);
				size_t currentIndx_ext = x_new(i_ext, j_ext, length_ext);

				error += (dataType)(pow((1 + coefEast[currentIndx] + coefNorth[currentIndx] + 
					coefWest[currentIndx] + coefSouth[currentIndx]) * gauss_Seidel_Sol[currentIndx_ext]
					- (coefEast[currentIndx] * gauss_Seidel_Sol[x_new(i_ext + 1, j_ext, length_ext)] + 
						coefNorth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext - 1, length_ext)] + 
						coefWest[currentIndx] * gauss_Seidel_Sol[x_new(i_ext - 1, j_ext, length_ext)] + 
						coefSouth[currentIndx] * gauss_Seidel_Sol[x_new(i_ext, j_ext + 1, length_ext)]) 
					- previous_Solution[currentIndx_ext], 2) * hx * hy);
			}
		}
	} while (cpt < maxIter && error > tol);
	printf("The number of iterations is %zd\n", cpt);

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

bool geodesicMeanCurvature(Image_Data inputImageData, const Filter_Parameters filterParameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL)
		return false;

	size_t k, i, j, x, steps = filterParameters.maxNumberOfSolverIteration;
	size_t k_ext, j_ext, i_ext, x_ext;

	dataType hx = inputImageData.spacing.sx;
	dataType hy = inputImageData.spacing.sy;
	dataType hz = inputImageData.spacing.sz;
	dataType hx_2 = hx * hx;
	dataType hy_2 = hy * hy;
	dataType hz_2 = hz * hz;
	dataType tau = filterParameters.timeStepSize;
	dataType coef_edge = filterParameters.edge_detector_coefficient;
	dataType eps2 = filterParameters.eps2;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType error, gauss_seidel;

	// Prepare variables toExplicitImage.height, toExplicitImage.length, toExplicitImage.width
	size_t height = inputImageData.height, length = inputImageData.length, width = inputImageData.width;
	size_t dim2D = length * width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	
	//Presmoothing step
	dataType** presmoothedImage = (dataType**)malloc(sizeof(dataType*) * height);
	if (presmoothedImage == NULL)
	{
		return false;
	}
	for(k = 0; k < height; k++)
	{
		presmoothedImage[k] = malloc(sizeof(dataType) * dim2D);
		if (presmoothedImage[k] == NULL)
		{
			return false;
		}
	}

	copyDataToAnotherArray(inputImageData.imageDataPtr, presmoothedImage, height, length, width);

	Image_Data presmoothingData = {height, length, width, presmoothedImage, inputImageData.origin, inputImageData.spacing, inputImageData.orientation};

	heatImplicitRectangularScheme(presmoothingData, filterParameters);

	////save smoothed image
	//char path_saving [] = "C:/Users/Konan Allaly/Documents/Tests/output/Data journal paper submission/smoothed_p4.raw";
	//Storage_Flags flags = { false,false };
	//store3dDataArrayD(presmoothedImage, length, width, height, path_saving, flags);

	// Compute edge detector coefficients on presmoothed image
	Pointers_Neighbours norm_of_gradient;
	norm_of_gradient.east = malloc(sizeof(dataType*) * height);
	norm_of_gradient.west = malloc(sizeof(dataType*) * height);
	norm_of_gradient.north = malloc(sizeof(dataType*) * height);
	norm_of_gradient.south = malloc(sizeof(dataType*) * height);
	norm_of_gradient.top = malloc(sizeof(dataType*) * height);
	norm_of_gradient.bottom = malloc(sizeof(dataType*) * height);
	for(k = 0; k < height; k++)
	{
		norm_of_gradient.east[k] = malloc(sizeof(dataType) * dim2D);
		norm_of_gradient.west[k] = malloc(sizeof(dataType) * dim2D);
		norm_of_gradient.north[k] = malloc(sizeof(dataType) * dim2D);
		norm_of_gradient.south[k] = malloc(sizeof(dataType) * dim2D);
		norm_of_gradient.top[k] = malloc(sizeof(dataType) * dim2D);
		norm_of_gradient.bottom[k] = malloc(sizeof(dataType) * dim2D);
		if (norm_of_gradient.east[k] == NULL || norm_of_gradient.west[k] == NULL || norm_of_gradient.north[k] == NULL ||
			norm_of_gradient.south[k] == NULL || norm_of_gradient.top[k] == NULL || norm_of_gradient.bottom[k] == NULL)
			return false;
	}
	if (norm_of_gradient.east == NULL || norm_of_gradient.west == NULL || norm_of_gradient.north == NULL ||
		norm_of_gradient.south == NULL || norm_of_gradient.top == NULL || norm_of_gradient.bottom == NULL)
		return false;

	//Copy to extended area
	// Create temporary Image Data holder for Previous time step data - with extended boundary because of boundary condition
	dataType** prevSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	// Create temporary Image Data holder for Current time step data - with extended boundary because of boundary condition
	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	//checks if the memory was allocated
	if (prevSolPtr == NULL || gauss_seidelPtr == NULL)
		return false;
	for (k = 0; k < height_ext; k++)
	{
		gauss_seidelPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		prevSolPtr[k] = malloc(sizeof(dataType) * length_ext * width_ext);
		//checks if the memory was allocated
		if (gauss_seidelPtr[k] == NULL || prevSolPtr[k] == NULL)
			return false;
	}

	//Copy the presmoothed image to the extended area
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				prevSolPtr[k_ext][x_ext] = presmoothedImage[k][x];
				gauss_seidelPtr[k_ext][x_ext] = presmoothedImage[k][x];
			}
		}
	}

	//perform reflection of the extended area to ensure zero Neumann boundary condition (for LHE)
	reflection3D(prevSolPtr, height_ext, length_ext, width_ext);
	reflection3D(gauss_seidelPtr, height_ext, length_ext, width_ext);

	//Compute the coefficients for the original image
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
		coefPtr_e[k] = malloc(sizeof(dataType) * dim2D);
		coefPtr_w[k] = malloc(sizeof(dataType) * dim2D);
		coefPtr_n[k] = malloc(sizeof(dataType) * dim2D);
		coefPtr_s[k] = malloc(sizeof(dataType) * dim2D);
		coefPtr_t[k] = malloc(sizeof(dataType) * dim2D);
		coefPtr_b[k] = malloc(sizeof(dataType) * dim2D);
		//checks if the memory was allocated
		if (coefPtr_e[k] == NULL || coefPtr_w[k] == NULL || coefPtr_n[k] == NULL || coefPtr_s[k] == NULL ||
			coefPtr_t[k] == NULL || coefPtr_b[k] == NULL)
			return false;
	}

	normOfGradientReducedDiamondCells(presmoothingData, norm_of_gradient);

	//calculation of coefficients
	dataType voxel_coef, average_face_coef;
	dataType g_east, g_west, g_north, g_south, g_top, g_bottom;
	dataType n_east, n_west, n_north, n_south, n_top, n_bottom;
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				
				//edge detector function at each voxel face
				g_east = gradientFunction(norm_of_gradient.east[k][x], coef_edge);
				g_west = gradientFunction(norm_of_gradient.west[k][x], coef_edge);
				g_north = gradientFunction(norm_of_gradient.north[k][x], coef_edge);
				g_south = gradientFunction(norm_of_gradient.south[k][x], coef_edge);
				g_top = gradientFunction(norm_of_gradient.top[k][x], coef_edge);
				g_bottom = gradientFunction(norm_of_gradient.bottom[k][x], coef_edge);

				//epsilon regularization
				n_east = sqrt(norm_of_gradient.west[k][x] + eps2);
				n_west = sqrt(norm_of_gradient.west[k][x] + eps2);
				n_north = sqrt(norm_of_gradient.north[k][x] + eps2);
				n_south = sqrt(norm_of_gradient.south[k][x] + eps2);
				n_top = sqrt(norm_of_gradient.top[k][x] + eps2);
				n_bottom = sqrt(norm_of_gradient.bottom[k][x] + eps2);

				//average of the norm of gradient at the voxel faces
				average_face_coef = (dataType)((n_east + n_west + n_north + n_south + n_top + n_bottom) / 6.0);
				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + eps2);

				//evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				//image at each voxel face and reciprocal of norm of gradient of image at each voxel face
				coefPtr_e[k][x] = (dataType)(tau * voxel_coef * g_east / (n_east * hx_2));
				coefPtr_w[k][x] = (dataType)(tau * voxel_coef * g_west / (n_west * hx_2));
				coefPtr_n[k][x] = (dataType)(tau * voxel_coef * g_north / (n_north * hy_2));
				coefPtr_s[k][x] = (dataType)(tau * voxel_coef * g_south / (n_south * hy_2));
				coefPtr_t[k][x] = (dataType)(tau * voxel_coef * g_top / (n_top * hz_2));
				coefPtr_b[k][x] = (dataType)(tau * voxel_coef * g_bottom / (n_bottom * hz_2));
			}
		}
	}

	// The Implicit Scheme Evaluation
	size_t count_iteration = 0; // Steps counter
	do
	{
		count_iteration++;
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
					gauss_seidel = (dataType)((prevSolPtr[k_ext][x_ext] + ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
						+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
						+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
						+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
						+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
						+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))) /
						(1 + coefPtr_e[k][x] + coefPtr_w[k][x] + coefPtr_n[k][x]
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

					error += (dataType)(pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coefPtr_e[k][x]
						+ coefPtr_w[k][x] + coefPtr_n[k][x] + coefPtr_s[k][x]
						+ coefPtr_t[k][x] + coefPtr_b[k][x])
						- ((coefPtr_e[k][x] * gauss_seidelPtr[k_ext][x_ext + 1])
							+ (coefPtr_w[k][x] * gauss_seidelPtr[k_ext][x_ext - 1])
							+ (coefPtr_s[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)])
							+ (coefPtr_n[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)])
							+ (coefPtr_b[k][x] * gauss_seidelPtr[k_ext + 1][x_ext])
							+ (coefPtr_t[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])) - prevSolPtr[k_ext][x_ext], 2));
				}
			}
		}
	} while (error > filterParameters.tolerance && count_iteration < filterParameters.maxNumberOfSolverIteration);
	printf("The number of iterations is %zd\n", count_iteration);
	printf("Error is %e\n", error);

	//Copy the current time step to original data holder after timeStepsNum
	copyDataToReducedArea(inputImageData.imageDataPtr, gauss_seidelPtr, height, length, width);
	
	for (k = 0; k < height_ext; k++) 
	{
		if(k < height)
		{
			free(presmoothedImage[k]);
			free(norm_of_gradient.east[k]);
			free(norm_of_gradient.west[k]);
			free(norm_of_gradient.north[k]);
			free(norm_of_gradient.south[k]);
			free(norm_of_gradient.top[k]);
			free(norm_of_gradient.bottom[k]);
			free(coefPtr_e[k]);
			free(coefPtr_w[k]);
			free(coefPtr_n[k]);
			free(coefPtr_s[k]);
			free(coefPtr_t[k]);
			free(coefPtr_b[k]);
		}
		free(prevSolPtr[k]);
		free(gauss_seidelPtr[k]);
	}
	free(presmoothedImage);
	free(norm_of_gradient.east);
	free(norm_of_gradient.west);
	free(norm_of_gradient.north);
	free(norm_of_gradient.south);
	free(norm_of_gradient.top);
	free(norm_of_gradient.bottom);
	free(coefPtr_e);
	free(coefPtr_w);
	free(coefPtr_n);
	free(coefPtr_s);
	free(coefPtr_t);
	free(coefPtr_b);
	free(prevSolPtr);
	free(gauss_seidelPtr);

	return true;
}