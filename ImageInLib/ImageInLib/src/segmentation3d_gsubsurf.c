#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <time.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include <string.h>

#include "segmentation3D_subsurf.h"
#include "segmentation3d_gsubsurf.h"

#include "file.h"
#include "heat_equation.h"
#include "non_linear_heat_equation.h"
#include "image_norm.h"
#include "data_initialization.h"
#include "edgedetection.h"
#include "data_storage.h"
#include "data_load.h"
#include "common_functions.h"
#include "setting_boundary_values.h"

bool generalizedSubsurfSegmentation(Image_Data inputImageData, dataType** initialSegment, Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters, unsigned char* outputPathPtr) {

	if (inputImageData.imageDataPtr == NULL || initialSegment == NULL || outputPathPtr == NULL)
		return false;

	size_t i, j, k, xd;
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t dim2D = length * width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t dim2D_ext = length_ext * width_ext;
	dataType h = segParameters.h;

	dataType coef_conv = segParameters.coef_conv;
	dataType coef_dif = segParameters.coef_dif;

	dataType difference_btw_current_and_previous_sol = 0.0;

	dataType** gauss_seidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** prevSol_extPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);

	dataType** imageToBeSegPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** segmFuntionPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeGradientPtr = (dataType**)malloc(sizeof(dataType*) * height);

	dataType** e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);

	dataType** VePtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VwPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VnPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VsPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VtPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** VbPtr = (dataType**)malloc(sizeof(dataType*) * height);

	//checks if the memory was allocated
	if (gauss_seidelPtr == NULL || prevSol_extPtr == NULL || imageToBeSegPtr == NULL || segmFuntionPtr == NULL || edgeGradientPtr == NULL ||
		e_Ptr == NULL || w_Ptr == NULL || n_Ptr == NULL || s_Ptr == NULL || t_Ptr == NULL || b_Ptr == NULL ||
		VePtr == NULL || VwPtr == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VbPtr == NULL)
		return false;
	for (i = 0; i < height; i++)
	{
		imageToBeSegPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segmFuntionPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		edgeGradientPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		e_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		w_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		n_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		s_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		t_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		b_Ptr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		VePtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VwPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VnPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VsPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VtPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);
		VbPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D);

		//checks if the memory was allocated
		if (imageToBeSegPtr[i] == NULL || segmFuntionPtr[i] == NULL || edgeGradientPtr[i] == NULL || e_Ptr[i] == NULL
			|| w_Ptr[i] == NULL || n_Ptr[i] == NULL || s_Ptr[i] == NULL || t_Ptr[i] == NULL || b_Ptr[i] == NULL
			|| VePtr[i] == NULL || VwPtr[i] == NULL || VnPtr == NULL || VsPtr == NULL || VtPtr == NULL || VbPtr == NULL)
			return false;
	}

	for (i = 0; i < height_ext; i++)
	{
		gauss_seidelPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		prevSol_extPtr[i] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		//checks if the memory was allocated
		if (gauss_seidelPtr[i] == NULL || prevSol_extPtr[i] == NULL)
			return false;
	}

	//Initialize the arrays to avoid unwanted values
	initialize3dArrayD(gauss_seidelPtr, length_ext, width_ext, height_ext, 0.0);
	initialize3dArrayD(prevSol_extPtr, length_ext, width_ext, height_ext, 0.0);

	Coefficient_Pointers CoefPtrs;
	CoefPtrs.e_Ptr = e_Ptr;
	CoefPtrs.w_Ptr = w_Ptr;
	CoefPtrs.n_Ptr = n_Ptr;
	CoefPtrs.s_Ptr = s_Ptr;
	CoefPtrs.t_Ptr = t_Ptr;
	CoefPtrs.b_Ptr = b_Ptr;

	Gradient_Pointers VPtrs;
	VPtrs.GePtr = VePtr;
	VPtrs.GwPtr = VwPtr;
	VPtrs.GnPtr = VnPtr;
	VPtrs.GsPtr = VsPtr;
	VPtrs.GtPtr = VtPtr;
	VPtrs.GbPtr = VbPtr;

	//Array for name construction
	unsigned char  name[500];
	unsigned char  name_ending[200];
	Storage_Flags flags = { false,false };
	
	size_t k_n, i_n, j_n, xd_n;
	for (k = 0, k_n = 1; k < height; k++, k_n++) {
		for (i = 0, i_n = 1; i < length; i++, i_n++) {
			for (j = 0, j_n = 1; j < width; j++, j_n++) {
				xd = x_new(i, j, length);
				xd_n = x_new(i_n, j_n, length_ext);
				segmFuntionPtr[k][xd] = initialSegment[k][xd];
				gauss_seidelPtr[k_n][xd_n] = initialSegment[k][xd];
				prevSol_extPtr[k_n][xd_n] = initialSegment[k][xd];
			}
		}
	}

	//Set the boundary values
	setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
	setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

	//compute coefficients from presmoothed image
	generalizedGFunctionForImageToBeSegmented(inputImageData, edgeGradientPtr, VPtrs, segParameters, explicit_lhe_Parameters);

	bool isFileSaved;
	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_smoothed.raw");
	strcat_s(name, sizeof(name), name_ending);
	isFileSaved = manageFile(inputImageData.imageDataPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
	if (isFileSaved == false) {
		printf("The file was not saved\n");
		return false;
	}

	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	isFileSaved = manageFile(edgeGradientPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
	if (isFileSaved == false) {
		printf("The file was not saved\n");
		return false;
	}
	
	//i = 0;
	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%05zd.raw", i);
	//strcat_s(name, sizeof(name), name_ending);
	//isFileSaved = manageFile(initialSegment, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
	//if (isFileSaved == false) {
	//	printf("The file was not saved\n");
	//	return false;
	//}

	//FILE* error_file;
	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_l2_norm.csv");
	//strcat_s(name, sizeof(name), name_ending);
	//if (fopen_s(&error_file, name, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(error_file, "ID,distance\n");

	Image_Data segmentationFunction;
	segmentationFunction.height = height;
	segmentationFunction.length = length;
	segmentationFunction.width = width;
	segmentationFunction.imageDataPtr = segmFuntionPtr;
	
	
	//loop for segmentation time steps	
	size_t number_time_step = 0;
	do
	{
		number_time_step++;
		setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
		setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

		////calcution of coefficients
		generalizedGaussSeidelCoefficients(segmentationFunction, edgeGradientPtr, CoefPtrs, VPtrs, segParameters);

		// Call to function that will evolve segmentation function in each discrete time step
		generalizedSubsurfSegmentationTimeStep(prevSol_extPtr, gauss_seidelPtr, segmentationFunction, segParameters, CoefPtrs);

		//Compute the L2 norm of the difference between the current and previous solutions
		difference_btw_current_and_previous_sol = l2normD(prevSol_extPtr, gauss_seidelPtr, length_ext, width_ext, height_ext, h);

		copyDataToAnotherArray(gauss_seidelPtr, prevSol_extPtr, height_ext, length_ext, width_ext);

		//writing density.
		if ((number_time_step % segParameters.mod) == 0)
		{
			strcpy_s(name, sizeof name, outputPathPtr);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%05zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			isFileSaved = manageFile(segmFuntionPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
			printf("Step : %zd\n", number_time_step);
			if (isFileSaved == false) {
				printf("The file was not saved\n");
				return false;
			}
			printf("Step %zd , residual = %e \n", number_time_step, difference_btw_current_and_previous_sol);
			//fprintf(error_file, "%d,%f\n", number_time_step, difference_btw_current_and_previous_sol);
		}
	} while ((number_time_step <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));
	

	//fclose(error_file);

	for (i = 0; i < height; i++)
	{
		free(imageToBeSegPtr[i]);
		free(segmFuntionPtr[i]);
		free(edgeGradientPtr[i]);

		free(e_Ptr[i]);
		free(w_Ptr[i]);
		free(n_Ptr[i]);
		free(s_Ptr[i]);
		free(t_Ptr[i]);
		free(b_Ptr[i]);

		free(VePtr[i]);
		free(VwPtr[i]);
		free(VnPtr[i]);
		free(VsPtr[i]);
		free(VtPtr[i]);
		free(VbPtr[i]);
	}
	free(imageToBeSegPtr);
	free(segmFuntionPtr);
	free(edgeGradientPtr);

	free(e_Ptr);
	free(w_Ptr);
	free(n_Ptr);
	free(s_Ptr);
	free(t_Ptr);
	free(b_Ptr);

	free(VePtr);
	free(VwPtr);
	free(VnPtr);
	free(VsPtr);
	free(VtPtr);
	free(VbPtr);

	for (i = 0; i < height_ext; i++)
	{
		free(prevSol_extPtr[i]);
		free(gauss_seidelPtr[i]);
	}
	free(prevSol_extPtr);
	free(gauss_seidelPtr);

	return true;
}

bool generalizedGFunctionForImageToBeSegmented(Image_Data inputImageData, dataType** edgeGradientPtr, Gradient_Pointers VPtrs,
	Segmentation_Parameters segParameters, Filter_Parameters explicit_lhe_Parameters)
{
	//checks if the memory was allocated
	if (inputImageData.imageDataPtr == NULL || edgeGradientPtr == NULL || VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL
		|| VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
		return false;

	size_t i, j, k, x, x_ext;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t height = inputImageData.height;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t dim2D = length * width;
	size_t k_ext, j_ext, i_ext;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	dataType h = segParameters.h;
	dataType quotient = 4.0 * h;
	dataType ux, uy, uz;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW, Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW, //current and surrounding voxel values
		Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType norm_image_smoothed_e, norm_image_smoothed_w, norm_image_smoothed_n, norm_image_smoothed_s, norm_image_smoothed_t, norm_image_smoothed_b;
	dataType norm_image_smoothed_average;

	dataType** gradient_coef_ext = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) {
		gradient_coef_ext[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
	}
	if (gradient_coef_ext == NULL || extendedCoefPtr == NULL) 
		return false;

	dataType coef_conv = segParameters.coef_conv;

	// Initialize array
	initialize3dArrayD(gradient_coef_ext, length_ext, width_ext, height_ext, 0.0);
	initialize3dArrayD(extendedCoefPtr, length_ext, width_ext, height_ext, 0.0);

	////perfom presmoothing
	heatImplicitScheme(inputImageData, explicit_lhe_Parameters); // unconditionnally stable

	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = inputImageData.imageDataPtr[k][x_new(i, j, length)];
			}
		}
	}
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext);

	//calculation of coefficients
	for (k = 0, k_ext = 1; k < inputImageData.height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < inputImageData.length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < inputImageData.width; j++, j_ext++)
			{
				// 2D to 1D representation for i, j
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, inputImageData.length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)];
				uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = extendedCoefPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)];
				uW = extendedCoefPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)];
				uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = extendedCoefPtr[kminus1][x_new(i_ext, j_ext, length_ext)];
				TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = extendedCoefPtr[kminus1][x_new(i_ext + 1, j_ext, length_ext)];
				TuW = extendedCoefPtr[kminus1][x_new(i_ext - 1, j_ext, length_ext)];
				TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = extendedCoefPtr[kplus1][x_new(i_ext, j_ext, length_ext)];
				BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = extendedCoefPtr[kplus1][x_new(i_ext + 1, j_ext, length_ext)];
				BuW = extendedCoefPtr[kplus1][x_new(i_ext - 1, j_ext, length_ext)];
				BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the presmooted image data

				// Calculation of coefficients in east direction
				ux = (uE - u) / segParameters.h;
				uy = ((uN + uNE) - (uS + uSE)) / quotient;
				uz = ((Tu + TuE) - (Bu + BuE)) / quotient;
				dataType val_east = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in west direction
				ux = (uW - u) / segParameters.h;
				uy = ((uNW + uN) - (uSW + uS)) / quotient;
				uz = ((TuW + Tu) - (BuW + Bu)) / quotient;
				dataType val_west = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in north direction
				ux = ((uNE + uE) - (uNW + uW)) / quotient;
				uy = (uN - u) / segParameters.h;
				uz = ((TuN + Tu) - (BuN + Bu)) / quotient;
				dataType val_north = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in south direction
				ux = ((uE + uSE) - (uW + uSW)) / quotient;
				uy = (uS - u) / segParameters.h;
				uz = ((TuS + Tu) - (BuS + Bu)) / quotient;
				dataType val_south = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in top direction
				ux = ((TuE + uE) - (TuW + uW)) / quotient;
				uy = ((TuN + uN) - (TuS + uS)) / quotient;
				uz = (Tu - u) / segParameters.h;
				dataType val_top = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				// Calculation of coefficients in bottom direction
				ux = ((BuW + uW) - (BuE + uE)) / quotient;
				uy = ((BuN + uN) - (BuS + uS)) / quotient;
				uz = (Bu - u) / segParameters.h;
				dataType val_bottom = gradientFunction((ux * ux) + (uy * uy) + (uz * uz), segParameters.coef);

				dataType val_average = (val_east + val_west + val_north + val_south + val_top + val_bottom) / 6.0;
				edgeGradientPtr[k][x] = gradientFunction(val_average * val_average, segParameters.coef);
				gradient_coef_ext[k_ext][x_new(i_ext, j_ext, length_ext)] = edgeGradientPtr[k][x];
			}
		}
	}
	reflection3D(gradient_coef_ext, height_ext, length_ext, width_ext);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				VPtrs.GePtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(iplus1, j_ext, length_ext)] - gradient_coef_ext[k_ext][x_new(iminus1, j_ext, length_ext)]) / (2 * h);
				VPtrs.GwPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(iminus1, j_ext, length_ext)] - gradient_coef_ext[k_ext][x_new(iplus1, j_ext, length_ext)]) / (2 * h);
				VPtrs.GnPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(i_ext, jminus1, length_ext)] - gradient_coef_ext[k_ext][x_new(i_ext, jplus1, length_ext)]) / (2 * h);
				VPtrs.GsPtr[k][x] = -coef_conv * (gradient_coef_ext[k_ext][x_new(i_ext, jplus1, length_ext)] - gradient_coef_ext[k_ext][x_new(i_ext, jminus1, length_ext)]) / (2 * h);
				VPtrs.GtPtr[k][x] = -coef_conv * (gradient_coef_ext[kminus1][x_ext] - gradient_coef_ext[kplus1][x_ext]) / (2 * h);
				VPtrs.GbPtr[k][x] = -coef_conv * (gradient_coef_ext[kplus1][x_ext] - gradient_coef_ext[kminus1][x_ext]) / (2 * h);
			}
		}
	}

	for (k = 0; k < height_ext; k++) {
		free(gradient_coef_ext[k]);
		free(extendedCoefPtr[k]);
	}
	free(gradient_coef_ext);
	free(extendedCoefPtr);

	return true;
}

bool generalizedGaussSeidelCoefficients(Image_Data segmentationData, dataType** edgeGradientPtr, Coefficient_Pointers CoefPtrs, Gradient_Pointers VPtrs, Segmentation_Parameters segParameters)
{
	//checks if the memory was allocated
	if (segmentationData.imageDataPtr == NULL || edgeGradientPtr == NULL
		|| CoefPtrs.w_Ptr == NULL || CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL
		|| VPtrs.GePtr == NULL || VPtrs.GwPtr == NULL || VPtrs.GnPtr == NULL || VPtrs.GsPtr == NULL || VPtrs.GtPtr == NULL || VPtrs.GbPtr == NULL)
		return false;

	size_t i, j, k, x, x_ext;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;
	size_t dim2D = segmentationData.length * segmentationData.width;
	dataType coef_dif = segParameters.coef_dif;

	size_t k_ext, j_ext, i_ext;
	size_t height = segmentationData.height, length = segmentationData.length, width = segmentationData.width;

	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	dataType h = segParameters.h, quotient = (dataType)(4.0 * segParameters.h);
	dataType orig_ux, orig_uy, orig_uz; //change in x, y and z respectively
	dataType orig_u, orig_uN, orig_uS, orig_uE, orig_uW, orig_uNW, orig_uNE, orig_uSE, orig_uSW, orig_Tu, orig_TuN, orig_TuS,
		orig_TuE, orig_TuW, orig_TuNW, orig_TuNE, orig_TuSE, orig_TuSW, //current and surrounding voxel values
		orig_Bu, orig_BuN, orig_BuS, orig_BuE, orig_BuW, orig_BuNW, orig_BuNE, orig_BuSE, orig_BuSW;
	dataType orig_e, orig_w, orig_n, orig_s, orig_t, orig_b;
	dataType voxel_coef, average_face_coef;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) {
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
	}
	if (extendedCoefPtr == NULL)
		return false;

	////copy data to extended area which will be used in each time step
	//copyDataToExtendedArea(inputImageData.segmentationFuntionPtr, extendedCoefPtr, height, length, width);
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = segmentationData.imageDataPtr[k][x_new(i, j, length)];
			}
		}
	}

	setBoundaryToZeroDirichletBC(extendedCoefPtr, length_ext, width_ext, height_ext); // In this case, It should work as initializing the array to 0
	
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

				//values of voxels in the extended data container for the original image
				orig_u = extendedCoefPtr[k_ext][x_ext];
				orig_uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				orig_uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				orig_uE = extendedCoefPtr[k_ext][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_uW = extendedCoefPtr[k_ext][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				orig_uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				orig_uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				orig_uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				orig_Tu = extendedCoefPtr[kminus1][x_ext];
				orig_TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				orig_TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				orig_TuE = extendedCoefPtr[kminus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_TuW = extendedCoefPtr[kminus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				orig_TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				orig_TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				orig_TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				orig_Bu = extendedCoefPtr[kplus1][x_ext];
				orig_BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				orig_BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				orig_BuE = extendedCoefPtr[kplus1][x_ext + 1];//x_new(i_ext + 1, j_ext, length_ext)
				orig_BuW = extendedCoefPtr[kplus1][x_ext - 1];//x_new(i_ext - 1, j_ext, length_ext)
				orig_BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				orig_BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				orig_BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				orig_BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				//calculation of coefficients in the original image data
				// Calculation of coefficients in east direction
				orig_ux = (orig_uE - orig_u) / h;
				orig_uy = ((orig_uN + orig_uNE) - (orig_uS + orig_uSE)) / quotient;
				orig_uz = ((orig_Tu + orig_TuE) - (orig_Bu + orig_BuE)) / quotient;
				orig_e = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in west direction
				orig_ux = (orig_uW - orig_u) / h;
				orig_uy = ((orig_uNW + orig_uN) - (orig_uSW + orig_uS)) / quotient;
				orig_uz = ((orig_TuW + orig_Tu) - (orig_BuW + orig_Bu)) / quotient;
				orig_w = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in north direction
				orig_ux = ((orig_uNE + orig_uE) - (orig_uNW + orig_uW)) / quotient;
				orig_uy = (orig_uN - orig_u) / h;
				orig_uz = ((orig_TuN + orig_Tu) - (orig_BuN + orig_Bu)) / quotient;
				orig_n = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in south direction
				orig_ux = ((orig_uE + orig_uSE) - (orig_uW + orig_uSW)) / quotient;
				orig_uy = (orig_uS - orig_u) / h;
				orig_uz = ((orig_TuS + orig_Tu) - (orig_BuS + orig_Bu)) / quotient;
				orig_s = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in top direction
				orig_ux = ((orig_TuE + orig_uE) - (orig_TuW + orig_uW)) / quotient;
				orig_uy = ((orig_TuN + orig_uN) - (orig_TuS + orig_uS)) / quotient;
				orig_uz = (orig_Tu - orig_u) / h;
				orig_t = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// Calculation of coefficients in bottom direction
				orig_ux = ((orig_BuW + orig_uW) - (orig_BuE + orig_uE)) / quotient;
				orig_uy = ((orig_BuN + orig_uN) - (orig_BuS + orig_uS)) / quotient;
				orig_uz = (orig_Bu - orig_u) / h;
				orig_b = (dataType)sqrt((orig_ux * orig_ux) + (orig_uy * orig_uy) + (orig_uz * orig_uz) + segParameters.eps2);

				// evaluation of norm of gradient of image at each voxel
				average_face_coef = (dataType)(((orig_e + orig_w + orig_n + orig_s + orig_t + orig_b) / 6.0));

				voxel_coef = (dataType)sqrt(pow(average_face_coef, 2) + segParameters.eps2);

				//evaluation of norm of gradient of image at each voxel, norm of gradient of presmoothed
				//image at each voxel face and reciprocal of norm of gradient of image at each voxel face
				CoefPtrs.e_Ptr[k][x] = (dataType)(-min(VPtrs.GePtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_e));
				CoefPtrs.w_Ptr[k][x] = (dataType)(-min(VPtrs.GwPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_w));
				CoefPtrs.n_Ptr[k][x] = (dataType)(-min(VPtrs.GnPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_n));
				CoefPtrs.s_Ptr[k][x] = (dataType)(-min(VPtrs.GsPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_s));
				CoefPtrs.t_Ptr[k][x] = (dataType)(-min(VPtrs.GtPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_t));
				CoefPtrs.b_Ptr[k][x] = (dataType)(-min(VPtrs.GbPtr[k][x],0) + coef_dif * voxel_coef * edgeGradientPtr[k][x] * (1.0 / orig_b));

			}
		}
	}
	
	for (i = 0; i < height_ext; i++) {
		free(extendedCoefPtr[i]);
	}
	free(extendedCoefPtr);

	return true;
}

bool generalizedSubsurfSegmentationTimeStep(dataType** prevSol_extPtr, dataType** gauss_seidelPtr, Image_Data segmentationData,
	Segmentation_Parameters segParameters, Coefficient_Pointers CoefPtrs)
{
	//check if the memory was allocated successfully
	if (segmentationData.imageDataPtr == NULL || prevSol_extPtr == NULL || gauss_seidelPtr == NULL || CoefPtrs.e_Ptr == NULL || CoefPtrs.w_Ptr == NULL
		|| CoefPtrs.n_Ptr == NULL || CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
		return false;

	size_t k, i, j;
	dataType hh = segParameters.h * segParameters.h * segParameters.h;
	dataType tau = segParameters.tau;

	// Error value used to check iteration
	// sor - successive over relation value, used in Gauss-Seidel formula
	dataType mean_square_residue = 0.0, gauss_seidel = 0.0;

	size_t height = segmentationData.height;
	size_t length = segmentationData.length;
	size_t width = segmentationData.width;
	size_t height_ext = height + 2;
	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t k_ext, j_ext, i_ext;
	size_t x;
	size_t x_ext;

	const dataType coef_tauh = tau / hh;
	dataType new_value = 0.0;

	// The Implicit Scheme Evaluation
	size_t count_step = 0;
	do
	{
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					// Gauss-Seidel Formula Evaluation
					//explicit for advection
					gauss_seidel = (dataType)(((prevSol_extPtr[k_ext][x_ext] + coef_tauh * (
					  CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]
					+ CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] 
					+ CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
					+ CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
					+ CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext] 
					+ CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext]))
					/ (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.s_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.b_Ptr[k][x] + CoefPtrs.t_Ptr[k][x]))));

					//SOR implementation using Gauss-Seidel
					gauss_seidelPtr[k_ext][x_ext] = gauss_seidelPtr[k_ext][x_ext] + segParameters.omega_c * (gauss_seidel - gauss_seidelPtr[k_ext][x_ext]);

				}
			}
		}
		
		// Error Evaluation
		mean_square_residue = 0.0; // Initialize
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					// 2D to 1D representation for i, j
					x_ext = x_new(i_ext, j_ext, length_ext);
					x = x_new(i, j, length);

					mean_square_residue += (dataType)(pow(gauss_seidelPtr[k_ext][x_ext] * (1 + coef_tauh * (CoefPtrs.e_Ptr[k][x] + CoefPtrs.w_Ptr[k][x] + CoefPtrs.s_Ptr[k][x] + CoefPtrs.n_Ptr[k][x] + CoefPtrs.b_Ptr[k][x] + CoefPtrs.t_Ptr[k][x]))
						- coef_tauh * (CoefPtrs.e_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext + 1] + CoefPtrs.w_Ptr[k][x] * gauss_seidelPtr[k_ext][x_ext - 1]
							+ CoefPtrs.s_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]
							+ CoefPtrs.n_Ptr[k][x] * gauss_seidelPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]
							+ CoefPtrs.b_Ptr[k][x] * gauss_seidelPtr[k_ext + 1][x_ext] + CoefPtrs.t_Ptr[k][x] * gauss_seidelPtr[k_ext - 1][x_ext])
						- prevSol_extPtr[k_ext][x_ext], 2));
				}
			}
		}
		
		count_step = count_step + 1;

	} while (mean_square_residue > segParameters.gauss_seidelTolerance && count_step < segParameters.maxNoGSIteration);

	//rescaling
	rescaleToIntervalZeroOne(gauss_seidelPtr, length_ext, width_ext, height_ext);

	////Copy the current time step to original data array after timeStepsNum
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				segmentationData.imageDataPtr[k][x_new(i, j, length)] = gauss_seidelPtr[k_ext][x_new(i_ext, j_ext, length_ext)];
			}
		}
	}

	return true;
}

//Functions for IIOE scheme

dataType getMinInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k)
{

	if (k == 0)
	{
		size_t kplus1 = k + 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uE, uSE)));
				dataType m2 = fmin(Bu, fmin(BuS, fmin(BuE, BuSE)));
				return fmin(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uE, uNE)));
				dataType m2 = fmin(Bu, fmin(BuN, fmin(BuE, BuNE)));
				return fmin(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, uE)));
				dataType m2 = fmin(uNE, fmin(uSE, fmin(Bu, BuN)));
				dataType m3 = fmin(BuS, fmin(BuE, fmin(BuNE, BuSE)));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uW, uSW)));
				dataType m2 = fmin(Bu, fmin(BuS, fmin(BuW, BuSW)));
				return fmin(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uW, uNW)));
				dataType m2 = fmin(Bu, fmin(BuN, fmin(BuW, BuNW)));
				return fmin(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, uW)));
				dataType m2 = fmin(uNW, fmin(uSW, fmin(Bu, BuN)));
				dataType m3 = fmin(BuS, fmin(BuW, fmin(BuNW, BuSW)));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uE, uW)));
				dataType m2 = fmin(uSE, fmin(uSW, fmin(Bu, BuS)));
				dataType m3 = fmin(BuE, fmin(BuW, fmin(BuSE, BuSW)));
				return fmin(m1, fmin(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uE, uW)));
				dataType m2 = fmin(uNW, fmin(uNE, fmin(Bu, BuN)));
				dataType m3 = fmin(BuE, fmin(BuW, fmin(BuNW, BuNE)));
				return fmin(m1, fmin(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, fmin(uE, fmin(uW, uNW)))));
				dataType m2 = fmin(uNE, fmin(uSE, fmin(uSW, fmin(Bu, fmin(BuN, BuS)))));
				dataType m3 = fmin(BuE, fmin(BuW, fmin(BuNW, fmin(BuNE, fmin(BuSE, BuSW)))));
				return fmin(m1, fmin(m2, m3));
			}
		}
	}
	else if (k == height - 1)
	{
		size_t kminus1 = k - 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uE, uSE)));
				dataType m2 = fmin(Tu, fmin(TuS, fmin(TuE, TuSE)));
				return fmin(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uE, uNE)));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuE, TuNE)));
				return fmin(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, uE)));
				dataType m2 = fmin(uNE, fmin(uSE, fmin(Tu, TuN)));
				dataType m3 = fmin(TuS, fmin(TuE, fmin(TuNE, TuSE)));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uW, uSW)));
				dataType m2 = fmin(Tu, fmin(TuS, fmin(TuW, TuSW)));
				return fmin(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uW, uNW)));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuW, TuNW)));
				return fmin(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, uW)));
				dataType m2 = fmin(uNW, fmin(uSW, fmin(Tu, TuN)));
				dataType m3 = fmin(TuS, fmin(TuW, fmin(TuNW, TuSW)));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uS, fmin(uE, uW)));
				dataType m2 = fmin(uSE, fmin(uSW, fmin(Tu, TuS)));
				dataType m3 = fmin(TuE, fmin(TuW, fmin(TuSE, TuSW)));
				return fmin(m1, fmin(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uE, uW)));
				dataType m2 = fmin(uNW, fmin(uNE, fmin(Tu, TuN)));
				dataType m3 = fmin(TuE, fmin(TuW, fmin(TuNW, TuNE)));
				return fmin(m1, fmin(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType m1 = fmin(u, fmin(uN, fmin(uS, fmin(uE, fmin(uW, uNW)))));
				dataType m2 = fmin(uNE, fmin(uSE, fmin(uSW, fmin(Tu, fmin(TuN, TuS)))));
				dataType m3 = fmin(TuE, fmin(TuW, fmin(TuNW, fmin(TuNE, fmin(TuSE, TuSW)))));
				return fmin(m1, fmin(m2, m3));
			}
		}
	}
	else
	{
		size_t kminus1 = k - 1;
		size_t kplus1 = k + 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uS, fmin(uE, uSE)));
				dataType m2 = fmin(Tu, fmin(TuS, fmin(TuE, TuSE)));
				dataType m3 = fmin(Bu, fmin(BuS, fmin(BuE, BuSE)));
				return fmin(m1, fmin(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uE, uNE)));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuE, TuNE)));
				dataType m3 = fmin(Bu, fmin(BuN, fmin(BuE, BuNE)));
				return fmin(m1, fmin(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uS, fmin(uE, fmin(uNE, uSE)))));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuS, fmin(TuE, fmin(TuNE, TuSE)))));
				dataType m3 = fmin(Bu, fmin(BuN, fmin(BuS, fmin(BuE, fmin(BuNE, BuSE)))));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;

			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uS, fmin(uW, uSW)));
				dataType m2 = fmin(Tu, fmin(TuS, fmin(TuW, TuSW)));
				dataType m3 = fmin(Bu, fmin(BuS, fmin(BuW, BuSW)));
				return fmin(m1, fmin(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uW, uNW)));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuW, TuNW)));
				dataType m3 = fmin(Bu, fmin(BuN, fmin(BuW, BuNW)));
				return fmin(m1, fmin(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uS, fmin(uW, fmin(uNW, uSW)))));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuS, fmin(TuW, fmin(TuNW, TuSW)))));
				dataType m3 = fmin(Bu, fmin(BuN, fmin(BuS, fmin(BuW, fmin(BuNW, BuSW)))));
				return fmin(m1, fmin(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;

			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uS, fmin(uE, fmin(uW, fmin(uSE, uSW)))));
				dataType m2 = fmin(Tu, fmin(TuS, fmin(TuE, fmin(TuW, fmin(TuSE, TuSW)))));
				dataType m3 = fmin(Bu, fmin(BuS, fmin(BuE, fmin(BuW, fmin(BuSE, BuSW)))));
				return fmin(m1, fmin(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uE, fmin(uW, fmin(uNW, uNE)))));
				dataType m2 = fmin(Tu, fmin(TuN, fmin(TuE, fmin(TuW, fmin(TuNW, TuNE)))));
				dataType m3 = fmin(Bu, fmin(BuN, fmin(BuE, fmin(BuW, fmin(BuNW, BuNE)))));
				return fmin(m1, fmin(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmin(u, fmin(uN, fmin(uS, fmin(uE, uW))));
				dataType m2 = fmin(uNW, fmin(uNE, fmin(uSE, uSW)));
				dataType m3 = fmin(Tu, fmin(TuN, fmin(TuS, fmin(TuE, TuW))));
				dataType m4 = fmin(TuNW, fmin(TuNE, fmin(TuSE, TuSW)));
				dataType m5 = fmin(Bu, fmin(BuN, fmin(BuS, fmin(BuE, BuW))));
				dataType m6 = fmin(BuNW, fmin(BuNE, fmin(BuSE, BuSW)));
				return fmin(m1, fmin(m2, fmin(m3, fmin(m4, fmin(m5, m6)))));
			}
		}
	}
}

dataType getMaxInNeighborhood3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k)
{

	if (k == 0)
	{
		size_t kplus1 = k + 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, uSE)));
				dataType m2 = fmax(Bu, fmax(BuS, fmax(BuE, BuSE)));
				return fmax(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, uNE)));
				dataType m2 = fmax(Bu, fmax(BuN, fmax(BuE, BuNE)));
				return fmax(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, uE)));
				dataType m2 = fmax(uNE, fmax(uSE, fmax(Bu, BuN)));
				dataType m3 = fmax(BuS, fmax(BuE, fmax(BuNE, BuSE)));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uW, uSW)));
				dataType m2 = fmax(Bu, fmax(BuS, fmax(BuW, BuSW)));
				return fmax(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uW, uNW)));
				dataType m2 = fmax(Bu, fmax(BuN, fmax(BuW, BuNW)));
				return fmax(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, uW)));
				dataType m2 = fmax(uNW, fmax(uSW, fmax(Bu, BuN)));
				dataType m3 = fmax(BuS, fmax(BuW, fmax(BuNW, BuSW)));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, uW)));
				dataType m2 = fmax(uSE, fmax(uSW, fmax(Bu, BuS)));
				dataType m3 = fmax(BuE, fmax(BuW, fmax(BuSE, BuSW)));
				return fmax(m1, fmax(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, uW)));
				dataType m2 = fmax(uNW, fmax(uNE, fmax(Bu, BuN)));
				dataType m3 = fmax(BuE, fmax(BuW, fmax(BuNW, BuNE)));
				return fmax(m1, fmax(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, fmax(uE, fmax(uW, uNW)))));
				dataType m2 = fmax(uNE, fmax(uSE, fmax(uSW, fmax(Bu, fmax(BuN, BuS)))));
				dataType m3 = fmax(BuE, fmax(BuW, fmax(BuNW, fmax(BuNE, fmax(BuSE, BuSW)))));
				return fmax(m1, fmax(m2, m3));
			}
		}
	}
	else if (k == height - 1)
	{
		size_t kminus1 = k - 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, uSE)));
				dataType m2 = fmax(Tu, fmax(TuS, fmax(TuE, TuSE)));
				return fmax(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;
				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, uNE)));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuE, TuNE)));
				return fmax(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, uE)));
				dataType m2 = fmax(uNE, fmax(uSE, fmax(Tu, TuN)));
				dataType m3 = fmax(TuS, fmax(TuE, fmax(TuNE, TuSE)));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uW, uSW)));
				dataType m2 = fmax(Tu, fmax(TuS, fmax(TuW, TuSW)));
				return fmax(m1, m2);
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uW, uNW)));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuW, TuNW)));
				return fmax(m1, m2);
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, uW)));
				dataType m2 = fmax(uNW, fmax(uSW, fmax(Tu, TuN)));
				dataType m3 = fmax(TuS, fmax(TuW, fmax(TuNW, TuSW)));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, uW)));
				dataType m2 = fmax(uSE, fmax(uSW, fmax(Tu, TuS)));
				dataType m3 = fmax(TuE, fmax(TuW, fmax(TuSE, TuSW)));
				return fmax(m1, fmax(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, uW)));
				dataType m2 = fmax(uNW, fmax(uNE, fmax(Tu, TuN)));
				dataType m3 = fmax(TuE, fmax(TuW, fmax(TuNW, TuNE)));
				return fmax(m1, fmax(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, fmax(uE, fmax(uW, uNW)))));
				dataType m2 = fmax(uNE, fmax(uSE, fmax(uSW, fmax(Tu, fmax(TuN, TuS)))));
				dataType m3 = fmax(TuE, fmax(TuW, fmax(TuNW, fmax(TuNE, fmax(TuSE, TuSW)))));
				return fmax(m1, fmax(m2, m3));
			}
		}
	}
	else
	{
		size_t kminus1 = k - 1;
		size_t kplus1 = k + 1;
		if (i == 0)
		{
			size_t iplus1 = i + 1;
			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, uSE)));
				dataType m2 = fmax(Tu, fmax(TuS, fmax(TuE, TuSE)));
				dataType m3 = fmax(Bu, fmax(BuS, fmax(BuE, BuSE)));
				return fmax(m1, fmax(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, uNE)));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuE, TuNE)));
				dataType m3 = fmax(Bu, fmax(BuN, fmax(BuE, BuNE)));
				return fmax(m1, fmax(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, fmax(uE, fmax(uNE, uSE)))));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuS, fmax(TuE, fmax(TuNE, TuSE)))));
				dataType m3 = fmax(Bu, fmax(BuN, fmax(BuS, fmax(BuE, fmax(BuNE, BuSE)))));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else if (i == length - 1)
		{
			size_t iminus1 = i - 1;

			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uW, uSW)));
				dataType m2 = fmax(Tu, fmax(TuS, fmax(TuW, TuSW)));
				dataType m3 = fmax(Bu, fmax(BuS, fmax(BuW, BuSW)));
				return fmax(m1, fmax(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uW, uNW)));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuW, TuNW)));
				dataType m3 = fmax(Bu, fmax(BuN, fmax(BuW, BuNW)));
				return fmax(m1, fmax(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, fmax(uW, fmax(uNW, uSW)))));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuS, fmax(TuW, fmax(TuNW, TuSW)))));
				dataType m3 = fmax(Bu, fmax(BuN, fmax(BuS, fmax(BuW, fmax(BuNW, BuSW)))));
				return fmax(m1, fmax(m2, m3));
			}
		}
		else
		{
			size_t iminus1 = i - 1;
			size_t iplus1 = i + 1;

			if (j == 0)
			{
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uS, fmax(uE, fmax(uW, fmax(uSE, uSW)))));
				dataType m2 = fmax(Tu, fmax(TuS, fmax(TuE, fmax(TuW, fmax(TuSE, TuSW)))));
				dataType m3 = fmax(Bu, fmax(BuS, fmax(BuE, fmax(BuW, fmax(BuSE, BuSW)))));
				return fmax(m1, fmax(m2, m3));
			}
			else if (j == width - 1)
			{
				size_t jminus1 = j - 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uE, fmax(uW, fmax(uNW, uNE)))));
				dataType m2 = fmax(Tu, fmax(TuN, fmax(TuE, fmax(TuW, fmax(TuNW, TuNE)))));
				dataType m3 = fmax(Bu, fmax(BuN, fmax(BuE, fmax(BuW, fmax(BuNW, BuNE)))));
				return fmax(m1, fmax(m2, m3));
			}
			else
			{
				size_t jminus1 = j - 1;
				size_t jplus1 = j + 1;

				dataType u = imageDataPtr[k][x_new(i, j, length)];
				dataType uN = imageDataPtr[k][x_new(i, jminus1, length)];
				dataType uS = imageDataPtr[k][x_new(i, jplus1, length)];
				dataType uE = imageDataPtr[k][x_new(iplus1, j, length)];
				dataType uW = imageDataPtr[k][x_new(iminus1, j, length)];
				dataType uNW = imageDataPtr[k][x_new(iminus1, jminus1, length)];
				dataType uNE = imageDataPtr[k][x_new(iplus1, jminus1, length)];
				dataType uSE = imageDataPtr[k][x_new(iplus1, jplus1, length)];
				dataType uSW = imageDataPtr[k][x_new(iminus1, jplus1, length)];
				dataType Tu = imageDataPtr[kminus1][x_new(i, j, length)];
				dataType TuN = imageDataPtr[kminus1][x_new(i, jminus1, length)];
				dataType TuS = imageDataPtr[kminus1][x_new(i, jplus1, length)];
				dataType TuE = imageDataPtr[kminus1][x_new(iplus1, j, length)];
				dataType TuW = imageDataPtr[kminus1][x_new(iminus1, j, length)];
				dataType TuNW = imageDataPtr[kminus1][x_new(iminus1, jminus1, length)];
				dataType TuNE = imageDataPtr[kminus1][x_new(iplus1, jminus1, length)];
				dataType TuSE = imageDataPtr[kminus1][x_new(iplus1, jplus1, length)];
				dataType TuSW = imageDataPtr[kminus1][x_new(iminus1, jplus1, length)];
				dataType Bu = imageDataPtr[kplus1][x_new(i, j, length)];
				dataType BuN = imageDataPtr[kplus1][x_new(i, jminus1, length)];
				dataType BuS = imageDataPtr[kplus1][x_new(i, jplus1, length)];
				dataType BuE = imageDataPtr[kplus1][x_new(iplus1, j, length)];
				dataType BuW = imageDataPtr[kplus1][x_new(iminus1, j, length)];
				dataType BuNW = imageDataPtr[kplus1][x_new(iminus1, jminus1, length)];
				dataType BuNE = imageDataPtr[kplus1][x_new(iplus1, jminus1, length)];
				dataType BuSE = imageDataPtr[kplus1][x_new(iplus1, jplus1, length)];
				dataType BuSW = imageDataPtr[kplus1][x_new(iminus1, jplus1, length)];

				dataType m1 = fmax(u, fmax(uN, fmax(uS, fmax(uE, uW))));
				dataType m2 = fmax(uNW, fmax(uNE, fmax(uSE, uSW)));
				dataType m3 = fmax(Tu, fmax(TuN, fmax(TuS, fmax(TuE, TuW))));
				dataType m4 = fmax(TuNW, fmax(TuNE, fmax(TuSE, TuSW)));
				dataType m5 = fmax(Bu, fmax(BuN, fmax(BuS, fmax(BuE, BuW))));
				dataType m6 = fmax(BuNW, fmax(BuNE, fmax(BuSE, BuSW)));
				return fmax(m1, fmax(m2, fmax(m3, fmax(m4, fmax(m5, m6)))));
			}
		}
	}
}

bool computeNormOfGradientDiamondCell3D(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Coefficient_Pointers nGrad)
{
	if(imageData == NULL || nGrad.e_Ptr == NULL || nGrad.w_Ptr == NULL ||
		nGrad.s_Ptr == NULL || nGrad.n_Ptr == NULL ||
		nGrad.t_Ptr == NULL || nGrad.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}

	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;

	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t height_ext = height + 2;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext );
	for( k = 0; k < height_ext; k++)
	{
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if(extendedCoefPtr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}
	if(extendedCoefPtr == NULL)
	{
		return false; // Memory allocation failed
	}

	//Copy to extended area
	for(k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for(i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for(j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = imageData[k][x_new(i, j, length)];
			}
		}
	}
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext); // Reflect the data in the extended area

	// Calculate the norm of gradient in diamond cell
	size_t iminus1, iplus1, jminus1, jplus1, kminus1, kplus1;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW;
	dataType Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW;
	dataType Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;
	dataType quotient = (dataType)(4.0 * h);
	for(k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x_ext = x_new(i_ext, j_ext, length_ext);
				x = x_new(i, j, length);
				iminus1 = i_ext - 1;
				iplus1 = i_ext + 1;
				jplus1 = j_ext + 1;
				jminus1 = j_ext - 1;
				kplus1 = k_ext + 1;
				kminus1 = k_ext - 1;

				//values of voxels in the extended data container for presmoothed image
				u = extendedCoefPtr[k_ext][x_ext];
				uN = extendedCoefPtr[k_ext][x_new(i_ext, jminus1, length_ext)];
				uS = extendedCoefPtr[k_ext][x_new(i_ext, jplus1, length_ext)];
				uE = extendedCoefPtr[k_ext][x_new(iplus1, j_ext, length_ext)];
				uW = extendedCoefPtr[k_ext][x_new(iminus1, j_ext, length_ext)];
				uNW = extendedCoefPtr[k_ext][x_new(iminus1, jminus1, length_ext)];
				uNE = extendedCoefPtr[k_ext][x_new(iplus1, jminus1, length_ext)];
				uSE = extendedCoefPtr[k_ext][x_new(iplus1, jplus1, length_ext)];
				uSW = extendedCoefPtr[k_ext][x_new(iminus1, jplus1, length_ext)];
				Tu = extendedCoefPtr[kminus1][x_ext];
				TuN = extendedCoefPtr[kminus1][x_new(i_ext, jminus1, length_ext)];
				TuS = extendedCoefPtr[kminus1][x_new(i_ext, jplus1, length_ext)];
				TuE = extendedCoefPtr[kminus1][x_new(iplus1, j_ext, length_ext)];
				TuW = extendedCoefPtr[kminus1][x_new(iminus1, j_ext, length_ext)];
				TuNW = extendedCoefPtr[kminus1][x_new(iminus1, jminus1, length_ext)];
				TuNE = extendedCoefPtr[kminus1][x_new(iplus1, jminus1, length_ext)];
				TuSE = extendedCoefPtr[kminus1][x_new(iplus1, jplus1, length_ext)];
				TuSW = extendedCoefPtr[kminus1][x_new(iminus1, jplus1, length_ext)];
				Bu = extendedCoefPtr[kplus1][x_ext];
				BuN = extendedCoefPtr[kplus1][x_new(i_ext, jminus1, length_ext)];
				BuS = extendedCoefPtr[kplus1][x_new(i_ext, jplus1, length_ext)];
				BuE = extendedCoefPtr[kplus1][x_new(iplus1, j_ext, length_ext)];
				BuW = extendedCoefPtr[kplus1][x_new(iminus1, j_ext, length_ext)];
				BuNW = extendedCoefPtr[kplus1][x_new(iminus1, jminus1, length_ext)];
				BuNE = extendedCoefPtr[kplus1][x_new(iplus1, jminus1, length_ext)];
				BuSE = extendedCoefPtr[kplus1][x_new(iplus1, jplus1, length_ext)];
				BuSW = extendedCoefPtr[kplus1][x_new(iminus1, jplus1, length_ext)];

				// Calculation of coefficients in East direction
				ux = (uE - u) / h;
				uy = ((uN + uNE) - (uS + uSE)) / quotient;
				uz = ((Tu + TuE) - (Bu + BuE)) / quotient;
				nGrad.e_Ptr[k][x] = ux * ux + uy * uy + uz * uz;

				// Calculation of coefficients in West direction
				ux = (uW - u) / h;
				uy = ((uNW + uN) - (uSW + uS)) / quotient;
				uz = ((TuW + Tu) - (BuW + Bu)) / quotient;
				nGrad.w_Ptr[k][x] = ux * ux + uy * uy + uz * uz;

				// Calculation of coefficients in North direction
				ux = ((uNE + uE) - (uNW + uW)) / quotient;
				uy = (uN - u) / h;
				uz = ((TuN + Tu) - (BuN + Bu)) / quotient;
				nGrad.n_Ptr[k][x] = ux * ux + uy * uy + uz * uz;

				// Calculation of coefficients in South direction
				ux = ((uE + uSE) - (uW + uSW)) / quotient;
				uy = (uS - u) / h;
				uz = ((TuS + Tu) - (BuS + Bu)) / quotient;
				nGrad.s_Ptr[k][x] = ux * ux + uy * uy + uz * uz;

				// Calculation of coefficients in Top direction
				ux = ((TuE + uE) - (TuW + uW)) / quotient;
				uy = ((TuN + uN) - (TuS + uS)) / quotient;
				uz = (Tu - u) / h;
				nGrad.t_Ptr[k][x] = ux * ux + uy * uy + uz * uz;

				// Calculation of coefficients in Bottom direction
				ux = ((BuW + uW) - (BuE + uE)) / quotient;
				uy = ((BuN + uN) - (BuS + uS)) / quotient;
				uz = (Bu - u) / h;
				nGrad.b_Ptr[k][x] = ux * ux + uy * uy + uz * uz;
			}
		}
	}

	for(k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);

	return true;
}

bool generalizedSubsurf_iioe(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;
	size_t kplus1, kminus1, iminus1, iplus1, jminus1, jplus1;

	VoxelSpacing spacing = imageData.spacing;
	
	const size_t height = imageData.height;
	const size_t height_ext = height + 2;

	const size_t length = imageData.length;
	const size_t length_ext = length + 2;
	
	const size_t width = imageData.width;
	const size_t width_ext = width + 2;
	
	size_t dim2D = height * width;
	size_t dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;

	dataType** segmentationPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeDetectorPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++) 
	{
		if(k < height) 
		{
			segmentationPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
			edgeDetectorPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
			if (segmentationPtr[k] == NULL || edgeDetectorPtr[k] == NULL) 
			{
				return false; // Memory allocation failed
			}
		}
		gaussSeidelPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		previousSolPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	}
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}

	Coefficient_Pointers gPtrs;
	gPtrs.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for (k = 0; k < height; k++)
	{
		gPtrs.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.e_Ptr[k] == NULL || gPtrs.w_Ptr[k] == NULL || gPtrs.n_Ptr[k] == NULL ||
			gPtrs.s_Ptr[k] == NULL || gPtrs.t_Ptr[k] == NULL || gPtrs.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers segPtrs;
	segPtrs.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for (k = 0; k < height; k++)
	{
		segPtrs.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segPtrs.e_Ptr[k] == NULL || segPtrs.w_Ptr[k] == NULL || segPtrs.n_Ptr[k] == NULL ||
			segPtrs.s_Ptr[k] == NULL || segPtrs.t_Ptr[k] == NULL || segPtrs.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	//Intitialization
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			segmentationPtr[k][i] = 0.0;
			edgeDetectorPtr[k][i] = 0.0;
			gPtrs.e_Ptr[k][i] = 0.0;
			gPtrs.w_Ptr[k][i] = 0.0;
			gPtrs.n_Ptr[k][i] = 0.0;
			gPtrs.s_Ptr[k][i] = 0.0;
			gPtrs.t_Ptr[k][i] = 0.0;
			gPtrs.b_Ptr[k][i] = 0.0;
			segPtrs.e_Ptr[k][i] = 0.0;
			segPtrs.w_Ptr[k][i] = 0.0;
			segPtrs.n_Ptr[k][i] = 0.0;
			segPtrs.s_Ptr[k][i] = 0.0;
			segPtrs.t_Ptr[k][i] = 0.0;
			segPtrs.b_Ptr[k][i] = 0.0;
		}
	}

	/*
	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);
	heatImplicitScheme(imageData, smooth_parms);

	//compute the morm of gradient for the edge detector
	computeNormOfGradientDiamondCell3D(imageData.imageDataPtr, length, width, height, h, gPtrs);

	//compute the edge detector
	dataType value_gF_e, value_gF_w, value_gF_n, value_gF_s, value_gF_t, value_gF_b, average_value;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				x = x_new(i, j, length);
				value_gF_e = gradientFunction(gPtrs.e_Ptr[k][x], coef_edge_detector);
				value_gF_w = gradientFunction(gPtrs.w_Ptr[k][x], coef_edge_detector);
				value_gF_n = gradientFunction(gPtrs.n_Ptr[k][x], coef_edge_detector);
				value_gF_s = gradientFunction(gPtrs.s_Ptr[k][x], coef_edge_detector);
				value_gF_t = gradientFunction(gPtrs.t_Ptr[k][x], coef_edge_detector);
				value_gF_b = gradientFunction(gPtrs.b_Ptr[k][x], coef_edge_detector);
				average_value = (value_gF_e + value_gF_w + value_gF_n + value_gF_s + value_gF_t + value_gF_b) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(average_value * average_value, coef_edge_detector);
			}
		}
	}
	*/

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_new.raw");
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_old.raw");
	strcat_s(name, sizeof(name), name_ending);
	//store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);
	load3dDataArrayD(edgeDetectorPtr, length, width, height, name);

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for (k = 0; k < height; k++)
	{
		uCoef.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		uCoef.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		uCoef.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		uCoef.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		uCoef.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		uCoef.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (uCoef.e_Ptr[k] == NULL || uCoef.w_Ptr[k] == NULL || uCoef.n_Ptr[k] == NULL ||
			uCoef.s_Ptr[k] == NULL || uCoef.t_Ptr[k] == NULL || uCoef.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}
	
	Coefficient_Pointers coefPtrs;
	coefPtrs.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPtrs.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPtrs.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPtrs.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPtrs.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPtrs.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for(k  = 0; k < height; k++) 
	{
		coefPtrs.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPtrs.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPtrs.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPtrs.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPtrs.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPtrs.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (coefPtrs.e_Ptr[k] == NULL || coefPtrs.w_Ptr[k] == NULL || coefPtrs.n_Ptr[k] == NULL ||
			coefPtrs.s_Ptr[k] == NULL || coefPtrs.t_Ptr[k] == NULL || coefPtrs.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers a_out;
	a_out.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_out.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_out.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_out.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_out.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_out.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for(k = 0; k < height; k++) 
	{
		a_out.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_out.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_out.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_out.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_out.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_out.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (a_out.e_Ptr[k] == NULL || a_out.w_Ptr[k] == NULL || a_out.n_Ptr[k] == NULL ||
			a_out.s_Ptr[k] == NULL || a_out.t_Ptr[k] == NULL || a_out.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers a_in;
	a_in.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_in.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_in.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_in.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_in.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	a_in.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for(k = 0; k < height; k++) 
	{
		a_in.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_in.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_in.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_in.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_in.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		a_in.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (a_in.e_Ptr[k] == NULL || a_in.w_Ptr[k] == NULL || a_in.n_Ptr[k] == NULL ||
			a_in.s_Ptr[k] == NULL || a_in.t_Ptr[k] == NULL || a_in.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers theta_out;
	theta_out.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_out.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_out.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_out.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_out.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_out.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for(k = 0; k < height; k++) 
	{
		theta_out.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_out.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_out.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_out.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_out.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_out.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (theta_out.e_Ptr[k] == NULL || theta_out.w_Ptr[k] == NULL || theta_out.n_Ptr[k] == NULL ||
			theta_out.s_Ptr[k] == NULL || theta_out.t_Ptr[k] == NULL || theta_out.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers theta_in;
	theta_in.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_in.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_in.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_in.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_in.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	theta_in.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	for (k = 0; k < height; k++)
	{
		theta_in.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_in.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_in.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_in.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_in.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		theta_in.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (theta_in.e_Ptr[k] == NULL || theta_in.w_Ptr[k] == NULL || theta_in.n_Ptr[k] == NULL ||
			theta_in.s_Ptr[k] == NULL || theta_in.t_Ptr[k] == NULL || theta_in.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	dataType** n_out_pq = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** n_out_qp = (dataType**)malloc(sizeof(dataType*) * height);
	for (k = 0; k < height; k++) 
	{
		n_out_pq[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		n_out_qp[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (n_out_pq[k] == NULL || n_out_qp[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}
	if(n_out_pq == NULL || n_out_qp == NULL)
	{
		return false;
	}

	//Initialization
	for(k = 0; k < height; k++)
	{
		for(i = 0; i < dim2D; i++)
		{
			uCoef.e_Ptr[k][i] = 0.0;
			uCoef.w_Ptr[k][i] = 0.0;
			uCoef.n_Ptr[k][i] = 0.0;
			uCoef.s_Ptr[k][i] = 0.0;
			uCoef.t_Ptr[k][i] = 0.0;
			uCoef.b_Ptr[k][i] = 0.0;
			coefPtrs.e_Ptr[k][i] = 0.0;
			coefPtrs.w_Ptr[k][i] = 0.0;
			coefPtrs.n_Ptr[k][i] = 0.0;
			coefPtrs.s_Ptr[k][i] = 0.0;
			coefPtrs.t_Ptr[k][i] = 0.0;
			coefPtrs.b_Ptr[k][i] = 0.0;
			a_out.e_Ptr[k][i] = 0.0;
			a_out.w_Ptr[k][i] = 0.0;
			a_out.n_Ptr[k][i] = 0.0;
			a_out.s_Ptr[k][i] = 0.0;
			a_out.t_Ptr[k][i] = 0.0;
			a_out.b_Ptr[k][i] = 0.0;
			a_in.e_Ptr[k][i] = 0.0;
			a_in.w_Ptr[k][i] = 0.0;
			a_in.n_Ptr[k][i] = 0.0;
			a_in.s_Ptr[k][i] = 0.0;
			a_in.t_Ptr[k][i] = 0.0;
			a_in.b_Ptr[k][i] = 0.0;
			theta_out.e_Ptr[k][i] = 0.0;
			theta_out.w_Ptr[k][i] = 0.0;
			theta_out.n_Ptr[k][i] = 0.0;
			theta_out.s_Ptr[k][i] = 0.0;
			theta_out.t_Ptr[k][i] = 0.0;
			theta_out.b_Ptr[k][i] = 0.0;
			theta_in.e_Ptr[k][i] = 0.0;
			theta_in.w_Ptr[k][i] = 0.0;
			theta_in.n_Ptr[k][i] = 0.0;
			theta_in.s_Ptr[k][i] = 0.0;
			theta_in.t_Ptr[k][i] = 0.0;
			theta_in.b_Ptr[k][i] = 0.0;
			n_out_pq[k][i] = 0.0;
			n_out_qp[k][i] = 0.0;
		}
	}

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	for (k = 0; k < height; k++) 
	{
		for (i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				x = x_new(i, j, length);

				if (i == 0) 
				{
					vpe = adv * (edgeDetectorPtr[k][x_new(i + 1, j, length)] - edgeDetectorPtr[k][x]);
					vpw = -vpe;
				}
				else if (i == length - 1) 
				{
					vpw = adv * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k][x_new(i - 1, j, length)]);
					vpe = -vpw;					
				}
				else {
					vpe = adv * 0.5 * (edgeDetectorPtr[k][x_new(i + 1, j, length)] - edgeDetectorPtr[k][x_new(i - 1, j, length)]);
					vpw = -vpe;
				}

				if (j == 0) 
				{
					vps = adv * (edgeDetectorPtr[k][x_new(i, j + 1, length)] - edgeDetectorPtr[k][x]);
					vpn = -vps;
				}
				else if (j == width - 1) 
				{
					vpn = adv * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k][x_new(i, j - 1, length)]);
					vps = -vpn;
				}
				else 
				{
					vps = adv * 0.5 * (edgeDetectorPtr[k][x_new(i, j + 1, length)] - edgeDetectorPtr[k][x_new(i, j - 1, length)]);
					vpn = -vps;
				}
				
				if (k == 0) 
				{
					vpb = adv * (edgeDetectorPtr[k + 1][x] - edgeDetectorPtr[k][x]);
					vpt = -vpb;
				}
				else if (k == height - 1)
				{
					vpb = adv * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k - 1][x]);
					vpt = -vpb;					
				}
				else {
					vpb = adv * 0.5 * (edgeDetectorPtr[k + 1][x] - edgeDetectorPtr[k - 1][x]);
					vpt = -vpb;
				}

				a_in.e_Ptr[k][x] = fmax(vpe, 0);
				a_in.w_Ptr[k][x] = fmax(vpw, 0);
				a_in.n_Ptr[k][x] = fmax(vpn, 0);
				a_in.s_Ptr[k][x] = fmax(vps, 0);
				a_in.t_Ptr[k][x] = fmax(vpt, 0);
				a_in.b_Ptr[k][x] = fmax(vpb, 0);

				a_out.e_Ptr[k][x] = fmin(vpe, 0);
				a_out.w_Ptr[k][x] = fmin(vpw, 0);
				a_out.n_Ptr[k][x] = fmin(vpn, 0);
				a_out.s_Ptr[k][x] = fmin(vps, 0);
				a_out.t_Ptr[k][x] = fmin(vpt, 0);
				a_out.b_Ptr[k][x] = fmin(vpb, 0);

				n_out_pq[k][x] = -(signum(a_out.e_Ptr[k][x]) + signum(a_out.w_Ptr[k][x]) + signum(a_out.n_Ptr[k][x]) + signum(a_out.s_Ptr[k][x]) + signum(a_out.t_Ptr[k][x]) + signum(a_out.b_Ptr[k][x]));

				//a_out_pq = -a_in_pq
				n_out_qp[k][x] = -(signum(-a_in.e_Ptr[k][x]) + signum(-a_in.w_Ptr[k][x]) + signum(-a_in.n_Ptr[k][x]) + signum(-a_in.s_Ptr[k][x]) + signum(-a_in.t_Ptr[k][x]) + signum(-a_in.b_Ptr[k][x]));
			}
		}
	}
	
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				segmentationPtr[k][x] = initialSegment[k][x];
				previousSolPtr[k_ext][x_ext] = initialSegment[k][x];
				gaussSeidelPtr[k_ext][x_ext] = initialSegment[k][x];
			}
		}
	}
	setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p_min = 0.0, u_p_max = 0.0;
	dataType u_p = 0.0;
	dataType mp = h * h * h;
	dataType coef_tau = tau / mp;
	dataType average_norm_gradient = 0.0, u_average = 0;
	dataType gauss_seidel_coef = 0.0;
	dataType norm_grad_e, norm_grad_w, norm_grad_n, norm_grad_s, norm_grad_t, norm_grad_b;
	size_t ind_east, ind_west, ind_north, ind_south;
	dataType numerator_max, numerator_min;
	size_t count_gauss_seidel_iteration = 0;
	dataType error_gauss_seidel = 0.0;
	dataType prod_pq = 0.0, prod_qp = 0.0;
	dataType value_pq = 0.0, value_qp = 0.0;
	
	do {
		number_time_step++;
		computeNormOfGradientDiamondCell3D(segmentationPtr, length, width, height, h, segPtrs);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) 
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) 
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					//epsilon regularization
					norm_grad_e = sqrt(segPtrs.e_Ptr[k][x] + eps);
					norm_grad_w = sqrt(segPtrs.w_Ptr[k][x] + eps);
					norm_grad_n = sqrt(segPtrs.n_Ptr[k][x] + eps);
					norm_grad_s = sqrt(segPtrs.s_Ptr[k][x] + eps);
					norm_grad_t = sqrt(segPtrs.t_Ptr[k][x] + eps);
					norm_grad_b = sqrt(segPtrs.b_Ptr[k][x] + eps);

					average_norm_gradient = (norm_grad_e + norm_grad_w + norm_grad_n + norm_grad_s + norm_grad_t + norm_grad_b) / 6.0;
					u_average = sqrt(average_norm_gradient * average_norm_gradient + eps);

					u_p = previousSolPtr[k_ext][x_ext];
					u_p_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					u_p_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					numerator_max = mp * (u_p_max - u_p);
					numerator_min = mp * (u_p_min - u_p);

					//============= Compute theta_out_pq ===========
					if (n_out_pq[k][x] == 0)
					{
						theta_out.e_Ptr[k][x] = 0.5;
						theta_out.w_Ptr[k][x] = 0.5;
						theta_out.s_Ptr[k][x] = 0.5;
						theta_out.n_Ptr[k][x] = 0.5;
						theta_out.t_Ptr[k][x] = 0.5;
						theta_out.b_Ptr[k][x] = 0.5;
					}
					else
					{
						//East
						prod_pq = a_out.e_Ptr[k][x] * (previousSolPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)] - u_p);
						if (prod_pq == 0)
						{
							theta_out.e_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.e_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.e_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//West
						prod_pq = a_out.w_Ptr[k][x] * (previousSolPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)] - u_p);
						if (prod_pq == 0)
						{
							theta_out.w_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.w_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.w_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//North
						prod_pq = a_out.n_Ptr[k][x] * (previousSolPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - u_p);
						if (prod_pq == 0)
						{
							theta_out.n_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.n_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.n_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//South
						prod_pq = a_out.s_Ptr[k][x] * (previousSolPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - u_p);
						if (prod_pq == 0)
						{
							theta_out.s_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.s_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.s_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//Top
						prod_pq = a_out.t_Ptr[k][x] * (previousSolPtr[k_ext - 1][x_ext] - u_p);
						if (prod_pq == 0)
						{
							theta_out.t_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.t_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.t_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//Bottom
						prod_pq = a_out.b_Ptr[k][x] * (previousSolPtr[k_ext + 1][x_ext] - u_p);
						if (prod_pq == 0)
						{
							theta_out.b_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.b_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.b_Ptr[k][x] = fmin(0.5, value_pq);
						}
					}

					//============= Compute theta_in_pq ===========
					//a_in_qp = - a_out_pq 
					//a_out_qp = - a_in_pq
					if (n_out_qp[k][x] == 0)
					{
						theta_in.e_Ptr[k][x] = 0.5;
						theta_in.w_Ptr[k][x] = 0.5;
						theta_in.s_Ptr[k][x] = 0.5;
						theta_in.n_Ptr[k][x] = 0.5;
						theta_in.t_Ptr[k][x] = 0.5;
						theta_in.b_Ptr[k][x] = 0.5;
					}
					else
					{
						//East
						prod_qp = -a_in.e_Ptr[k][x] * (u_p - previousSolPtr[k_ext][x_new(i_ext + 1, j_ext, length_ext)]);
						if (prod_qp == 0)
						{
							theta_in.e_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.e_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.e_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
						//West
						prod_qp = -a_in.w_Ptr[k][x] * (u_p - previousSolPtr[k_ext][x_new(i_ext - 1, j_ext, length_ext)]);
						if (prod_qp == 0)
						{
							theta_in.w_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.w_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.w_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
						//North
						prod_qp = -a_in.n_Ptr[k][x] * (u_p - previousSolPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)]);
						if (prod_qp == 0)
						{
							theta_in.n_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.n_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.n_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
						//South
						prod_qp = -a_in.s_Ptr[k][x] * (u_p - previousSolPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)]);
						if (prod_qp == 0)
						{
							theta_in.s_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.s_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.s_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
						//Top
						prod_qp = -a_in.t_Ptr[k][x] * (u_p - previousSolPtr[k_ext - 1][x_ext]);
						if (prod_qp == 0)
						{
							theta_in.t_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.t_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.t_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
						//Bottom
						prod_qp = -a_in.b_Ptr[k][x] * (u_p - previousSolPtr[k_ext + 1][x_ext]);
						if (prod_qp == 0)
						{
							theta_in.b_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = numerator_max / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.b_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = numerator_min / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.b_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						
					}

					uCoef.e_Ptr[k][x] = (dataType)(theta_in.e_Ptr[k][x] * a_in.e_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_e);
					uCoef.w_Ptr[k][x] = (dataType)(theta_in.w_Ptr[k][x] * a_in.w_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_w);
					uCoef.n_Ptr[k][x] = (dataType)(theta_in.n_Ptr[k][x] * a_in.n_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_n);
					uCoef.s_Ptr[k][x] = (dataType)(theta_in.s_Ptr[k][x] * a_in.s_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_s);
					uCoef.t_Ptr[k][x] = (dataType)(theta_in.t_Ptr[k][x] * a_in.t_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_t);
					uCoef.b_Ptr[k][x] = (dataType)(theta_in.b_Ptr[k][x] * a_in.b_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / norm_grad_b);
				}
			}
		}

		//Copy to extended area
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);
					previousSolPtr[k_ext][x_ext] = segmentationPtr[k][x];
					gaussSeidelPtr[k_ext][x_ext] = segmentationPtr[k][x];
				}
			}
		}
		setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);
		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		
		//gauss seidel for segmentation function
		count_gauss_seidel_iteration = 0;
		
		do {
			count_gauss_seidel_iteration++;
			for (k = 0, k_ext = 1; k < height; k++, k_ext++) 
			{
				for (i = 0, i_ext = 1; i < length; i++, i_ext++) 
				{
					for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
					{
						ind_east = x_new(i_ext + 1, j_ext, length_ext);
						ind_west = x_new(i_ext - 1, j_ext, length_ext);
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);

						gauss_seidel_coef = (dataType)(((1 - coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] + theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] + theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] + theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext]
							+ coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
								+ theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
								+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
								+ theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
								+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
								+ theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext])
							+ coef_tau * (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] + uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north]
								+ uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] + uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
							/ (1.0 + coef_tau * (uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x])));

						gaussSeidelPtr[k_ext][x_ext] = gaussSeidelPtr[k_ext][x_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[k_ext][x_ext]);
					}
				}
			}
		
			error_gauss_seidel = 0.0;
			for (k = 0, k_ext = 1; k < height; k++, k_ext++)
			{
				for (i = 0, i_ext = 1; i < length; i++, i_ext++)
				{
					for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
					{
						ind_east = x_new(i_ext + 1, j_ext, length_ext);
						ind_west = x_new(i_ext - 1, j_ext, length_ext);
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);

						u1 = (1.0 + coef_tau * (uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x])) * gaussSeidelPtr[k_ext][x_ext];
						u2 = (1 - coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] + theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] + theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] + theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext];
						u3 = coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
							+ theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
							+ theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
							+ theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext]);
						u4 = coef_tau * (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]);
						error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
					}
				}
			}
		} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

		//rescall to data range 0-1
		rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);
		
		//compute L2-norm
		error_segmentation = l2normD(previousSolPtr, gaussSeidelPtr, length_ext, width_ext, height_ext, h);

		//copy to reduce array
		copyDataToReducedArea(segmentationPtr, gaussSeidelPtr, height, length, width);
		
		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store3dDataArrayD(segmentationPtr, length, width, height, name, flags);
			printf("Step %zd , residual = %e \n", number_time_step, error_segmentation);
		}
	
	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > seg_parms.segTolerance);
	
	for (k = 0; k < height_ext; k++) 
	{
		if(k < height) 
		{
			free(segmentationPtr[k]);
			free(edgeDetectorPtr[k]);
			free(gPtrs.e_Ptr[k]);
			free(gPtrs.w_Ptr[k]);
			free(gPtrs.n_Ptr[k]);
			free(gPtrs.s_Ptr[k]);
			free(gPtrs.t_Ptr[k]);
			free(gPtrs.b_Ptr[k]);
			free(segPtrs.e_Ptr[k]);
			free(segPtrs.w_Ptr[k]);
			free(segPtrs.n_Ptr[k]);
			free(segPtrs.s_Ptr[k]);
			free(segPtrs.t_Ptr[k]);
			free(segPtrs.b_Ptr[k]);
			free(uCoef.e_Ptr[k]);
			free(uCoef.w_Ptr[k]);
			free(uCoef.n_Ptr[k]);
			free(uCoef.s_Ptr[k]);
			free(uCoef.t_Ptr[k]);
			free(uCoef.b_Ptr[k]);
			free(coefPtrs.e_Ptr[k]);
			free(coefPtrs.w_Ptr[k]);
			free(coefPtrs.n_Ptr[k]);
			free(coefPtrs.s_Ptr[k]);
			free(coefPtrs.t_Ptr[k]);
			free(coefPtrs.b_Ptr[k]);
			free(a_out.e_Ptr[k]);
			free(a_out.w_Ptr[k]);
			free(a_out.n_Ptr[k]);
			free(a_out.s_Ptr[k]);
			free(a_out.t_Ptr[k]);
			free(a_out.b_Ptr[k]);
			free(a_in.e_Ptr[k]);
			free(a_in.w_Ptr[k]);
			free(a_in.n_Ptr[k]);
			free(a_in.s_Ptr[k]);
			free(a_in.t_Ptr[k]);
			free(a_in.b_Ptr[k]);
			free(theta_out.e_Ptr[k]);
			free(theta_out.w_Ptr[k]);
			free(theta_out.n_Ptr[k]);
			free(theta_out.s_Ptr[k]);
			free(theta_out.t_Ptr[k]);
			free(theta_out.b_Ptr[k]);
			free(theta_in.e_Ptr[k]);
			free(theta_in.w_Ptr[k]);
			free(theta_in.n_Ptr[k]);
			free(theta_in.s_Ptr[k]);
			free(theta_in.t_Ptr[k]);
			free(theta_in.b_Ptr[k]);
			free(n_out_pq[k]);	
			free(n_out_qp[k]);
		}
		free(previousSolPtr[k]);
		free(gaussSeidelPtr[k]);
	}
	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(previousSolPtr);
	free(gaussSeidelPtr);
	free(gPtrs.e_Ptr);
	free(gPtrs.w_Ptr);
	free(gPtrs.n_Ptr);
	free(gPtrs.s_Ptr);
	free(gPtrs.t_Ptr);
	free(gPtrs.b_Ptr);
	free(segPtrs.e_Ptr);
	free(segPtrs.w_Ptr);
	free(segPtrs.n_Ptr);
	free(segPtrs.s_Ptr);
	free(segPtrs.t_Ptr);
	free(segPtrs.b_Ptr);
	free(uCoef.e_Ptr);
	free(uCoef.w_Ptr);
	free(uCoef.n_Ptr);
	free(uCoef.s_Ptr);
	free(uCoef.t_Ptr);
	free(uCoef.b_Ptr);
	free(coefPtrs.e_Ptr);
	free(coefPtrs.w_Ptr);
	free(coefPtrs.n_Ptr);
	free(coefPtrs.s_Ptr);
	free(coefPtrs.t_Ptr);
	free(coefPtrs.b_Ptr);
	free(a_out.e_Ptr);
	free(a_out.w_Ptr);
	free(a_out.n_Ptr);
	free(a_out.s_Ptr);
	free(a_out.t_Ptr);
	free(a_out.b_Ptr);
	free(a_in.e_Ptr);
	free(a_in.w_Ptr);
	free(a_in.n_Ptr);
	free(a_in.s_Ptr);
	free(a_in.t_Ptr);
	free(a_in.b_Ptr);
	free(theta_out.e_Ptr);
	free(theta_out.w_Ptr);
	free(theta_out.n_Ptr);
	free(theta_out.s_Ptr);
	free(theta_out.t_Ptr);
	free(theta_out.b_Ptr);
	free(theta_in.e_Ptr);
	free(theta_in.w_Ptr);
	free(theta_in.n_Ptr);
	free(theta_in.s_Ptr);
	free(theta_in.t_Ptr);
	free(theta_in.b_Ptr);
	free(n_out_pq);
	free(n_out_qp);

	return true;
}

void computeNormOfGradientReducedDiamondCells(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Coefficient_Pointers nGrad)
{
	if (imageData == NULL || nGrad.e_Ptr == NULL || nGrad.w_Ptr == NULL ||
		nGrad.s_Ptr == NULL || nGrad.n_Ptr == NULL ||
		nGrad.t_Ptr == NULL || nGrad.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}

	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;

	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t height_ext = height + 2;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++)
	{
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if (extendedCoefPtr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}
	if (extendedCoefPtr == NULL)
	{
		return false; // Memory allocation failed
	}

	//Copy to extended area
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = imageData[k][x_new(i, j, length)];
			}
		}
	}
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext); // Reflect the data in the extended area

	// Calculate the norm of gradient in diamond cell
	size_t iminus1, iplus1, jminus1, jplus1, kminus1, kplus1;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW;
	dataType Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW;
	dataType Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType ux, uy, uz;
	dataType quotient = (dataType)(4.0 * h);

	for (k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);
}

void computeNormOfGradientSplitDiamondCells(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Coefficient_Pointers nGrad)
{
	if (imageData == NULL || nGrad.e_Ptr == NULL || nGrad.w_Ptr == NULL ||
		nGrad.s_Ptr == NULL || nGrad.n_Ptr == NULL ||
		nGrad.t_Ptr == NULL || nGrad.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}

	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;

	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t height_ext = height + 2;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	for (k = 0; k < height_ext; k++)
	{
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if (extendedCoefPtr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}
	if (extendedCoefPtr == NULL)
	{
		return false; // Memory allocation failed
	}

	//Copy to extended area
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = imageData[k][x_new(i, j, length)];
			}
		}
	}
	reflection3D(extendedCoefPtr, height_ext, length_ext, width_ext); // Reflect the data in the extended area

	// Calculate the norm of gradient in diamond cell
	size_t iminus1, iplus1, jminus1, jplus1, kminus1, kplus1;
	dataType u, uN, uS, uE, uW, uNW, uNE, uSE, uSW;
	dataType Tu, TuN, TuS, TuE, TuW, TuNW, TuNE, TuSE, TuSW;
	dataType Bu, BuN, BuS, BuE, BuW, BuNW, BuNE, BuSE, BuSW;
	dataType ux, uy, uz;
	dataType quotient = (dataType)(4.0 * h);

	for (k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);
}
