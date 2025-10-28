#include <stdio.h> // Standard lib for input and output functions
#include <stdlib.h>
#include <time.h>
#include <math.h> // Maths functions i.e. pow, sin, cos
#include <stdbool.h> // Boolean function bool
#include <string.h>

#include "conmon_filtering.h"
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
	if (gauss_seidelPtr == NULL || prevSol_extPtr == NULL)
		return false;
	for (k = 0; k < height_ext; k++)
	{
		gauss_seidelPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		prevSol_extPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		if (gauss_seidelPtr[k] == NULL || prevSol_extPtr[k] == NULL)
			return false;
	}
	
	dataType** segmFuntionPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeGradientPtr = (dataType**)malloc(sizeof(dataType*) * height);
	if (segmFuntionPtr == NULL || edgeGradientPtr == NULL)
		return false;
	for (k = 0; k < height; k++)
	{
		segmFuntionPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		edgeGradientPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segmFuntionPtr[k] == NULL || edgeGradientPtr[k] == NULL)
			return false;
	}
	
	Coefficient_Pointers CoefPtrs;
	CoefPtrs.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	CoefPtrs.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	CoefPtrs.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	CoefPtrs.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	CoefPtrs.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	CoefPtrs.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (CoefPtrs.e_Ptr == NULL || CoefPtrs.w_Ptr == NULL || CoefPtrs.n_Ptr == NULL
		|| CoefPtrs.s_Ptr == NULL || CoefPtrs.t_Ptr == NULL || CoefPtrs.b_Ptr == NULL)
	{
		return false;
	}
	for (k = 0; k < height; k++)
	{
		CoefPtrs.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		CoefPtrs.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		CoefPtrs.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		CoefPtrs.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		CoefPtrs.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		CoefPtrs.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (CoefPtrs.e_Ptr[k] == NULL || CoefPtrs.w_Ptr[k] == NULL || CoefPtrs.n_Ptr[k] == NULL
			|| CoefPtrs.s_Ptr[k] == NULL || CoefPtrs.t_Ptr[k] == NULL || CoefPtrs.b_Ptr[k] == NULL)
		{
			return false;
		}
	}

	Gradient_Pointers vPtrs;
	vPtrs.GePtr = (dataType**)malloc(sizeof(dataType*) * height);
	vPtrs.GwPtr = (dataType**)malloc(sizeof(dataType*) * height);
	vPtrs.GnPtr = (dataType**)malloc(sizeof(dataType*) * height);
	vPtrs.GsPtr = (dataType**)malloc(sizeof(dataType*) * height);
	vPtrs.GtPtr = (dataType**)malloc(sizeof(dataType*) * height);
	vPtrs.GbPtr = (dataType**)malloc(sizeof(dataType*) * height);
	if (vPtrs.GePtr == NULL || vPtrs.GwPtr == NULL || vPtrs.GnPtr == NULL
		|| vPtrs.GsPtr == NULL || vPtrs.GtPtr == NULL || vPtrs.GbPtr == NULL)
	{
		return false;
	}
	for (k = 0; k < height; k++)
	{
		vPtrs.GePtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		vPtrs.GwPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		vPtrs.GnPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		vPtrs.GsPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		vPtrs.GtPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		vPtrs.GbPtr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (vPtrs.GePtr[k] == NULL || vPtrs.GwPtr[k] == NULL || vPtrs.GnPtr[k] == NULL
			|| vPtrs.GsPtr[k] == NULL || vPtrs.GtPtr[k] == NULL || vPtrs.GbPtr[k] == NULL)
		{
			return false;
		}
	}

	Coefficient_Pointers gPtrs;
	gPtrs.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.e_Ptr == NULL || gPtrs.w_Ptr == NULL || gPtrs.n_Ptr == NULL
		|| gPtrs.s_Ptr == NULL || gPtrs.t_Ptr == NULL || gPtrs.b_Ptr == NULL)
	{
		return false;
	}
	for(k = 0; k < height; k++)
	{
		gPtrs.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if(gPtrs.e_Ptr[k] == NULL || gPtrs.w_Ptr[k] == NULL || gPtrs.n_Ptr[k] == NULL
			|| gPtrs.s_Ptr[k] == NULL || gPtrs.t_Ptr[k] == NULL || gPtrs.b_Ptr[k] == NULL)
		{
			return false;
		}
	}

	//Initialization
	for(k = 0; k < height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			segmFuntionPtr[k][i] = 0.0;
			edgeGradientPtr[k][i] = 0.0;
			CoefPtrs.e_Ptr[k][i] = 0.0;
			CoefPtrs.w_Ptr[k][i] = 0.0;
			CoefPtrs.n_Ptr[k][i] = 0.0;
			CoefPtrs.s_Ptr[k][i] = 0.0;
			CoefPtrs.t_Ptr[k][i] = 0.0;
			CoefPtrs.b_Ptr[k][i] = 0.0;
			vPtrs.GePtr[k][i] = 0.0;
			vPtrs.GwPtr[k][i] = 0.0;
			vPtrs.GnPtr[k][i] = 0.0;
			vPtrs.GsPtr[k][i] = 0.0;
			vPtrs.GtPtr[k][i] = 0.0;
			vPtrs.GbPtr[k][i] = 0.0;
			gPtrs.e_Ptr[k][i] = 0.0;
			gPtrs.w_Ptr[k][i] = 0.0;
			gPtrs.n_Ptr[k][i] = 0.0;
			gPtrs.s_Ptr[k][i] = 0.0;
			gPtrs.t_Ptr[k][i] = 0.0;
			gPtrs.b_Ptr[k][i] = 0.0;
		}
	}

	//Array for name construction
	unsigned char  name[500];
	unsigned char  name_ending[200];
	Storage_Flags flags = { false,false };
	
	size_t k_ext, i_ext, j_ext, xd_ext;
	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				xd = x_new(i, j, length);
				xd_ext = x_new(i_ext, j_ext, length_ext);
				segmFuntionPtr[k][xd] = initialSegment[k][xd];
				gauss_seidelPtr[k_ext][xd_ext] = initialSegment[k][xd];
				prevSol_extPtr[k_ext][xd_ext] = initialSegment[k][xd];
			}
		}
	}

	copyDataToAnotherArray(initialSegment, segmFuntionPtr, height, length, width);
	
	copyDataToExtendedArea(initialSegment, gauss_seidelPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gauss_seidelPtr, length_ext, width_ext, height_ext);
	
	copyDataToExtendedArea(initialSegment, prevSol_extPtr, height, length, width);
	setBoundaryToZeroDirichletBC(prevSol_extPtr, length_ext, width_ext, height_ext);

	////smoothing
	heatImplicitScheme(inputImageData, explicit_lhe_Parameters);

	//compute the morm of gradient for the edge detector
	computeNormOfGradientDiamondCell3D(inputImageData.imageDataPtr, length, width, height, h, gPtrs);

	//edge detection
	dataType average_g_grad = 0.0;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	for(k = 0; k < height; k++)
	{
		for(i = 0; i < length; i++)
		{
			for(j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				average_g_grad = (sqrt(gPtrs.e_Ptr[k][xd]) + sqrt(gPtrs.w_Ptr[k][xd]) + sqrt(gPtrs.n_Ptr[k][xd]) + 
					sqrt(gPtrs.s_Ptr[k][xd]) + sqrt(gPtrs.t_Ptr[k][xd]) + sqrt(gPtrs.b_Ptr[k][xd])) / 6.0;
				edgeGradientPtr[k][xd] = gradientFunction(average_g_grad * average_g_grad, segParameters.coef);

				if (i == 0)
				{
					vpe = -coef_conv * h * (edgeGradientPtr[k][x_new(i + 1, j, length)] - edgeGradientPtr[k][xd]);
					vpw = -vpe;
				}
				else if (i == length - 1)
				{
					vpw = -coef_conv * h * (edgeGradientPtr[k][xd] - edgeGradientPtr[k][x_new(i - 1, j, length)]);
					vpe = -vpw;
				}
				else {
					vpe = -coef_conv * h * 0.5 * (edgeGradientPtr[k][x_new(i + 1, j, length)] - edgeGradientPtr[k][x_new(i - 1, j, length)]);
					vpw = -vpe;
				}

				if (j == 0)
				{
					vps = -coef_conv * h * (edgeGradientPtr[k][x_new(i, j + 1, length)] - edgeGradientPtr[k][xd]);
					vpn = -vps;
				}
				else if (j == width - 1)
				{
					vpn = -coef_conv * h * (edgeGradientPtr[k][xd] - edgeGradientPtr[k][x_new(i, j - 1, length)]);
					vps = -vpn;
				}
				else
				{
					vps = -coef_conv * 0.5 * h * (edgeGradientPtr[k][x_new(i, j + 1, length)] - edgeGradientPtr[k][x_new(i, j - 1, length)]);
					vpn = -vps;
				}

				if (k == 0)
				{
					vpb = -coef_conv * h * (edgeGradientPtr[k + 1][xd] - edgeGradientPtr[k][xd]);
					vpt = -vpb;
				}
				else if (k == height - 1)
				{
					vpt = -coef_conv * h * (edgeGradientPtr[k][xd] - edgeGradientPtr[k - 1][xd]);
					vpb = -vpt;
				}
				else 
				{
					vpb = -coef_conv * 0.5 * h * (edgeGradientPtr[k + 1][xd] - edgeGradientPtr[k - 1][xd]);
					vpt = -vpb;
				}

				vPtrs.GePtr[k][xd] = fmin(vpe, 0.0);
				vPtrs.GwPtr[k][xd] = fmin(vpw, 0.0);
				vPtrs.GnPtr[k][xd] = fmin(vpn, 0.0);
				vPtrs.GsPtr[k][xd] = fmin(vps, 0.0);
				vPtrs.GtPtr[k][xd] = fmin(vpt, 0.0);
				vPtrs.GbPtr[k][xd] = fmin(vpb, 0.0);
			}
		}
	}

	bool isFileSaved;
	//strcpy_s(name, sizeof name, outputPathPtr);
	//sprintf_s(name_ending, sizeof(name_ending), "_smoothed.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//isFileSaved = manageFile(inputImageData.imageDataPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
	//if (isFileSaved == false) {
	//	printf("The file was not saved\n");
	//	return false;
	//}

	strcpy_s(name, sizeof name, outputPathPtr);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	isFileSaved = manageFile(edgeGradientPtr, length, width, height, name, STORE_DATA_RAW, BINARY_DATA, flags);
	if (isFileSaved == false) {
		printf("The file was not saved\n");
		return false;
	}

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

		//calcution of coefficients
		generalizedGaussSeidelCoefficients(segmentationFunction, edgeGradientPtr, CoefPtrs, vPtrs, segParameters);

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
			//printf("Step : %zd\n", number_time_step);
			if (isFileSaved == false) {
				printf("The file was not saved\n");
				return false;
			}
			printf("Step %zd , residual = %e \n", number_time_step, difference_btw_current_and_previous_sol);
			//fprintf(error_file, "%d,%f\n", number_time_step, difference_btw_current_and_previous_sol);
		}

	} while ((number_time_step <= segParameters.maxNoOfTimeSteps) && (difference_btw_current_and_previous_sol > segParameters.segTolerance));
	
	for (i = 0; i < height; i++)
	{
		free(segmFuntionPtr[i]);
		free(edgeGradientPtr[i]);
		free(gPtrs.e_Ptr[i]);
		free(gPtrs.w_Ptr[i]);
		free(gPtrs.n_Ptr[i]);
		free(gPtrs.s_Ptr[i]);
		free(gPtrs.t_Ptr[i]);
		free(gPtrs.b_Ptr[i]);
		free(vPtrs.GePtr[i]);
		free(vPtrs.GwPtr[i]);
		free(vPtrs.GnPtr[i]);
		free(vPtrs.GsPtr[i]);
		free(vPtrs.GtPtr[i]);
		free(vPtrs.GbPtr[i]);
	}
	free(segmFuntionPtr);
	free(edgeGradientPtr);
	free(gPtrs.e_Ptr);
	free(gPtrs.w_Ptr);
	free(gPtrs.n_Ptr);
	free(gPtrs.s_Ptr);
	free(gPtrs.t_Ptr);
	free(gPtrs.b_Ptr);
	free(vPtrs.GePtr);
	free(vPtrs.GwPtr);
	free(vPtrs.GnPtr);
	free(vPtrs.GsPtr);
	free(vPtrs.GtPtr);
	free(vPtrs.GbPtr);

	for (i = 0; i < height_ext; i++)
	{
		free(prevSol_extPtr[i]);
		free(gauss_seidelPtr[i]);
	}
	free(prevSol_extPtr);
	free(gauss_seidelPtr);

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

	if(height_ext < 3 || length_ext < 3 || width_ext < 3)
	{
		return false;
	}	
	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (extendedCoefPtr == NULL)
		return false;
	for (k = 0; k < height_ext; k++) {
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if (extendedCoefPtr[k] == NULL)
			return false;
	}
	
	copyDataToExtendedArea(segmentationData.imageDataPtr, extendedCoefPtr, length, width, height);
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
				CoefPtrs.e_Ptr[k][x] = (dataType)(-VPtrs.GePtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_e);
				CoefPtrs.w_Ptr[k][x] = (dataType)(-VPtrs.GwPtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_w);
				CoefPtrs.n_Ptr[k][x] = (dataType)(-VPtrs.GnPtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_n);
				CoefPtrs.s_Ptr[k][x] = (dataType)(-VPtrs.GsPtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_s);
				CoefPtrs.t_Ptr[k][x] = (dataType)(-VPtrs.GtPtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_t);
				CoefPtrs.b_Ptr[k][x] = (dataType)(-VPtrs.GbPtr[k][x] + coef_dif * voxel_coef * h * edgeGradientPtr[k][x] / orig_b);

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
	dataType hhh = segParameters.h * segParameters.h * segParameters.h;
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

	const dataType coef_tauh = tau / hhh;
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

	copyDataToReducedArea(segmentationData.imageDataPtr, gauss_seidelPtr, height, length, width);

	return true;
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
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL ||
		gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}
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
		if(gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL)
		{
			return false;
		}
	}

	Pointers_Neighbours gPtrs;
	gPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.east == NULL || gPtrs.west == NULL || gPtrs.north == NULL ||
		gPtrs.south == NULL || gPtrs.top == NULL || gPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		gPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.east[k] == NULL || gPtrs.west[k] == NULL || gPtrs.north[k] == NULL ||
			gPtrs.south[k] == NULL || gPtrs.top[k] == NULL || gPtrs.bottom[k] == NULL)
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
	if (segPtrs.e_Ptr == NULL || segPtrs.w_Ptr == NULL || segPtrs.n_Ptr == NULL ||
		segPtrs.s_Ptr == NULL || segPtrs.t_Ptr == NULL || segPtrs.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (uCoef.e_Ptr == NULL || uCoef.w_Ptr == NULL || uCoef.n_Ptr == NULL ||
		uCoef.s_Ptr == NULL || uCoef.t_Ptr == NULL || uCoef.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (coefPtrs.e_Ptr == NULL || coefPtrs.w_Ptr == NULL || coefPtrs.n_Ptr == NULL ||
		coefPtrs.s_Ptr == NULL || coefPtrs.t_Ptr == NULL || coefPtrs.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (a_out.e_Ptr == NULL || a_out.w_Ptr == NULL || a_out.n_Ptr == NULL ||
		a_out.s_Ptr == NULL || a_out.t_Ptr == NULL || a_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (a_in.e_Ptr == NULL || a_in.w_Ptr == NULL || a_in.n_Ptr == NULL ||
		a_in.s_Ptr == NULL || a_in.t_Ptr == NULL || a_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (theta_out.e_Ptr == NULL || theta_out.w_Ptr == NULL || theta_out.n_Ptr == NULL ||
		theta_out.s_Ptr == NULL || theta_out.t_Ptr == NULL || theta_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (theta_in.e_Ptr == NULL || theta_in.w_Ptr == NULL || theta_in.n_Ptr == NULL ||
		theta_in.s_Ptr == NULL || theta_in.t_Ptr == NULL || theta_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (n_out_pq == NULL || n_out_qp == NULL)
	{
		return false;
	}
	for (k = 0; k < height; k++) 
	{
		n_out_pq[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		n_out_qp[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (n_out_pq[k] == NULL || n_out_qp[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	//Initialization
	for(k = 0; k < height; k++)
	{
		for(i = 0; i < dim2D; i++)
		{
			segmentationPtr[k][i] = 0.0;
			edgeDetectorPtr[k][i] = 0.0;
			gPtrs.east[k][i] = 0.0;
			gPtrs.west[k][i] = 0.0;
			gPtrs.north[k][i] = 0.0;
			gPtrs.south[k][i] = 0.0;
			gPtrs.top[k][i] = 0.0;
			gPtrs.bottom[k][i] = 0.0;
			segPtrs.e_Ptr[k][i] = 0.0;
			segPtrs.w_Ptr[k][i] = 0.0;
			segPtrs.n_Ptr[k][i] = 0.0;
			segPtrs.s_Ptr[k][i] = 0.0;
			segPtrs.t_Ptr[k][i] = 0.0;
			segPtrs.b_Ptr[k][i] = 0.0;
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

	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);
	heatImplicitScheme(imageData, smooth_parms);

	////compute the morm of gradient for the edge detector
	//computeNormOfGradientDiamondCell3D(imageData.imageDataPtr, length, width, height, h, gPtrs);
	normOfGradientReducedDiamondCells(imageData, gPtrs);

	//compute the edge detector : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	dataType average_value = 0.0;
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	for (k = 0; k < height; k++) 
	{
		for (i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				x = x_new(i, j, length);

				average_value = (sqrt(gPtrs.east[k][x]) + sqrt(gPtrs.west[k][x]) + sqrt(gPtrs.north[k][x]) +
					sqrt(gPtrs.south[k][x]) + sqrt(gPtrs.top[k][x]) + sqrt(gPtrs.bottom[k][x])) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(average_value * average_value, coef_edge_detector);

				if (i == 0) 
				{
					vpe = adv * h * (edgeDetectorPtr[k][x_new(i + 1, j, length)] - edgeDetectorPtr[k][x]);
					vpw = -vpe;
				}
				else if (i == length - 1) 
				{
					vpw = adv * h * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k][x_new(i - 1, j, length)]);
					vpe = -vpw;					
				}
				else {
					vpe = adv * 0.5 * h * (edgeDetectorPtr[k][x_new(i + 1, j, length)] - edgeDetectorPtr[k][x_new(i - 1, j, length)]);
					vpw = -vpe;
				}

				if (j == 0) 
				{
					vps = adv * h * (edgeDetectorPtr[k][x_new(i, j + 1, length)] - edgeDetectorPtr[k][x]);
					vpn = -vps;
				}
				else if (j == width - 1) 
				{
					vpn = adv * h * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k][x_new(i, j - 1, length)]);
					vps = -vpn;
				}
				else 
				{
					vps = adv * 0.5 * h * (edgeDetectorPtr[k][x_new(i, j + 1, length)] - edgeDetectorPtr[k][x_new(i, j - 1, length)]);
					vpn = -vps;
				}
				
				if (k == 0) 
				{
					vpb = adv * h * (edgeDetectorPtr[k + 1][x] - edgeDetectorPtr[k][x]);
					vpt = -vpb;
				}
				else if (k == height - 1)
				{
					vpt = adv * h * (edgeDetectorPtr[k][x] - edgeDetectorPtr[k - 1][x]);
					vpb = -vpt;					
				}
				else {
					vpb = adv * 0.5 * h * (edgeDetectorPtr[k + 1][x] - edgeDetectorPtr[k - 1][x]);
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

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_iioe.raw");
	strcat_s(name, sizeof(name), name_ending);
	store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);
	
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
	dataType numerator_max = 0.0, numerator_min = 0.0;
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
					segPtrs.e_Ptr[k][x] = sqrt(segPtrs.e_Ptr[k][x] + eps);
					segPtrs.w_Ptr[k][x] = sqrt(segPtrs.w_Ptr[k][x] + eps);
					segPtrs.n_Ptr[k][x] = sqrt(segPtrs.n_Ptr[k][x] + eps);
					segPtrs.s_Ptr[k][x] = sqrt(segPtrs.s_Ptr[k][x] + eps);
					segPtrs.t_Ptr[k][x] = sqrt(segPtrs.t_Ptr[k][x] + eps);
					segPtrs.b_Ptr[k][x] = sqrt(segPtrs.b_Ptr[k][x] + eps);

					//average norm of gradient
					average_norm_gradient = (segPtrs.e_Ptr[k][x] + segPtrs.w_Ptr[k][x] + segPtrs.n_Ptr[k][x] + 
						segPtrs.s_Ptr[k][x] + segPtrs.t_Ptr[k][x] + segPtrs.b_Ptr[k][x]) / 6.0;
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

					uCoef.e_Ptr[k][x] = theta_in.e_Ptr[k][x] * a_in.e_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.e_Ptr[k][x];
					uCoef.w_Ptr[k][x] = theta_in.w_Ptr[k][x] * a_in.w_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.w_Ptr[k][x];
					uCoef.n_Ptr[k][x] = theta_in.n_Ptr[k][x] * a_in.n_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.n_Ptr[k][x];
					uCoef.s_Ptr[k][x] = theta_in.s_Ptr[k][x] * a_in.s_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.s_Ptr[k][x];
					uCoef.t_Ptr[k][x] = theta_in.t_Ptr[k][x] * a_in.t_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.t_Ptr[k][x];
					uCoef.b_Ptr[k][x] = theta_in.b_Ptr[k][x] * a_in.b_Ptr[k][x] + diff * h * u_average * edgeDetectorPtr[k][x] / segPtrs.b_Ptr[k][x];
				}
			}
		}
		
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

		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		copyDataToAnotherArray(gaussSeidelPtr, previousSolPtr, height_ext, length_ext, width_ext);

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
			free(gPtrs.east[k]);
			free(gPtrs.west[k]);
			free(gPtrs.north[k]);
			free(gPtrs.south[k]);
			free(gPtrs.top[k]);
			free(gPtrs.bottom[k]);
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
	free(gPtrs.east);
	free(gPtrs.west);
	free(gPtrs.north);
	free(gPtrs.south);
	free(gPtrs.top);
	free(gPtrs.bottom);
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

bool GSUBSURF(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
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

	size_t dim2D = length * width;
	size_t dim2D_ext = length_ext * width_ext;

	dataType hx = imageData.spacing.sx;
	dataType hy = imageData.spacing.sy;
	dataType hz = imageData.spacing.sz;
	dataType hx2 = hx * hx, hy2 = hy * hy, hz2 = hz * hz;
	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps2 = seg_parms.eps2;
	dataType coef_diff = seg_parms.coef_dif, coef_adv = seg_parms.coef_conv;

	dataType** segmentationPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeDetectorPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedEdge = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL ||
		gaussSeidelPtr == NULL || previousSolPtr == NULL || extendedEdge == NULL)
	{
		return false;
	}
	for (k = 0; k < height_ext; k++)
	{
		if (k < height)
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
		extendedEdge[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		if (gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL || extendedEdge[k] == NULL)
		{
			return false;
		}
	}

	Pointers_Neighbours gPtrs;
	gPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.east == NULL || gPtrs.west == NULL || gPtrs.north == NULL ||
		gPtrs.south == NULL || gPtrs.top == NULL || gPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		gPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.east[k] == NULL || gPtrs.west[k] == NULL || gPtrs.north[k] == NULL ||
			gPtrs.south[k] == NULL || gPtrs.top[k] == NULL || gPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Pointers_Neighbours segPtrs;
	segPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (segPtrs.east == NULL || segPtrs.west == NULL || segPtrs.north == NULL ||
		segPtrs.south == NULL || segPtrs.top == NULL || segPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		segPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segPtrs.east[k] == NULL || segPtrs.west[k] == NULL || segPtrs.north[k] == NULL ||
			segPtrs.south[k] == NULL || segPtrs.top[k] == NULL || segPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (uCoef.e_Ptr == NULL || uCoef.w_Ptr == NULL || uCoef.n_Ptr == NULL ||
		uCoef.s_Ptr == NULL || uCoef.t_Ptr == NULL || uCoef.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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

	Coefficient_Pointers v_in;
	v_in.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	v_in.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	v_in.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	v_in.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	v_in.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	v_in.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (v_in.e_Ptr == NULL || v_in.w_Ptr == NULL || v_in.n_Ptr == NULL ||
		v_in.s_Ptr == NULL || v_in.t_Ptr == NULL || v_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		v_in.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		v_in.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		v_in.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		v_in.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		v_in.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		v_in.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (v_in.e_Ptr[k] == NULL || v_in.w_Ptr[k] == NULL || v_in.n_Ptr[k] == NULL ||
			v_in.s_Ptr[k] == NULL || v_in.t_Ptr[k] == NULL || v_in.b_Ptr[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);

	//compute the morm of gradient for the edge detector
	normOfGradientReducedDiamondCells(imageData, gPtrs);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				dataType avg_value = (sqrt(gPtrs.east[k][x]) + sqrt(gPtrs.west[k][x]) + sqrt(gPtrs.north[k][x]) +
					sqrt(gPtrs.south[k][x]) + sqrt(gPtrs.top[k][x]) + sqrt(gPtrs.bottom[k][x])) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(avg_value * avg_value, coef_edge_detector);
				extendedEdge[k_ext][x_ext] = edgeDetectorPtr[k][x];
			}
		}
	}
	reflection3D(extendedEdge, height_ext, length_ext, width_ext);

	//compute the edge detector : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	dataType average_value = 0.0;
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	dataType mp = hx * hy * hz;
	dataType coef_tau = tau / mp;

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				vpe = -coef_adv * hy * hz * (extendedEdge[k_ext][x_ext + 1] - extendedEdge[k_ext][x_ext - 1]) / (2 * hx);
				vpw = -coef_adv * hy * hz * (extendedEdge[k_ext][x_ext - 1] - extendedEdge[k_ext][x_ext + 1]) / (2 * hx);
				vps = -coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)]) / (2 * hy);
				vpn = -coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)]) / (2 * hy);
				vpt = -coef_adv * hx * hy * (extendedEdge[k_ext - 1][x_ext] - extendedEdge[k_ext + 1][x_ext]) / (2 * hz);
				vpb = -coef_adv * hx * hy * (extendedEdge[k_ext + 1][x_ext] - extendedEdge[k_ext - 1][x_ext]) / (2 * hz);

				v_in.e_Ptr[k][x] = fmin(vpe, 0);
				v_in.w_Ptr[k][x] = fmin(vpw, 0);
				v_in.n_Ptr[k][x] = fmin(vpn, 0);
				v_in.s_Ptr[k][x] = fmin(vps, 0);
				v_in.t_Ptr[k][x] = fmin(vpt, 0);
				v_in.b_Ptr[k][x] = fmin(vpb, 0);
			}
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_iioe.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);

	copyDataToAnotherArray(initialSegment, segmentationPtr, height, length, width);

	copyDataToExtendedArea(initialSegment, previousSolPtr, height, length, width);
	setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);

	copyDataToExtendedArea(initialSegment, gaussSeidelPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType u1 = 0.0, u2 = 0.0;
	dataType average_norm_gradient = 0.0, u_average = 0;
	dataType gauss_seidel_coef = 0.0;
	size_t ind_east, ind_west, ind_north, ind_south;
	size_t count_gauss_seidel_iteration = 0;
	dataType error_gauss_seidel = 0.0;

	Image_Data segmentationData = { height, length, width, segmentationPtr, imageData.origin, imageData.spacing, imageData.orientation };

	do {
		number_time_step++;

		normOfGradientReducedDiamondCells(segmentationData, segPtrs);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					//epsilon regularization
					segPtrs.east[k][x] = sqrt(segPtrs.east[k][x] + eps2);
					segPtrs.west[k][x] = sqrt(segPtrs.west[k][x] + eps2);
					segPtrs.north[k][x] = sqrt(segPtrs.north[k][x] + eps2);
					segPtrs.south[k][x] = sqrt(segPtrs.south[k][x] + eps2);
					segPtrs.top[k][x] = sqrt(segPtrs.top[k][x] + eps2);
					segPtrs.bottom[k][x] = sqrt(segPtrs.bottom[k][x] + eps2);

					//average norm of gradient
					average_norm_gradient = (segPtrs.east[k][x] + segPtrs.west[k][x] + segPtrs.north[k][x] +
						segPtrs.south[k][x] + segPtrs.top[k][x] + segPtrs.bottom[k][x]) / 6.0;
					u_average = sqrt(average_norm_gradient * average_norm_gradient + eps2);

					uCoef.e_Ptr[k][x] = - coef_tau * v_in.e_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.east[k][x]);
					uCoef.w_Ptr[k][x] = - coef_tau * v_in.w_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.west[k][x]);
					uCoef.n_Ptr[k][x] = - coef_tau * v_in.n_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.north[k][x]);
					uCoef.s_Ptr[k][x] = - coef_tau * v_in.s_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.south[k][x]);
					uCoef.t_Ptr[k][x] = - coef_tau * v_in.t_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.top[k][x]);
					uCoef.b_Ptr[k][x] = - coef_tau * v_in.b_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.bottom[k][x]);
				}
			}
		}

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
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);

						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						gauss_seidel_coef = (dataType)((previousSolPtr[k_ext][x_ext] 
							+ (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] 
							+ uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] 
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] 
							+ uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] 
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] 
							+ uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
							/ (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] 
								+ uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] 
								+ uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]));

						gaussSeidelPtr[k_ext][x_ext] = gaussSeidelPtr[k_ext][x_ext] 
							+ omega * (gauss_seidel_coef - gaussSeidelPtr[k_ext][x_ext]);
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

						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);
						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						u1 = (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]) * gaussSeidelPtr[k_ext][x_ext];
						u2 = uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext];
						error_gauss_seidel += pow(u1 - u2 - previousSolPtr[k_ext][x_ext], 2);
					}
				}
			}
		} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

		//rescall to data range 0-1
		rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);

		//compute L2-norm
		error_segmentation = l2normRectangularGrid(previousSolPtr, gaussSeidelPtr, length_ext, width_ext, height_ext, spacing);

		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		copyDataToAnotherArray(gaussSeidelPtr, previousSolPtr, height_ext, length_ext, width_ext);

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
		if (k < height)
		{
			free(segmentationPtr[k]);
			free(edgeDetectorPtr[k]);
			free(gPtrs.east[k]);
			free(gPtrs.west[k]);
			free(gPtrs.north[k]);
			free(gPtrs.south[k]);
			free(gPtrs.top[k]);
			free(gPtrs.bottom[k]);
			free(segPtrs.east[k]);
			free(segPtrs.west[k]);
			free(segPtrs.north[k]);
			free(segPtrs.south[k]);
			free(segPtrs.top[k]);
			free(segPtrs.bottom[k]);
			free(uCoef.e_Ptr[k]);
			free(uCoef.w_Ptr[k]);
			free(uCoef.n_Ptr[k]);
			free(uCoef.s_Ptr[k]);
			free(uCoef.t_Ptr[k]);
			free(uCoef.b_Ptr[k]);
			free(v_in.e_Ptr[k]);
			free(v_in.w_Ptr[k]);
			free(v_in.n_Ptr[k]);
			free(v_in.s_Ptr[k]);
			free(v_in.t_Ptr[k]);
			free(v_in.b_Ptr[k]);
		}
		free(previousSolPtr[k]);
		free(gaussSeidelPtr[k]);
		free(extendedEdge[k]);
	}
	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(previousSolPtr);
	free(gaussSeidelPtr);
	free(extendedEdge);
	free(gPtrs.east);
	free(gPtrs.west);
	free(gPtrs.north);
	free(gPtrs.south);
	free(gPtrs.top);
	free(gPtrs.bottom);
	free(segPtrs.east);
	free(segPtrs.west);
	free(segPtrs.north);
	free(segPtrs.south);
	free(segPtrs.top);
	free(segPtrs.bottom);
	free(uCoef.e_Ptr);
	free(uCoef.w_Ptr);
	free(uCoef.n_Ptr);
	free(uCoef.s_Ptr);
	free(uCoef.t_Ptr);
	free(uCoef.b_Ptr);
	free(v_in.e_Ptr);
	free(v_in.w_Ptr);
	free(v_in.n_Ptr);
	free(v_in.s_Ptr);
	free(v_in.t_Ptr);
	free(v_in.b_Ptr);

	return true;
}

bool GSUBSURF_IIOE(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
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

	size_t dim2D = length * width;
	size_t dim2D_ext = length_ext * width_ext;

	dataType hx = imageData.spacing.sx;
	dataType hy = imageData.spacing.sy;
	dataType hz = imageData.spacing.sz;
	dataType hx2 = hx * hx, hy2 = hy * hy, hz2 = hz * hz;
	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps2 = seg_parms.eps2;
	dataType coef_diff = seg_parms.coef_dif, coef_adv = seg_parms.coef_conv;

	dataType** segmentationPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeDetectorPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedEdge = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL ||
		gaussSeidelPtr == NULL || previousSolPtr == NULL || extendedEdge == NULL)
	{
		return false;
	}
	for (k = 0; k < height_ext; k++)
	{
		if (k < height)
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
		extendedEdge[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		if (gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL || extendedEdge[k] == NULL)
		{
			return false;
		}
	}

	Pointers_Neighbours gPtrs;
	gPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.east == NULL || gPtrs.west == NULL || gPtrs.north == NULL ||
		gPtrs.south == NULL || gPtrs.top == NULL || gPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		gPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.east[k] == NULL || gPtrs.west[k] == NULL || gPtrs.north[k] == NULL ||
			gPtrs.south[k] == NULL || gPtrs.top[k] == NULL || gPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Pointers_Neighbours segPtrs;
	segPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (segPtrs.east == NULL || segPtrs.west == NULL || segPtrs.north == NULL ||
		segPtrs.south == NULL || segPtrs.top == NULL || segPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		segPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segPtrs.east[k] == NULL || segPtrs.west[k] == NULL || segPtrs.north[k] == NULL ||
			segPtrs.south[k] == NULL || segPtrs.top[k] == NULL || segPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (uCoef.e_Ptr == NULL || uCoef.w_Ptr == NULL || uCoef.n_Ptr == NULL ||
		uCoef.s_Ptr == NULL || uCoef.t_Ptr == NULL || uCoef.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (coefPtrs.e_Ptr == NULL || coefPtrs.w_Ptr == NULL || coefPtrs.n_Ptr == NULL ||
		coefPtrs.s_Ptr == NULL || coefPtrs.t_Ptr == NULL || coefPtrs.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (a_out.e_Ptr == NULL || a_out.w_Ptr == NULL || a_out.n_Ptr == NULL ||
		a_out.s_Ptr == NULL || a_out.t_Ptr == NULL || a_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (a_in.e_Ptr == NULL || a_in.w_Ptr == NULL || a_in.n_Ptr == NULL ||
		a_in.s_Ptr == NULL || a_in.t_Ptr == NULL || a_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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

	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);

	//compute the morm of gradient for the edge detector
	normOfGradientReducedDiamondCells(imageData, gPtrs);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				dataType avg_value = (sqrt(gPtrs.east[k][x]) + sqrt(gPtrs.west[k][x]) + sqrt(gPtrs.north[k][x]) +
					sqrt(gPtrs.south[k][x]) + sqrt(gPtrs.top[k][x]) + sqrt(gPtrs.bottom[k][x])) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(avg_value * avg_value, coef_edge_detector);
				extendedEdge[k_ext][x_ext] = edgeDetectorPtr[k][x];
			}
		}
	}
	reflection3D(extendedEdge, height_ext, length_ext, width_ext);

	//compute the edge detector : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	dataType average_value = 0.0;
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	dataType mp = hx * hy * hz;
	dataType coef_tau = tau / mp;

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				vpe = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext + 1] - extendedEdge[k_ext][x_ext - 1]) / (2 * hx);
				vpw = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext - 1] - extendedEdge[k_ext][x_ext + 1]) / (2 * hx);
				vps = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)]) / (2 * hy);
				vpn = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)]) / (2 * hy);
				vpt = coef_adv * hx * hy * (extendedEdge[k_ext - 1][x_ext] - extendedEdge[k_ext + 1][x_ext]) / (2 * hz);
				vpb = coef_adv * hx * hy * (extendedEdge[k_ext + 1][x_ext] - extendedEdge[k_ext - 1][x_ext]) / (2 * hz);

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
			}
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_iioe.raw");
	strcat_s(name, sizeof(name), name_ending);
	store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);

	copyDataToAnotherArray(initialSegment, segmentationPtr, height, length, width);
	copyDataToExtendedArea(initialSegment, previousSolPtr, height, length, width);
	setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);
	copyDataToExtendedArea(initialSegment, gaussSeidelPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType average_norm_gradient = 0.0, u_average = 0;
	dataType gauss_seidel_coef = 0.0;
	size_t ind_east, ind_west, ind_north, ind_south;
	size_t count_gauss_seidel_iteration = 0;
	dataType error_gauss_seidel = 0.0;

	Image_Data segmentationData = { height, length, width, segmentationPtr, imageData.origin, imageData.spacing, imageData.orientation };

	do {
		number_time_step++;

		normOfGradientReducedDiamondCells(segmentationData, segPtrs);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					//epsilon regularization
					segPtrs.east[k][x] = sqrt(segPtrs.east[k][x] + eps2);
					segPtrs.west[k][x] = sqrt(segPtrs.west[k][x] + eps2);
					segPtrs.north[k][x] = sqrt(segPtrs.north[k][x] + eps2);
					segPtrs.south[k][x] = sqrt(segPtrs.south[k][x] + eps2);
					segPtrs.top[k][x] = sqrt(segPtrs.top[k][x] + eps2);
					segPtrs.bottom[k][x] = sqrt(segPtrs.bottom[k][x] + eps2);

					//average norm of gradient
					average_norm_gradient = (segPtrs.east[k][x] + segPtrs.west[k][x] + segPtrs.north[k][x] +
						segPtrs.south[k][x] + segPtrs.top[k][x] + segPtrs.bottom[k][x]) / 6.0;
					u_average = sqrt(average_norm_gradient * average_norm_gradient + eps2);
	
					uCoef.e_Ptr[k][x] = 0.5 * coef_tau * a_in.e_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.east[k][x]);
					uCoef.w_Ptr[k][x] = 0.5 * coef_tau * a_in.w_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.west[k][x]);
					uCoef.n_Ptr[k][x] = 0.5 * coef_tau * a_in.n_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.north[k][x]);
					uCoef.s_Ptr[k][x] = 0.5 * coef_tau * a_in.s_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.south[k][x]);
					uCoef.t_Ptr[k][x] = 0.5 * coef_tau * a_in.t_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.top[k][x]);
					uCoef.b_Ptr[k][x] = 0.5 * coef_tau * a_in.b_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.bottom[k][x]);
				}
			}
		}

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
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);
						
						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						gauss_seidel_coef = (dataType)(((1 - 0.5 * coef_tau * (a_out.e_Ptr[k][x] + a_out.w_Ptr[k][x]
							+ a_out.n_Ptr[k][x] + a_out.s_Ptr[k][x] + a_out.t_Ptr[k][x] + a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext]
							+ 0.5 * coef_tau * (a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
								+ a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
								+ a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
								+ a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
								+ a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
								+ a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext])
							+ (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] + uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north]
								+ uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] + uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
							/ (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]));

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

						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);
						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						u1 = (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]) * gaussSeidelPtr[k_ext][x_ext];
						u2 = (1 - 0.5 * coef_tau * (a_out.e_Ptr[k][x] + a_out.w_Ptr[k][x]
							+ a_out.n_Ptr[k][x] + a_out.s_Ptr[k][x]
							+ a_out.t_Ptr[k][x] + a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext];
						u3 = 0.5 * coef_tau * (a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
							+ a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
							+ a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
							+ a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
							+ a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
							+ a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext]);
						u4 = uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext];
						error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
					}
				}
			}

		} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

		//rescall to data range 0-1
		rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);

		//compute L2-norm
		error_segmentation = l2normRectangularGrid(previousSolPtr, gaussSeidelPtr, length_ext, width_ext, height_ext, spacing);

		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		copyDataToAnotherArray(gaussSeidelPtr, previousSolPtr, height_ext, length_ext, width_ext);

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
		if (k < height)
		{
			free(segmentationPtr[k]);
			free(edgeDetectorPtr[k]);
			free(gPtrs.east[k]);
			free(gPtrs.west[k]);
			free(gPtrs.north[k]);
			free(gPtrs.south[k]);
			free(gPtrs.top[k]);
			free(gPtrs.bottom[k]);
			free(segPtrs.east[k]);
			free(segPtrs.west[k]);
			free(segPtrs.north[k]);
			free(segPtrs.south[k]);
			free(segPtrs.top[k]);
			free(segPtrs.bottom[k]);
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
		}
		free(previousSolPtr[k]);
		free(gaussSeidelPtr[k]);
		free(extendedEdge[k]);
	}
	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(previousSolPtr);
	free(gaussSeidelPtr);
	free(extendedEdge);
	free(gPtrs.east);
	free(gPtrs.west);
	free(gPtrs.north);
	free(gPtrs.south);
	free(gPtrs.top);
	free(gPtrs.bottom);
	free(segPtrs.east);
	free(segPtrs.west);
	free(segPtrs.north);
	free(segPtrs.south);
	free(segPtrs.top);
	free(segPtrs.bottom);
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

	return true;
}

bool GSUBSURF_S_ONE_IIOE(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
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

	size_t dim2D = length * width;
	size_t dim2D_ext = length_ext * width_ext;

	dataType hx = imageData.spacing.sx;
	dataType hy = imageData.spacing.sy;
	dataType hz = imageData.spacing.sz;
	dataType hx2 = hx * hx, hy2 = hy * hy, hz2 = hz * hz;
	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps2 = seg_parms.eps2;
	dataType coef_diff = seg_parms.coef_dif, coef_adv = seg_parms.coef_conv;

	dataType** segmentationPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeDetectorPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedEdge = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL ||
		gaussSeidelPtr == NULL || previousSolPtr == NULL || extendedEdge == NULL)
	{
		return false;
	}
	for (k = 0; k < height_ext; k++)
	{
		if (k < height)
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
		extendedEdge[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		if (gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL || extendedEdge[k] == NULL)
		{
			return false;
		}
	}

	Pointers_Neighbours gPtrs;
	gPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.east == NULL || gPtrs.west == NULL || gPtrs.north == NULL ||
		gPtrs.south == NULL || gPtrs.top == NULL || gPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		gPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.east[k] == NULL || gPtrs.west[k] == NULL || gPtrs.north[k] == NULL ||
			gPtrs.south[k] == NULL || gPtrs.top[k] == NULL || gPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Pointers_Neighbours segPtrs;
	segPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (segPtrs.east == NULL || segPtrs.west == NULL || segPtrs.north == NULL ||
		segPtrs.south == NULL || segPtrs.top == NULL || segPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		segPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segPtrs.east[k] == NULL || segPtrs.west[k] == NULL || segPtrs.north[k] == NULL ||
			segPtrs.south[k] == NULL || segPtrs.top[k] == NULL || segPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (uCoef.e_Ptr == NULL || uCoef.w_Ptr == NULL || uCoef.n_Ptr == NULL ||
		uCoef.s_Ptr == NULL || uCoef.t_Ptr == NULL || uCoef.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (coefPtrs.e_Ptr == NULL || coefPtrs.w_Ptr == NULL || coefPtrs.n_Ptr == NULL ||
		coefPtrs.s_Ptr == NULL || coefPtrs.t_Ptr == NULL || coefPtrs.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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

	Coefficient_Pointers coefPrev;
	coefPrev.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPrev.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPrev.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPrev.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPrev.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	coefPrev.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (coefPrev.e_Ptr == NULL || coefPrev.w_Ptr == NULL || coefPrev.n_Ptr == NULL ||
		coefPrev.s_Ptr == NULL || coefPrev.t_Ptr == NULL || coefPrev.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		coefPrev.e_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPrev.w_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPrev.n_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPrev.s_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPrev.t_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		coefPrev.b_Ptr[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (coefPrev.e_Ptr[k] == NULL || coefPrev.w_Ptr[k] == NULL || coefPrev.n_Ptr[k] == NULL ||
			coefPrev.s_Ptr[k] == NULL || coefPrev.t_Ptr[k] == NULL || coefPrev.b_Ptr[k] == NULL)
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
	if (a_out.e_Ptr == NULL || a_out.w_Ptr == NULL || a_out.n_Ptr == NULL ||
		a_out.s_Ptr == NULL || a_out.t_Ptr == NULL || a_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (a_in.e_Ptr == NULL || a_in.w_Ptr == NULL || a_in.n_Ptr == NULL ||
		a_in.s_Ptr == NULL || a_in.t_Ptr == NULL || a_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (theta_out.e_Ptr == NULL || theta_out.w_Ptr == NULL || theta_out.n_Ptr == NULL ||
		theta_out.s_Ptr == NULL || theta_out.t_Ptr == NULL || theta_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (theta_in.e_Ptr == NULL || theta_in.w_Ptr == NULL || theta_in.n_Ptr == NULL ||
		theta_in.s_Ptr == NULL || theta_in.t_Ptr == NULL || theta_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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

	dataType** n_out = (dataType**)malloc(sizeof(dataType*) * height);
	if (n_out == NULL)
	{
		return false;
	}
	for (k = 0; k < height; k++)
	{
		n_out[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (n_out[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);

	//compute the morm of gradient for the edge detector
	normOfGradientReducedDiamondCells(imageData, gPtrs);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				dataType avg_value = (sqrt(gPtrs.east[k][x]) + sqrt(gPtrs.west[k][x]) + sqrt(gPtrs.north[k][x]) +
					sqrt(gPtrs.south[k][x]) + sqrt(gPtrs.top[k][x]) + sqrt(gPtrs.bottom[k][x])) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(avg_value * avg_value, coef_edge_detector);
				extendedEdge[k_ext][x_ext] = edgeDetectorPtr[k][x];
			}
		}
	}
	reflection3D(extendedEdge, height_ext, length_ext, width_ext);

	//compute the edge detector : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	dataType average_value = 0.0;
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	dataType mp = hx * hy * hz;
	dataType coef_tau = tau / mp;

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				vpe = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext + 1] - extendedEdge[k_ext][x_ext - 1]) / (2 * hx);
				vpw = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext - 1] - extendedEdge[k_ext][x_ext + 1]) / (2 * hx);
				vps = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)]) / (2 * hy);
				vpn = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)]) / (2 * hy);
				vpt = coef_adv * hx * hy * (extendedEdge[k_ext - 1][x_ext] - extendedEdge[k_ext + 1][x_ext]) / (2 * hz);
				vpb = coef_adv * hx * hy * (extendedEdge[k_ext + 1][x_ext] - extendedEdge[k_ext - 1][x_ext]) / (2 * hz);

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

				if(a_out.e_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				if (a_out.w_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				if (a_out.n_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				if (a_out.s_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				if (a_out.t_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				if (a_out.b_Ptr[k][x] != 0)
				{
					n_out[k][x] += 1.0;
				}
				//n_out[k][x] = -(signum(a_out.e_Ptr[k][x]) + signum(a_out.w_Ptr[k][x]) + signum(a_out.n_Ptr[k][x]) + signum(a_out.s_Ptr[k][x]) + signum(a_out.t_Ptr[k][x]) + signum(a_out.b_Ptr[k][x]));

			}
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_iioe.raw");
	strcat_s(name, sizeof(name), name_ending);
	store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);

	copyDataToAnotherArray(initialSegment, segmentationPtr, height, length, width);

	copyDataToExtendedArea(initialSegment, previousSolPtr, height, length, width);
	setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);

	copyDataToExtendedArea(initialSegment, gaussSeidelPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p = 0.0, u_p_min = 0.0, u_p_max = 0.0;
	dataType average_norm_gradient = 0.0, u_average = 0;
	dataType gauss_seidel_coef = 0.0;
	size_t ind_east, ind_west, ind_north, ind_south;
	dataType numerator_max_p = 0.0, numerator_min_p = 0.0;
	size_t count_gauss_seidel_iteration = 0;
	dataType error_gauss_seidel = 0.0;

	Image_Data segmentationData = { height, length, width, segmentationPtr, imageData.origin, imageData.spacing, imageData.orientation };

	dataType u_e_min, u_w_min, u_n_min, u_s_min, u_t_min, u_b_min;
	dataType u_e_max, u_w_max, u_n_max, u_s_max, u_t_max, u_b_max;
	dataType u_east, u_west, u_north, u_south, u_top, u_bottom;

	dataType theta_out_east, theta_out_west, theta_out_north, theta_out_south, theta_out_top, theta_out_bottom;
	dataType theta_in_east, theta_in_west, theta_in_north, theta_in_south, theta_in_top, theta_in_bottom;
	dataType prod_east, prod_west, prod_north, prod_south, prod_top, prod_bottom;
	
	dataType u_east_min, u_west_min, u_north_min, u_south_min, u_top_min, u_bottom_min;
	dataType u_east_max, u_west_max, u_north_max, u_south_max, u_top_max, u_bottom_max;
	
	do {
		number_time_step++;

		normOfGradientReducedDiamondCells(segmentationData, segPtrs);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					//average norm of gradient
					average_norm_gradient = (dataType)((segPtrs.east[k][x] + segPtrs.west[k][x] + segPtrs.north[k][x] +
						segPtrs.south[k][x] + segPtrs.top[k][x] + segPtrs.bottom[k][x]) / 6.0);
					u_average = sqrt(average_norm_gradient + eps2);

					//epsilon regularization
					segPtrs.east[k][x] = sqrt(segPtrs.east[k][x] + eps2);
					segPtrs.west[k][x] = sqrt(segPtrs.west[k][x] + eps2);
					segPtrs.north[k][x] = sqrt(segPtrs.north[k][x] + eps2);
					segPtrs.south[k][x] = sqrt(segPtrs.south[k][x] + eps2);
					segPtrs.top[k][x] = sqrt(segPtrs.top[k][x] + eps2);
					segPtrs.bottom[k][x] = sqrt(segPtrs.bottom[k][x] + eps2);

					u_p = previousSolPtr[k_ext][x_ext];
					u_east = previousSolPtr[k_ext][x_ext + 1];
					u_west = previousSolPtr[k_ext][x_ext - 1];
					u_north = previousSolPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)];
					u_south = previousSolPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)];
					u_top = previousSolPtr[k_ext - 1][x_ext];
					u_bottom = previousSolPtr[k_ext + 1][x_ext];

					//============= Compute theta_out_pq ===========
					u_p_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					u_p_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					numerator_max_p = mp * (u_p_max - u_p);
					numerator_min_p = mp * (u_p_min - u_p);

					if (n_out[k][x] == 0)
					{
						theta_out_east = 0.5;
						theta_out_west = 0.5;
						theta_out_south = 0.5;
						theta_out_north = 0.5;
						theta_out_top = 0.5;
						theta_out_bottom = 0.5;
					}
					else
					{
						//East
						prod_east = a_out.e_Ptr[k][x] * (u_east - u_p);
						if (prod_east == 0)
						{
							theta_out_east = 0.5;
						}
						else if (prod_east > 0)
						{
							theta_out_east = numerator_max_p / (tau * n_out[k][x] * prod_east);
						}
						else
						{
							theta_out_east = numerator_min_p / (tau * n_out[k][x] * prod_east);
						}

						//West
						prod_west = a_out.w_Ptr[k][x] * (u_west - u_p);
						if (prod_west == 0)
						{
							theta_out_west = 0.5;
						}
						else if (prod_west > 0)
						{
							theta_out_west = numerator_max_p / (tau * n_out[k][x] * prod_west);
						}
						else
						{
							theta_out_west = numerator_min_p / (tau * n_out[k][x] * prod_west);
						}

						//North
						prod_north = a_out.n_Ptr[k][x] * (u_north - u_p);
						if (prod_north == 0)
						{
							theta_out_north = 0.5;
						}
						else if (prod_north > 0)
						{
							theta_out_north = numerator_max_p / (tau * n_out[k][x] * prod_north);
						}
						else
						{
							theta_out_north = numerator_min_p / (tau * n_out[k][x] * prod_north);
						}

						//South
						prod_south = a_out.s_Ptr[k][x] * (u_south - u_p);
						if (prod_south == 0)
						{
							theta_out_south = 0.5;
						}
						else if (prod_south > 0)
						{
							theta_out_south = numerator_max_p / (tau * n_out[k][x] * prod_south);
						}
						else
						{
							theta_out_south = numerator_min_p / (tau * n_out[k][x] * prod_south);
						}

						//Top
						prod_top = a_out.t_Ptr[k][x] * (u_top - u_p);
						if (prod_top == 0)
						{
							theta_out_top = 0.5;
						}
						else if (prod_top > 0)
						{
							theta_out_top = numerator_max_p / (tau * n_out[k][x] * prod_top);
						}
						else
						{
							theta_out_top = numerator_min_p / (tau * n_out[k][x] * prod_top);
						}

						//Bottom
						prod_bottom = a_out.b_Ptr[k][x] * (u_bottom - u_p);
						if (prod_bottom == 0)
						{
							theta_out_bottom = 0.5;
						}
						else if (prod_bottom > 0)
						{
							theta_out_bottom = numerator_max_p / (tau * n_out[k][x] * prod_bottom);
						}
						else
						{
							theta_out_bottom = numerator_min_p / (tau * n_out[k][x] * prod_bottom);
						}
					}

					theta_out.e_Ptr[k][x] = fmin(0.5, theta_out_east);
					theta_out.w_Ptr[k][x] = fmin(0.5, theta_out_west);
					theta_out.n_Ptr[k][x] = fmin(0.5, theta_out_north);
					theta_out.s_Ptr[k][x] = fmin(0.5, theta_out_south);
					theta_out.t_Ptr[k][x] = fmin(0.5, theta_out_top);
					theta_out.b_Ptr[k][x] = fmin(0.5, theta_out_bottom);

					//============= Compute theta_in_pq ===========
					//a_in_qp = - a_out_pq 
					//a_out_qp = - a_in_pq
					
					if (n_out[k][x] == 0)
					{
						theta_in_east = 0.5;
						theta_in_west = 0.5;
						theta_in_south = 0.5;
						theta_in_north = 0.5;
						theta_in_top = 0.5;
						theta_in_bottom = 0.5;
					}
					else
					{
						//East
						u_e_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext + 1, j_ext, k_ext);
						u_e_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext + 1, j_ext, k_ext);
						prod_east = -a_in.e_Ptr[k][x] * (u_p - u_east);
						if (prod_east == 0)
						{
							theta_in_east = 0.5;
						}
						else if (prod_east > 0)
						{
							theta_in_east = mp * (u_e_max - u_east) / (tau * n_out[k][x] * prod_east);
						}
						else
						{
							theta_in_east = mp * (u_e_min - u_east) / (tau * n_out[k][x] * prod_east);
						}

						//West
						u_w_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext - 1, j_ext, k_ext);
						u_w_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext - 1, j_ext, k_ext);
						prod_west = -a_in.w_Ptr[k][x] * (u_p - u_west);
						if (prod_west == 0)
						{
							theta_in_west = 0.5;
						}
						else if (prod_west > 0)
						{
							theta_in_west = mp * (u_w_max - u_west) / (tau * n_out[k][x] * prod_west);
						}
						else
						{
							theta_in_west = mp * (u_w_min - u_west) / (tau * n_out[k][x] * prod_west);
						}

						//North
						u_n_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext - 1, k_ext);
						u_n_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext - 1, k_ext);
						prod_north = -a_in.n_Ptr[k][x] * (u_p - u_north);
						if (prod_north == 0)
						{
							theta_in_north = 0.5;
						}
						else if (prod_north > 0)
						{
							theta_in_north = mp * (u_n_max - u_north) / (tau * n_out[k][x] * prod_north);
						}
						else
						{
							theta_in_north = mp * (u_n_min - u_north) / (tau * n_out[k][x] * prod_north);
						}

						//South
						u_s_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext + 1, k_ext);
						u_s_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext + 1, k_ext);
						prod_south = -a_in.s_Ptr[k][x] * (u_p - u_south);
						if (prod_south == 0)
						{
							theta_in_south = 0.5;
						}
						else if (prod_south > 0)
						{
							theta_in_south = mp * (u_s_max - u_south) / (tau * n_out[k][x] * prod_south);
						}
						else
						{
							theta_in_south = mp * (u_s_min - u_south) / (tau * n_out[k][x] * prod_south);
						}

						//Top
						u_t_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext - 1);
						u_t_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext - 1);
						prod_top = -a_in.t_Ptr[k][x] * (u_p - u_top);
						if (prod_top == 0)
						{
							theta_in_top = 0.5;
						}
						else if (prod_top > 0)
						{
							theta_in_top = mp * (u_t_max - u_top) / (tau * n_out[k][x] * prod_top);
						}
						else
						{
							theta_in_top = mp * (u_t_min - u_top) / (tau * n_out[k][x] * prod_top);
						}

						//Bottom
						u_b_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext + 1);
						u_b_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext + 1);
						prod_bottom = -a_in.b_Ptr[k][x] * (u_p - u_bottom);
						if (prod_bottom == 0)
						{
							theta_in_bottom = 0.5;
						}
						else if (prod_bottom > 0)
						{
							theta_in_bottom = mp * (u_b_max - u_bottom) / (tau * n_out[k][x] * prod_bottom);
						}
						else
						{
							theta_in_bottom = mp * (u_b_min - u_bottom) / (tau * n_out[k][x] * prod_bottom);
						}
					}

					theta_in.e_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_east);
					theta_in.w_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_west);
					theta_in.n_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_north);
					theta_in.s_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_south);
					theta_in.t_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_top);
					theta_in.b_Ptr[k][x] = 1.0 - fmin(0.5, theta_in_bottom);

					//============= Compute coefficients ===========

					coefPrev.e_Ptr[k][x] = (dataType)(coef_tau * theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x]);
					coefPrev.w_Ptr[k][x] = (dataType)(coef_tau * theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x]);
					coefPrev.n_Ptr[k][x] = (dataType)(coef_tau * theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x]);
					coefPrev.s_Ptr[k][x] = (dataType)(coef_tau * theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x]);
					coefPrev.t_Ptr[k][x] = (dataType)(coef_tau * theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x]);
					coefPrev.b_Ptr[k][x] = (dataType)(coef_tau * theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x]);

					uCoef.e_Ptr[k][x] = (dataType)(coef_tau * theta_in.e_Ptr[k][x] * a_in.e_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.east[k][x]));
					uCoef.w_Ptr[k][x] = (dataType)(coef_tau * theta_in.w_Ptr[k][x] * a_in.w_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.west[k][x]));
					uCoef.n_Ptr[k][x] = (dataType)(coef_tau * theta_in.n_Ptr[k][x] * a_in.n_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.north[k][x]));
					uCoef.s_Ptr[k][x] = (dataType)(coef_tau * theta_in.s_Ptr[k][x] * a_in.s_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.south[k][x]));
					uCoef.t_Ptr[k][x] = (dataType)(coef_tau * theta_in.t_Ptr[k][x] * a_in.t_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.top[k][x]));
					uCoef.b_Ptr[k][x] = (dataType)(coef_tau * theta_in.b_Ptr[k][x] * a_in.b_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.bottom[k][x]));
				}
			}
		}

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
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);

						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						gauss_seidel_coef = (dataType)(((1 - (coefPrev.e_Ptr[k][x] + coefPrev.w_Ptr[k][x] + coefPrev.n_Ptr[k][x] 
							+ coefPrev.s_Ptr[k][x] + coefPrev.t_Ptr[k][x] + coefPrev.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext]
							+ (coefPrev.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
								+ coefPrev.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
								+ coefPrev.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
								+ coefPrev.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
								+ coefPrev.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
								+ coefPrev.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext])
							+ (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] 
								+ uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] 
								+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north]
								+ uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] 
								+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] 
								+ uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
							/ (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] 
								+ uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]));

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

						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);
						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						u1 = (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]) * gaussSeidelPtr[k_ext][x_ext];
						u2 = (1 - (coefPrev.e_Ptr[k][x] + coefPrev.w_Ptr[k][x] + coefPrev.n_Ptr[k][x] 
							+ coefPrev.s_Ptr[k][x] + coefPrev.t_Ptr[k][x] + coefPrev.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext];
						u3 = coefPrev.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
							+ coefPrev.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
							+ coefPrev.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
							+ coefPrev.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
							+ coefPrev.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
							+ coefPrev.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext];
						u4 = uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext];
						error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
					}
				}
			}
		} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

		//rescall to data range 0-1
		rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);

		//compute L2-norm
		error_segmentation = l2normRectangularGrid(previousSolPtr, gaussSeidelPtr, length_ext, width_ext, height_ext, spacing);

		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		copyDataToAnotherArray(gaussSeidelPtr, previousSolPtr, height_ext, length_ext, width_ext);

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
		if (k < height)
		{
			free(segmentationPtr[k]);
			free(edgeDetectorPtr[k]);
			free(gPtrs.east[k]);
			free(gPtrs.west[k]);
			free(gPtrs.north[k]);
			free(gPtrs.south[k]);
			free(gPtrs.top[k]);
			free(gPtrs.bottom[k]);
			free(segPtrs.east[k]);
			free(segPtrs.west[k]);
			free(segPtrs.north[k]);
			free(segPtrs.south[k]);
			free(segPtrs.top[k]);
			free(segPtrs.bottom[k]);
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
			free(n_out[k]);
			free(coefPrev.e_Ptr[k]);
			free(coefPrev.w_Ptr[k]);
			free(coefPrev.n_Ptr[k]);
			free(coefPrev.s_Ptr[k]);
			free(coefPrev.t_Ptr[k]);
			free(coefPrev.b_Ptr[k]);
		}
		free(previousSolPtr[k]);
		free(gaussSeidelPtr[k]);
		free(extendedEdge[k]);
	}
	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(previousSolPtr);
	free(gaussSeidelPtr);
	free(extendedEdge);
	free(gPtrs.east);
	free(gPtrs.west);
	free(gPtrs.north);
	free(gPtrs.south);
	free(gPtrs.top);
	free(gPtrs.bottom);
	free(segPtrs.east);
	free(segPtrs.west);
	free(segPtrs.north);
	free(segPtrs.south);
	free(segPtrs.top);
	free(segPtrs.bottom);
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
	free(n_out);
	free(coefPrev.e_Ptr);
	free(coefPrev.w_Ptr);
	free(coefPrev.n_Ptr);
	free(coefPrev.s_Ptr);
	free(coefPrev.t_Ptr);
	free(coefPrev.b_Ptr);

	return true;
}

bool GSUBSURF_S_TWO_IIOE(Image_Data imageData, dataType** initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
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

	size_t dim2D = length * width;
	size_t dim2D_ext = length_ext * width_ext;

	dataType hx = imageData.spacing.sx;
	dataType hy = imageData.spacing.sy;
	dataType hz = imageData.spacing.sz;
	dataType hx2 = hx * hx, hy2 = hy * hy, hz2 = hz * hz;
	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps2 = seg_parms.eps2;
	dataType coef_diff = seg_parms.coef_dif, coef_adv = seg_parms.coef_conv;

	dataType** segmentationPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** edgeDetectorPtr = (dataType**)malloc(sizeof(dataType*) * height);
	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** extendedEdge = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL ||
		gaussSeidelPtr == NULL || previousSolPtr == NULL || extendedEdge == NULL)
	{
		return false;
	}
	for (k = 0; k < height_ext; k++)
	{
		if (k < height)
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
		extendedEdge[k] = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
		if (gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL || extendedEdge[k] == NULL)
		{
			return false;
		}
	}

	Pointers_Neighbours gPtrs;
	gPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	gPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (gPtrs.east == NULL || gPtrs.west == NULL || gPtrs.north == NULL ||
		gPtrs.south == NULL || gPtrs.top == NULL || gPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		gPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		gPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (gPtrs.east[k] == NULL || gPtrs.west[k] == NULL || gPtrs.north[k] == NULL ||
			gPtrs.south[k] == NULL || gPtrs.top[k] == NULL || gPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Pointers_Neighbours segPtrs;
	segPtrs.east = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.west = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.north = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.south = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.top = (dataType**)malloc(sizeof(dataType*) * height);
	segPtrs.bottom = (dataType**)malloc(sizeof(dataType*) * height);
	if (segPtrs.east == NULL || segPtrs.west == NULL || segPtrs.north == NULL ||
		segPtrs.south == NULL || segPtrs.top == NULL || segPtrs.bottom == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
	{
		segPtrs.east[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.west[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.north[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.south[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.top[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		segPtrs.bottom[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (segPtrs.east[k] == NULL || segPtrs.west[k] == NULL || segPtrs.north[k] == NULL ||
			segPtrs.south[k] == NULL || segPtrs.top[k] == NULL || segPtrs.bottom[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	Coefficient_Pointers uCoef;
	uCoef.e_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.w_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.n_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.s_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.t_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	uCoef.b_Ptr = (dataType**)malloc(sizeof(dataType*) * height);
	if (uCoef.e_Ptr == NULL || uCoef.w_Ptr == NULL || uCoef.n_Ptr == NULL ||
		uCoef.s_Ptr == NULL || uCoef.t_Ptr == NULL || uCoef.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (coefPtrs.e_Ptr == NULL || coefPtrs.w_Ptr == NULL || coefPtrs.n_Ptr == NULL ||
		coefPtrs.s_Ptr == NULL || coefPtrs.t_Ptr == NULL || coefPtrs.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (a_out.e_Ptr == NULL || a_out.w_Ptr == NULL || a_out.n_Ptr == NULL ||
		a_out.s_Ptr == NULL || a_out.t_Ptr == NULL || a_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (a_in.e_Ptr == NULL || a_in.w_Ptr == NULL || a_in.n_Ptr == NULL ||
		a_in.s_Ptr == NULL || a_in.t_Ptr == NULL || a_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (theta_out.e_Ptr == NULL || theta_out.w_Ptr == NULL || theta_out.n_Ptr == NULL ||
		theta_out.s_Ptr == NULL || theta_out.t_Ptr == NULL || theta_out.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++)
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
	if (theta_in.e_Ptr == NULL || theta_in.w_Ptr == NULL || theta_in.n_Ptr == NULL ||
		theta_in.s_Ptr == NULL || theta_in.t_Ptr == NULL || theta_in.b_Ptr == NULL)
	{
		return false; // Memory allocation failed
	}
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
	if (n_out_pq == NULL || n_out_qp == NULL)
	{
		return false;
	}
	for (k = 0; k < height; k++)
	{
		n_out_pq[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		n_out_qp[k] = (dataType*)malloc(sizeof(dataType) * dim2D);
		if (n_out_pq[k] == NULL || n_out_qp[k] == NULL)
		{
			return false; // Memory allocation failed
		}
	}

	////smoothing
	//heatImplicitRectangularScheme(imageData, smooth_parms);

	//compute the morm of gradient for the edge detector
	normOfGradientReducedDiamondCells(imageData, gPtrs);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				dataType avg_value = (sqrt(gPtrs.east[k][x]) + sqrt(gPtrs.west[k][x]) + sqrt(gPtrs.north[k][x]) +
					sqrt(gPtrs.south[k][x]) + sqrt(gPtrs.top[k][x]) + sqrt(gPtrs.bottom[k][x])) / 6.0;
				edgeDetectorPtr[k][x] = gradientFunction(avg_value * avg_value, coef_edge_detector);
				extendedEdge[k_ext][x_ext] = edgeDetectorPtr[k][x];
			}
		}
	}
	reflection3D(extendedEdge, height_ext, length_ext, width_ext);

	//compute the edge detector : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	dataType average_value = 0.0;
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps, vpt, vpb;
	dataType mp = hx * hy * hz;
	dataType mq = hx * hy * hz;
	dataType coef_tau = tau / mp;

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				vpe = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext + 1] - extendedEdge[k_ext][x_ext - 1]) / (2 * hx);
				vpw = coef_adv * hy * hz * (extendedEdge[k_ext][x_ext - 1] - extendedEdge[k_ext][x_ext + 1]) / (2 * hx);
				vps = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)]) / (2 * hy);
				vpn = coef_adv * hx * hz * (extendedEdge[k_ext][x_new(i_ext, j_ext - 1, length_ext)] - extendedEdge[k_ext][x_new(i_ext, j_ext + 1, length_ext)]) / (2 * hy);
				vpt = coef_adv * hx * hy * (extendedEdge[k_ext - 1][x_ext] - extendedEdge[k_ext + 1][x_ext]) / (2 * hz);
				vpb = coef_adv * hx * hy * (extendedEdge[k_ext + 1][x_ext] - extendedEdge[k_ext - 1][x_ext]) / (2 * hz);

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

				////a_out_pq = -a_in_pq
				n_out_qp[k][x] = -(signum(-a_in.e_Ptr[k][x]) + signum(-a_in.w_Ptr[k][x]) + signum(-a_in.n_Ptr[k][x]) + signum(-a_in.s_Ptr[k][x]) + signum(-a_in.t_Ptr[k][x]) + signum(-a_in.b_Ptr[k][x]));
			}
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_iioe.raw");
	strcat_s(name, sizeof(name), name_ending);
	store3dDataArrayD(edgeDetectorPtr, length, width, height, name, flags);

	copyDataToAnotherArray(initialSegment, segmentationPtr, height, length, width);

	copyDataToExtendedArea(initialSegment, previousSolPtr, height, length, width);
	setBoundaryToZeroDirichletBC(previousSolPtr, length_ext, width_ext, height_ext);

	copyDataToExtendedArea(initialSegment, gaussSeidelPtr, height, length, width);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p = 0.0, u_p_min = 0.0, u_p_max = 0.0;
	dataType average_norm_gradient = 0.0, u_average = 0;
	dataType gauss_seidel_coef = 0.0;
	size_t ind_east, ind_west, ind_north, ind_south;
	dataType numerator_max_p = 0.0, numerator_min_p = 0.0;
	size_t count_gauss_seidel_iteration = 0;
	dataType error_gauss_seidel = 0.0;
	dataType prod_pq = 0.0, prod_qp = 0.0;
	dataType value_pq = 0.0, value_qp = 0.0;

	Image_Data segmentationData = { height, length, width, segmentationPtr, imageData.origin, imageData.spacing, imageData.orientation };

	dataType u_e_min, u_w_min, u_n_min, u_s_min, u_t_min, u_b_min;
	dataType u_e_max, u_w_max, u_n_max, u_s_max, u_t_max, u_b_max;
	dataType u_east, u_west, u_north, u_south, u_top, u_bottom;
	do {
		number_time_step++;

		normOfGradientReducedDiamondCells(segmentationData, segPtrs);
		for (k = 0, k_ext = 1; k < height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++)
				{
					x = x_new(i, j, length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					//epsilon regularization
					segPtrs.east[k][x] = sqrt(segPtrs.east[k][x] + eps2);
					segPtrs.west[k][x] = sqrt(segPtrs.west[k][x] + eps2);
					segPtrs.north[k][x] = sqrt(segPtrs.north[k][x] + eps2);
					segPtrs.south[k][x] = sqrt(segPtrs.south[k][x] + eps2);
					segPtrs.top[k][x] = sqrt(segPtrs.top[k][x] + eps2);
					segPtrs.bottom[k][x] = sqrt(segPtrs.bottom[k][x] + eps2);

					//average norm of gradient
					average_norm_gradient = (segPtrs.east[k][x] + segPtrs.west[k][x] + segPtrs.north[k][x] +
						segPtrs.south[k][x] + segPtrs.top[k][x] + segPtrs.bottom[k][x]) / 6.0;
					u_average = sqrt(average_norm_gradient * average_norm_gradient + eps2);

					u_p = previousSolPtr[k_ext][x_ext];
					u_east = previousSolPtr[k_ext][x_ext + 1];
					u_west = previousSolPtr[k_ext][x_ext - 1];
					u_north = previousSolPtr[k_ext][x_new(i_ext, j_ext - 1, length_ext)];
					u_south = previousSolPtr[k_ext][x_new(i_ext, j_ext + 1, length_ext)];
					u_top = previousSolPtr[k_ext - 1][x_ext];
					u_bottom = previousSolPtr[k_ext + 1][x_ext];

					//============= Compute theta_out_pq ===========
					u_p_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					u_p_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext);
					numerator_max_p = mp * (u_p_max - u_p);
					numerator_min_p = mp * (u_p_min - u_p);

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
						prod_pq = a_out.e_Ptr[k][x] * (u_east - u_p);
						if (prod_pq == 0)
						{
							theta_out.e_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.e_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.e_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//West
						prod_pq = a_out.w_Ptr[k][x] * (u_west - u_p);
						if (prod_pq == 0)
						{
							theta_out.w_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.w_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.w_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//North
						prod_pq = a_out.n_Ptr[k][x] * (u_north - u_p);
						if (prod_pq == 0)
						{
							theta_out.n_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.n_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.n_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//South
						prod_pq = a_out.s_Ptr[k][x] * (u_south - u_p);
						if (prod_pq == 0)
						{
							theta_out.s_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.s_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.s_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//Top
						prod_pq = a_out.t_Ptr[k][x] * (u_top - u_p);
						if (prod_pq == 0)
						{
							theta_out.t_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.t_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.t_Ptr[k][x] = fmin(0.5, value_pq);
						}

						//Bottom
						prod_pq = a_out.b_Ptr[k][x] * (u_bottom - u_p);
						if (prod_pq == 0)
						{
							theta_out.b_Ptr[k][x] = 0.5;
						}
						else if (prod_pq > 0)
						{
							value_pq = numerator_max_p / (tau * n_out_pq[k][x] * prod_pq);
							theta_out.b_Ptr[k][x] = fmin(0.5, value_pq);
						}
						else
						{
							value_pq = numerator_min_p / (tau * n_out_pq[k][x] * prod_pq);
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
						u_e_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext + 1, j_ext, k_ext);
						u_e_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext + 1, j_ext, k_ext);
						prod_qp = -a_in.e_Ptr[k][x] * (u_p - u_east);
						if (prod_qp == 0)
						{
							theta_in.e_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_e_max - u_east) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.e_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_e_min - u_east) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.e_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}

						//West
						u_w_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext - 1, j_ext, k_ext);
						u_w_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext - 1, j_ext, k_ext);
						prod_qp = -a_in.w_Ptr[k][x] * (u_p - u_west);
						if (prod_qp == 0)
						{
							theta_in.w_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_w_max - u_west) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.w_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_w_min - u_west) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.w_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}

						//North
						u_n_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext - 1, k_ext);
						u_n_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext - 1, k_ext);
						prod_qp = -a_in.n_Ptr[k][x] * (u_p - u_north);
						if (prod_qp == 0)
						{
							theta_in.n_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_n_max - u_north) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.n_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_n_min - u_north) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.n_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}

						//South
						u_s_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext + 1, k_ext);
						u_s_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext + 1, k_ext);
						prod_qp = -a_in.s_Ptr[k][x] * (u_p - u_south);
						if (prod_qp == 0)
						{
							theta_in.s_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_s_max - u_south) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.s_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_s_min - u_south) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.s_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}

						//Top
						u_t_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext - 1);
						u_t_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext - 1);
						prod_qp = -a_in.t_Ptr[k][x] * (u_p - u_top);
						if (prod_qp == 0)
						{
							theta_in.t_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_t_max - u_top) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.t_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_t_min - u_top) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.t_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}

						//Bottom
						u_b_min = getMinInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext + 1);
						u_b_max = getMaxInNeighborhood3D(previousSolPtr, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext + 1);
						prod_qp = -a_in.b_Ptr[k][x] * (u_p - u_bottom);
						if (prod_qp == 0)
						{
							theta_in.b_Ptr[k][x] = 0.5;
						}
						else if (prod_qp > 0)
						{
							value_qp = mq * (u_b_max - u_bottom) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.b_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
						else
						{
							value_qp = mq * (u_b_min - u_bottom) / (tau * n_out_qp[k][x] * prod_qp);
							theta_in.b_Ptr[k][x] = 1 - fmin(0.5, value_qp);
						}
					}

					uCoef.e_Ptr[k][x] = coef_tau * theta_in.e_Ptr[k][x] * a_in.e_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.east[k][x]);
					uCoef.w_Ptr[k][x] = coef_tau * theta_in.w_Ptr[k][x] * a_in.w_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hx2 * segPtrs.west[k][x]);
					uCoef.n_Ptr[k][x] = coef_tau * theta_in.n_Ptr[k][x] * a_in.n_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.north[k][x]);
					uCoef.s_Ptr[k][x] = coef_tau * theta_in.s_Ptr[k][x] * a_in.s_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hy2 * segPtrs.south[k][x]);
					uCoef.t_Ptr[k][x] = coef_tau * theta_in.t_Ptr[k][x] * a_in.t_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.top[k][x]);
					uCoef.b_Ptr[k][x] = coef_tau * theta_in.b_Ptr[k][x] * a_in.b_Ptr[k][x] + tau * coef_diff * u_average * edgeDetectorPtr[k][x] / (hz2 * segPtrs.bottom[k][x]);
				}
			}
		}

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
						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);

						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						gauss_seidel_coef = (dataType)(((1 - coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] + theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] + theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] + theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext]
							+ coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
								+ theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
								+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
								+ theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
								+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
								+ theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext])
							+ (uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] + uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north]
								+ uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] + uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
							/ (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]));

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

						x = x_new(i, j, length);
						x_ext = x_new(i_ext, j_ext, length_ext);
						ind_east = x_ext + 1;
						ind_west = x_ext - 1;
						ind_north = x_new(i_ext, j_ext - 1, length_ext);
						ind_south = x_new(i_ext, j_ext + 1, length_ext);

						u1 = (1.0 + uCoef.e_Ptr[k][x] + uCoef.w_Ptr[k][x] + uCoef.n_Ptr[k][x] + uCoef.s_Ptr[k][x] + uCoef.t_Ptr[k][x] + uCoef.b_Ptr[k][x]) * gaussSeidelPtr[k_ext][x_ext];
						u2 = (1 - coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] + theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] + theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] + theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext];
						u3 = coef_tau * (theta_out.e_Ptr[k][x] * a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
							+ theta_out.w_Ptr[k][x] * a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
							+ theta_out.n_Ptr[k][x] * a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
							+ theta_out.s_Ptr[k][x] * a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
							+ theta_out.t_Ptr[k][x] * a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
							+ theta_out.b_Ptr[k][x] * a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext]);
						u4 = uCoef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + uCoef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
							+ uCoef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + uCoef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
							+ uCoef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + uCoef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext];
						error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
					}
				}
			}
		} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

		//rescall to data range 0-1
		rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);

		//compute L2-norm
		error_segmentation = l2normRectangularGrid(previousSolPtr, gaussSeidelPtr, length_ext, width_ext, height_ext, spacing);

		setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);
		copyDataToAnotherArray(gaussSeidelPtr, previousSolPtr, height_ext, length_ext, width_ext);

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
		if (k < height)
		{
			free(segmentationPtr[k]);
			free(edgeDetectorPtr[k]);
			free(gPtrs.east[k]);
			free(gPtrs.west[k]);
			free(gPtrs.north[k]);
			free(gPtrs.south[k]);
			free(gPtrs.top[k]);
			free(gPtrs.bottom[k]);
			free(segPtrs.east[k]);
			free(segPtrs.west[k]);
			free(segPtrs.north[k]);
			free(segPtrs.south[k]);
			free(segPtrs.top[k]);
			free(segPtrs.bottom[k]);
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
		free(extendedEdge[k]);
	}
	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(previousSolPtr);
	free(gaussSeidelPtr);
	free(extendedEdge);
	free(gPtrs.east);
	free(gPtrs.west);
	free(gPtrs.north);
	free(gPtrs.south);
	free(gPtrs.top);
	free(gPtrs.bottom);
	free(segPtrs.east);
	free(segPtrs.west);
	free(segPtrs.north);
	free(segPtrs.south);
	free(segPtrs.top);
	free(segPtrs.bottom);
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

bool gsubsurf_s_one_iioe_time_step(Image_Data segmentationData, Coefficient_Pointers a_out, Coefficient_Pointers coef, Segmentation_Parameters seg_parms)
{
	if(segmentationData.imageDataPtr == NULL || coef.e_Ptr == NULL || 
		coef.w_Ptr == NULL || coef.n_Ptr == NULL || coef.s_Ptr ||
		coef.t_Ptr == NULL || coef.b_Ptr == NULL)
	{
		return true;
	}

	size_t length_ext = segmentationData.length + 2;
	size_t width_ext = segmentationData.width + 2;
	size_t height_ext = segmentationData.height + 2;

	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;

	dataType** gaussSeidelPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	dataType** previousSolPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if(gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}
	for(k = 0; k < height_ext; k++)
	{
		gaussSeidelPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		previousSolPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if(gaussSeidelPtr[k] == NULL || previousSolPtr[k] == NULL)
		{
			return false;
		}
	}

	copyDataToExtendedArea(segmentationData.imageDataPtr, gaussSeidelPtr, segmentationData.height, segmentationData.length, segmentationData.width);
	setBoundaryToZeroDirichletBC(gaussSeidelPtr, length_ext, width_ext, height_ext);

	//gauss seidel for segmentation function
	size_t count_gauss_seidel_iteration = 0;
	size_t ind_east, ind_west, ind_north, ind_south;
	dataType gauss_seidel_coef = 0;
	dataType tau = seg_parms.tau;
	dataType mp = segmentationData.spacing.sx * segmentationData.spacing.sy * segmentationData.spacing.sz;
	dataType coef_tau = seg_parms.tau;
	dataType omega = seg_parms.omega_c;
	dataType error_gauss_seidel = 0.0;

	do {
		count_gauss_seidel_iteration++;
		for (k = 0, k_ext = 1; k < segmentationData.height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < segmentationData.length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < segmentationData.width; j++, j_ext++)
				{
					x = x_new(i, j, segmentationData.length);
					x_ext = x_new(i_ext, j_ext, length_ext);

					ind_east = x_ext + 1;
					ind_west = x_ext - 1;
					ind_north = x_new(i_ext, j_ext - 1, length_ext);
					ind_south = x_new(i_ext, j_ext + 1, length_ext);

					gauss_seidel_coef = (dataType)(((1 - 0.5 * coef_tau * (a_out.e_Ptr[k][x] + a_out.w_Ptr[k][x]
						+ a_out.n_Ptr[k][x] + a_out.s_Ptr[k][x] + a_out.t_Ptr[k][x] + a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext]
						+ 0.5 * coef_tau * (a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
							+ a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
							+ a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
							+ a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
							+ a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
							+ a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext])
						+ (coef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + coef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west] + coef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north]
							+ coef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south] + coef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + coef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext]))
						/ (1.0 + coef.e_Ptr[k][x] + coef.w_Ptr[k][x] + coef.n_Ptr[k][x] + coef.s_Ptr[k][x] + coef.t_Ptr[k][x] + coef.b_Ptr[k][x]));

					gaussSeidelPtr[k_ext][x_ext] = gaussSeidelPtr[k_ext][x_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[k_ext][x_ext]);
				}
			}
		}

		error_gauss_seidel = 0.0;
		for (k = 0, k_ext = 1; k < segmentationData.height; k++, k_ext++)
		{
			for (i = 0, i_ext = 1; i < segmentationData.length; i++, i_ext++)
			{
				for (j = 0, j_ext = 1; j < segmentationData.width; j++, j_ext++)
				{

					x = x_new(i, j, segmentationData.length);
					x_ext = x_new(i_ext, j_ext, length_ext);
					ind_east = x_ext + 1;
					ind_west = x_ext - 1;
					ind_north = x_new(i_ext, j_ext - 1, length_ext);
					ind_south = x_new(i_ext, j_ext + 1, length_ext);

					dataType u1 = (1.0 + coef.e_Ptr[k][x] + coef.w_Ptr[k][x] + coef.n_Ptr[k][x] + coef.s_Ptr[k][x] + coef.t_Ptr[k][x] + coef.b_Ptr[k][x]) * gaussSeidelPtr[k_ext][x_ext];
					dataType u2 = (1 - 0.5 * coef_tau * (a_out.e_Ptr[k][x] + a_out.w_Ptr[k][x]
						+ a_out.n_Ptr[k][x] + a_out.s_Ptr[k][x]
						+ a_out.t_Ptr[k][x] + a_out.b_Ptr[k][x])) * previousSolPtr[k_ext][x_ext];
					dataType u3 = 0.5 * coef_tau * (a_out.e_Ptr[k][x] * previousSolPtr[k_ext][ind_east]
						+ a_out.w_Ptr[k][x] * previousSolPtr[k_ext][ind_west]
						+ a_out.n_Ptr[k][x] * previousSolPtr[k_ext][ind_north]
						+ a_out.s_Ptr[k][x] * previousSolPtr[k_ext][ind_south]
						+ a_out.t_Ptr[k][x] * previousSolPtr[k_ext - 1][x_ext]
						+ a_out.b_Ptr[k][x] * previousSolPtr[k_ext + 1][x_ext]);
					dataType u4 = coef.e_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_east] + coef.w_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_west]
						+ coef.n_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_north] + coef.s_Ptr[k][x] * gaussSeidelPtr[k_ext][ind_south]
						+ coef.t_Ptr[k][x] * gaussSeidelPtr[k_ext - 1][x_ext] + coef.b_Ptr[k][x] * gaussSeidelPtr[k_ext + 1][x_ext];
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}
		}
	} while (count_gauss_seidel_iteration < seg_parms.maxNoGSIteration && error_gauss_seidel > seg_parms.gauss_seidelTolerance);

	rescaleToIntervalZeroOne(gaussSeidelPtr, length_ext, width_ext, height_ext);
	copyDataToReducedArea(segmentationData.imageDataPtr, gaussSeidelPtr, segmentationData.height, segmentationData.length, segmentationData.width);

	for(k = 0; k < height_ext; k++)
	{
		free(gaussSeidelPtr[k]);
		free(previousSolPtr[k]);
	}
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}