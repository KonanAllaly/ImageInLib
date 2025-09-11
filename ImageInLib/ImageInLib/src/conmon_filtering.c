#pragma warning(disable:6385) // the compiler doesn't understand the indexing

#include <stdio.h>      // Standard lib for input and output functions
#include <stdlib.h>
#include <time.h>
#include <math.h>       // Maths functions i.e. pow, sin, cos
#include <stdbool.h>    // Boolean function bool
#include <string.h>

#include "conmon_filtering.h"

bool normOfGradientReducedDiamondCells(Image_Data inputImageData, Pointers_Neighbours vGrad)
{
	if (inputImageData.imageDataPtr == NULL || vGrad.east == NULL || vGrad.west == NULL ||
		vGrad.south == NULL || vGrad.north == NULL ||
		vGrad.top == NULL || vGrad.bottom == NULL)
	{
		return false; // Memory allocation failed
	}

	size_t length = inputImageData.length;
	size_t length_ext = length + 2;

	size_t width = inputImageData.width;
	size_t width_ext = width + 2;

	size_t height = inputImageData.height;
	size_t height_ext = height + 2;

	dataType hx = inputImageData.spacing.sx;
	dataType hy = inputImageData.spacing.sy;
	dataType hz = inputImageData.spacing.sz;

	size_t i, j, k, x;
	size_t i_ext, j_ext, k_ext, x_ext;

	dataType** extendedCoefPtr = (dataType**)malloc(sizeof(dataType*) * height_ext);
	if (extendedCoefPtr == NULL)
	{
		return false;
	}
	for (k = 0; k < height_ext; k++)
	{
		extendedCoefPtr[k] = (dataType*)malloc(sizeof(dataType) * length_ext * width_ext);
		if (extendedCoefPtr[k] == NULL)
		{
			return false;
		}
	}

	//Copy to extended area
	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = 1; i < length; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				extendedCoefPtr[k_ext][x_new(i_ext, j_ext, length_ext)] = inputImageData.imageDataPtr[k][x_new(i, j, length)];
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
	dataType quotient_x = (dataType)(4.0 * hx);
	dataType quotient_y = (dataType)(4.0 * hy);
	dataType quotient_z = (dataType)(4.0 * hz);

	for (k = 0, k_ext = 1; k < height; k++, k_ext++)
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
				ux = (dataType)((uE - u) / hx);
				uy = (dataType)(((uN + uNE) - (uS + uSE)) / quotient_y);
				uz = (dataType)(((Tu + TuE) - (Bu + BuE)) / quotient_z);
				vGrad.east[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);

				// Calculation of coefficients in West direction
				ux = (dataType)((uW - u) / hx);
				uy = (dataType)(((uNW + uN) - (uSW + uS)) / quotient_y);
				uz = (dataType)(((TuW + Tu) - (BuW + Bu)) / quotient_z);
				vGrad.west[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);

				// Calculation of coefficients in North direction
				ux = (dataType)(((uNE + uE) - (uNW + uW)) / quotient_x);
				uy = (dataType)((uN - u) / hy);
				uz = (dataType)(((TuN + Tu) - (BuN + Bu)) / quotient_z);
				vGrad.north[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);

				// Calculation of coefficients in South direction
				ux = (dataType)(((uE + uSE) - (uW + uSW)) / quotient_x);
				uy = (dataType)((uS - u) / hy);
				uz = (dataType)(((TuS + Tu) - (BuS + Bu)) / quotient_z);
				vGrad.south[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);

				// Calculation of coefficients in Top direction
				ux = (dataType)(((TuE + uE) - (TuW + uW)) / quotient_x);
				uy = (dataType)(((TuN + uN) - (TuS + uS)) / quotient_y);
				uz = (dataType)((Tu - u) / hz);
				vGrad.top[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);

				// Calculation of coefficients in Bottom direction
				ux = (dataType)(((BuW + uW) - (BuE + uE)) / quotient_x);
				uy = (dataType)(((BuN + uN) - (BuS + uS)) / quotient_y);
				uz = (dataType)((Bu - u) / hz);
				vGrad.bottom[k][x] = (dataType)(ux * ux + uy * uy + uz * uz);
			}
		}
	}

	for (k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);

	return true;
}

bool normOfGradientDiamondCells(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Pointers_Neighbours nGrad)
{
	if (imageData == NULL || nGrad.east == NULL || nGrad.west == NULL ||
		nGrad.south == NULL || nGrad.north == NULL ||
		nGrad.top == NULL || nGrad.bottom == NULL)
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

	//TODO : write the code

	for (k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);
	
	return true;
}

bool normOfGradientSplitDiamondCells(dataType** imageData, const size_t length, const size_t width, const size_t height, const dataType h, Pointers_Neighbours nGrad)
{
	if (imageData == NULL || nGrad.east == NULL || nGrad.west == NULL ||
		nGrad.south == NULL || nGrad.north == NULL ||
		nGrad.top == NULL || nGrad.bottom == NULL)
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

	//TODO : write the code

	for (k = 0; k < height_ext; k++)
	{
		free(extendedCoefPtr[k]);
	}
	free(extendedCoefPtr);

	return true;
}
