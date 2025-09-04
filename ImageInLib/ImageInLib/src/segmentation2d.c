#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include <math.h>
#include "segmentation2d.h"

void initialize2dArrayWithZero(dataType* arrayPtr, const size_t height, const size_t width)
{
	for (size_t i = 0; i < height * width; i++) 
	{
		arrayPtr[i] = 0.0;
	}
}

dataType getMinInNeighborhood(dataType* imageDataPtr, const size_t height, const size_t width, const size_t i, const size_t j)
{
	if(i == 0)
	{
		if(j == 0)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i + 1, j + 1, height)]);
			return fmin(m1, m2);
		}
		else if(j == width - 1)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i + 1, j - 1, height)]);
			return fmin(m1, m2);
		}
		else
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i + 1, j - 1, height)]);
			dataType m3 = fmin(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i + 1, j + 1, height)]);
			return fmin(m1, fmin(m2, m3));
		}
	}
	else if (i == height - 1)
	{
		if (j == 0)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]);
			return fmin(m1, m2);
		}
		else if (j == width - 1)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]);
			return fmin(m1, m2);
		}
		else
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]);
			dataType m3 = fmin(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]);
			return fmin(m1, fmin(m2, m3));
		}
	}
	else 
	{
		if (j == 0)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], fmin(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmin(imageDataPtr[x_new(i, j + 1, height)], fmin(imageDataPtr[x_new(i + 1, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]));
			return fmin(m1, m2);
		}
		else if (j == width - 1)
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], fmin(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], fmin(imageDataPtr[x_new(i + 1, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]));
			return fmin(m1, m2);
		}
		else
		{
			dataType m1 = fmin(imageDataPtr[x_new(i, j, height)], fmin(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmin(imageDataPtr[x_new(i, j - 1, height)], fmin(imageDataPtr[x_new(i + 1, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]));
			dataType m3 = fmin(imageDataPtr[x_new(i, j + 1, height)], fmin(imageDataPtr[x_new(i + 1, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]));
			return fmin(m1, fmin(m2, m3));
		}
	}
}

dataType getMaxInNeighborhood(dataType* imageDataPtr, const size_t height, const size_t width, const size_t i, const size_t j)
{
	if (i == 0)
	{
		if (j == 0)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i + 1, j + 1, height)]);
			return fmax(m1, m2);
		}
		else if (j == width - 1)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i + 1, j - 1, height)]);
			return fmax(m1, m2);
		}
		else
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i + 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i + 1, j - 1, height)]);
			dataType m3 = fmax(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i + 1, j + 1, height)]);
			return fmax(m1, fmax(m2, m3));
		}
	}
	else if (i == height - 1)
	{
		if (j == 0)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]);
			return fmax(m1, m2);
		}
		else if (j == width - 1)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]);
			return fmax(m1, m2);
		}
		else
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], imageDataPtr[x_new(i - 1, j, height)]);
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]);
			dataType m3 = fmax(imageDataPtr[x_new(i, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]);
			return fmax(m1, fmax(m2, m3));
		}
	}
	else
	{
		if (j == 0)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], fmax(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmax(imageDataPtr[x_new(i, j + 1, height)], fmax(imageDataPtr[x_new(i + 1, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]));
			return fmax(m1, m2);
		}
		else if (j == width - 1)
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], fmax(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], fmax(imageDataPtr[x_new(i + 1, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]));
			return fmax(m1, m2);
		}
		else
		{
			dataType m1 = fmax(imageDataPtr[x_new(i, j, height)], fmax(imageDataPtr[x_new(i + 1, j, height)], imageDataPtr[x_new(i - 1, j, height)]));
			dataType m2 = fmax(imageDataPtr[x_new(i, j - 1, height)], fmax(imageDataPtr[x_new(i + 1, j - 1, height)], imageDataPtr[x_new(i - 1, j - 1, height)]));
			dataType m3 = fmax(imageDataPtr[x_new(i, j + 1, height)], fmax(imageDataPtr[x_new(i + 1, j + 1, height)], imageDataPtr[x_new(i - 1, j + 1, height)]));
			return fmax(m1, fmax(m2, m3));
		}
	}
}

dataType l2norm(dataType* arrayPtr1, dataType* arrayPtr2, const size_t height, const size_t width, dataType h) {
	size_t i;

	dataType sumPower = 0.0, norm = 0.0;
	dataType hh = h * h;

	for (i = 0; i < height * width; i++) {
		//sumPower += (dataType)((pow(arrayPtr1[i] - arrayPtr2[i], 2) * hh));
		sumPower += (dataType)(pow(arrayPtr1[i] - arrayPtr2[i], 2));
	}
	norm = sqrt(sumPower);

	return norm;
}

bool rescaleToZeroOne2d(dataType* imageDataPtr, const size_t height, const size_t width)
{
	//check if the memory was allocated successfully
	if (imageDataPtr == NULL)
		return false;

	size_t i, dim2D = height * width;
	dataType max = 0, min = 100000, quotient, offset;

	//Determine minimum and maximum value
	for (i = 0; i < dim2D; i++) {
		if (imageDataPtr[i] < min)
		{
			min = imageDataPtr[i];
		}
			
		if (imageDataPtr[i] > max)
		{
			max = imageDataPtr[i];
		}
	}

	quotient = (dataType)(1.0 / (max - min));
	offset = min * quotient;
	//Rescale values to interval (0, 1)
	for (i = 0; i < dim2D; i++) {
		imageDataPtr[i] = (quotient * imageDataPtr[i] - offset);
	}

	return true;
}

bool generateInitialSegmentationFunction(dataType* imageDataPtr, const size_t height, const size_t width, Point2D* center, dataType v, dataType R)
{
	size_t i, j;
	int dx, dy;
	dataType distance_to_center = 0.0, new_value = 0.0;

	if (imageDataPtr == NULL)
		return false;

	for (i = 0; i < height; i++) {
		dx = i - center->x;
		for (j = 0; j < width; j++) {
			dy = j - center->y;
			distance_to_center = sqrt(dx * dx + dy * dy);
			new_value = (dataType)((1.0 / (distance_to_center + v)) - (1.0 / (R + v)));
			if (distance_to_center > R) {
				imageDataPtr[x_new(i, j, height)] = 0;
			}
			else {
				imageDataPtr[x_new(i, j, height)] = new_value;
			}
		}
	}

	rescaleToZeroOne2d(imageDataPtr, height, width);

	return true;
}

bool set2dDirichletBoundaryCondition(dataType* imageDataPtr, const size_t height, const size_t width) {
	size_t i, j;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			if (i == 0 || i == height - 1 || j == 0 || j == width - 1) {
				imageDataPtr[x_new(i, j, height)] = 0.0;
			}
		}
	}
	return true;
}

bool computeNormOfGradientDiamondCells(dataType* imageDataPtr, neighPtrs neigbours, const size_t height, const size_t width, dataType h) {

	size_t i, j, i_ext, j_ext, xd;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D_ext = height_ext * width_ext;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (extendedArray == NULL)
		return false;

	copyDataTo2dExtendedArea(imageDataPtr, extendedArray, height, width);
	//reflection2D(extendedArray, height_ext, width_ext);
	set2dDirichletBoundaryCondition(extendedArray, height_ext, width_ext);

	dataType uP, uN, uNW, uNE, uS, uSW, uSE, uW, uE;
	dataType ux, uy;

	for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

			size_t iplus = i_ext + 1;
			size_t iminus = i_ext - 1;
			size_t jplus = j_ext + 1;
			size_t jminus = j_ext - 1;

			xd = x_new(i, j, height);
			uP = extendedArray[x_new(i_ext, j_ext, height_ext)];
			uE = extendedArray[x_new(iplus, j_ext, height_ext)];
			uW = extendedArray[x_new(iminus, j_ext, height_ext)];
			uN = extendedArray[x_new(i_ext, jminus, height_ext)];
			uS = extendedArray[x_new(i_ext, jplus, height_ext)];
			uNE = extendedArray[x_new(iplus, jminus, height_ext)];
			uNW = extendedArray[x_new(iminus, jminus, height_ext)];
			uSE = extendedArray[x_new(iplus, jplus, height_ext)];
			uSW = extendedArray[x_new(iminus, jplus, height_ext)];

			//East
			ux = (uE - uP) / h;
			uy = (uNE + uN - uS - uSE) / (4.0 * h);
			neigbours.East[xd] = sqrt(ux * ux + uy * uy);

			//West
			ux = (uP - uW) / h;
			uy = (uNW + uN - uSW - uS) / (4.0 * h);
			neigbours.West[xd] = sqrt(ux * ux + uy * uy);

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * h);
			uy = (uN - uP) / h;
			neigbours.North[xd] = sqrt(ux * ux + uy * uy);

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * h);
			uy = (uP - uS) / h;
			neigbours.South[xd] = sqrt(ux * ux + uy * uy);
		}
	}

	free(extendedArray);

	return true;
}

bool epsilonRegularization(neighPtrs neighbours, const size_t height, const size_t width, dataType epsilon) {
	size_t i, j, currentIndx;
	dataType current = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			currentIndx = x_new(i, j, height);

			current = neighbours.East[currentIndx];
			neighbours.East[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.West[currentIndx];
			neighbours.West[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.North[currentIndx];
			neighbours.North[currentIndx] = (dataType)(sqrt(current * current + epsilon));

			current = neighbours.South[currentIndx];
			neighbours.South[currentIndx] = (dataType)(sqrt(current * current + epsilon));
		}
	}
	return true;
}

bool subsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;

	size_t maxIter = seg_parms.maxNoGSIteration;

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);

	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	if (segmentationPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
		return false;

	heatImplicit2dScheme(imageData, smooth_parms);

	////Save filtered image
	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_filtered.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//store2dRawData(imageData.imageDataPtr, height, width, name, flags);

	dataType* uNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uNorth == NULL || uSouth == NULL || uEast == NULL || uWest == NULL)
		return false;

	neighPtrs U;
	U.West = uWest;
	U.East = uEast;
	U.North = uNorth;
	U.South = uSouth;

	dataType* gNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (gNorth == NULL || gSouth == NULL || gEast == NULL || gWest == NULL)
		return false;

	//Norm of gradient computed on input image for edge detector
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, U, height, width, h);

	////visualize edge detector
	//dataType* edgeAverage = (dataType*)malloc(sizeof(dataType) * dim2D);

	dataType current = 0.0;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t currentIndx = x_new(i, j, height);
			gEast[currentIndx] = gradientFunction(pow(U.East[currentIndx], 2), coef_edge_detector);
			gWest[currentIndx] = gradientFunction(pow(U.West[currentIndx], 2), coef_edge_detector);
			gNorth[currentIndx] = gradientFunction(pow(U.North[currentIndx], 2), coef_edge_detector);
			gSouth[currentIndx] = gradientFunction(pow(U.South[currentIndx], 2), coef_edge_detector);
			//edgeAverage[currentIndx] = (dataType)((gEast[currentIndx] + gWest[currentIndx] + gNorth[currentIndx] + gSouth[currentIndx]) / 4.0);
		}
	}

	////Save edge detector
	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//store2dRawData(edgeAverage, height, width, name, flags);
	//free(edgeAverage);


	dataType* coefNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	dataType average_norm_gradient, u_average;

	copyDataToAnother2dArray(initialSegment, segmentationPtr, height, width);

	copyDataTo2dExtendedArea(initialSegment, previousSolPtr, height, width);
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);
	copyDataTo2dExtendedArea(initialSegment, gaussSeidelPtr, height, width);
	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;

	do {
		number_time_step++;

		//compute the coefficents
		computeNormOfGradientDiamondCells(segmentationPtr, U, height, width, h);
		epsilonRegularization(U, height, width, eps);
		
		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				size_t currentIndx = x_new(i, j, height);

				average_norm_gradient = (dataType)((U.East[currentIndx] + U.West[currentIndx] + U.North[currentIndx] + U.South[currentIndx]) / 4.0);
				u_average = sqrt(average_norm_gradient * average_norm_gradient + eps * eps);

				coefEast[currentIndx] = coef_tau * u_average * gEast[currentIndx] * (1.0 / U.East[currentIndx]);
				coefNorth[currentIndx] = coef_tau * u_average * gNorth[currentIndx] * (1.0 / U.North[currentIndx]);
				coefWest[currentIndx] = coef_tau * u_average * gWest[currentIndx] * (1.0 / U.West[currentIndx]);
				coefSouth[currentIndx] = coef_tau * u_average * gSouth[currentIndx] * (1.0 / U.South[currentIndx]);
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;
		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t xd = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)((previousSolPtr[currentIndx_ext] + coefEast[xd] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[xd] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
						+ coefWest[xd] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[xd] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)])
						/ (1 + coefEast[xd] + coefNorth[xd] + coefWest[xd] + coefSouth[xd]));

					gaussSeidelPtr[currentIndx_ext] = gaussSeidelPtr[currentIndx_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[currentIndx_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t xd = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					error_gauss_seidel += (dataType)(pow((1 + coefEast[xd] + coefNorth[xd] + coefWest[xd] + coefSouth[xd]) * gaussSeidelPtr[currentIndx_ext]
						- (coefEast[xd] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[xd] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
							+ coefWest[xd] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[xd] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)]) - previousSolPtr[currentIndx_ext], 2));
				}
			}
		} while (cpt < maxIter && error_gauss_seidel > tol);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, h);

		//Dirichlet Boundary condition
		set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

		//copy
		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
			printf("Step  %zd : residu = %e \n", number_time_step, error_segmentation);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > tol);

	//FILE* file_peak;
	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "final_segment.csv");
	//strcat_s(name, sizeof(name), name_ending);
	//if (fopen_s(&file_peak, name, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_peak, "x,y\n");
	//for (i = 0; i < dim2D; i++) {
	//	fprintf(file_peak, "%d,%f\n", i, segmentationPtr[i]);
	//}
	//fclose(file_peak);

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

	free(segmentationPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}

bool gsubsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	size_t maxIter = seg_parms.maxNoGSIteration;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);

	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	if (segmentationPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
		return false;

	dataType* uNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* uWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uNorth == NULL || uSouth == NULL || uEast == NULL || uWest == NULL)
		return false;

	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (edgeDetectorPtr == NULL)
		return false;

	neighPtrs uCoef;
	uCoef.West = uWest;
	uCoef.East = uEast;
	uCoef.North = uNorth;
	uCoef.South = uSouth;

	dataType* vNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* vWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (vNorth == NULL || vSouth == NULL || vEast == NULL || vWest == NULL)
		return false;

	dataType* coefNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	dataType current, average_norm_gradient, average_gFunction;
	dataType u_average;

	heatImplicit2dScheme(imageData, smooth_parms);

	//compute g function
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uCoef, height, width, h);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t currentIndx = x_new(i, j, height);
			average_gFunction = (dataType)((uCoef.East[currentIndx] + uCoef.West[currentIndx] + uCoef.North[currentIndx] + uCoef.South[currentIndx]) / 4.0);
			edgeDetectorPtr[currentIndx] = gradientFunction(pow(average_gFunction, 2), coef_edge_detector);
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edgeDetector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t currentIndx = x_new(i, j, height);

			if (i == 0) {
				vEast[currentIndx] = -adv * (edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[currentIndx]);
				vWest[currentIndx] = -vEast[currentIndx];
			}
			else {
				if (i == height - 1) {
					vEast[currentIndx] = -adv * (edgeDetectorPtr[currentIndx] - edgeDetectorPtr[x_new(i - 1, j, height)]);
					vWest[currentIndx] = -vEast[currentIndx];
				}
				else {
					vEast[currentIndx] = -adv * 0.5 * (edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[x_new(i - 1, j, height)]);
					vWest[currentIndx] = -vEast[currentIndx];
				}
			}
			if (j == 0) {
				vSouth[currentIndx] = -adv * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[currentIndx]);
				vNorth[currentIndx] = -vSouth[currentIndx];
			}
			else {
				if (j == width - 1) {
					vSouth[currentIndx] = -adv * (edgeDetectorPtr[currentIndx] - edgeDetectorPtr[x_new(i, j - 1, height)]);
					vNorth[currentIndx] = -vSouth[currentIndx];
				}
				else {
					vSouth[currentIndx] = -adv * 0.5 * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[x_new(i, j - 1, height)]);
					vNorth[currentIndx] = -vSouth[currentIndx];
				}
			}
		}
	}

	copyDataToAnother2dArray(initialSegment, segmentationPtr, height, width);
	copyDataTo2dExtendedArea(initialSegment, gaussSeidelPtr, height, width);
	copyDataTo2dExtendedArea(initialSegment, previousSolPtr, height, width);

	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	do {
		number_time_step++;

		computeNormOfGradientDiamondCells(segmentationPtr, uCoef, height, width, h);
		epsilonRegularization(uCoef, height, width, eps);

		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				size_t currentIndx = x_new(i, j, height);
				average_norm_gradient = (dataType)((uCoef.East[currentIndx] + uCoef.West[currentIndx] + uCoef.North[currentIndx] + uCoef.South[currentIndx]) / 4.0);
				u_average = sqrt(average_norm_gradient * average_norm_gradient + eps);
				coefEast[currentIndx] = (dataType)(coef_tau * (-fmin(vEast[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.East[currentIndx])));
				coefNorth[currentIndx] = (dataType)(coef_tau * (-fmin(vNorth[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.North[currentIndx])));
				coefWest[currentIndx] = (dataType)(coef_tau * (-fmin(vWest[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.West[currentIndx])));
				coefSouth[currentIndx] = (dataType)(coef_tau * (-fmin(vSouth[currentIndx], 0) + diff * edgeDetectorPtr[currentIndx] * u_average * (1.0 / uCoef.South[currentIndx])));
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;

		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)((previousSolPtr[currentIndx_ext] + coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
						+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
						+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) / (1 + coefEast[currentIndx] + coefWest[currentIndx] +
							coefNorth[currentIndx] + coefSouth[currentIndx]));
					gaussSeidelPtr[currentIndx_ext] = gaussSeidelPtr[currentIndx_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[currentIndx_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					error_gauss_seidel += (dataType)(pow((1 + coefEast[currentIndx] + coefWest[currentIndx] + coefNorth[currentIndx] + coefSouth[currentIndx]) * gaussSeidelPtr[currentIndx_ext]
						- (coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
							+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
							+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) - previousSolPtr[currentIndx_ext], 2));
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, h);

		//set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
			printf("Step %zd , residual = %e \n", number_time_step, error_segmentation);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > seg_parms.segTolerance);

	//FILE* file_peak;
	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "final_segment_v1.csv");
	//strcat_s(name, sizeof(name), name_ending);
	//if (fopen_s(&file_peak, name, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_peak, "x,y\n");
	//for (i = 0; i < dim2D; i++) {
	//	fprintf(file_peak, "%d,%f\n", i, segmentationPtr[i]);
	//}
	//fclose(file_peak);

	free(uNorth);
	free(uSouth);
	free(uEast);
	free(uWest);

	free(edgeDetectorPtr);

	free(vNorth);
	free(vSouth);
	free(vEast);
	free(vWest);

	free(coefNorth);
	free(coefSouth);
	free(coefEast);
	free(coefWest);

	free(segmentationPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}

bool gsubsurf_iioe(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;

	dataType tau = seg_parms.tau, h = seg_parms.h;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (h * h), gauss_seidel_coef = 0.0;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	size_t maxIter = seg_parms.maxNoGSIteration;
	
	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}

	neighPtrs uCoef;
	uCoef.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uCoef.North == NULL || uCoef.East == NULL || uCoef.West == NULL || uCoef.South == NULL)
	{
		return false;
	}
		
	neighPtrs uGrad;
	uGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uGrad.North == NULL || uGrad.East == NULL || 
		uGrad.West == NULL || uGrad.South == NULL)
	{
		return false;
	}

	neighPtrs normGrad;
	normGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	normGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	normGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	normGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (normGrad.North == NULL || normGrad.East == NULL || 
		normGrad.West == NULL || normGrad.South == NULL)
	{
		return false;
	}

	neighPtrs a_in;
	a_in.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_in.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_in.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_in.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (a_in.East == NULL || a_in.West == NULL ||
		a_in.North == NULL || a_in.South == NULL)
	{
		return false;
	}

	neighPtrs a_out_pq;
	a_out_pq.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_pq.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_pq.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_pq.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (a_out_pq.East == NULL || a_out_pq.West == NULL ||
		a_out_pq.North == NULL || a_out_pq.South == NULL)
	{
		return false;
	}

	neighPtrs a_out_qp;
	a_out_qp.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_qp.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_qp.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out_qp.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (a_out_qp.East == NULL || a_out_qp.West == NULL ||
		a_out_qp.North == NULL || a_out_qp.South == NULL)
	{
		return false;
	}

	neighPtrs theta_out;
	theta_out.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_out.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_out.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_out.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (theta_out.East == NULL || theta_out.West == NULL || 
		theta_out.North == NULL || theta_out.South == NULL)
	{
		return false;
	}

	neighPtrs theta_in;
	theta_in.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_in.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_in.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	theta_in.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if(theta_in.East == NULL || theta_in.West == NULL || 
		theta_in.North == NULL || theta_in.South == NULL)
	{
		return false;
	}

	dataType* n_out_pq = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* n_out_qp = (dataType*)malloc(sizeof(dataType) * dim2D);

	for (i = 0; i < dim2D; i++) 
	{
		segmentationPtr[i] = 0.0;
		edgeDetectorPtr[i] = 0.0;
		uCoef.North[i] = 0.0;
		uCoef.South[i] = 0.0;
		uCoef.East[i] = 0.0;
		uCoef.West[i] = 0.0;
		uGrad.North[i] = 0.0;
		uGrad.South[i] = 0.0;
		uGrad.East[i] = 0.0;
		uGrad.West[i] = 0.0;
		normGrad.North[i] = 0.0;
		normGrad.South[i] = 0.0;
		normGrad.East[i] = 0.0;
		normGrad.West[i] = 0.0;
		a_in.East[i] = 0.0;
		a_in.West[i] = 0.0;
		a_in.North[i] = 0.0;
		a_in.South[i] = 0.0;
		a_out_pq.East[i] = 0.0;
		a_out_pq.West[i] = 0.0;
		a_out_pq.North[i] = 0.0;
		a_out_pq.South[i] = 0.0;
		a_out_qp.East[i] = 0.0;
		a_out_qp.West[i] = 0.0;
		a_out_qp.North[i] = 0.0;
		a_out_qp.South[i] = 0.0;
		theta_out.East[i] = 0.0;
		theta_out.West[i] = 0.0;
		theta_out.North[i] = 0.0;
		theta_out.South[i] = 0.0;
		n_out_pq[i] = 0.0;
		n_out_qp[i] = 0.0;
	}

	dataType current = 0.0, average_norm_gradient = 0.0, average_gFunction = 0.0;
	dataType u_average = 0.0;

	heatImplicit2dScheme(imageData, smooth_parms);

	//compute g function
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uGrad, height, width, h);
	for (i = 0; i < dim2D; i++) 
	{
		average_gFunction = (dataType)((uGrad.East[i] + uGrad.West[i] + uGrad.North[i] + uGrad.South[i]) / 4.0);
		edgeDetectorPtr[i] = gradientFunction(average_gFunction * average_gFunction, coef_edge_detector);
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edgeDetector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType vpe, vpw, vpn, vps;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			size_t xd = x_new(i, j, height);

			if (i == 0) 
			{
				vpe = adv * (edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[xd]);
				vpw = -vpe;
			}
			else if (i == height - 1)
			{
				vpw = adv * (edgeDetectorPtr[xd] - edgeDetectorPtr[x_new(i - 1, j, height)]);
				vpe = -vpw;
			}
			else 
			{
				vpe = adv * 0.5 * (edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[x_new(i - 1, j, height)]);
				vpw = -vpe;
			}

			if (j == 0) 
			{
				vps = adv * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[xd]);
				vpn = -vps;
			}
			else if (j == width - 1) 
			{
				vpn = adv * (edgeDetectorPtr[xd] - edgeDetectorPtr[x_new(i, j - 1, height)]);
				vps = -vpn;
			}
			else 
			{
				vps = adv * 0.5 * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[x_new(i, j - 1, height)]);
				vpn = -vps;
			}

			a_in.East[xd] = fmax(vpe, 0);
			a_in.West[xd] = fmax(vpw, 0);
			a_in.North[xd] = fmax(vpn, 0);
			a_in.South[xd] = fmax(vps, 0);

			a_out_pq.East[xd] = fmin(vpe, 0);
			a_out_pq.West[xd] = fmin(vpw, 0);
			a_out_pq.North[xd] = fmin(vpn, 0);
			a_out_pq.South[xd] = fmin(vps, 0);

			a_out_qp.East[xd] = -a_in.East[xd];
			a_out_qp.West[xd] = -a_in.West[xd];
			a_out_qp.North[xd] = -a_in.North[xd];
			a_out_qp.South[xd] = -a_in.South[xd];

			n_out_pq[xd] = -(signum(a_out_pq.East[xd]) + signum(a_out_pq.West[xd]) + signum(a_out_pq.North[xd]) + signum(a_out_pq.South[xd]));
			n_out_qp[xd] = -(signum(a_out_qp.East[xd]) + signum(a_out_qp.West[xd]) + signum(a_out_qp.North[xd]) + signum(a_out_qp.South[xd]));
		}
	}

	copyDataToAnother2dArray(initialSegment, segmentationPtr, height, width);
	
	copyDataTo2dExtendedArea(initialSegment, previousSolPtr, height, width);
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);

	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);
	copyDataTo2dExtendedArea(initialSegment, gaussSeidelPtr, height, width);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;
	dataType mp = h * h;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p_min = 0.0, u_p_max = 0.0;
	dataType u_p = 0.0, u_east = 0.0, u_west = 0.0, u_north = 0.0, u_south = 0.0;
	dataType value = 0.0, prod_pq = 0.0, prod_qp = 0.0;
	dataType numerator_max = 0.0, numerator_min = 0.0;
	do {
		number_time_step++;

		computeNormOfGradientDiamondCells(segmentationPtr, normGrad, height, width, h);
		epsilonRegularization(normGrad, height, width, eps);

		for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				size_t xd = x_new(i, j, height);
				
				average_norm_gradient = (dataType)((normGrad.East[xd] + normGrad.West[xd] + normGrad.North[xd] + normGrad.South[xd]) / 4.0);
				u_average = sqrt(pow(average_norm_gradient,2) + eps);

				u_p = previousSolPtr[x_new(i_ext, j_ext, height_ext)];
				u_p_min = getMinInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
				u_p_max = getMaxInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
				numerator_max = mp * (u_p_max - u_p);
				numerator_min = mp * (u_p_min - u_p);
				
				//Compute theta_out_pq
				if (n_out_pq[xd] == 0)
				{
					theta_out.East[xd] = 0.5;
					theta_out.West[xd] = 0.5;
					theta_out.North[xd] = 0.5;
					theta_out.South[xd] = 0.5;
				}
				else 
				{
					//East
					prod_pq = a_out_pq.East[xd] * (previousSolPtr[x_new(i_ext + 1, j_ext, height_ext)] - u_p);
					if (prod_pq == 0)
					{
						theta_out.East[xd] = 0.5;
					}
					else if (prod_pq > 0)
					{
						value = numerator_max / (tau * n_out_pq[xd] * prod_pq);
						theta_out.East[xd] = fmin(0.5, value);
					}
					else
					{
						value = numerator_min / (tau * n_out_pq[xd] * prod_pq);
						theta_out.East[xd] = fmin(0.5, value);
					}

					//West
					prod_pq = a_out_pq.West[xd] * (previousSolPtr[x_new(i_ext - 1, j_ext, height_ext)] - u_p);
					if (prod_pq == 0)
					{
						theta_out.West[xd] = 0.5;
					}
					else if (prod_pq > 0)
					{
						value = numerator_max / (tau * n_out_pq[xd] * prod_pq);
						theta_out.West[xd] = fmin(0.5, value);
					}
					else 
					{
						value = numerator_min / (tau * n_out_pq[xd] * prod_pq);
						theta_out.West[xd] = fmin(0.5, value);
					}

					//North
					prod_pq = a_out_pq.North[xd] * (previousSolPtr[x_new(i_ext, j_ext - 1, height_ext)] - u_p);
					if (prod_pq == 0)
					{
						theta_out.North[xd] = 0.5;
					}
					else if (prod_pq > 0)
					{
						value = numerator_max / (tau * n_out_pq[xd] * prod_pq);
						theta_out.North[xd] = fmin(0.5, value);
					}
					else 
					{
						value = numerator_min / (tau * n_out_pq[xd] * prod_pq);
						theta_out.North[xd] = fmin(0.5, value);
					}

					//South
					prod_pq = a_out_pq.South[xd] * (previousSolPtr[x_new(i_ext, j_ext + 1, height_ext)] - u_p);
					if (prod_pq == 0)
					{
						theta_out.South[xd] = 0.5;
					}
					else if (prod_pq > 0)
					{
						value = numerator_max / (tau * n_out_pq[xd] * prod_pq);
						theta_out.South[xd] = fmin(0.5, value);
					}
					else 
					{
						value = numerator_min / (tau * n_out_pq[xd] * prod_pq);
						theta_out.South[xd] = fmin(0.5, value);
					}
				}

				//Compute thata_in_pq = 1 - theta_out_pq
				if (n_out_qp[xd] == 0)
				{
					theta_in.East[xd] = 0.5;
					theta_in.West[xd] = 0.5;
					theta_in.North[xd] = 0.5;
					theta_in.South[xd] = 0.5;
				}
				else
				{
					//East
					prod_qp = a_out_qp.East[xd] * (u_p - previousSolPtr[x_new(i_ext + 1, j_ext, height_ext)]);
					if (prod_qp == 0)
					{
						theta_in.East[xd] = 0.5;
					}
					else if (prod_qp > 0)
					{
						value = numerator_max / (tau * n_out_qp[xd] * prod_qp);
						theta_in.East[xd] = 1 - fmin(0.5, value);
					}
					else
					{
						value = numerator_min / (tau * n_out_qp[xd] * prod_qp);
						theta_in.East[xd] = 1 - fmin(0.5, value);
					}

					//West
					prod_qp = a_out_qp.West[xd] * (u_p - previousSolPtr[x_new(i_ext - 1, j_ext, height_ext)]);
					if (prod_qp == 0)
					{
						theta_in.West[xd] = 0.5;
					}
					else if (prod_qp > 0)
					{
						value = numerator_max / (tau * n_out_qp[xd] * prod_qp);
						theta_in.West[xd] = 1 - fmin(0.5, value);
					}
					else
					{
						value = numerator_min / (tau * n_out_qp[xd] * prod_qp);
						theta_in.West[xd] = 1 - fmin(0.5, value);
					}

					//North
					prod_qp = a_out_qp.North[xd] * (u_p - previousSolPtr[x_new(i_ext, j_ext - 1, height_ext)]);
					if (prod_qp == 0)
					{
						theta_in.North[xd] = 0.5;
					}
					else if (prod_qp > 0)
					{
						value = numerator_max / (tau * n_out_qp[xd] * prod_qp);
						theta_in.North[xd] = 1 - fmin(0.5, value);
					}
					else
					{
						value = numerator_min / (tau * n_out_qp[xd] * prod_qp);
						theta_in.North[xd] = 1 - fmin(0.5, value);
					}

					//South
					prod_qp = a_out_qp.South[xd] * (u_p - previousSolPtr[x_new(i_ext, j_ext + 1, height_ext)]);
					if (prod_qp == 0)
					{
						theta_in.South[xd] = 0.5;
					}
					else if (prod_qp > 0)
					{
						value = numerator_max / (tau * n_out_qp[xd] * prod_qp);
						theta_in.South[xd] = 1 - fmin(0.5, value);
					}
					else
					{
						value = numerator_min / (tau * n_out_qp[xd] * prod_qp);
						theta_in.South[xd] = 1 - fmin(0.5, value);
					}
				}
				
				uCoef.East[xd] = (dataType)(theta_in.East[xd] * a_in.East[xd] + diff * u_average * edgeDetectorPtr[xd] / normGrad.East[xd]);
				uCoef.West[xd] = (dataType)(theta_in.West[xd] * a_in.West[xd] + diff * u_average * edgeDetectorPtr[xd] / normGrad.West[xd]);
				uCoef.North[xd] = (dataType)(theta_in.North[xd] * a_in.North[xd] + diff * u_average * edgeDetectorPtr[xd] / normGrad.North[xd]);
				uCoef.South[xd] = (dataType)(theta_in.South[xd] * a_in.South[xd] + diff * u_average * edgeDetectorPtr[xd] / normGrad.South[xd]);
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;

		dataType error_gauss_seidel = 0.0;
		do 
		{
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t ind_east = x_new(i_ext + 1, j_ext, height_ext);
					size_t ind_west = x_new(i_ext - 1, j_ext, height_ext);
					size_t ind_north = x_new(i_ext, j_ext - 1, height_ext);
					size_t ind_south = x_new(i_ext, j_ext + 1, height_ext);
					size_t xd = x_new(i, j, height);
					size_t xd_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)(((1 - coef_tau * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]
						+ coef_tau * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
							+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
							+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
							+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south])
						+ coef_tau * (uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
							+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]))
						/ (1.0 + coef_tau * (uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd])));
					gaussSeidelPtr[xd_ext] = gaussSeidelPtr[xd_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[xd_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					size_t ind_east = x_new(i_ext + 1, j_ext, height_ext);
					size_t ind_west = x_new(i_ext - 1, j_ext, height_ext);
					size_t ind_north = x_new(i_ext, j_ext - 1, height_ext);
					size_t ind_south = x_new(i_ext, j_ext + 1, height_ext);
					size_t xd = x_new(i, j, height);
					size_t xd_ext = x_new(i_ext, j_ext, height_ext);

					u1 = (1.0 + coef_tau * (uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd])) * gaussSeidelPtr[xd_ext];
					u2 = (1 - coef_tau * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext];
					u3 = coef_tau * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
						+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
						+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
						+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south]);
					u4 = coef_tau * (uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, h);

		set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) {
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
			printf("Step %zd , residual = %e \n", number_time_step, error_segmentation);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > seg_parms.segTolerance);

	free(edgeDetectorPtr);

	free(uCoef.North);
	free(uCoef.South);
	free(uCoef.East);
	free(uCoef.West);

	free(uGrad.North);
	free(uGrad.South);
	free(uGrad.East);
	free(uGrad.West);

	free(normGrad.North);
	free(normGrad.South);
	free(normGrad.East);
	free(normGrad.West);

	free(segmentationPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	free(a_in.East);
	free(a_in.West);
	free(a_in.North);
	free(a_in.South);

	free(a_out_pq.East);
	free(a_out_pq.West);
	free(a_out_pq.North);
	free(a_out_pq.South);

	free(a_out_qp.East);
	free(a_out_qp.West);
	free(a_out_qp.North);
	free(a_out_qp.South);

	free(theta_out.East);
	free(theta_out.West);
	free(theta_out.North);
	free(theta_out.South);

	free(n_out_pq);
	free(n_out_qp);

	return true;
}
