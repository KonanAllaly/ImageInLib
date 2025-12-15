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

dataType getMinInNeighborhood(dataType* imageDataPtr, const size_t length, const size_t width, const size_t x, const size_t y)
{
	dataType min_value = 10000.0;
	size_t i_min, i_max, j_min, j_max;

	if (x == 0)
	{
		i_min = x;
	}
	else
	{
		i_min = x - 1;
	}
	if (x == length - 1)
	{
		i_max = x;
	}
	else
	{
		i_max = x + 1;
	}

	if (y == 0)
	{
		j_min = y;
	}
	else
	{
		j_min = y - 1;
	}
	if (y == width - 1)
	{
		j_max = y;
	}
	else
	{
		j_max = y + 1;
	}

	for (size_t i = i_min; i <= i_max; i++)
	{
		for (size_t j = j_min; j <= j_max; j++)
		{
			if (imageDataPtr[x_new(i, j, length)] < min_value)
			{
				min_value = imageDataPtr[x_new(i, j, length)];
			}
		}
	}
	return min_value;
}

dataType getMaxInNeighborhood(dataType* imageDataPtr, const size_t length, const size_t width, const size_t x, const size_t y)
{
	dataType max_value = 0.0;
	size_t i_min, i_max, j_min, j_max;

	if (x == 0)
	{
		i_min = x;
	}
	else
	{
		i_min = x - 1;
	}

	if (x == length - 1)
	{
		i_max = x;
	}
	else
	{
		i_max = x + 1;
	}

	if (y == 0)
	{
		j_min = y;
	}
	else
	{
		j_min = y - 1;
	}
	if (y == width - 1)
	{
		j_max = y;
	}
	else
	{
		j_max = y + 1;
	}

	for (size_t i = i_min; i <= i_max; i++)
	{
		for (size_t j = j_min; j <= j_max; j++)
		{
			if (imageDataPtr[x_new(i, j, length)] > max_value)
			{
				max_value = imageDataPtr[x_new(i, j, length)];
			}
		}
	}

	return max_value;
}

dataType l2norm(dataType* arrayPtr1, dataType* arrayPtr2, const size_t height, const size_t width, dataType h) {

	dataType sumPower = 0.0, norm = 0.0;
	dataType hh = h * h;

	for (size_t i = 0; i < height * width; i++) {
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

bool computeNormOfGradientDiamondCells(dataType* imageDataPtr, neighPtrs neigbours, const size_t height, const size_t width, PixelSpacing spacing) {

	size_t i, j, i_ext, j_ext, xd;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D_ext = height_ext * width_ext;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (extendedArray == NULL)
		return false;

	copyDataTo2dExtendedArea(imageDataPtr, extendedArray, height, width);
	reflection2D(extendedArray, height_ext, width_ext);

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
			ux = (uE - uP) / spacing.sx;
			uy = (uNE + uN - uS - uSE) / (4.0 * spacing.sy);
			neigbours.East[xd] = ux * ux + uy * uy;

			//West
			ux = (uP - uW) / spacing.sx;
			uy = (uNW + uN - uSW - uS) / (4.0 * spacing.sy);
			neigbours.West[xd] = ux * ux + uy * uy;

			//North
			ux = (uNE + uE - uNW - uW) / (4.0 * spacing.sx);
			uy = (uN - uP) / spacing.sy;
			neigbours.North[xd] = ux * ux + uy * uy;

			//South
			ux = (uSE + uE - uSW - uW) / (4.0 * spacing.sx);
			uy = (uP - uS) / spacing.sy;
			neigbours.South[xd] = ux * ux + uy * uy;
		}
	}

	free(extendedArray);

	return true;
}

bool subsurf(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	size_t x, x_ext;
	size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;
	dataType hx = spacing.sx * spacing.sx, hy = spacing.sy * spacing.sy;

	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType coef_tau = tau / (spacing.sx * spacing.sy), gauss_seidel_coef = 0.0;

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

	//heatImplicit2dScheme(imageData, smooth_parms);
	////Save filtered image
	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_smoothed.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//store2dRawData(imageData.imageDataPtr, height, width, name, flags);

	neighPtrs uGrad;
	uGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uGrad.West == NULL || uGrad.East == NULL || uGrad.North == NULL || uGrad.South == NULL)
		return false;

	neighPtrs gGrad;
	gGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (gGrad.East == NULL || gGrad.West == NULL || gGrad.North == NULL || gGrad.South == NULL)
		return false;

	dataType* coefNorth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefSouth = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefEast = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* coefWest = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (coefNorth == NULL || coefSouth == NULL || coefEast == NULL || coefWest == NULL)
		return false;

	//Initialize arrays
	for(i = 0; i < dim2D; i++) {
		segmentationPtr[i] = 0.0;
		uGrad.East[i] = 0.0;
		uGrad.West[i] = 0.0;
		uGrad.North[i] = 0.0;
		uGrad.South[i] = 0.0;
		gGrad.East[i] = 0.0;
		gGrad.West[i] = 0.0;
		gGrad.North[i] = 0.0;
		gGrad.South[i] = 0.0;
		coefNorth[i] = 0.0;
		coefSouth[i] = 0.0;
		coefEast[i] = 0.0;
		coefWest[i] = 0.0;
	}

	for(i = 0; i < dim2D_ext; i++) {
		gaussSeidelPtr[i] = 0.0;
		previousSolPtr[i] = 0.0;
	}

	//Norm of gradient computed on input image for edge detector
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, gGrad, height, width, spacing);

	dataType current = 0.0;
	for (i = 0; i < dim2D; i++) 
	{
		gGrad.East[i] = gradientFunction(gGrad.East[i], coef_edge_detector);
		gGrad.West[i] = gradientFunction(gGrad.West[i], coef_edge_detector);
		gGrad.North[i] = gradientFunction(gGrad.North[i], coef_edge_detector);
		gGrad.South[i] = gradientFunction(gGrad.South[i], coef_edge_detector);
	}

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_east.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(gGrad.East, height, width, name, flags);

	dataType average_norm_gradient = 0.0, u_average = 0.0;

	//Copy to arrays
	for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
	{
		for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
		{
			x = x_new(i, j, height);
			x_ext = x_new(i_ext, j_ext, height_ext);
			previousSolPtr[x_ext] = initialSegment[x];
			gaussSeidelPtr[x_ext] = initialSegment[x];
			segmentationPtr[x] = initialSegment[x];
		}
	}
	set2dDirichletBoundaryCondition(previousSolPtr, height_ext, width_ext);
	set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);

	//segmentation loop
	size_t number_time_step = 0;
	dataType error_segmentation = 0.0;

	do {
		number_time_step++;

		//compute the coefficents
		computeNormOfGradientDiamondCells(segmentationPtr, uGrad, height, width, spacing);
		for (i = 0; i < height; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				x = x_new(i, j, height);

				average_norm_gradient = (dataType)((uGrad.East[x] + uGrad.West[x] + uGrad.North[x] + uGrad.South[x]) / 4.0);
				u_average = sqrt(average_norm_gradient * average_norm_gradient + eps);

				//epsilon regularization
				uGrad.East[x] = sqrt(uGrad.East[x] + eps);
				uGrad.West[x] = sqrt(uGrad.West[x] + eps);
				uGrad.North[x] = sqrt(uGrad.North[x] + eps);
				uGrad.South[x] = sqrt(uGrad.South[x] + eps);

				coefEast[x] = (dataType)(tau * u_average * gGrad.East[x] / (hx * uGrad.East[x]));
				coefNorth[x] = (dataType)(tau * u_average * gGrad.West[x] / (hx * uGrad.West[x]));
				coefWest[x] = (dataType)(tau * u_average * gGrad.North[x] / (hy * uGrad.North[x]));
				coefSouth[x] = (dataType)(tau * u_average * gGrad.South[x] / (hy * uGrad.South[x]));
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;
		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
				{
					x = x_new(i, j, height);
					x_ext = x_new(i_ext, j_ext, height_ext);
					gauss_seidel_coef = (dataType)((previousSolPtr[x_ext] + coefEast[x] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[x] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
						+ coefWest[x] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[x] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)])
						/ (1 + coefEast[x] + coefNorth[x] + coefWest[x] + coefSouth[x]));
					gaussSeidelPtr[x_ext] = gaussSeidelPtr[x_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[x_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
				{
					x = x_new(i, j, height);
					x_ext = x_new(i_ext, j_ext, height_ext);
					error_gauss_seidel += (dataType)(pow((1 + coefEast[x] + coefNorth[x] + coefWest[x] + coefSouth[x]) * gaussSeidelPtr[x_ext]
						- (coefEast[x] * gaussSeidelPtr[x_new(i_ext + 1, j_ext, height_ext)] + coefNorth[x] * gaussSeidelPtr[x_new(i_ext, j_ext - 1, height_ext)]
							+ coefWest[x] * gaussSeidelPtr[x_new(i_ext - 1, j_ext, height_ext)] + coefSouth[x] * gaussSeidelPtr[x_new(i_ext, j_ext + 1, height_ext)]) - previousSolPtr[x_ext], 2));
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > tol);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

		set2dDirichletBoundaryCondition(gaussSeidelPtr, height_ext, width_ext);
		copyDataToAnother2dArray(gaussSeidelPtr, previousSolPtr, height_ext, width_ext);

		//copy to reduce array
		copyDataTo2dReducedArea(segmentationPtr, gaussSeidelPtr, height, width);

		//save the solution
		if (number_time_step % seg_parms.mod == 0) 
		{
			strcpy_s(name, sizeof name, segmentPath);
			sprintf_s(name_ending, sizeof(name_ending), "_seg_func_%03zd.raw", number_time_step);
			strcat_s(name, sizeof(name), name_ending);
			store2dRawData(segmentationPtr, height, width, name, flags);
			printf("Step  %zd : residu = %e \n", number_time_step, error_segmentation);
		}

	} while (number_time_step <= seg_parms.maxNoOfTimeSteps && error_segmentation > tol);

	free(uGrad.East);
	free(uGrad.West);
	free(uGrad.North);
	free(uGrad.South);

	free(gGrad.East);
	free(gGrad.South);
	free(gGrad.North);
	free(gGrad.West);

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
	size_t i, j, x, i_ext, j_ext, x_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;

	dataType tau = seg_parms.tau;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	dataType gauss_seidel_coef = 0.0;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	size_t maxIter = seg_parms.maxNoGSIteration;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL)
	{
		return false;
	}
		
	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	if (segmentationPtr == NULL || gaussSeidelPtr == NULL || previousSolPtr == NULL)
		return false;

	neighPtrs gGrad;
	gGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	gGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);

	neighPtrs uGrad;
	uGrad.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uGrad.South = (dataType*)malloc(sizeof(dataType) * dim2D);

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

	//heatImplicit2dScheme(imageData, smooth_parms);

	//compute g function
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, gGrad, height, width, spacing);
	for (i = 0; i < height; i++) 
	{
		for (j = 0; j < width; j++) 
		{
			x = x_new(i, j, height);
			average_gFunction = (sqrt(gGrad.East[x]) + sqrt(gGrad.West[x]) + sqrt(gGrad.North[x]) + sqrt(gGrad.South[x])) / 4.0;
			edgeDetectorPtr[x] = gradientFunction(pow(average_gFunction, 2), coef_edge_detector);
		}
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	//strcpy_s(name, sizeof name, segmentPath);
	//sprintf_s(name_ending, sizeof(name_ending), "_edge_detector_gsubsurf.raw");
	//strcat_s(name, sizeof(name), name_ending);
	//store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	dataType vpe = 0.0, vpw = 0.0, vpn = 0.0, vps = 0.0;
	for (i = 0; i < height; i++) 
	{
		for (j = 0; j < width; j++) 
		{
			
			x = x_new(i, j, height);
			dataType edge_value = edgeDetectorPtr[x];

			////Smooth edge
			//if (i == 0) 
			//{
			//	vpw = -adv * (edgeDetectorPtr[x_new(i + 1, j, height)] - edge_value);
			//	vpe = -vpw;
			//}
			//else {
			//	if (i == height - 1) 
			//	{
			//		vpe = -adv * (edge_value - edgeDetectorPtr[x_new(i - 1, j, height)]);
			//		vpw = -vpe;
			//	}
			//	else 
			//	{
			//		vpe = -adv * 0.5 * (edgeDetectorPtr[x_new(i + 1, j, height)] - edgeDetectorPtr[x_new(i - 1, j, height)]);
			//		vpw = -vpe;
			//	}
			//}
			//if (j == 0) {
			//	vps = -adv * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[x]);
			//	vpn = -vps;
			//}
			//else {
			//	if (j == width - 1) {
			//		vpn = -adv * (edgeDetectorPtr[x] - edgeDetectorPtr[x_new(i, j - 1, height)]);
			//		vps = -vpn;
			//	}
			//	else {
			//		vps = -adv * 0.5 * (edgeDetectorPtr[x_new(i, j + 1, height)] - edgeDetectorPtr[x_new(i, j - 1, height)]);
			//		vpn = -vps;
			//	}
			//}

			//sharp edge
			if (i == 0) 
			{
				vpe = -adv * (edgeDetectorPtr[x_new(i + 1, j, height)] - edge_value);
				vpw = -vpe;
			}
			else if (i == height - 1) 
			{
				vpw = -adv * (edgeDetectorPtr[x_new(i - 1, j, height)] - edge_value);
				vpe = -vpw;
			}
			else 
			{
				vpe = -adv * (edgeDetectorPtr[x_new(i + 1, j, height)] - edge_value);
				vpw = -adv * (edgeDetectorPtr[x_new(i - 1, j, height)] - edge_value);
			}
			
			if (j == 0)
			{
				vpn = -adv * (edgeDetectorPtr[x_new(i, j + 1, height)] - edge_value);
				vps = -vpn;
			}
			else if (j == width - 1)
			{
				vps = -adv * (edgeDetectorPtr[x_new(i, j - 1, height)] - edge_value);
				vpn = -vps;
			}
			else
			{
				vps = -adv * (edgeDetectorPtr[x_new(i, j + 1, height)] - edge_value);
				vpn = -adv * (edgeDetectorPtr[x_new(i, j - 1, height)] - edge_value);
			}
			
			vEast[x] = fmin(vpe, 0.0);
			vWest[x] = fmin(vpw, 0.0);
			vNorth[x] = fmin(vpn, 0.0);
			vSouth[x] = fmin(vps, 0.0);
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
	dataType tau_mp = tau / (spacing.sx * spacing.sy);
	dataType hx = spacing.sx * spacing.sx, hy = spacing.sy * spacing.sy;
	do {
		number_time_step++;

		computeNormOfGradientDiamondCells(segmentationPtr, uGrad, height, width, spacing);
		for (i = 0; i < height; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				x = x_new(i, j, height);
				average_norm_gradient = (dataType)((uGrad.East[x] + uGrad.West[x] + uGrad.North[x] + uGrad.South[x]) / 4.0);
				u_average = sqrt(average_norm_gradient + eps);

				//epsilon regularization
				uGrad.East[x] = sqrt(uGrad.East[x] + eps);
				uGrad.West[x] = sqrt(uGrad.West[x] + eps);
				uGrad.North[x] = sqrt(uGrad.North[x] + eps);
				uGrad.South[x] = sqrt(uGrad.South[x] + eps);

				coefEast[x] = (dataType)(- tau_mp * vEast[x] + tau * diff * edgeDetectorPtr[x] * u_average / (hx * uGrad.East[x]));
				coefWest[x] = (dataType)(- tau_mp * vWest[x] + tau * diff * edgeDetectorPtr[x] * u_average / (hx * uGrad.West[x]));
				coefNorth[x] = (dataType)(- tau_mp * vNorth[x] + tau * diff * edgeDetectorPtr[x] * u_average / (hy * uGrad.North[x]));
				coefSouth[x] = (dataType)(- tau_mp * vSouth[x] + tau * diff * edgeDetectorPtr[x] * u_average / (hy * uGrad.South[x]));
			}
		}

		//gauss seidel for segmentation function
		size_t cpt = 0;

		dataType error_gauss_seidel = 0.0;
		do {
			cpt++;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
				{

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)((previousSolPtr[currentIndx_ext] + coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
						+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
						+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) 
						/ (1 + coefEast[currentIndx] + coefWest[currentIndx] +
							coefNorth[currentIndx] + coefSouth[currentIndx]));
					gaussSeidelPtr[currentIndx_ext] = gaussSeidelPtr[currentIndx_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[currentIndx_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
			{
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
				{

					size_t iplus = i_ext + 1;
					size_t iminus = i_ext - 1;
					size_t jplus = j_ext + 1;
					size_t jminus = j_ext - 1;

					size_t currentIndx = x_new(i, j, height);
					size_t currentIndx_ext = x_new(i_ext, j_ext, height_ext);

					error_gauss_seidel += (dataType)(pow((1 + coefEast[currentIndx] + coefWest[currentIndx] + 
						coefNorth[currentIndx] + coefSouth[currentIndx]) * gaussSeidelPtr[currentIndx_ext]
						- (coefEast[currentIndx] * gaussSeidelPtr[x_new(iplus, j_ext, height_ext)]
							+ coefWest[currentIndx] * gaussSeidelPtr[x_new(iminus, j_ext, height_ext)] + 
							coefNorth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jminus, height_ext)]
							+ coefSouth[currentIndx] * gaussSeidelPtr[x_new(i_ext, jplus, height_ext)]) -
						previousSolPtr[currentIndx_ext], 2));
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

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

	free(uGrad.East);
	free(uGrad.West);
	free(uGrad.North);
	free(uGrad.South);

	free(gGrad.East);
	free(gGrad.West);
	free(gGrad.North);
	free(gGrad.South);

	free(vNorth);
	free(vSouth);
	free(vEast);
	free(vWest);

	free(coefNorth);
	free(coefSouth);
	free(coefEast);
	free(coefWest);

	free(segmentationPtr);
	free(edgeDetectorPtr);
	free(gaussSeidelPtr);
	free(previousSolPtr);

	return true;
}

bool gsubsurf_s_one_iioe(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;
	dataType hx = spacing.sx, hy = spacing.sy;
	dataType hx2 = hx * hx, hy2 = hy * hy;
	
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

	dataType* n_out = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* n_out_qp = (dataType*)malloc(sizeof(dataType) * dim2D);

	//Initialization
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
		n_out[i] = 0.0;
	}

	dataType current = 0.0, average_norm_gradient = 0.0, average_g = 0.0;
	dataType u_average = 0.0;

	heatImplicit2dScheme(imageData, smooth_parms);

	//Compute g function : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uGrad, height, width, spacing);
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	for (i = 0; i < dim2D; i++) 
	{
		average_g = (dataType)((sqrt(uGrad.East[i]) + sqrt(uGrad.West[i]) + sqrt(uGrad.North[i]) + sqrt(uGrad.South[i])) / 4.0);
		edgeDetectorPtr[i] = gradientFunction(average_g * average_g, coef_edge_detector);
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	dataType vpe, vpw, vpn, vps;
	for (i = 0; i < height; i++) 
	{
		for (j = 0; j < width; j++) 
		{
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

			//a_in_pq = max(v_pq, 0)
			a_in.East[xd] = fmax(vpe, 0);
			a_in.West[xd] = fmax(vpw, 0);
			a_in.North[xd] = fmax(vpn, 0);
			a_in.South[xd] = fmax(vps, 0);

			//a_out_pq = min(v_pq, 0)
			a_out_pq.East[xd] = fmin(vpe, 0);
			a_out_pq.West[xd] = fmin(vpw, 0);
			a_out_pq.North[xd] = fmin(vpn, 0);
			a_out_pq.South[xd] = fmin(vps, 0);

			//a_out_qp = -a_in_pq
			a_out_qp.East[xd] = -a_in.East[xd];
			a_out_qp.West[xd] = -a_in.West[xd];
			a_out_qp.North[xd] = -a_in.North[xd];
			a_out_qp.South[xd] = -a_in.South[xd];

			//n_out_pq = -sum(sign(a_out_pq))
			n_out[xd] = 0.0;
			if (a_out_pq.East[xd] != 0.0) 
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.West[xd] != 0.0) 
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.North[xd] != 0.0) 
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.South[xd] != 0.0) 
			{
				n_out[xd] += 1.0;
			}

			//n_out_qp = -sum(sign(a_out_qp))
			n_out_qp[xd] = 0.0;
			if (a_out_qp.East[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.West[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.North[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.South[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
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
	dataType mp = spacing.sx * spacing.sy;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p_min = 0.0, u_p_max = 0.0;
	dataType value = 0.0, prod_pq = 0.0, prod_qp = 0.0;
	dataType numerator_max = 0.0, numerator_min = 0.0;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType gauss_seidel_coef = 0.0;
	size_t maxIter = seg_parms.maxNoGSIteration;
	dataType tau = seg_parms.tau;
	dataType tau_mp = tau / mp;
	
	dataType prod_east = 0.0, prod_west = 0.0, prod_north = 0.0, prod_south = 0.0;
	dataType prod_east_p = 0.0, prod_west_p = 0.0, prod_north_p = 0.0, prod_south_p = 0.0;
	
	dataType theta_east = 0.0, theta_west = 0.0, theta_north = 0.0, theta_south = 0.0;
	dataType theta_east_p = 0.0, theta_west_p = 0.0, theta_north_p = 0.0, theta_south_p = 0.0;
	
	dataType u_p = 0.0, u_east = 0.0, u_west = 0.0, u_north = 0.0, u_south = 0.0;
	
	//dataType u_east_min = 0.0, u_east_max = 0.0;
	//dataType u_west_min = 0.0, u_west_max = 0.0;
	//dataType u_north_min = 0.0, u_north_max = 0.0;
	//dataType u_south_min = 0.0, u_south_max = 0.0;
	
	size_t iplus, iminus, jplus, jminus;
	do {
		number_time_step++;

		//The return norm of gradient is |grad(u)|^2
		computeNormOfGradientDiamondCells(segmentationPtr, normGrad, height, width, spacing);

		for (i = 0, i_ext = 1; i < height; i++, i_ext++) 
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) 
			{
				size_t xd = x_new(i, j, height);
				iplus = i_ext + 1;
				iminus = i_ext - 1;
				jplus = j_ext + 1;
				jminus = j_ext - 1;

				//Compute average of norm of gradient
				average_norm_gradient = (dataType)((normGrad.East[xd] + normGrad.West[xd] + normGrad.North[xd] + normGrad.South[xd]) / 4.0);
				u_average = sqrt(average_norm_gradient + eps);

				//Epsilon regularization
				normGrad.East[xd] = sqrt(normGrad.East[xd] + eps);
				normGrad.West[xd] = sqrt(normGrad.West[xd] + eps);
				normGrad.North[xd] = sqrt(normGrad.North[xd] + eps);
				normGrad.South[xd] = sqrt(normGrad.South[xd] + eps);

				u_p = previousSolPtr[x_new(i_ext, j_ext, height_ext)];
				u_east = previousSolPtr[x_new(iplus, j_ext, height_ext)];
				u_west = previousSolPtr[x_new(iminus, j_ext, height_ext)];
				u_north = previousSolPtr[x_new(i_ext, jminus, height_ext)];
				u_south = previousSolPtr[x_new(i_ext, jplus, height_ext)];
				
				//Compute theta_out_pq
				u_p_min = getMinInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
				u_p_max = getMaxInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
				numerator_max = mp * (u_p_max - u_p);
				numerator_min = mp * (u_p_min - u_p);
				
				if (n_out[xd] == 0)
				{
					theta_east = 0.5;
					theta_west = 0.5;
					theta_north = 0.5;
					theta_south = 0.5;
				}
				else 
				{
					//East
					prod_east = a_out_pq.East[xd] * (u_east - u_p);
					if (prod_east == 0)
					{
						theta_east = 0.5;
					}
					else if (prod_east > 0)
					{
						theta_east = numerator_max / (tau * n_out[xd] * prod_east);
					}
					else
					{
						theta_east = numerator_min / (tau * n_out[xd] * prod_east);
					}
					
					//West
					prod_west = a_out_pq.West[xd] * (u_west - u_p);
					if (prod_west == 0)
					{
						theta_west = 0.5;
					}
					else if (prod_west > 0)
					{
						theta_west = numerator_max / (tau * n_out[xd] * prod_west);
					}
					else 
					{
						theta_west = numerator_min / (tau * n_out[xd] * prod_west);
					}

					//North
					prod_north = a_out_pq.North[xd] * (u_north - u_p);
					if (prod_north == 0)
					{
						theta_north = 0.5;
					}
					else if (prod_north > 0)
					{
						theta_north = numerator_max / (tau * n_out[xd] * prod_north);
					}
					else 
					{
						theta_north = numerator_min / (tau * n_out[xd] * prod_north);
					}

					//South
					prod_south = a_out_pq.South[xd] * (u_south - u_p);
					if (prod_south == 0)
					{
						theta_south = 0.5;
					}
					else if (prod_south > 0)
					{
						theta_south = numerator_max / (tau * n_out[xd] * prod_south);
					}
					else 
					{
						theta_south = numerator_min / (tau * n_out[xd] * prod_south);
					}
				}

				theta_out.East[xd] = fmin(0.5, theta_east);
				theta_out.West[xd] = fmin(0.5, theta_west);
				theta_out.North[xd] = fmin(0.5, theta_north);
				theta_out.South[xd] = fmin(0.5, theta_south);

				//Compute thata_in_pq = 1 - theta_out_pq

				if (n_out_qp[xd] == 0)
				{
					theta_east_p = 0.5;
					theta_west_p = 0.5;
					theta_north_p = 0.5;
					theta_south_p = 0.5;
				}
				else
				{
					//East
					prod_east_p = a_out_qp.East[xd] * (u_east - u_p);
					if (prod_east_p == 0)
					{
						theta_east_p = 0.5;
					}
					else if (prod_east_p > 0)
					{
						theta_east_p = numerator_max / (tau * n_out_qp[xd] * prod_east_p);
					}
					else
					{
						theta_east_p = numerator_min / (tau * n_out_qp[xd] * prod_east_p);
					}

					//West
					prod_west_p = a_out_qp.West[xd] * (u_west - u_p);
					if (prod_west_p == 0)
					{
						theta_west_p = 0.5;
					}
					else if (prod_west_p > 0)
					{
						theta_west_p = numerator_max / (tau * n_out_qp[xd] * prod_west_p);
					}
					else
					{
						theta_west_p = numerator_min / (tau * n_out_qp[xd] * prod_west_p);
					}

					//North
					prod_north_p = a_out_qp.North[xd] * (u_north - u_p);
					if (prod_north_p == 0)
					{
						theta_north_p = 0.5;
					}
					else if (prod_north_p > 0)
					{
						theta_north_p = numerator_max / (tau * n_out_qp[xd] * prod_north_p);
					}
					else
					{
						theta_north_p = numerator_min / (tau * n_out_qp[xd] * prod_north_p);
					}

					//South
					prod_south_p = a_out_qp.South[xd] * (u_south - u_p);
					if (prod_south_p == 0)
					{
						theta_south_p = 0.5;
					}
					else if (prod_south_p > 0)
					{
						theta_south_p = numerator_max / (tau * n_out_qp[xd] * prod_south_p);
					}
					else
					{
						theta_south_p = numerator_min / (tau * n_out_qp[xd] * prod_south_p);
					}
				}

				theta_in.East[xd] = 1.0 - fmin(0.5, theta_east_p);
				theta_in.West[xd] = 1.0 - fmin(0.5, theta_west_p);
				theta_in.North[xd] = 1.0 - fmin(0.5, theta_north_p);
				theta_in.South[xd] = 1.0 - fmin(0.5, theta_south_p);
				
				uCoef.East[xd] = (dataType)(tau_mp * theta_in.East[xd] * a_in.East[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.East[xd]));
				uCoef.West[xd] = (dataType)(tau_mp * theta_in.West[xd] * a_in.West[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.West[xd]));
				uCoef.North[xd] = (dataType)(tau_mp * theta_in.North[xd] * a_in.North[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.North[xd]));
				uCoef.South[xd] = (dataType)(tau_mp * theta_in.South[xd] * a_in.South[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.South[xd]));
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

					gauss_seidel_coef = (dataType)(((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]
						+ tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
							+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
							+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
							+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south])
						+ (uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
							+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]))
						/ (1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]));
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

					u1 = (dataType)((1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]) * gaussSeidelPtr[xd_ext]);
					u2 = (dataType)((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]);
					u3 = (dataType)(tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
						+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
						+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
						+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south]));
					u4 = (dataType)(uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

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

	free(n_out);
	free(n_out_qp);

	return true;
}

bool gsubsurf_implicit(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;
	dataType hx = spacing.sx, hy = spacing.sy;
	dataType hx2 = hx * hx, hy2 = hy * hy;

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

	//Initialization
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
	}

	dataType current = 0.0, average_norm_gradient = 0.0, average_g = 0.0;
	dataType u_average = 0.0;

	heatImplicit2dScheme(imageData, smooth_parms);

	//Compute g function : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uGrad, height, width, spacing);
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	for (i = 0; i < dim2D; i++)
	{
		average_g = (dataType)((sqrt(uGrad.East[i]) + sqrt(uGrad.West[i]) + sqrt(uGrad.North[i]) + sqrt(uGrad.South[i])) / 4.0);
		edgeDetectorPtr[i] = gradientFunction(average_g * average_g, coef_edge_detector);
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	dataType vpe, vpw, vpn, vps;
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
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

			//a_in_pq = max(v_pq, 0)
			a_in.East[xd] = fmax(vpe, 0);
			a_in.West[xd] = fmax(vpw, 0);
			a_in.North[xd] = fmax(vpn, 0);
			a_in.South[xd] = fmax(vps, 0);
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
	dataType mp = spacing.sx * spacing.sy;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType gauss_seidel_coef = 0.0;
	size_t maxIter = seg_parms.maxNoGSIteration;
	dataType tau = seg_parms.tau;
	dataType tau_mp = tau / mp;

	size_t iplus, iminus, jplus, jminus;
	do {
		number_time_step++;

		//The return norm of gradient is |grad(u)|^2
		computeNormOfGradientDiamondCells(segmentationPtr, normGrad, height, width, spacing);

		for (i = 0, i_ext = 1; i < height; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				size_t xd = x_new(i, j, height);

				//Compute average of norm of gradient
				average_norm_gradient = (dataType)((normGrad.East[xd] + normGrad.West[xd] + normGrad.North[xd] + normGrad.South[xd]) / 4.0);
				u_average = sqrt(average_norm_gradient + eps);

				//Epsilon regularization
				normGrad.East[xd] = sqrt(normGrad.East[xd] + eps);
				normGrad.West[xd] = sqrt(normGrad.West[xd] + eps);
				normGrad.North[xd] = sqrt(normGrad.North[xd] + eps);
				normGrad.South[xd] = sqrt(normGrad.South[xd] + eps);

				uCoef.East[xd] = (dataType)(tau_mp * (a_in.East[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.East[xd]));
				uCoef.West[xd] = (dataType)(tau_mp * (a_in.West[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.West[xd]));
				uCoef.North[xd] = (dataType)(tau_mp * (a_in.North[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.North[xd]));
				uCoef.South[xd] = (dataType)(tau_mp * (a_in.South[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.South[xd]));
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

					gauss_seidel_coef = (dataType)((previousSolPtr[xd_ext] + (uCoef.East[xd] * gaussSeidelPtr[ind_east] 
						+ uCoef.West[xd] * gaussSeidelPtr[ind_west] 
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] 
						+ uCoef.South[xd] * gaussSeidelPtr[ind_south]))
						/ (1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]));
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

					u1 = (dataType)((1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]) * gaussSeidelPtr[xd_ext]);
					u2 = (dataType)(previousSolPtr[xd_ext]);
					u3 = (dataType)(uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

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

	return true;
}

bool gsubsurf_iioe(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL || 
		gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}

	neighPtrs uCoef;
	uCoef.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uCoef.North == NULL || uCoef.East == NULL || 
		uCoef.West == NULL || uCoef.South == NULL)
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

	neighPtrs a_out;
	a_out.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	a_out.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (a_out.East == NULL || a_out.West == NULL ||
		a_out.North == NULL || a_out.South == NULL)
	{
		return false;
	}

	//Initialization
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
		a_out.East[i] = 0.0;
		a_out.West[i] = 0.0;
		a_out.North[i] = 0.0;
		a_out.South[i] = 0.0;
	}

	dataType current = 0.0, average_norm_gradient = 0.0, average_g = 0.0;
	dataType u_average = 0.0;

	heatImplicit2dScheme(imageData, smooth_parms);

	//Compute g function : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uGrad, height, width, spacing);
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	for (i = 0; i < dim2D; i++)
	{
		average_g = (dataType)((sqrt(uGrad.East[i]) + sqrt(uGrad.West[i]) + sqrt(uGrad.North[i]) + sqrt(uGrad.South[i])) / 4.0);
		edgeDetectorPtr[i] = gradientFunction(average_g * average_g, coef_edge_detector);
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	dataType vpe, vpw, vpn, vps;
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
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

			//a_in_pq = max(v_pq, 0)
			a_in.East[xd] = fmax(vpe, 0);
			a_in.West[xd] = fmax(vpw, 0);
			a_in.North[xd] = fmax(vpn, 0);
			a_in.South[xd] = fmax(vps, 0);

			//a_out_pq = min(v_pq, 0)
			a_out.East[xd] = fmin(vpe, 0);
			a_out.West[xd] = fmin(vpw, 0);
			a_out.North[xd] = fmin(vpn, 0);
			a_out.South[xd] = fmin(vps, 0);
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
	dataType mp = spacing.sx * spacing.sy;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType gauss_seidel_coef = 0.0;
	size_t maxIter = seg_parms.maxNoGSIteration;
	dataType tau = seg_parms.tau;
	dataType tau_mp = tau / mp;

	size_t xd, xd_ext, ind_east, ind_west, ind_north, ind_south;
	do {
		number_time_step++;

		//The return norm of gradient is |grad(u)|^2
		computeNormOfGradientDiamondCells(segmentationPtr, normGrad, height, width, spacing);

		for (i = 0, i_ext = 1; i < height; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				xd = x_new(i, j, height);

				//Compute average of norm of gradient
				average_norm_gradient = (dataType)((normGrad.East[xd] + normGrad.West[xd] + normGrad.North[xd] + normGrad.South[xd]) / 4.0);
				u_average = sqrt(average_norm_gradient + eps);

				//Epsilon regularization
				normGrad.East[xd] = sqrt(normGrad.East[xd] + eps);
				normGrad.West[xd] = sqrt(normGrad.West[xd] + eps);
				normGrad.North[xd] = sqrt(normGrad.North[xd] + eps);
				normGrad.South[xd] = sqrt(normGrad.South[xd] + eps);

				uCoef.East[xd] = (dataType)(tau_mp * (0.5 * a_in.East[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.East[xd]));
				uCoef.West[xd] = (dataType)(tau_mp * (0.5 * a_in.West[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.West[xd]));
				uCoef.North[xd] = (dataType)(tau_mp * (0.5 * a_in.North[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.North[xd]));
				uCoef.South[xd] = (dataType)(tau_mp * (0.5 * a_in.South[xd] + (diff * u_average * edgeDetectorPtr[xd]) / normGrad.South[xd]));
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

					ind_east = x_new(i_ext + 1, j_ext, height_ext);
					ind_west = x_new(i_ext - 1, j_ext, height_ext);
					ind_north = x_new(i_ext, j_ext - 1, height_ext);
					ind_south = x_new(i_ext, j_ext + 1, height_ext);
					xd = x_new(i, j, height);
					xd_ext = x_new(i_ext, j_ext, height_ext);

					gauss_seidel_coef = (dataType)(((1.0 - 0.5 * tau_mp * (a_out.East[xd] + a_out.West[xd]
						+ a_out.North[xd] + a_out.South[xd])) * previousSolPtr[xd_ext]
						+ 0.5 * tau_mp * (a_out.East[xd] * previousSolPtr[ind_east]
							+ a_out.West[xd] * previousSolPtr[ind_west]
							+ a_out.North[xd] * previousSolPtr[ind_north]
							+ a_out.South[xd] * previousSolPtr[ind_south])
						+ uCoef.East[xd] * gaussSeidelPtr[ind_east]
						+ uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north]
						+ uCoef.South[xd] * gaussSeidelPtr[ind_south])
						/ (1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]));
					gaussSeidelPtr[xd_ext] = gaussSeidelPtr[xd_ext] + omega * (gauss_seidel_coef - gaussSeidelPtr[xd_ext]);
				}
			}

			error_gauss_seidel = 0.0;
			for (i = 0, i_ext = 1; i < height; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {

					ind_east = x_new(i_ext + 1, j_ext, height_ext);
					ind_west = x_new(i_ext - 1, j_ext, height_ext);
					ind_north = x_new(i_ext, j_ext - 1, height_ext);
					ind_south = x_new(i_ext, j_ext + 1, height_ext);
					xd = x_new(i, j, height);
					xd_ext = x_new(i_ext, j_ext, height_ext);

					u1 = (dataType)((1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]) * gaussSeidelPtr[xd_ext]);
					u2 = (dataType)((1.0 - 0.5 * tau_mp * (a_out.East[xd] + a_out.West[xd]
						+ a_out.North[xd] + a_out.South[xd])) * previousSolPtr[xd_ext]);
					u3 = (dataType)(0.5 * tau_mp * (a_out.East[xd] * previousSolPtr[ind_east]
						+ a_out.West[xd] * previousSolPtr[ind_west]
						+ a_out.North[xd] * previousSolPtr[ind_north]
						+ a_out.South[xd] * previousSolPtr[ind_south]));
					u4 = (dataType)(uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

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

	free(a_out.East);
	free(a_out.West);
	free(a_out.North);
	free(a_out.South);

	return true;
}

bool gsubsurf_s_two_iioe(Image_Data2D imageData, dataType* initialSegment, const char* segmentPath, const Filter_Parameters smooth_parms, Segmentation_Parameters seg_parms)
{
	size_t i, j, i_ext, j_ext;
	const size_t height = imageData.height, width = imageData.width;
	const size_t height_ext = height + 2, width_ext = width + 2;
	size_t dim2D = height * width, dim2D_ext = height_ext * width_ext;
	PixelSpacing spacing = imageData.spacing;
	dataType hx = spacing.sx, hy = spacing.sy;
	dataType hx2 = hx * hx, hy2 = hy * hy;

	dataType* segmentationPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* edgeDetectorPtr = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gaussSeidelPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* previousSolPtr = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	if (segmentationPtr == NULL || edgeDetectorPtr == NULL || 
		gaussSeidelPtr == NULL || previousSolPtr == NULL)
	{
		return false;
	}

	neighPtrs uCoef;
	uCoef.North = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.South = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.East = (dataType*)malloc(sizeof(dataType) * dim2D);
	uCoef.West = (dataType*)malloc(sizeof(dataType) * dim2D);
	if (uCoef.North == NULL || uCoef.East == NULL || 
		uCoef.West == NULL || uCoef.South == NULL)
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
	if (theta_in.East == NULL || theta_in.West == NULL ||
		theta_in.North == NULL || theta_in.South == NULL)
	{
		return false;
	}

	dataType* n_out = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* n_out_qp = (dataType*)malloc(sizeof(dataType) * dim2D);

	//Initialization
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
		n_out[i] = 0.0;
	}

	dataType current = 0.0, average_norm_gradient = 0.0, average_g = 0.0;
	dataType u_average = 0.0;

	heatImplicit2dScheme(imageData, smooth_parms);

	//Compute g function : 1 / (1 + s^2), s = (1 / card(N_p)) * sum(|I_smooth_q|)
	computeNormOfGradientDiamondCells(imageData.imageDataPtr, uGrad, height, width, spacing);
	dataType coef_edge_detector = seg_parms.coef, eps = seg_parms.eps2;
	for (i = 0; i < dim2D; i++)
	{
		average_g = (dataType)((sqrt(uGrad.East[i]) + sqrt(uGrad.West[i]) + sqrt(uGrad.North[i]) + sqrt(uGrad.South[i])) / 4.0);
		edgeDetectorPtr[i] = gradientFunction(average_g * average_g, coef_edge_detector);
	}

	//Array for name construction
	char name[350];
	char name_ending[100];
	Storage_Flags flags = { false,false };

	strcpy_s(name, sizeof name, segmentPath);
	sprintf_s(name_ending, sizeof(name_ending), "_edge_detector.raw");
	strcat_s(name, sizeof(name), name_ending);
	store2dRawData(edgeDetectorPtr, height, width, name, flags);

	//compute gradient of edge detector function
	// v_pq = -w_a * m(e_pq) * G_pq;
	dataType diff = seg_parms.coef_dif, adv = seg_parms.coef_conv;
	dataType vpe, vpw, vpn, vps;
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
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

			//a_in_pq = max(v_pq, 0)
			a_in.East[xd] = fmax(vpe, 0);
			a_in.West[xd] = fmax(vpw, 0);
			a_in.North[xd] = fmax(vpn, 0);
			a_in.South[xd] = fmax(vps, 0);

			//a_out_pq = min(v_pq, 0)
			a_out_pq.East[xd] = fmin(vpe, 0);
			a_out_pq.West[xd] = fmin(vpw, 0);
			a_out_pq.North[xd] = fmin(vpn, 0);
			a_out_pq.South[xd] = fmin(vps, 0);

			//a_out_qp = -a_in_pq
			a_out_qp.East[xd] = -a_in.East[xd];
			a_out_qp.West[xd] = -a_in.West[xd];
			a_out_qp.North[xd] = -a_in.North[xd];
			a_out_qp.South[xd] = -a_in.South[xd];

			//n_out_pq = -sum(sign(a_out_pq))
			n_out[xd] = 0.0;
			if (a_out_pq.East[xd] != 0.0)
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.West[xd] != 0.0)
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.North[xd] != 0.0)
			{
				n_out[xd] += 1.0;
			}
			if (a_out_pq.South[xd] != 0.0)
			{
				n_out[xd] += 1.0;
			}

			//n_out_qp = -sum(sign(a_out_qp))
			n_out_qp[xd] = 0.0;
			if (a_out_qp.East[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.West[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.North[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
			if (a_out_qp.South[xd] != 0.0)
			{
				n_out_qp[xd] += 1.0;
			}
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
	dataType mp = spacing.sx * spacing.sy;
	dataType u1 = 0.0, u2 = 0.0, u3 = 0.0, u4 = 0.0;
	dataType u_p_min = 0.0, u_p_max = 0.0;
	dataType value = 0.0, prod_pq = 0.0, prod_qp = 0.0;
	dataType numerator_max = 0.0, numerator_min = 0.0;
	dataType tol = seg_parms.segTolerance, omega = seg_parms.omega_c;
	dataType gauss_seidel_coef = 0.0;
	size_t maxIter = seg_parms.maxNoGSIteration;
	dataType tau = seg_parms.tau;
	dataType tau_mp = tau / mp;

	dataType prod_east = 0.0, prod_west = 0.0, prod_north = 0.0, prod_south = 0.0;
	dataType prod_east_p = 0.0, prod_west_p = 0.0, prod_north_p = 0.0, prod_south_p = 0.0;

	dataType theta_east = 0.0, theta_west = 0.0, theta_north = 0.0, theta_south = 0.0;
	dataType theta_east_p = 0.0, theta_west_p = 0.0, theta_north_p = 0.0, theta_south_p = 0.0;

	dataType u_p_prev = 0.0, u_p_cur = 0.0, u_east = 0.0, u_west = 0.0, u_north = 0.0, u_south = 0.0, error_gauss_seidel = 0.0;

	size_t cpt, xd, iplus, iminus, jplus, jminus;
	do {
		number_time_step++;

		//The return norm of gradient is |grad(u)|^2
		computeNormOfGradientDiamondCells(segmentationPtr, normGrad, height, width, spacing);

		//=================================================

		//============ First setp ========================
		for (i = 0, i_ext = 1; i < height; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				xd = x_new(i, j, height);
				iplus = i_ext + 1;
				iminus = i_ext - 1;
				jplus = j_ext + 1;
				jminus = j_ext - 1;

				//Compute average of norm of gradient
				average_norm_gradient = (dataType)((normGrad.East[xd] + normGrad.West[xd] + normGrad.North[xd] + normGrad.South[xd]) / 4.0);
				u_average = sqrt(average_norm_gradient + eps);

				//Epsilon regularization
				normGrad.East[xd] = sqrt(normGrad.East[xd] + eps);
				normGrad.West[xd] = sqrt(normGrad.West[xd] + eps);
				normGrad.North[xd] = sqrt(normGrad.North[xd] + eps);
				normGrad.South[xd] = sqrt(normGrad.South[xd] + eps);

				theta_out.East[xd] = 0.5;
				theta_out.West[xd] = 0.5;
				theta_out.North[xd] = 0.5;
				theta_out.South[xd] = 0.5;

				theta_in.East[xd] = 0.5;
				theta_in.West[xd] = 0.5;
				theta_in.North[xd] = 0.5;
				theta_in.South[xd] = 0.5;

				uCoef.East[xd] = (dataType)(tau_mp * theta_in.East[xd] * a_in.East[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.East[xd]));
				uCoef.West[xd] = (dataType)(tau_mp * theta_in.West[xd] * a_in.West[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.West[xd]));
				uCoef.North[xd] = (dataType)(tau_mp * theta_in.North[xd] * a_in.North[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.North[xd]));
				uCoef.South[xd] = (dataType)(tau_mp * theta_in.South[xd] * a_in.South[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.South[xd]));
			}
		}

		//gauss seidel for segmentation function
		cpt = 0;

		error_gauss_seidel = 0.0;
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

					gauss_seidel_coef = (dataType)(((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]
						+ tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
							+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
							+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
							+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south])
						+ (uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
							+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]))
						/ (1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]));
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

					u1 = (dataType)((1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]) * gaussSeidelPtr[xd_ext]);
					u2 = (dataType)((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]);
					u3 = (dataType)(tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
						+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
						+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
						+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south]));
					u4 = (dataType)(uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//============ Second step ========================

		for (i = 0, i_ext = 1; i < height; i++, i_ext++)
		{
			for (j = 0, j_ext = 1; j < width; j++, j_ext++)
			{
				u_p_cur = gaussSeidelPtr[x_new(i_ext, j_ext, height_ext)];
				u_p_min = getMinInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
				u_p_max = getMaxInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);

				if (u_p_cur < u_p_min || u_p_cur > u_p_max) 
				{
					xd = x_new(i, j, height);
					iplus = i_ext + 1;
					iminus = i_ext - 1;
					jplus = j_ext + 1;
					jminus = j_ext - 1;

					u_p_prev = previousSolPtr[x_new(i_ext, j_ext, height_ext)];
					u_east = previousSolPtr[x_new(iplus, j_ext, height_ext)];
					u_west = previousSolPtr[x_new(iminus, j_ext, height_ext)];
					u_north = previousSolPtr[x_new(i_ext, jminus, height_ext)];
					u_south = previousSolPtr[x_new(i_ext, jplus, height_ext)];

					//Compute theta_out_pq
					u_p_min = getMinInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
					u_p_max = getMaxInNeighborhood(previousSolPtr, height_ext, width_ext, i_ext, j_ext);
					numerator_max = mp * (u_p_max - u_p_prev);
					numerator_min = mp * (u_p_min - u_p_prev);

					if (n_out[xd] == 0)
					{
						theta_east = 0.5;
						theta_west = 0.5;
						theta_north = 0.5;
						theta_south = 0.5;
					}
					else
					{
						//East
						prod_east = a_out_pq.East[xd] * (u_east - u_p_prev);
						if (prod_east == 0)
						{
							theta_east = 0.5;
						}
						else if (prod_east > 0)
						{
							theta_east = numerator_max / (tau * n_out[xd] * prod_east);
						}
						else
						{
							theta_east = numerator_min / (tau * n_out[xd] * prod_east);
						}

						//West
						prod_west = a_out_pq.West[xd] * (u_west - u_p_prev);
						if (prod_west == 0)
						{
							theta_west = 0.5;
						}
						else if (prod_west > 0)
						{
							theta_west = numerator_max / (tau * n_out[xd] * prod_west);
						}
						else
						{
							theta_west = numerator_min / (tau * n_out[xd] * prod_west);
						}

						//North
						prod_north = a_out_pq.North[xd] * (u_north - u_p_prev);
						if (prod_north == 0)
						{
							theta_north = 0.5;
						}
						else if (prod_north > 0)
						{
							theta_north = numerator_max / (tau * n_out[xd] * prod_north);
						}
						else
						{
							theta_north = numerator_min / (tau * n_out[xd] * prod_north);
						}

						//South
						prod_south = a_out_pq.South[xd] * (u_south - u_p_prev);
						if (prod_south == 0)
						{
							theta_south = 0.5;
						}
						else if (prod_south > 0)
						{
							theta_south = numerator_max / (tau * n_out[xd] * prod_south);
						}
						else
						{
							theta_south = numerator_min / (tau * n_out[xd] * prod_south);
						}
					}

					theta_out.East[xd] = fmin(0.5, theta_east);
					theta_out.West[xd] = fmin(0.5, theta_west);
					theta_out.North[xd] = fmin(0.5, theta_north);
					theta_out.South[xd] = fmin(0.5, theta_south);

					//Compute thata_in_pq = 1 - theta_out_pq

					if (n_out_qp[xd] == 0)
					{
						theta_east_p = 0.5;
						theta_west_p = 0.5;
						theta_north_p = 0.5;
						theta_south_p = 0.5;
					}
					else
					{
						//East
						prod_east_p = a_out_qp.East[xd] * (u_east - u_p_prev);
						if (prod_east_p == 0)
						{
							theta_east_p = 0.5;
						}
						else if (prod_east_p > 0)
						{
							theta_east_p = numerator_max / (tau * n_out_qp[xd] * prod_east_p);
						}
						else
						{
							theta_east_p = numerator_min / (tau * n_out_qp[xd] * prod_east_p);
						}

						//West
						prod_west_p = a_out_qp.West[xd] * (u_west - u_p_prev);
						if (prod_west_p == 0)
						{
							theta_west_p = 0.5;
						}
						else if (prod_west_p > 0)
						{
							theta_west_p = numerator_max / (tau * n_out_qp[xd] * prod_west_p);
						}
						else
						{
							theta_west_p = numerator_min / (tau * n_out_qp[xd] * prod_west_p);
						}

						//North
						prod_north_p = a_out_qp.North[xd] * (u_north - u_p_prev);
						if (prod_north_p == 0)
						{
							theta_north_p = 0.5;
						}
						else if (prod_north_p > 0)
						{
							theta_north_p = numerator_max / (tau * n_out_qp[xd] * prod_north_p);
						}
						else
						{
							theta_north_p = numerator_min / (tau * n_out_qp[xd] * prod_north_p);
						}

						//South
						prod_south_p = a_out_qp.South[xd] * (u_south - u_p_prev);
						if (prod_south_p == 0)
						{
							theta_south_p = 0.5;
						}
						else if (prod_south_p > 0)
						{
							theta_south_p = numerator_max / (tau * n_out_qp[xd] * prod_south_p);
						}
						else
						{
							theta_south_p = numerator_min / (tau * n_out_qp[xd] * prod_south_p);
						}
					}

					theta_in.East[xd] = 1.0 - fmin(0.5, theta_east_p);
					theta_in.West[xd] = 1.0 - fmin(0.5, theta_west_p);
					theta_in.North[xd] = 1.0 - fmin(0.5, theta_north_p);
					theta_in.South[xd] = 1.0 - fmin(0.5, theta_south_p);

					uCoef.East[xd] = (dataType)(tau_mp * theta_in.East[xd] * a_in.East[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.East[xd]));
					uCoef.West[xd] = (dataType)(tau_mp * theta_in.West[xd] * a_in.West[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hx2 * normGrad.West[xd]));
					uCoef.North[xd] = (dataType)(tau_mp * theta_in.North[xd] * a_in.North[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.North[xd]));
					uCoef.South[xd] = (dataType)(tau_mp * theta_in.South[xd] * a_in.South[xd] + tau * (diff * u_average * edgeDetectorPtr[xd]) / (hy2 * normGrad.South[xd]));
				}
			}
		}

		//gauss seidel for segmentation function
		cpt = 0;

		error_gauss_seidel = 0.0;
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

					gauss_seidel_coef = (dataType)(((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]
						+ tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
							+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
							+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
							+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south])
						+ (uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
							+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]))
						/ (1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]));
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

					u1 = (dataType)((1.0 + uCoef.East[xd] + uCoef.West[xd] + uCoef.North[xd] + uCoef.South[xd]) * gaussSeidelPtr[xd_ext]);
					u2 = (dataType)((1 - tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] + theta_out.West[xd] * a_out_pq.West[xd]
						+ theta_out.North[xd] * a_out_pq.North[xd] + theta_out.South[xd] * a_out_pq.South[xd])) * previousSolPtr[xd_ext]);
					u3 = (dataType)(tau_mp * (theta_out.East[xd] * a_out_pq.East[xd] * previousSolPtr[ind_east]
						+ theta_out.West[xd] * a_out_pq.West[xd] * previousSolPtr[ind_west]
						+ theta_out.North[xd] * a_out_pq.North[xd] * previousSolPtr[ind_north]
						+ theta_out.South[xd] * a_out_pq.South[xd] * previousSolPtr[ind_south]));
					u4 = (dataType)(uCoef.East[xd] * gaussSeidelPtr[ind_east] + uCoef.West[xd] * gaussSeidelPtr[ind_west]
						+ uCoef.North[xd] * gaussSeidelPtr[ind_north] + uCoef.South[xd] * gaussSeidelPtr[ind_south]);
					error_gauss_seidel += pow(u1 - u2 - u3 - u4, 2);
				}
			}

		} while (cpt < maxIter && error_gauss_seidel > 0.001);

		//=================================================
		
		//rescall to data range 0-1
		rescaleToZeroOne2d(gaussSeidelPtr, height_ext, width_ext);

		//compute L2-norm
		error_segmentation = l2norm(gaussSeidelPtr, previousSolPtr, height_ext, width_ext, 1.0);

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

	free(n_out);
	free(n_out_qp);

	return true;
}

