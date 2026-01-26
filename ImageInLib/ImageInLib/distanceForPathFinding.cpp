#include <iostream>
#include <sstream>  
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include <cmath>
#include <omp.h>
#include <vector>
#include <tuple>
#include "distanceForPathFinding.h"
#include <template_functions.h>
#include "filtering.h"
#include "../src/heat_equation.h"
#include "hough_transform.h"
#include "../src/distance_function.h"
#include "percentile.h"

#define BIG_VALUE INFINITY

//J.A Sethian, A Fast Marching Level Set method for Monotonically advancing fronts, 1995, page 8 and 10.
//link to article ---> http://ugweb.cs.ualberta.ca/~vis/courses/CompVis/readings/modelrec/sethian95fastlev.pdf
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//=================== Functions for the 2D Fast Marching and Path Tracking ==============================
/////////////////////////////////////////////////////////////////////////////////////////////////////////

dataType min0(dataType x, dataType y) {
	if (y - x > 0)
		return pow(x - y, 2);
	else
		return 0;
}

/*
dataType upwindFiniteDifference2dX(dataType* actionPtr, const size_t length, const size_t width, size_t i, size_t j) {
	
	dataType i_minus, i_plus;

	if (i == 0) {
		i_minus = INFINITY;
	}
	else {
		i_minus = actionPtr[x_new(i - 1, j, length)];
	}

	if (i >= length - 1) {
		i_plus = INFINITY;
	}
	else {
		i_plus = actionPtr[x_new(i + 1, j, length)];
	}

	return min(i_minus, i_plus);
}

dataType upwindFiniteDifference2dY(dataType* actionPtr, const size_t length, const size_t width, size_t i, size_t j) {
	
	dataType j_minus, j_plus;

	if (j == 0) {
		j_minus = INFINITY;
	}
	else {
		j_minus = actionPtr[x_new(i, j - 1, length)];
	}

	if (j >= width - 1) {
		j_plus = INFINITY;
	}
	else {
		j_plus = actionPtr[x_new(i, j + 1, length)];
	}

	return min(j_minus, j_plus);
}

dataType solve2dQuadratic(dataType X, dataType Y, dataType P, PixelSpacing h, size_t indx, size_t indy, FILE* pFile) {

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P_2 = P * P;

	dataType hx = h.sx;
	dataType hx_2 = h.sx * h.sx;

	dataType hy = h.sy;
	dataType hy_2 = h.sy * h.sy;

	if (P <= 0.0) {
		std::cerr << "Error: P must be positive." << std::endl;
		return solution; // Return 0 if P is not positive
	}

	if (X == INFINITY && Y != INFINITY)
	{
		a = 1.0;
		b = -2 * Y;
		c = (dataType)(Y * Y - hy_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > Y) {
				return solution;
			}
			else {
				return (dataType)(Y + hy * P);
			}
		}
		else {
			return (dataType)(Y + hy * P);
		}
	}

	if (X != INFINITY && Y == INFINITY)
	{
		a = 1.0;
		b = -2 * X;
		c = (dataType)(X * X - hx_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > X)
			{
				return solution;
			}
			else {
				return (dataType)(X + hx * P);
			}
		}
		else {
			return (dataType)(X + hx * P);
		}
	}

	if (X != INFINITY && Y != INFINITY)
	{
		a = hx_2 + hy_2;
		b = -2 * (hy_2 * X + hx_2 * Y);
		c = (dataType)(hy_2 * X * X + hx_2 * Y * Y - hx_2 * hy_2 * P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > max(X, Y)) {
				return solution;
			}
			else {
				return (dataType)(min(X + hx * P, Y + hy * P));
			}
		}
		else {
			fprintf(pFile, "%d,%d\n", indx, indy);
			return (dataType)(min(X + hx * P, Y + hy * P));
		}
	}

}

bool computePotential(Image_Data2D imageDataStr, dataType* potentialFuncPtr, Point2D * seedPoints, Potential_Parameters parameters)
{

	if (imageDataStr.imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	const size_t length = imageDataStr.height;
	const size_t width = imageDataStr.width;
	size_t i, j, dim2D = length * width;

	dataType seedVal1 = (dataType)(imageDataStr.imageDataPtr[x_new((size_t)seedPoints[0].x, (size_t)seedPoints[0].y, length)]);
	//dataType seedVal2 = (dataType)(imageDataStr.imageDataPtr[x_new((size_t)seedPoints[1].x, (size_t)seedPoints[1].y, length)]);
	//dataType seedVal = (dataType)(seedVal1 + seedVal2) / 2.0; //Average of the two seed points

	//dataType* distanceMap = new dataType[dim2D]{ 0 };
	//dataType* edgeImage = new dataType[dim2D]{ 0 };	

	dataType norm_of_gradient = 0.0, edgeValue = 0.0;
	bool isGradientComputed = false;
	Point2D grad_vector;
	PixelSpacing fVolume = imageDataStr.spacing;

	//dataType* arraySorted = new dataType[dim2D];
	//copyDataToAnother2dArray(imageDataStr.imageDataPtr, arraySorted, length, width);

	//for (i = 0; i < length; i++) 
	//{
	//	for (j = 0; j < width; j++) 
	//	{
	//		isGradientComputed = getGradient2D(imageDataStr.imageDataPtr, length, width, i, j, fVolume, &grad_vector);
	//		if (isGradientComputed == true) {
	//			norm_of_gradient = sqrt(grad_vector.x * grad_vector.x + grad_vector.y * grad_vector.y);
	//			edgeValue = gradientFunction(norm_of_gradient, parameters.K);
	//		}
	//		else {
	//			std::cout << "Error in computing gradient at point (" << i << ", " << j << ")" << std::endl;
	//			return false;
	//		}
	//		//Threshold
	//		if( edgeValue <= parameters.thres) 
	//		{
	//			edgeImage[x_new(i, j, length)] = 0.0;
	//		}
	//		else 
	//		{
	//			edgeImage[x_new(i, j, length)] = 1.0;
	//		}
	//	}
	//}

	//string path_file = "C:/Users/Konan Allaly/Documents/Tests/output/edge_image.raw";
	//manageRAWFile2D<dataType>(edgeImage, length, width, path_file.c_str(), STORE_DATA, false);

	//////fastSweepingFunction_2D(distanceMap, edgeImage, length, width, 1.0, 1000000.0, 1.0);
	//Image_Data2D todistanceMapStr = { length, width, imageDataStr.imageDataPtr, imageDataStr.origin, imageDataStr.spacing, imageDataStr.orientation };
	////fastMarchingForDistanceMap(todistanceMapStr, distanceMap, 0.0);
	//bruteForceDistanceMap2D(todistanceMapStr, distanceMap, 1.0);
	//std::string path_file = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_2d.raw";
	//manageRAWFile2D<dataType>(distanceMap, length, width, path_file.c_str(), STORE_DATA, false);

	for (i = 0; i < dim2D; i++) {
		potentialFuncPtr[i] = fabs(imageDataStr.imageDataPtr[i] - seedVal1);
		//potentialFuncPtr[i] = fabs(imageDataStr.imageDataPtr[i] - seedVal1);
	}

	////Find max difference
	//dataType maxDiff = 0.0;
	//for (i = 0; i < dim2D; i++) {
	//	if (potentialFuncPtr[i] > maxDiff) {
	//		maxDiff = potentialFuncPtr[i];
	//	}
	//}
	//std::cout << "MaxDiff = " << maxDiff << std::endl;

	//quickSort(arraySorted, 0, dim2D - 1);
	//size_t index = (size_t)round((95 / 100.0) * (dim2D - 1));
	//dataType factor = arraySorted[index];
	//double percentile = 95.0;
	//dataType factor = 0.596723;//computePercentile(arraySorted, length, width, percentile);
	//std::cout << "Percentile 95 = " << factor << std::endl;

	//Normalization
	dataType weight = 0.0;
	for (i = 0; i < dim2D; i++) {
		potentialFuncPtr[i] = parameters.eps + potentialFuncPtr[i];// / maxDiff;
		//if(factor != 0.0)
		//{
		//	potentialFuncPtr[i] = parameters.eps + potentialFuncPtr[i] / factor;
		//}
		//else {
		//	potentialFuncPtr[i] = parameters.eps + potentialFuncPtr[i];
		//}
	}

	//delete[] distanceMap;
	//delete[] edgeImage;
	//delete[] arraySorted;

	return true;
}

void heapifyDown2D(vector<pointFastMarching2D>& in_Process, int pos) {
	
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int current = pos;
	int left_child = 2 * pos + 1;
	int right_child = 2 * pos + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}

	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left <= val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right <= val_current) {
			current = right_child;
		}
	}

	if (current != pos) {
		swap_elts(&in_Process[pos], &in_Process[current], sizeof(pointFastMarching2D));
		heapifyDown2D(in_Process, current);
	}
}

void heapifyUp2D(vector<pointFastMarching2D>& in_Process, int i) {

	int current = i;

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current <= val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		swap_elts(&in_Process[current], &in_Process[i], sizeof(pointFastMarching2D));
		heapifyUp2D(in_Process, current);
	}

}

void heapifyVector2D(vector<pointFastMarching2D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int indx, start = length_array / 2 - 1;
	for (indx = start; indx >= 0; indx--) {
		heapifyDown2D(in_Process, indx);
	}
}

void deleteRootHeap2D(vector<pointFastMarching2D>& in_Process) {
	int l = in_Process.size();
	if(l > 1)
	{
		swap_elts(&in_Process[0], &in_Process[l - 1], sizeof(pointFastMarching2D));
		in_Process.pop_back();
		heapifyDown2D(in_Process, 0);
	}
	else if(l == 1) {
		in_Process.pop_back();
	}
}

void addPointHeap2D(vector<pointFastMarching2D>& in_Process, pointFastMarching2D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp2D(in_Process, l - 1);
}

int getIndexFromHeap2D(vector<pointFastMarching2D>& in_Process, size_t i, size_t j) {
	for (int ind = 0; ind < in_Process.size(); ind++) {
		if (in_Process[ind].x == i && in_Process[ind].y == j) {
			return ind;
		}
	}
	return -1; //not found
}

void updateNeighbor2D(size_t ind_x, size_t ind_y, size_t length, size_t width,
	dataType* action, dataType* potential, short* labelArray,
	PixelSpacing spacing, vector<pointFastMarching2D>& narrowBand, FILE* pFile)
{

	if (ind_x >= length || ind_y >= width)
		return;

	size_t neighborIndx = x_new(ind_x, ind_y, length);
	dataType ux = upwindFiniteDifference2dX(action, length, width, ind_x, ind_y);
	dataType uy = upwindFiniteDifference2dY(action, length, width, ind_x, ind_y);
	dataType coefSpeed = potential[neighborIndx];
	dataType solution = solve2dQuadratic(ux, uy, coefSpeed, spacing, ind_x, ind_y, pFile);
	pointFastMarching2D neighbor = { ind_x, ind_y, solution };
	if (labelArray[neighborIndx] == 3) {
		addPointHeap2D(narrowBand, neighbor);
		action[neighborIndx] = solution;
		labelArray[neighborIndx] = 2;
	}
	else if (labelArray[neighborIndx] == 2 && solution < action[neighborIndx]) {
		action[neighborIndx] = solution;
		int pIndex = getIndexFromHeap2D(narrowBand, ind_x, ind_y);
		if (pIndex != -1) {
			heapifyUp2D(narrowBand, pIndex);
		}
	}

}

bool fastMarching2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* seedPoints) {

	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;
	size_t dim2D = length * width;

	vector<pointFastMarching2D> inProcess;
	short* labelArray = new short[dim2D] { 0 };
	if (imageData.imageDataPtr == NULL || distancePtr == NULL || potentialPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}
	
	//1 ---> already processed, 
	//2 ---> in process and 
	//3 ---> not processed
	for (size_t k = 0; k < dim2D; k++) {
		distancePtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	std::string path_discriminant = "C:/Users/Konan Allaly/Documents/Tests/output/negative_discrinant_01.csv";
	FILE* dFile;
	if (fopen_s(&dFile, path_discriminant.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(dFile, "x,y\n");

	size_t i = (size_t)seedPoints[0].x;
	size_t j = (size_t)seedPoints[0].y;
	size_t currentIndx = x_new(i, j, length);
	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	if (i > 0)
	{
		if( labelArray[x_new(i - 1, j, length)] != 1 ) {
			updateNeighbor2D(i - 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (i < length_minus)
	{
		if (labelArray[x_new(i + 1, j, length)] != 1) {
			updateNeighbor2D(i + 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (j > 0)
	{
		if (labelArray[x_new(i, j - 1, length)] != 1) {
			updateNeighbor2D(i, j - 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (j < width_minus)
	{
		if (labelArray[x_new(i, j + 1, length)] != 1) {
			updateNeighbor2D(i, j + 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	
	size_t nb_pt_processed = 0;
	while (inProcess.size() > 0) {

		pointFastMarching2D current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;
		deleteRootHeap2D(inProcess);
		nb_pt_processed++;

		if (i > 0)
		{
			if (labelArray[x_new(i - 1, j, length)] != 1) {
				updateNeighbor2D(i - 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (i < length_minus)
		{
			if (labelArray[x_new(i + 1, j, length)] != 1) {
				updateNeighbor2D(i + 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (j > 0)
		{
			if (labelArray[x_new(i, j - 1, length)] != 1) {
				updateNeighbor2D(i, j - 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (j < width_minus)
		{
			if (labelArray[x_new(i, j + 1, length)] != 1) {
				updateNeighbor2D(i, j + 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}

	}

	fclose(dFile);

	delete[] labelArray;
}

bool partialFrontPropagation2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* endPoints) {

	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	short* labelArray = new short[dim2D] {0};

	if (imageData.imageDataPtr == NULL || distancePtr == NULL || potentialPtr == NULL || endPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0;
	vector<pointFastMarching2D> inProcess;

	i = (size_t)endPoints[0].x;
	j = (size_t)endPoints[0].y;
	size_t currentIndx = x_new(i, j, length);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distancePtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	std::string path_discriminant = "C:/Users/Konan Allaly/Documents/Tests/output/negative_discrinant.csv";
	FILE* dFile;
	if (fopen_s(&dFile, path_discriminant.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(dFile, "x,y\n");

	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	if(i > 0)
	{
		if(labelArray[x_new(i - 1, j, length)] != 1) {
			updateNeighbor2D(i - 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (i < length_minus)
	{
		if(labelArray[x_new(i + 1, j, length)] != 1) {
			updateNeighbor2D(i + 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (j > 0)
	{
		if(labelArray[x_new(i, j - 1, length)] != 1) {
			updateNeighbor2D(i, j - 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}
	if (j < width_minus)
	{
		if(labelArray[x_new(i, j + 1, length)] != 1) {
			updateNeighbor2D(i, j + 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
		}
	}

	size_t seedI = endPoints[1].x, seedJ = endPoints[1].y, seedIndex = x_new(seedI, seedJ, length);

	dataType max_save_action = 0.0;
	while (labelArray[seedIndex] != 1) {

		pointFastMarching2D current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		if(distancePtr[currentIndx] > max_save_action) {
			max_save_action = distancePtr[currentIndx];
		}
		
		deleteRootHeap2D(inProcess);

		if (i > 0)
		{
			if( labelArray[x_new(i - 1, j, length)] != 1) {
				updateNeighbor2D(i - 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (i < length_minus)
		{
			if (labelArray[x_new(i + 1, j, length)] != 1) {
				updateNeighbor2D(i + 1, j, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (j > 0)
		{
			if (labelArray[x_new(i, j - 1, length)] != 1) {
				updateNeighbor2D(i, j - 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
		if (j < width_minus)
		{
			if (labelArray[x_new(i, j + 1, length)] != 1) {
				updateNeighbor2D(i, j + 1, length, width, distancePtr, potentialPtr, labelArray, spacing, inProcess, dFile);
			}
		}
	}

	fclose(dFile);

	//Set the distance of the end point to the maximum value
	for(i = 0; i < dim2D; i++) {
		if (labelArray[i] == 3) {
			distancePtr[i] = max_save_action + 0.1;
		}
	}
	inProcess.clear();
	
	delete[] labelArray;
}
*/

/*
bool doubleFrontPropagation2D(Image_Data2D imageData, dataType* actionFirstFront, dataType* actionSecondFront, dataType* potentialPtr, Point2D* endPoints, string savingPath) {

	if (imageData.imageDataPtr == NULL || actionFirstFront == NULL || actionSecondFront == NULL || potentialPtr == NULL || endPoints == NULL) {
		return false;
	}
	const size_t length = imageData.height;
	const size_t width = imageData.width;
	PixelSpacing spacing = imageData.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	short* firstLabelArray = new short[dim2D];
	short* secondLabelArray = new short[dim2D];
	if(firstLabelArray == NULL || secondLabelArray == NULL || actionFirstFront == NULL || actionSecondFront == NULL) {
		return false;
	}

	//Initialize action
	for (size_t k = 0; k < dim2D; k++) {
		actionFirstFront[k] = INFINITY;
		actionSecondFront[k] = INFINITY;
		firstLabelArray[k] = 3; //3 ---> not processed
		secondLabelArray[k] = 3; //3 ---> not processed
	}

	size_t x1 = (size_t)endPoints[0].x;
	size_t y1 = (size_t)endPoints[0].y;
	if( x1 >= length || y1 >= width) {
		delete[] firstLabelArray;
		delete[] secondLabelArray;
		return false; //Invalid end point
	}
	firstLabelArray[x_new(x1, y1, length)] = 1; //1 ---> already processed
	actionFirstFront[x_new(x1, y1, length)] = 0.0;

	size_t x2 = (size_t)endPoints[1].x;
	size_t y2 = (size_t)endPoints[1].y;
	if( x2 >= length || y2 >= width) {
		delete[] secondLabelArray;
		delete[] firstLabelArray;
		return false; //Invalid end point
	}
	secondLabelArray[x_new(x2, y2, length)] = 1; //1 ---> already processed
	actionSecondFront[x_new(x2, y2, length)] = 0.0;
	
	vector<pointFastMarching2D> narrowBandFirstFront, narrowBandSecondFront;

	//Initialize neighbors for the first front

	if (x1 > 0)
	{
		updateNeighbor2D(x1 - 1, y1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (x1 < length_minus)
	{
		updateNeighbor2D(x1 + 1, y1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (y1 > 0)
	{
		updateNeighbor2D(x1, y1 - 1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (y1 < width_minus)
	{
		updateNeighbor2D(x1, y1 + 1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}

	//Initialize neighbors for the second front

	if (x2 > 0)
	{
		updateNeighbor2D(x2 - 1, y2, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (x2 < length_minus)
	{
		updateNeighbor2D(x2 + 1, y2, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (y2 > 0)
	{
		updateNeighbor2D(x2, y2 - 1, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (y2 < width_minus)
	{
		updateNeighbor2D(x2, y2 + 1, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}

	dataType max_save_action = 0.0;
	size_t nb_computed_points = 0;

	while (narrowBandFirstFront.size() > 0 && narrowBandSecondFront.size() > 0) {

		pointFastMarching2D firstFrontPoint = narrowBandFirstFront[0];
		x1 = firstFrontPoint.x;
		y1 = firstFrontPoint.y;
		size_t indexFirst = x_new(x1, y1, length);
		if (secondLabelArray[indexFirst] == 1) {
			//The fronts have met
			endPoints[2].x = x1;
			endPoints[2].y = y1;
			break;
		}
		else {
			firstLabelArray[indexFirst] = 1;
		}

		if(firstFrontPoint.arrival > max_save_action) {
			max_save_action = firstFrontPoint.arrival;
		}	
		deleteRootHeap2D(narrowBandFirstFront);
		nb_computed_points++;

		pointFastMarching2D secondFrontPoint = narrowBandSecondFront[0];
		x2 = secondFrontPoint.x;
		y2 = secondFrontPoint.y;
		size_t indexSecond = x_new(x2, y2, length);
		if(firstLabelArray[indexSecond] == 1) {
			//The fronts have met
			endPoints[2].x = x2;
			endPoints[2].y = y2;
			break;
		}
		else {
			secondLabelArray[indexSecond] = 1;
		}

		if (secondFrontPoint.arrival > max_save_action) {
			max_save_action = secondFrontPoint.arrival;
		}
		deleteRootHeap2D(narrowBandSecondFront);

		//update neighbors for the first front

		if (x1 > 0)
		{
			updateNeighbor2D(x1 - 1, y1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (x1 < length_minus)
		{
			updateNeighbor2D(x1 + 1, y1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (y1 > 0)
		{
			updateNeighbor2D(x1, y1 - 1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (y1 < width_minus)
		{
			updateNeighbor2D(x1, y1 + 1, length, width, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}

		//update  neighbors for the second front

		if (x2 > 0)
		{
			updateNeighbor2D(x2 - 1, y2, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (x2 < length_minus)
		{
			updateNeighbor2D(x2 + 1, y2, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (y2 > 0)
		{
			updateNeighbor2D(x2, y2 - 1, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (y2 < width_minus)
		{
			updateNeighbor2D(x2, y2 + 1, length, width, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}

	}

	for (size_t k = 0; k < dim2D; k++) {
		if(actionFirstFront[k] == INFINITY) {
			actionFirstFront[k] = max_save_action + 1;
		}
		if(actionSecondFront[k] == INFINITY) {
			actionSecondFront[k] = max_save_action + 1;
		}
	}
	narrowBandFirstFront.clear();
	narrowBandSecondFront.clear();

	delete[] firstLabelArray;
	delete[] secondLabelArray;
	
	return true;
}

bool frontPropagationWithKeyPointDetection2D(Image_Data2D actionMapStr, dataType* potentialFuncPtr, Point2D* seedPoint, const double LengthKeyPoints, vector<Point2D>& key_points, std::string path_saving) {

	if (actionMapStr.imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoint == NULL) {
		return false;
	}

	const size_t length = actionMapStr.height;
	const size_t width = actionMapStr.width;

	PixelSpacing spacing = actionMapStr.spacing;
	if (seedPoint[0].x < 0 || seedPoint[0].x > length || seedPoint[0].y < 0 || seedPoint[0].y > width) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	if (seedPoint[1].x < 0 || seedPoint[1].x > length || seedPoint[1].y < 0 || seedPoint[1].y > width) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}

	vector <pointFastMarching2D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short* labelArray = new short[dim2D];

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < dim2D; k++) {
		actionMapStr.imageDataPtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	pointFastMarching2D current;
	i = (size_t)seedPoint[0].x;
	j = (size_t)seedPoint[0].y;
	size_t currentIndx = x_new(i, j, length);
	actionMapStr.imageDataPtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProcess
	size_t length_minus = length - 1, width_minus = width - 1;

	if (i > 0)
	{
		updateNeighbor2D(i - 1, j, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
	}
	if (i < length_minus)
	{
		updateNeighbor2D(i + 1, j, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
	}
	if (j > 0)
	{
		updateNeighbor2D(i, j - 1, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
	}
	if (j < width_minus)
	{
		updateNeighbor2D(i, j + 1, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
	}

	size_t iEnd = (size_t)seedPoint[1].x;
	size_t jEnd = (size_t)seedPoint[1].y;

	//Visualize the front propagation
	size_t id_keyPoint = 1;

	double distanceToCurrentSourcePoint = 0.0;
	//Set the starting point as initial source point
	Point2D currentSourcePoint = { (dataType)i, (dataType)j};
	key_points.push_back(currentSourcePoint);

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		//Exit of the while loop if the end point is reached
		if (labelArray[x_new(iEnd, jEnd, length)] == 1) {
			break;
		}

		//check if the current point is a key point
		Point2D pSource = { i, j};
		Point2D pSourceReal = getRealCoordFromImageCoord2D(pSource, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		Point2D pCurrent = getRealCoordFromImageCoord2D(currentSourcePoint, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		distanceToCurrentSourcePoint = getPoint2DDistance(pCurrent, pSourceReal);
		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			//If the condition is true ---> new key point is found so we need to initilize it neighbors
			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionMapStr.imageDataPtr[currentIndx] = 0;

			//Initialize all the points inside the narrow band as processed ---> label = 1
			for (size_t it = 0; it < inProcess.size(); it++) {
				labelArray[x_new(inProcess[it].x, inProcess[it].y, length)] = 1;
			}

			//West
			if (i > 0 && i < length && j >= 0 && j < width)
			{
				size_t iminus = i - 1;
				size_t indxWest = x_new(iminus, j, length);
				actionMapStr.imageDataPtr[indxWest] = INFINITY;
				labelArray[indxWest] = 3;
			}

			//East
			if (i < length_minus && i >= 0 && j >= 0 && j < width)
			{
				size_t iplus = i + 1;
				size_t indxEast = x_new(iplus, j, length);
				actionMapStr.imageDataPtr[indxEast] = INFINITY;
				labelArray[indxEast] = 3;
			}

			//North
			if (j > 0 && j < width && i >= 0 && i < length)
			{
				size_t jminus = j - 1;
				size_t indxNorth = x_new(i, jminus, length);
				actionMapStr.imageDataPtr[indxNorth] = INFINITY;
				labelArray[indxNorth] = 3;
			}

			//South
			if (j >= 0 && j < width_minus && i >= 0 && i < length)
			{
				size_t jplus = j + 1;
				size_t indxSouth = x_new(i, jplus, length);
				actionMapStr.imageDataPtr[indxSouth] = INFINITY;
				labelArray[indxSouth] = 3;
			}

			inProcess.clear();
			id_keyPoint++;
		}
		else {
			//actionMapStr.imageDataPtr[currentIndx] = current.arrival;
			deleteRootHeap2D(inProcess);
		}

		//====================
		
		//update neighbors of the minimum in the narrow band

		if (i > 0)
		{
			updateNeighbor2D(i - 1, j, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
		}
		if (i < length_minus)
		{
			updateNeighbor2D(i + 1, j, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
		}
		if (j > 0)
		{
			updateNeighbor2D(i, j - 1, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
		}
		if (j < width_minus)
		{
			updateNeighbor2D(i, j + 1, length, width, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, inProcess);
		}

	}

	key_points.push_back(seedPoint[1]);

	delete[] labelArray;

	return true;
}
*/

/*
bool shortestPath2d(Image_Data2D distanceFuncPtr, Point2D* seedPoints, vector<Point2D>& path_points, Path_Parameters parameters) {

	if (distanceFuncPtr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	const size_t height = distanceFuncPtr.height;
	const size_t width = distanceFuncPtr.width;
	dataType dist_min = 0.0;

	size_t i = (size_t)seedPoints[1].x;
	size_t j = (size_t)seedPoints[1].y;
	dataType x = seedPoints[1].x;
	dataType y = seedPoints[1].y;

	Point2D grad;
	dataType tau = parameters.tau;
	size_t count_iter = 0;

	do{

		getGradient2D(distanceFuncPtr.imageDataPtr, height, width, i, j, distanceFuncPtr.spacing, &grad);
		dataType gradNorm = sqrt(grad.x * grad.x + grad.y * grad.y);
		x -= tau * grad.x / gradNorm;
		y -= tau * grad.y / gradNorm;
		
		Point2D current = { x, y };
		dist_min = getPoint2DDistance(current, seedPoints[0]);
		path_points.push_back(current);

		i = (size_t)round(x); 
		j = (size_t)round(y);

		count_iter++;
	}
	while(dist_min > parameters.tolerance && count_iter < parameters.max_iteration);

	return true;
}
*/

/*
bool bruteForceDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distancePtr == NULL) {
		return false;
	}
	const size_t length = ctImageData.height;
	const size_t width = ctImageData.width;
	PixelSpacing spacing = ctImageData.spacing;

	double min_distance = 0.0;
	for (size_t i = 0; i < length; i++) {
		for (size_t j = 0; j < width; j++) {
			size_t xd = x_new(i, j, length);
			Point2D cPoint = { i, j};
			cPoint = getRealCoordFromImageCoord2D(cPoint, ctImageData.origin, spacing, ctImageData.orientation);

			min_distance = INFINITY;

			for (size_t ti = 0; ti < length; ti++) {
				for (size_t tj = 0; tj < width; tj++) {
					size_t txd = x_new(ti, tj, length);
					if (ctImageData.imageDataPtr[txd] == foregroundValue) {
						Point2D tPoint = { ti, tj};
						tPoint = getRealCoordFromImageCoord2D(tPoint, ctImageData.origin, spacing, ctImageData.orientation);
						double pDistance = getPoint2DDistance(cPoint, tPoint);
						if (pDistance < min_distance) {
							min_distance = pDistance;
						}
					}
				}
			}
			if (min_distance == INFINITY) {
				distancePtr[xd] = 0.0;
			}
			else {
				distancePtr[xd] = min_distance;
			}
		}
	}
	return true;
}

bool fastMarchingForDistanceMap(Image_Data2D ctImageData, dataType* distanceFuncPtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distanceFuncPtr == NULL) {
		return false;
	}

	const size_t length = ctImageData.height;
	const size_t width = ctImageData.width;
	PixelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching2D> inProcess;
	size_t i, j, xd;
	size_t dim2D = length * width;

	//Define and array to follow up the label of each point
	short* labelArray = new short[dim2D] { 0 };
	dataType* potentialPtr = new dataType[dim2D] { 0 };

	//Initialization : Proceed the sources points
	//All the points are notProcessed ---> label = 3
	//Proceed points ---> label = 1
	//Points in narrow band ---> label = 2

	for(i = 0; i < dim2D; i++) {
		if(ctImageData.imageDataPtr[i] == foregroundValue) {
			distanceFuncPtr[i] = 0.0;
			labelArray[i] = 1;
		}
		else {
			distanceFuncPtr[i] = INFINITY;
			labelArray[i] = 3;
		}
		potentialPtr[i] = 1.0; //for distance map the potential is set equal to 1.0
	}

	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	dataType x = 0, y = 0, z = 0, coefSpeed = 0;
	dataType dWest = 0, dEast = 0, dNorth = 0, dSouth = 0;

	//Initialize the narrow band
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {

			xd = x_new(i, j, length);

			if (ctImageData.imageDataPtr[xd] == foregroundValue) {

				if (i > 0)
				{
					updateNeighbor2D(i - 1, j, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
				}
				if (i < length_minus)
				{
					updateNeighbor2D(i + 1, j, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
				}
				if (j > 0)
				{
					updateNeighbor2D(i, j - 1, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
				}
				if (j < width_minus)
				{
					updateNeighbor2D(i, j + 1, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
				}

			}
		}
	}

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching2D current = inProcess[0];
		i = current.x;
		j = current.y;
		xd = x_new(i, j, length);
		if (i >= 0 && i < length && j >= 0 && j < width) {
			labelArray[xd] = 1;
		}
		else {
			return false; //out of bounds check
		}

		deleteRootHeap2D(inProcess);
		if (i > 0)
		{
			updateNeighbor2D(i - 1, j, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
		}
		if (i < length_minus)
		{
			updateNeighbor2D(i + 1, j, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
		}
		if (j > 0)
		{
			updateNeighbor2D(i, j - 1, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
		}
		if (j < width_minus)
		{
			updateNeighbor2D(i, j + 1, length, width, distanceFuncPtr, potentialPtr, labelArray, spacing, inProcess);
		}

	}

	delete[] labelArray;

	return true;
}

bool fastSweepingDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType foregroundValue)
{

	const size_t length = ctImageData.height;
	const size_t width = ctImageData.width;
	PixelSpacing spacing = ctImageData.spacing;
	const size_t dim2D = length * width;
	
	dataType coefSpeed = 1.0;
	dataType hx = ctImageData.spacing.sx;
	dataType hy = ctImageData.spacing.sy;
	dataType hx_2 = pow(hx, 2);
	dataType hy_2 = pow(hy, 2);

	dataType x = 0.0, y = 0.0, new_val = 0.0;
	dataType a = 0.0, b = 0.0, c = 0.0, delta = 0.0;

	size_t length_minus = length - 1;
	size_t width_minus = width - 1;
	
	//Initialization
	for (size_t ij = 0; ij < dim2D; ij++)
	{
		if (ctImageData.imageDataPtr[ij] == foregroundValue)
		{
			distancePtr[ij] = 0.0;
		}
		else
		{
			distancePtr[ij] = INFINITY;
		}
	}
	
	size_t xd = 0;
	dataType grad_x = 0.0, grad_y = 0.0;

	//sweep 1
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < width; j++)
		{
			xd = x_new(i, j, length);
			if (i == 0) {
				x = min(distancePtr[xd], distancePtr[x_new(i + 1, j, length)]);
			}
			else {
				if (i == length_minus) {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[xd]);
				}
				else {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[x_new(i + 1, j, length)]);
				}
			}
			if (j == 0) {
				y = min(distancePtr[xd], distancePtr[x_new(i, j + 1, length)]);
			}
			else {
				if (j == width_minus) {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[xd]);
				}
				else {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[x_new(i, j + 1, length)]);
				}
			}
			a = hx_2 + hy_2;
			b = -2 * (hy_2 * x + hx_2 * y);
			c = hy_2 * pow(x, 2) + hx_2 * pow(y, 2) - hx_2 * hy_2 * coefSpeed;
			delta = pow(b, 2) - 4 * a * c;
			grad_x = x / hx;
			grad_y = y / hy;
			if (fabs(grad_x - grad_y) >= coefSpeed) {
				new_val = min(x + hx * coefSpeed, y + hy * coefSpeed);
			}
			else {
				new_val = ( -b + sqrt(delta) ) / (2 * a);
			}
			if (new_val < distancePtr[xd])
			{
				distancePtr[xd] = new_val;
			}
		}
	}

	//sweep 2
	for (int i = length_minus; i > -1; i--)
	{
		for (int j = 0; j < width; j++)
		{
			xd = x_new(i, j, length);
			if (i == 0) {
				x = min(distancePtr[xd], distancePtr[x_new(i + 1, j, length)]);
			}
			else {
				if (i == length_minus) {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[xd]);
				}
				else {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[x_new(i + 1, j, length)]);
				}
			}
			if (j == 0) {
				y = min(distancePtr[xd], distancePtr[x_new(i, j + 1, length)]);
			}
			else {
				if (j == width_minus) {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[xd]);
				}
				else {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[x_new(i, j + 1, length)]);
				}
			}
			a = hx_2 + hy_2;
			b = -2 * (hy_2 * x + hx_2 * y);
			c = hy_2 * pow(x, 2) + hx_2 * pow(y, 2) - hx_2 * hy_2 * coefSpeed;
			delta = pow(b, 2) - 4 * a * c;
			grad_x = x / hx;
			grad_y = y / hy;
			if (fabs(grad_x - grad_y) >= coefSpeed) {
				new_val = min(x + hx * coefSpeed, y + hy * coefSpeed);
			}
			else {
				new_val = ( -b + sqrt(delta) ) / (2 * a);
			}
			if (new_val < distancePtr[xd])
			{
				distancePtr[xd] = new_val;
			}
		}
	}

	//sweep 3
	for (int i = 0; i < length; i++)
	{
		for (int j = width_minus; j > -1; j--)
		{
			xd = x_new(i, j, length);
			if (i == 0) {
				x = min(distancePtr[xd], distancePtr[x_new(i + 1, j, length)]);
			}
			else {
				if (i == length_minus) {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[xd]);
				}
				else {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[x_new(i + 1, j, length)]);
				}
			}
			if (j == 0) {
				y = min(distancePtr[xd], distancePtr[x_new(i, j + 1, length)]);
			}
			else {
				if (j == width_minus) {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[xd]);
				}
				else {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[x_new(i, j + 1, length)]);
				}
			}
			a = hx_2 + hy_2;
			b = -2 * (hy_2 * x + hx_2 * y);
			c = hy_2 * pow(x, 2) + hx_2 * pow(y, 2) - hx_2 * hy_2 * coefSpeed;
			delta = pow(b, 2) - 4 * a * c;
			grad_x = x / hx;
			grad_y = y / hy;
			if (fabs(grad_x - grad_y) >= coefSpeed) {
				new_val = min(x + hx * coefSpeed, y + hy * coefSpeed);
			}
			else {
				new_val = ( -b + sqrt(delta) ) / (2 * a);
			}
			if (new_val < distancePtr[xd])
			{
				distancePtr[xd] = new_val;
			}
		}
	}

	//sweep 4
	for (int i = length_minus; i > -1; i--)
	{
		for (int j = width_minus; j > -1; j--)
		{
			xd = x_new(i, j, length);
			if (i == 0) {
				x = min(distancePtr[xd], distancePtr[x_new(i + 1, j, length)]);
			}
			else {
				if (i == length_minus) {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[xd]);
				}
				else {
					x = min(distancePtr[x_new(i - 1, j, length)], distancePtr[x_new(i + 1, j, length)]);
				}
			}
			if (j == 0) {
				y = min(distancePtr[xd], distancePtr[x_new(i, j + 1, length)]);
			}
			else {
				if (j == width_minus) {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[xd]);
				}
				else {
					y = min(distancePtr[x_new(i, j - 1, length)], distancePtr[x_new(i, j + 1, length)]);
				}
			}
			a = hx_2 + hy_2;
			b = -2 * (hy_2 * x + hx_2 * y);
			c = hy_2 * pow(x, 2) + hx_2 * pow(y, 2) - hx_2 * hy_2 * coefSpeed;
			delta = pow(b, 2) - 4 * a * c;
			grad_x = x / hx;
			grad_y = y / hy;
			if (fabs(grad_x - grad_y) >= coefSpeed) {
				new_val = min(x + hx * coefSpeed, y + hy * coefSpeed);
			}
			else {
				new_val = ( -b + sqrt(delta) ) / (2 * a);
			}
			if (new_val < distancePtr[xd])
			{
				distancePtr[xd] = new_val;
			}
		}
	}

	return true;
}

bool rouyTourinDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType tolerance, size_t max_iteration, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distancePtr == NULL)
		return false;

	size_t length = ctImageData.height;
	size_t width = ctImageData.width;
	size_t i, j, x;

	size_t length_ext = length + 2;
	size_t width_ext = length + 2;
	size_t i_ext, j_ext, x_ext;
	size_t dim2D = length_ext * width_ext;

	dataType* previousSolution = new dataType[dim2D]{ 0 };
	if (previousSolution == NULL) {
		return false;
	}

	double mass = 10.0;
	dataType hx = ctImageData.spacing.sx;
	dataType hy = ctImageData.spacing.sy;
	
	dataType hx_2 = 1.0 / (hx * hx);
	dataType hy_2 = 1.0 / (hy * hy);
	dataType value = 0.0;

	dataType tau = (hx * hy) / (2 * sqrt(hx * hx + hy * hy));
	std::cout << "Tau = " << tau << std::endl;

	size_t count_iteration = 0;

	while (mass > tolerance && count_iteration < max_iteration) {

		copyDataTo2dExtendedArea(distancePtr, previousSolution, width, length);
		reflection2D(previousSolution, length_ext, width_ext);

		count_iteration++;
		mass = 0.0;
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				x = x_new(i, j, length);
				x_ext = x_new(i_ext, j_ext, length_ext);
				if (ctImageData.imageDataPtr[x] != foregroundValue) {
					value = previousSolution[x_ext];
					distancePtr[x] = value + tau - tau * sqrt(hx_2 * max(min0(previousSolution[x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[x_new(i_ext + 1, j_ext, length_ext)], value))
						+ hy_2 * max(min0(previousSolution[x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[x_new(i_ext, j_ext + 1, length_ext)], value)));

					//Compute the mass
					mass += pow(previousSolution[x_ext] - distancePtr[x], 2);
				}
			}
		}
		mass = sqrt(mass);
	}
	std::cout << "Convergence is reached after : " << count_iteration << " iterations and the mass is : " << mass << std::endl;

	delete[] previousSolution;

	return true;
}
*/

/*
bool rouyTourinFrontPropagation2D(Image_Data2D ctImageData, dataType* distancePtr, dataType* potential, dataType tolerance, size_t max_iteration) {

	if (ctImageData.imageDataPtr == NULL || distancePtr == NULL)
		return false;

	size_t length = ctImageData.height;
	size_t width = ctImageData.width;
	size_t i, j, k, xd, xd_ext;
	size_t dim2D = length * width;

	size_t length_ext = length + 2;
	size_t width_ext = length + 2;
	size_t i_ext, j_ext, k_ext, x_ext;

	bool* status = new bool[dim2D]{ false };
	dataType* previousSolution = new dataType[length_ext * width_ext] {0};
	if(previousSolution == NULL || status == NULL)
	{
		return false;
	}

	//initialization
	//endPoints[0] = { 170.0, 12.0 };
	size_t indx_seed = x_new(175, 310, length);
	distancePtr[indx_seed] = 0.0;
	status[indx_seed] = true;

	double mass = 100.0;
	dataType hx = ctImageData.spacing.sx;
	dataType hy = ctImageData.spacing.sy;
	dataType value = 0.0;

	dataType hx_2 = 1.0 / (hx * hx);
	dataType hy_2 = 1.0 / (hy * hy);

	//dataType tau = hx * hy / sqrt(hx * hx + hy * hy);
	dataType tau = 0.5 * min(hx, hy);
	//std::cout << "tau = " << tau << std::endl;

	//Read the file of points with negative discriminant
	vector<Point2D> nd;
	dataType x, y;
	string path_discriminant = "C:/Users/Konan Allaly/Documents/Tests/output/negative_discrinant_simple_file.csv";
	FILE* pFile;
	if (fopen_s(&pFile, path_discriminant.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	while (feof(pFile) == 0) {
		fscanf_s(pFile, "%f", &x);
		fscanf_s(pFile, ",");
		fscanf_s(pFile, "%f", &y);
		fscanf_s(pFile, "\n");
		Point2D pt = { x, y };
		nd.push_back(pt);
	}
	//std::cout << "Number of points with negative discriminant: " << nd.size() << std::endl;
	fclose(pFile);
	
	string path_evolution = "C:/Users/Konan Allaly/Documents/Tests/output/action_evolution.csv";
	string header_name;
	FILE* dFile;
	if (fopen_s(&dFile, path_evolution.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(dFile, "ID,");
	for (i = 0; i < nd.size(); i++) 
	{
		if(i == nd.size() - 1) {
			header_name = "pt" + to_string(i+1);
			header_name += "\n";
		}
		else {
			header_name = "pt" + to_string(i+1) + ",";
		}
		fprintf(dFile, header_name.c_str());
	}
	
	size_t count_iteration = 0;
	size_t nb_pt_processed = 1;
	size_t xd_current;

	
	//while (nb_pt_processed > 0 && count_iteration < max_iteration) {
	//	copyDataTo2dExtendedArea(distancePtr, previousSolution, length, width);
	//	reflection2D(previousSolution, length_ext, width_ext);
	//	count_iteration++;
	//	mass = 0.0;
	//	nb_pt_processed = 0;
	//	for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
	//		for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
	//			xd = x_new(i_ext, j_ext, length_ext);
	//			x = x_new(i, j, length);
	//			if (status[x] == false)
	//			{
	//				nb_pt_processed++;
	//				value = previousSolution[xd];
	//				w = potential[x];
	//				distancePtr[x] = value + tau * (w - sqrt(hx_2 * max(min0(previousSolution[x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[x_new(i_ext + 1, j_ext, length_ext)], value))
	//					+ hy_2 * max(min0(previousSolution[x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[x_new(i_ext, j_ext + 1, length_ext)], value))));
	//
	//				ind_res = fabs(previousSolution[xd] - distancePtr[x]);
	//				if (ind_res <= tolerance) {
	//					status[x] = true;
	//				}
	//			}
	//		}
	//	}	
	//}

	while (mass > tolerance && count_iteration < max_iteration) {
		copyDataTo2dExtendedArea(distancePtr, previousSolution, length, width);
		reflection2D(previousSolution, length_ext, width_ext);
		count_iteration++;
		fprintf(dFile, "%d,", count_iteration);
		mass = 0.0;
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				xd_ext = x_new(i_ext, j_ext, length_ext);
				xd = x_new(i, j, length);
				if(xd != indx_seed)
				{
					value = previousSolution[xd_ext];
					distancePtr[xd] = value + tau * (potential[xd] - sqrt(hx_2 * max(min0(previousSolution[x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[x_new(i_ext + 1, j_ext, length_ext)], value))
						+ hy_2 * max(min0(previousSolution[x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[x_new(i_ext, j_ext + 1, length_ext)], value))));
				}
				////Compute the mass
				mass += pow(previousSolution[xd_ext] - distancePtr[xd], 2);
				//mass += abs(previousSolution[xd] - distancePtr[x]);
			}
		}
		mass = sqrt(mass);
		for (k = 0; k < nd.size(); k++)
		{
			xd_current = x_new((size_t)nd[k].x, (size_t)nd[k].y, length);
			if (k == nd.size() - 1) {
				fprintf(dFile, "%f\n", distancePtr[xd_current]);
			}
			else {
				fprintf(dFile, "%f,", distancePtr[xd_current]);
			}
		}
	}
	
	fclose(dFile);
	
	std::cout << "Iteration: " << count_iteration << ", Mass: " << mass << std::endl;

	delete[] previousSolution;
	delete[] status;

	return true;
}
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//=================== Functions for the 3D Fast Marching and Path Tracking ==============================
/////////////////////////////////////////////////////////////////////////////////////////////////////////

dataType upwindFiniteDifferenceX(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t x, const size_t y, const size_t z) {

	dataType x_minus, x_plus;
	if(x == 0) {
		x_minus = INFINITY;
	}
	else {
		x_minus = actionPtr[z][x_new(x - 1, y, length)];
	}
	if(x >= length - 1) {
		x_plus = INFINITY;
	}
	else {
		x_plus = actionPtr[z][x_new(x + 1, y, length)];
	}
	return min(x_minus, x_plus);
}

dataType upwindFiniteDifferenceY(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t x, const size_t y, const size_t z) {

	dataType y_minus, y_plus;
	if(y == 0) {
		y_minus = INFINITY;
	}
	else {
		y_minus = actionPtr[z][x_new(x, y - 1, length)];
	}
	if(y >= width - 1) {
		y_plus = INFINITY;
	}
	else {
		y_plus = actionPtr[z][x_new(x, y + 1, length)];
	}
	return min(y_minus, y_plus);
}

dataType upwindFiniteDifferenceZ(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t x, const size_t y, const size_t z) {
	
	size_t xd = x_new(x, y, length);
	dataType z_minus, z_plus;
	if(z == 0) {
		z_minus = INFINITY;
	}
	else {
		z_minus = actionPtr[z - 1][xd];
	}
	if(z >= height - 1) {
		z_plus = INFINITY;
	}
	else {
		z_plus = actionPtr[z + 1][xd];
	}
	return min(z_minus, z_plus);
}

dataType solve3dQuadraticEikonalEquation(dataType X, dataType Y, dataType Z, dataType P, VoxelSpacing h) {

	if (h.sx <= 0 || h.sy <= 0 || h.sz <= 0) {
		std::cout << "Error: Voxel spacing must be positive." << std::endl;
		return INFINITY; // Return a large value or handle the error appropriately
	}
	if (P <= 0) {
		std::cout << "Error: Propagation speed must be positive." << std::endl;
		return INFINITY; // Return a large value or handle the error appropriately
	}

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P_2 = P * P;
	dataType hx_2 = 1.0 / (h.sx * h.sx);
	dataType hy_2 = 1.0 / (h.sy * h.sy);
	dataType hz_2 = 1.0 / (h.sz * h.sz);

	if (X == INFINITY && Y == INFINITY && Z == INFINITY) {
		return INFINITY; // No solution if all coordinates are infinite
		std::cout << "Error: All coordinates are infinite." << std::endl;
	}

	if (X != INFINITY && Y == INFINITY && Z == INFINITY) {
		a = hx_2;
		b = (dataType)(-2 * X * hx_2);
		c = (dataType)(X * X * hx_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > X) {
				return solution;
			}
			else {
				return (dataType)(X + h.sx * P);
			}
		}
		else {
			return (dataType)(X + h.sx * P);
		}
	}

	if (Y != INFINITY && X == INFINITY && Z == INFINITY) {
		a = hy_2;
		b = (dataType)(-2 * Y * hy_2);
		c = (dataType)(Y * Y * hy_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > Y) {
				return solution;
			}
			else {
				return (dataType)(Y + h.sy * P);
			}
		}
		else {
			return (dataType)(Y + h.sy * P);
		}
	}

	if (Z != INFINITY && X == INFINITY && Y == INFINITY) {
		a = hz_2;
		b = (dataType)(-2 * Z * hz_2);
		c = (dataType)(Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > Z) {
				return solution;
			}
			else {
				return (dataType)(Z + h.sz * P);
			}
		}
		else {
			return (dataType)(Z + h.sz * P);
		}
	}

	if (X != INFINITY && Y != INFINITY && Z == INFINITY) {
		a = (dataType)(hx_2 + hy_2);
		b = (dataType)(-2 * (X * hx_2 + Y * hy_2));
		c = (dataType)(X * X * hx_2 + Y * Y * hy_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > max(X, Y)) {
				return solution;
			}
			else {
				return (dataType)(min(X + h.sx * P, Y + h.sy * P));
			}
		}
		else {
			return (dataType)(min(X + h.sx * P, Y + h.sy * P));
		}
	}

	if (X != INFINITY && Z != INFINITY && Y == INFINITY) {
		a = (dataType)(hx_2 + hz_2);
		b = (dataType)(-2 * (X * hx_2 + Z * hz_2));
		c = (dataType)(X * X * hx_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > max(X, Z)) {
				return solution;
			}
			else {
				return (dataType)(min(X + h.sx * P, Z + h.sz * P));
			}
		}
		else {
			return (dataType)(min(X + h.sx * P, Z + h.sz * P));
		}
	}

	if (Y != INFINITY && Z != INFINITY && X == INFINITY) {
		a = (dataType)(hy_2 + hz_2);
		b = (dataType)(-2 * (Y * hy_2 + Z * hz_2));
		c = (dataType)(Y * Y * hy_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > max(Y, Z)) {
				return solution;
			}
			else {
				return (dataType)(min(Y + h.sy * P, Z + h.sz * P));
			}
		}
		else {
			return (dataType)(min(Y + h.sy * P, Z + h.sz * P));
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		a = (dataType)(hx_2 + hy_2 + hz_2);
		b = (dataType)(-2 * (X * hx_2 + Y * hy_2 + Z * hz_2));
		c = (dataType)(X * X * hx_2 + Y * Y * hy_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution > max(X, max(Y, Z))) {
				return solution;
			}
			else {
				return (dataType)(min(X + h.sx * P, min(Y + h.sy * P, Z + h.sz * P)));
			}
		}
		else {
			return (dataType)(min(X + h.sx * P, min(Y + h.sy * P, Z + h.sz * P)));
		}
	}

}

void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b) {
	pointFastMarching3D temp = *a;
	*a = *b;
	*b = temp;
}

void heapifyDown3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex, int i) {

	int length_array = in_Process.size();
	if(length_array == 0) {
		return; //nothing to heapify
	}
	int current = i;
	int left_child = 2 * i + 1;
	int right_child = 2 * i + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}

	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}

	if (current != i) {
		//swap3dPoints(&in_Process[i], &in_Process[current]);
		//heapifyDown3D(in_Process, current);

		size_t idx1 = in_Process[i].index;
		size_t idx2 = in_Process[current].index;

		swap_elts(&heapIndex[idx1], &heapIndex[idx2], sizeof(int));
		swap_elts(&in_Process[i], &in_Process[current], sizeof(pointFastMarching3D));
		
		heapifyDown3D(in_Process, heapIndex, current);
	}

}

void heapifyUp3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex, int i) {

	if(i <= 0 || i >= in_Process.size()) {
		return; //nothing to heapify
	}
	int current = i;

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		//swap3dPoints(&in_Process[current], &in_Process[i]);
		//heapifyUp3D(in_Process, current);

		size_t ind1 = in_Process[i].index;
		size_t ind2 = in_Process[current].index;

		swap_elts(&heapIndex[ind1], &heapIndex[ind2], sizeof(int));
		swap_elts(&in_Process[current], &in_Process[i], sizeof(pointFastMarching3D));

		heapifyUp3D(in_Process, heapIndex, current);
	}
}

void heapifyVector3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex) {
	int length_array = in_Process.size();
	if(length_array < 2) {
		return; //nothing to heapify
	}
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		//heapifyDown3D(in_Process, ind);
		heapifyDown3D(in_Process, heapIndex, ind);
	}
}

void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	if (l > 1) {
		//swap3dPoints(&in_Process[0], &in_Process[l - 1]);

		size_t idx1 = in_Process[0].index;
		size_t idx2 = in_Process[l - 1].index;

		swap_elts(&heapIndex[idx1], &heapIndex[idx2], sizeof(int));
		swap_elts(&in_Process[0], &in_Process[l - 1], sizeof(pointFastMarching3D));

		heapIndex[idx2] = -1; //the point is removed from the heap
		in_Process.pop_back();

		//heapifyDown3D(in_Process, 0);
		heapifyDown3D(in_Process, heapIndex, 0);
	}
	else {
		in_Process.pop_back();
	}
}

void addPointHeap3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex, pointFastMarching3D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	size_t heapIndex_point = in_Process[l - 1].index;
	heapIndex[heapIndex_point] = l - 1;
	
	if (l < 2) 
	{
		return; //nothing to heapify
	}
	else 
	{
		//heapifyUp3D(in_Process, l - 1);
		heapifyUp3D(in_Process, heapIndex, l - 1);
	}
}

int getIndexFromHeap3D(vector<pointFastMarching3D>& in_Process, size_t i, size_t j, size_t k) {
	//we use type int for indexes because we do operations like pos--
	for (int ind = 0; ind < in_Process.size(); ind++) {
		if (in_Process[ind].x == i && in_Process[ind].y == j && in_Process[ind].z == k) {
			return ind;
		}
	}
	return -1; //not found
}

void updateNeighbor3D(size_t ind_x, size_t ind_y, size_t ind_z, size_t length, size_t width, size_t height,
	dataType** action, dataType** potential, short** labelArray,
	VoxelSpacing spacing, vector<pointFastMarching3D>& narrowBand, vector<int>& heapIndex)
{

	if (ind_x >= length || ind_y >= width || ind_z >= height)
		return;

	size_t xd = x_new(ind_x, ind_y, length);
	dataType ux = upwindFiniteDifferenceX(action, length, width, height, ind_x, ind_y, ind_z);
	dataType uy = upwindFiniteDifferenceY(action, length, width, height, ind_x, ind_y, ind_z);
	dataType uz = upwindFiniteDifferenceZ(action, length, width, height, ind_x, ind_y, ind_z);
	dataType coefSpeed = potential[ind_z][xd];
	dataType solution = solve3dQuadraticEikonalEquation(ux, uy, uz, coefSpeed, spacing);
	size_t pos = x_flat(ind_x, ind_y, ind_z, length, width);
	pointFastMarching3D neighbor = { ind_x, ind_y, ind_z, pos, solution };
	if (labelArray[ind_z][xd] == 3) 
	{
		//addPointHeap3D(narrowBand, neighbor);
		addPointHeap3D(narrowBand, heapIndex, neighbor);
		action[ind_z][xd] = solution;
		labelArray[ind_z][xd] = 2;
	}
	else if (labelArray[ind_z][xd] == 2 && solution < action[ind_z][xd]) 
	{
		action[ind_z][xd] = solution;
		//int pIndex = getIndexFromHeap3D(narrowBand, ind_x, ind_y, ind_z);
		int pIndex = heapIndex[pos];
		if (pIndex != -1) 
		{
			//heapifyUp3D(narrowBand, pIndex);
			heapifyUp3D(narrowBand, heapIndex, pIndex);
		}
	}
}

bool compute3DPotential(Image_Data ctImageData, dataType** potential, Point3D* seedPoint, Potential_Parameters parameters) {

	if (ctImageData.imageDataPtr == NULL || potential == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, xd = 0;
	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t dim2D = length * width;

	dataType** maskThreshold = new dataType * [height];
	dataType** distance = new dataType * [height];
	for (k = 0; k < height; k++) {
		distance[k] = new dataType[dim2D]{ 0 };
		maskThreshold[k] = new dataType[dim2D]{ 0 };
	}
	if (distance == NULL || maskThreshold == NULL)
		return false;

	dataType norm_of_gradient = 0.0, edgeValue = 0.0;
	bool isGradientComputed = false;
	Point3D grad_vector;
	
	for (k = 0; k < height; k++) 
	{
		for (i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				xd = x_new(i, j, length);
				isGradientComputed = getGradient3D(ctImageData, i, j, k, &grad_vector);
				if (isGradientComputed == true) {
					norm_of_gradient = sqrt(grad_vector.x * grad_vector.x + grad_vector.y * grad_vector.y + grad_vector.z * grad_vector.z);
				}
				else {
					std::cout << "Error in computing gradient at point (" << i << ", " << j << ", " << k << ")" << std::endl;
					return false;
				}
				
				dataType edgeValue = gradientFunction(norm_of_gradient, parameters.K);
				//maskThreshold[k][xd] = edgeValue;
				//threshold : real image
				if (edgeValue <= parameters.thres) {
					maskThreshold[k][xd] = 1.0;
				}
				else 
				{
					maskThreshold[k][xd] = 0.0;
				}
			}
		}
	}
	
	////Real image
	Image_Data toDistanceMap = { height, length, width, maskThreshold, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/Data journal paper submission/edge_image_p4.raw";
	manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), STORE_DATA, false);
	
	fastMarchingDistanceMap(toDistanceMap, distance, 1.0);
	storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/Data journal paper submission/distance_map_p4.raw";
	manageRAWFile3D<dataType>(distance, length, width, height, storing_path.c_str(), STORE_DATA, false);

	////////Artificial image : no need to compute the edge image when empty inside
	//Image_Data toDistanceMap = { height, length, width, ctImageData.imageDataPtr, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	//std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/edge_image.raw";
	////manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), STORE_DATA, false);
	//fastMarchingDistanceMap(toDistanceMap, distance, 1.0);
	//storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_fm.raw";
	//manageRAWFile3D<dataType>(distance, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };
	seedStats = getPointNeighborhoodStats(ctImageData, seedPoint[0], parameters.radius);
	dataType value_first_pt = seedStats.mean_data;
	std::cout << "Seed 1 value : " << value_first_pt << std::endl;
	dataType seedValCT = value_first_pt;
	
	dataType var_epsilon = 0;//0.01;
	if (seedStats.sd_data != 0) 
	{
		var_epsilon = seedStats.sd_data;
	}
	else {
		var_epsilon = parameters.eps;
	}
	
	std::cout << "Epsilon to be used : " << var_epsilon << std::endl;

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = fabs(seedValCT - ctImageData.imageDataPtr[k][i]);
		}
	}

	//Find the max of the difference
	dataType maxImage = 0.0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (potential[k][i] > maxImage) 
			{
				maxImage = potential[k][i];
			}
		}
	}
	std::cout << "Max pot: " << maxImage << std::endl;

	for (k = 0; k < height; k++) 
	{
		for (i = 0; i < dim2D; i++) 
		{
			potential[k][i] = (dataType)(var_epsilon + (potential[k][i] / maxImage)) * (1.0 / (1.0 + distance[k][i]));
		}
	}

	for (k = 0; k < height; k++) {
		delete[] maskThreshold[k];
		delete[] distance[k];
	}
	delete[] maskThreshold;
	delete[] distance;

	return true;
}

bool shortestPath3D(Image_Data actionMapStr, Point3D* seedPoints, vector<Point3D>& path_points, Path_Parameters parameters) {

	if (actionMapStr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	const size_t length = actionMapStr.length;
	const size_t width = actionMapStr.width;
	const size_t height = actionMapStr.height;
	VoxelSpacing spacing = actionMapStr.spacing;

	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	//Find the closest point till the last point
	i = (size_t)seedPoints[1].x;
	j = (size_t)seedPoints[1].y;
	k = (size_t)seedPoints[1].z;
	size_t currentIndx = x_new(i, j, length);

	dataType x = seedPoints[1].x;
	dataType y = seedPoints[1].y;
	dataType z = seedPoints[1].z;
	double dist_to_end = 0.0;

	size_t count_iter = 1;
	Point3D final_point = getRealCoordFromImageCoord3D(seedPoints[0], actionMapStr.origin, spacing, actionMapStr.orientation);

	bool isGradientComputed = false;
	Point3D grad_vector;
	dataType norm_of_gradient = 0.0;
	do {

		isGradientComputed = getGradient3D(actionMapStr, i, j, k, &grad_vector);
		if (isGradientComputed == true) {
			norm_of_gradient = sqrt(grad_vector.x * grad_vector.x + grad_vector.y * grad_vector.y + grad_vector.z * grad_vector.z);
		}
		else {
			std::cout << "Error in computing gradient at point (" << i << ", " << j << ", " << k << ")" << std::endl;
			return false;
		}

		x -= parameters.tau * (grad_vector.x / norm_of_gradient);
		y -= parameters.tau * (grad_vector.y / norm_of_gradient);
		z -= parameters.tau * (grad_vector.z / norm_of_gradient);

		if(x < 0.0 || x >= length || y < 0.0 || y >= width || z < 0.0 || z >= height) 
		{
			std::cout << "Error: Point out of bounds (" << x << ", " << y << ", " << z << ")" << std::endl;
			return false;
		}
		
		Point3D point_current = { x, y, z };
		Point3D pDistance = getRealCoordFromImageCoord3D(point_current, actionMapStr.origin, spacing, actionMapStr.orientation);
		
		//compute distance current Point - last point
		dist_to_end = getPoint3DDistance(pDistance, final_point);
		
		i = (size_t)(round(x));
		j = (size_t)(round(y));
		k = (size_t)(round(z));
		currentIndx = x_new(i, j, length);
		Point3D point_save = { i, j, k };
		path_points.push_back(point_current);
		
		count_iter++;

	} while (dist_to_end > parameters.tolerance && count_iter < parameters.max_iteration);

	return true;
}

bool partialFrontPropagation(Image_Data actionPtr, dataType** potentialFuncPtr, Point3D* endPoints) {

	if (actionPtr.imageDataPtr == NULL || potentialFuncPtr == NULL || endPoints == NULL) {
		return false;
	}

	const size_t height = actionPtr.height;
	const size_t length = actionPtr.length;
	const size_t width = actionPtr.width;
	VoxelSpacing spacing = actionPtr.spacing;

	vector <pointFastMarching3D> narrowBand;
	size_t i = 0, j = 0, k = 0, dim2D = length * width, dim3D = length * width * height;

	short** labelArray = new short* [height];
	if (labelArray == NULL) {
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
		if (labelArray[k] == NULL) {
			return false; // Memory allocation failed
		}
	}
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			actionPtr.imageDataPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	vector<int> heapIndex(dim3D);
	for (i = 0; i < dim3D; i++) {
		heapIndex[i] = -1;
	}

	//Processed the starting point
	if (endPoints[0].x < 0 || endPoints[0].x > length || 
		endPoints[0].y < 0 || endPoints[0].y > width || 
		endPoints[0].z < 0 || endPoints[0].z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	i = (size_t)endPoints[0].x;
	j = (size_t)endPoints[0].y;
	k = (size_t)endPoints[0].z;
	size_t currentIndx = x_new(i, j, length);
	actionPtr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	if (k > 0) 
	{
		if(labelArray[k - 1][currentIndx] != 1) 
		{
			updateNeighbor3D(i, j, k - 1, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (k < (height - 1))
	{
		if(labelArray[k + 1][currentIndx] != 1) {
			updateNeighbor3D(i, j, k + 1, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (i > 0)
	{
		if(labelArray[k][x_new(i - 1, j, length)] != 1) {
			updateNeighbor3D(i - 1, j, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (i < (length - 1))
	{
		if(labelArray[k][x_new(i + 1, j, length)] != 1) {
			updateNeighbor3D(i + 1, j, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (j > 0)
	{
		if (labelArray[k][x_new(i, j - 1, length)] != 1) {
			updateNeighbor3D(i, j - 1, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (j < (width - 1))
	{
		if (labelArray[k][x_new(i, j + 1, length)] != 1) {
			updateNeighbor3D(i, j + 1, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}

	if (endPoints[1].x < 0 || endPoints[1].x > length || 
		endPoints[1].y < 0 || endPoints[1].y > width || 
		endPoints[1].z < 0 || endPoints[1].z > height) 
	{
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	size_t seedI = (size_t)endPoints[1].x;
	size_t seedJ = (size_t)endPoints[1].y;
	size_t seedK = (size_t)endPoints[1].z;
	size_t seedIndex = x_new(seedI, seedJ, length);

	dataType max_action = 0;
	size_t num_proceed = 1;

	//FILE* processed_points;
	//string save_points = "C:/Users/Konan Allaly/Documents/Tests/output/processed_partial.csv";
	//if (fopen_s(&processed_points, save_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(processed_points, "x,y,z\n");
	//Point3D processed = getRealCoordFromImageCoord3D(endPoints[0], actionPtr.origin, actionPtr.spacing, actionPtr.orientation);
	//fprintf(processed_points, "%f,%f,%f\n", processed.x, processed.y, processed.z);
	
	//string save_front = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//dataType** saveAction = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	saveAction[k] = new dataType[dim2D]{ 0 };
	//}

	while (narrowBand.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching3D current = narrowBand[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		num_proceed++;

		////store the processed points
		//processed = { (dataType)i, (dataType)j, (dataType)k };
		//processed = getRealCoordFromImageCoord3D(processed, actionPtr.origin, actionPtr.spacing, actionPtr.orientation);
		//fprintf(processed_points, "%f,%f,%f\n", processed.x, processed.y, processed.z);

		//Find the maximum action value
		if (actionPtr.imageDataPtr[k][currentIndx] > max_action) {
			max_action = actionPtr.imageDataPtr[k][currentIndx];
		}
		
		//Exit the while loop when the final point is reached/computed
		if (labelArray[seedK][seedIndex] == 1) {
			break;
		}

		deleteRootHeap3D(narrowBand, heapIndex);

		if (k > 0)
		{
			if (labelArray[k - 1][currentIndx] != 1) 
			{
				updateNeighbor3D(i, j, k - 1, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (k < (height - 1))
		{
			if (labelArray[k + 1][currentIndx] != 1) 
			{
				updateNeighbor3D(i, j, k + 1, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (i > 0)
		{
			if (labelArray[k][x_new(i - 1, j, length)] != 1) 
			{
				updateNeighbor3D(i - 1, j, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (i < (length - 1))
		{
			if (labelArray[k][x_new(i + 1, j, length)] != 1) 
			{
				updateNeighbor3D(i + 1, j, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (j > 0)
		{
			if (labelArray[k][x_new(i, j - 1, length)] != 1) 
			{
				updateNeighbor3D(i, j - 1, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (j < (width - 1))
		{
			if (labelArray[k][x_new(i, j + 1, length)] != 1) 
			{
				updateNeighbor3D(i, j + 1, k, length, width, height, actionPtr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}

		//if(num_proceed == 30000)
		//{
		//	string save_first = save_front + "front1.raw";
		//	for(int ik = 0; ik < height; ik++)
		//	{
		//		for(int ij = 0; ij < dim2D; ij++)
		//		{
		//			if(labelArray[ik][ij] == 1)
		//			{
		//				saveAction[ik][ij] = actionPtr.imageDataPtr[ik][ij];
		//			}
		//			else {
		//				saveAction[ik][ij] = max_action + 1;
		//			}
		//		}
		//	}
		//	manageRAWFile3D<dataType>(saveAction, length, width, height, save_first.c_str(), STORE_DATA, false);
		//}
		//else if (num_proceed == 300000)
		//{
		//	string save_second = save_front + "front2.raw";
		//	for (int ik = 0; ik < height; ik++)
		//	{
		//		for (int ij = 0; ij < dim2D; ij++)
		//		{
		//			if (labelArray[ik][ij] == 1)
		//			{
		//				saveAction[ik][ij] = actionPtr.imageDataPtr[ik][ij];
		//			}
		//			else {
		//				saveAction[ik][ij] = max_action + 1;
		//			}
		//		}
		//	}
		//	manageRAWFile3D<dataType>(saveAction, length, width, height, save_second.c_str(), STORE_DATA, false);
		//}

	}

	//for(k = 0; k < height; k++) {
	//	delete[] saveAction[k];
	//}
	//delete[] saveAction;

	//fclose(processed_points);

	//std::cout << num_proceed << " have been processed before the front reaches final points." << std::endl;
	
	for (k = 0; k < height; k++) {
		for(i = 0; i < dim2D; i++) {
			if (actionPtr.imageDataPtr[k][i] == INFINITY) {
				actionPtr.imageDataPtr[k][i] = max_action + 1;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;

}

bool frontPropagation(Image_Data ctImageData, dataType** actionPtr, dataType** potentialFuncPtr, Point3D seedPoint) 
{

	if (actionPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> narrowBand;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	short** labelArray = new short* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
		if (labelArray[k] == NULL) {
			return false; // Memory allocation failed
		}
	}

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			actionPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the starting point
	if (seedPoint.x < 0 || seedPoint.x > length || seedPoint.y < 0 || seedPoint.y > width || seedPoint.z < 0 || seedPoint.z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	i = (size_t)seedPoint.x;
	j = (size_t)seedPoint.y;
	k = (size_t)seedPoint.z;
	size_t xd = x_new(i, j, length);
	actionPtr[k][xd] = 0.0;
	labelArray[k][xd] = 1;

	size_t dim3D = length * width * height;
	vector<int> heapIndex(dim3D);
	for (size_t n = 0; n < dim3D; n++) {
		heapIndex[n] = -1;
	}

	//find the neighbours of the initial point add add them to inProces

	if (k > 0)
	{
		if(labelArray[k - 1][xd] != 1)
		{
			updateNeighbor3D(i, j, k - 1, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (k < (height - 1))
	{
		if (labelArray[k + 1][xd] != 1)
		{
			updateNeighbor3D(i, j, k + 1, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (i > 0)
	{
		if (labelArray[k][x_new(i - 1, j, length)] != 1)
		{
			updateNeighbor3D(i - 1, j, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (i < (length - 1))
	{
		if (labelArray[k][x_new(i + 1, j, length)] != 1)
		{
			updateNeighbor3D(i + 1, j, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (j > 0)
	{
		if (labelArray[k][x_new(i, j - 1, length)] != 1)
		{
			updateNeighbor3D(i, j - 1, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
	}
	if (j < (width - 1))
	{
		if (labelArray[k][x_new(i, j + 1, length)] != 1)
		{
			updateNeighbor3D(i, j + 1, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}	
	}

	while (narrowBand.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching3D current = narrowBand[0];
		i = current.x;
		j = current.y;
		k = current.z;
		size_t currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		deleteRootHeap3D(narrowBand, heapIndex);

		if (k > 0)
		{
			if (labelArray[k - 1][xd] != 1)
			{
				updateNeighbor3D(i, j, k - 1, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (k < (height - 1))
		{
			if (labelArray[k + 1][xd] != 1)
			{
				updateNeighbor3D(i, j, k + 1, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (i > 0)
		{
			if (labelArray[k][x_new(i - 1, j, length)] != 1)
			{
				updateNeighbor3D(i - 1, j, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (i < (length - 1))
		{
			if (labelArray[k][x_new(i + 1, j, length)] != 1)
			{
				updateNeighbor3D(i + 1, j, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
			
		}
		if (j > 0)
		{
			if (labelArray[k][x_new(i, j - 1, length)] != 1)
			{
				updateNeighbor3D(i, j - 1, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
		if (j < (width - 1))
		{
			if (labelArray[k][x_new(i, j + 1, length)] != 1)
			{
				updateNeighbor3D(i, j + 1, k, length, width, height, actionPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

bool frontPropagationWithKeyPointDetection(Image_Data actionMapStr, dataType** potentialFuncPtr, Point3D* seedPoint, const double LengthKeyPoints, vector<Point3D>& key_points) {

	if (actionMapStr.imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoint == NULL) {
		return false;
	}
	
	const size_t length = actionMapStr.length;
	const size_t width = actionMapStr.width;
	const size_t height = actionMapStr.height;
	VoxelSpacing spacing = actionMapStr.spacing;
	if (seedPoint[0].x < 0 || seedPoint[0].x > length || seedPoint[0].y < 0 || seedPoint[0].y > width || seedPoint[0].z < 0 || seedPoint[0].z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	if (seedPoint[1].x < 0 || seedPoint[1].x > length || seedPoint[1].y < 0 || seedPoint[1].y > width || seedPoint[1].z < 0 || seedPoint[1].z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}

	vector <pointFastMarching3D> narrowBand;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
	}
	if (labelArray == NULL)
		return false;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			actionMapStr.imageDataPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	pointFastMarching3D current;
	i = (size_t)seedPoint[0].x;
	j = (size_t)seedPoint[0].y;
	k = (size_t)seedPoint[0].z;
	size_t currentIndx = x_new(i, j, length);
	actionMapStr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	size_t dim3D = length * width * height;
	vector<int> heapIndex(dim3D);
	for (size_t n = 0; n < dim3D; n++) {
		heapIndex[n] = -1;
	}

	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	if (k > 0)
	{
		updateNeighbor3D(i, j, k - 1, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}
	if (k < height_minus)
	{
		updateNeighbor3D(i, j, k + 1, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}
	if (i > 0)
	{
		updateNeighbor3D(i - 1, j, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}
	if (i < length_minus)
	{
		updateNeighbor3D(i + 1, j, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}
	if (j > 0)
	{
		updateNeighbor3D(i, j - 1, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}
	if (j < width_minus)
	{
		updateNeighbor3D(i, j + 1, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
	}

	size_t iEnd = (size_t)seedPoint[1].x;
	size_t jEnd = (size_t)seedPoint[1].y;
	size_t kEnd = (size_t)seedPoint[1].z;

	//Visualize the front propagation
	size_t id_keyPoint = 1;
	std::string storing_path;
	vector<Point3D> savingList;

	double distanceToCurrentSourcePoint = 0.0;
	//Set the starting point as initial source point
	Point3D currentSourcePoint = { (dataType)i, (dataType)j, (dataType)k };
	key_points.push_back(currentSourcePoint);

	dataType max_action = 0;
	size_t id_key = 0;

	//FILE* processed_points;
	//string save_points = "C:/Users/Konan Allaly/Documents/Tests/output/processed_keyp.csv";
	//if (fopen_s(&processed_points, save_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(processed_points, "x,y,z\n");
	//Point3D processed = getRealCoordFromImageCoord3D(seedPoint[0], actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
	//fprintf(processed_points, "%f,%f,%f\n", processed.x, processed.y, processed.z);

	//string save_front = "C:/Users/Konan Allaly/Documents/Tests/output/front_";
	//string save_front_key;
	//dataType** saveAction = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	saveAction[k] = new dataType[dim2D]{ 0 };
	//}
	
	while (narrowBand.size() > 0) {

		//processed the point with minimum distance
		current = narrowBand[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;

		//processed = { (dataType)i, (dataType)j, (dataType)k };
		//processed = getRealCoordFromImageCoord3D(processed, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		//fprintf(processed_points, "%f,%f,%f\n", processed.x, processed.y, processed.z);

		if(max_action < actionMapStr.imageDataPtr[k][currentIndx]) {
			max_action = actionMapStr.imageDataPtr[k][currentIndx];
		}

		//Exit of the while loop if the end point is reached
		if (labelArray[kEnd][x_new(iEnd, jEnd, length)] == 1) {
			break;
		}

		//check if the current point is a key point
		Point3D pSource = { i, j, k };
		Point3D pSourceReal = getRealCoordFromImageCoord3D(pSource, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		Point3D pCurrent = getRealCoordFromImageCoord3D(currentSourcePoint, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		distanceToCurrentSourcePoint = getPoint3DDistance(pCurrent, pSourceReal);
		
		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			id_key++;

			//If the condition is true ---> new key point is found so we need to initilize it neighbors
			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionMapStr.imageDataPtr[k][currentIndx] = 0;

			//Top
			if (k > 0)
			{
				size_t kminus = k - 1;
				actionMapStr.imageDataPtr[kminus][currentIndx] = INFINITY;
				labelArray[kminus][currentIndx] = 3;
			}

			//Bottom
			if (k < height_minus)
			{
				size_t kplus = k + 1;
				actionMapStr.imageDataPtr[kplus][currentIndx] = INFINITY;
				labelArray[kplus][currentIndx] = 3;
			}

			//East
			if (i < length_minus)
			{
				size_t iplus = i + 1;
				size_t indxEast = x_new(iplus, j, length);
				actionMapStr.imageDataPtr[k][indxEast] = INFINITY;
				labelArray[k][indxEast] = 3;
			}

			//North
			if (j > 0)
			{
				size_t jminus = j - 1;
				size_t indxNorth = x_new(i, jminus, length);
				actionMapStr.imageDataPtr[k][indxNorth] = INFINITY;
				labelArray[k][indxNorth] = 3;
			}

			//West
			if (i > 0)
			{
				size_t iminus = i - 1;
				size_t indxWest = x_new(iminus, j, length);
				actionMapStr.imageDataPtr[k][indxWest] = INFINITY;
				labelArray[k][indxWest] = 3;
			}

			//South
			if (j < width_minus)
			{
				size_t jplus = j + 1;
				size_t indxSouth = x_new(i, jplus, length);
				actionMapStr.imageDataPtr[k][indxSouth] = INFINITY;
				labelArray[k][indxSouth] = 3;
			}

			narrowBand.clear();

			//save_front_key = save_front + to_string(id_key) + ".raw";
			//for (int ik = 0; ik < height; ik++)
			//{
			//	for (int ij = 0; ij < dim2D; ij++)
			//	{
			//		if (labelArray[ik][ij] == 1)
			//		{
			//			saveAction[ik][ij] = actionMapStr.imageDataPtr[ik][ij];
			//		}
			//		else {
			//			saveAction[ik][ij] = max_action + 1;
			//		}
			//	}
			//}
			//manageRAWFile3D<dataType>(saveAction, length, width, height, save_front_key.c_str(), STORE_DATA, false);

		}
		else {
			actionMapStr.imageDataPtr[k][currentIndx] = current.arrival;
			deleteRootHeap3D(narrowBand, heapIndex);
		}

		//====================
		//processed neighbors of the minimum in the narrow band
		if (k > 0)
		{
			updateNeighbor3D(i, j, k - 1, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
		if (k < height_minus)
		{
			updateNeighbor3D(i, j, k + 1, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
		if (i > 0)
		{
			updateNeighbor3D(i - 1, j, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
		if (i < length_minus)
		{
			updateNeighbor3D(i + 1, j, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
		if (j > 0)
		{
			updateNeighbor3D(i, j - 1, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}
		if (j < width_minus)
		{
			updateNeighbor3D(i, j + 1, k, length, width, height, actionMapStr.imageDataPtr, potentialFuncPtr, labelArray, spacing, narrowBand, heapIndex);
		}

	}

	key_points.push_back(seedPoint[1]);

	//set all the unvisited points to a maximum value
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (actionMapStr.imageDataPtr[k][i] == INFINITY) {
				actionMapStr.imageDataPtr[k][i] = max_action + 1;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

/*
bool doubleFrontPropagation(Image_Data imageData, dataType** actionFirstFront, dataType** actionSecondFront, dataType** potentialPtr, Point3D* endPoints, string savingPath) {

	if (imageData.imageDataPtr == NULL || actionFirstFront == NULL || actionSecondFront == NULL || potentialPtr == NULL || endPoints == NULL) {
		return false;
	}

	const size_t length = imageData.length;
	const size_t width = imageData.width;
	const size_t height = imageData.height;
	VoxelSpacing spacing = imageData.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;
	size_t height_minus = height - 1;

	short** firstLabelArray = new short*[height];
	short** secondLabelArray = new short*[height];
	for (size_t k = 0; k < height; k++) {
		firstLabelArray[k] = new short[dim2D];
		secondLabelArray[k] = new short[dim2D];
		if (firstLabelArray[k] == NULL || secondLabelArray[k] == NULL) {
			return false;
		}
	}
	
	//Initialize action
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < dim2D; i++) {
			actionFirstFront[k][i] = INFINITY;
			actionSecondFront[k][i] = INFINITY;
			firstLabelArray[k][i] = 3; //3 ---> not processed
			secondLabelArray[k][i] = 3; //3 ---> not processed
		}
	}

	size_t x1 = (size_t)endPoints[0].x;
	size_t y1 = (size_t)endPoints[0].y;
	size_t z1 = (size_t)endPoints[0].z;
	if (x1 >= length || y1 >= width || z1 >= height) {
		delete[] firstLabelArray;
		delete[] secondLabelArray;
		return false; //Invalid end point
	}
	firstLabelArray[z1][x_new(x1, y1, length)] = 1; //1 ---> already processed
	actionFirstFront[z1][x_new(x1, y1, length)] = 0.0;

	size_t x2 = (size_t)endPoints[1].x;
	size_t y2 = (size_t)endPoints[1].y;
	size_t z2 = (size_t)endPoints[1].z;
	if (x2 >= length || y2 >= width || z2 >= height) {
		delete[] secondLabelArray;
		delete[] firstLabelArray;
		return false; //Invalid end point
	}
	secondLabelArray[z2][x_new(x2, y2, length)] = 1; //1 ---> already processed
	actionSecondFront[z2][x_new(x2, y2, length)] = 0.0;

	vector<pointFastMarching3D> narrowBandFirstFront, narrowBandSecondFront;

	//Initialize neighbors for the first front

	if(x1 > 0)
	{
		updateNeighbor3D(x1 - 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (x1 < length_minus)
	{
		updateNeighbor3D(x1 + 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (y1 > 0)
	{
		updateNeighbor3D(x1, y1 - 1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (y1 < width_minus)
	{
		updateNeighbor3D(x1, y1 + 1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (z1 > 0)
	{
		updateNeighbor3D(x1, y1, z1 - 1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}
	if (z1 < height_minus)
	{
		updateNeighbor3D(x1 - 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
	}

	//Initialize neighbors for the second front

	if (x2 > 0)
	{
		updateNeighbor3D(x2 - 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (x2 < length_minus)
	{
		updateNeighbor3D(x2 + 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (y2 > 0)
	{
		updateNeighbor3D(x2, y2 - 1, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (y2 < width_minus)
	{
		updateNeighbor3D(x2, y2 + 1, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (z2 > 0)
	{
		updateNeighbor3D(x2, y2, z2 - 1, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	if (z2 < height_minus)
	{
		updateNeighbor3D(x2 - 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
	}
	

	dataType max_save_action = 0.0;

	while (narrowBandFirstFront.size() > 0 && narrowBandSecondFront.size() > 0) {

		pointFastMarching3D firstFrontPoint = narrowBandFirstFront[0];
		x1 = firstFrontPoint.x;
		y1 = firstFrontPoint.y;
		z1 = firstFrontPoint.z;
		size_t indexFirst = x_new(x1, y1, length);
		if (secondLabelArray[z1][indexFirst] == 1) {
			//The fronts have met
			endPoints[2].x = x1;
			endPoints[2].y = y1;
			endPoints[2].z = z1;
			break;
		}
		else {
			firstLabelArray[z1][indexFirst] = 1;
		}

		if (firstFrontPoint.arrival > max_save_action) {
			max_save_action = firstFrontPoint.arrival;
		}
		deleteRootHeap3D(narrowBandFirstFront);

		pointFastMarching3D secondFrontPoint = narrowBandSecondFront[0];
		x2 = secondFrontPoint.x;
		y2 = secondFrontPoint.y;
		z2 = secondFrontPoint.z;
		size_t indexSecond = x_new(x2, y2, length);
		if (firstLabelArray[z2][indexSecond] == 1) {
			//The fronts have met
			endPoints[2].x = x2;
			endPoints[2].y = y2;
			endPoints[2].z = z2;
			break;
		}
		else {
			secondLabelArray[z2][indexSecond] = 1;
		}

		if (secondFrontPoint.arrival > max_save_action) {
			max_save_action = secondFrontPoint.arrival;
		}
		deleteRootHeap3D(narrowBandSecondFront);
		

		//Update neighbors for the first front

		if (x1 > 0)
		{
			updateNeighbor3D(x1 - 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (x1 < length_minus)
		{
			updateNeighbor3D(x1 + 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (y1 > 0)
		{
			updateNeighbor3D(x1, y1 - 1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (y1 < width_minus)
		{
			updateNeighbor3D(x1, y1 + 1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (z1 > 0)
		{
			updateNeighbor3D(x1, y1, z1 - 1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}
		if (z1 < height_minus)
		{
			updateNeighbor3D(x1 - 1, y1, z1, length, width, height, actionFirstFront, potentialPtr, firstLabelArray, spacing, narrowBandFirstFront);
		}

		//Update neighbors for the second front

		if (x2 > 0)
		{
			updateNeighbor3D(x2 - 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (x2 < length_minus)
		{
			updateNeighbor3D(x2 + 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (y2 > 0)
		{
			updateNeighbor3D(x2, y2 - 1, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (y2 < width_minus)
		{
			updateNeighbor3D(x2, y2 + 1, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (z2 > 0)
		{
			updateNeighbor3D(x2, y2, z2 - 1, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}
		if (z2 < height_minus)
		{
			updateNeighbor3D(x2 - 1, y2, z2, length, width, height, actionSecondFront, potentialPtr, secondLabelArray, spacing, narrowBandSecondFront);
		}


	}

	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < dim2D; i++) {
			if (actionFirstFront[k][i] == INFINITY) {
				actionFirstFront[k][i] = max_save_action + 1;
			}
			if (actionSecondFront[k][i] == INFINITY) {
				actionSecondFront[k][i] = max_save_action + 1;
			}
		}
	}
	narrowBandFirstFront.clear();
	narrowBandSecondFront.clear();

	delete[] firstLabelArray;
	delete[] secondLabelArray;

	return true;
}
*/

bool rouyTourinDistanceMap(Image_Data ctImageData, dataType** distancePtr, dataType tolerance, size_t max_iteration, dataType foregroundValue) {

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

	dataType** previousSolution = new dataType * [height_ext];
	for (k = 0; k < height_ext; k++) {
		previousSolution[k] = new dataType[length_ext * width_ext]{ 0 };
		if (previousSolution[k] == NULL) {
			return false;
		}
	}
	if(previousSolution == NULL) {
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

	dataType tau = hx * hy * hz / ( 2.0 * sqrt(hx * hx + hy * hy + hz * hz) );
	std::cout << "tau = " << tau << std::endl;

	size_t count_iteration = 0;

	while (mass > tolerance && count_iteration < max_iteration) {
		copyDataToExtendedArea(distancePtr, previousSolution, height, length, width);
		reflection3D(previousSolution, height_ext, length_ext, width_ext);
		count_iteration++;
		mass = 0.0;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					if (ctImageData.imageDataPtr[k][x_new(i, j, length)] != foregroundValue) {
						value = previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)];
						distancePtr[k][x_new(i, j, length)] = value + tau - tau * sqrt(hx_2 * max(min0(previousSolution[k_ext][x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext + 1, j_ext, length_ext)], value))
							+ hy_2 * max(min0(previousSolution[k_ext][x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext, j_ext + 1, length_ext)], value))
							+ hz_2 * max(min0(previousSolution[k_ext - 1][x_new(i_ext, j_ext, length_ext)], value), min0(previousSolution[k_ext + 1][x_new(i_ext, j_ext, length_ext)], value)));
						//Compute the mass
						mass += pow(previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)] - distancePtr[k][x_new(i, j, length)], 2);
					}
				}
			}
		}
		mass = sqrt(mass);
	}
	std::cout << "Iteration: " << count_iteration << ", Mass: " << mass << std::endl;

	for (k = 0; k < height_ext; k++) {
		delete[] previousSolution[k];
	}
	delete[] previousSolution;

	return true;
}

bool fastMarchingDistanceMap(Image_Data ctImageData, dataType** distanceFuncPtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distanceFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short* [height];
	dataType** potentialFuncPtr = new dataType * [height];
	if(labelArray == NULL || potentialFuncPtr == NULL) {
		return false; // Memory allocation failed
	}
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
		potentialFuncPtr[k] = new dataType[dim2D];
		if (labelArray[k] == NULL || potentialFuncPtr[k] == NULL) {
			return false; // Memory allocation failed
		}
	}

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (ctImageData.imageDataPtr[k][i] == foregroundValue) 
			{
				distanceFuncPtr[k][i] = 0;
				labelArray[k][i] = 1;
			}
			else 
			{
				distanceFuncPtr[k][i] = INFINITY;
				labelArray[k][i] = 3;
			}
			potentialFuncPtr[k][i] = 1.0;//For distance map, the potential is set equal to 1.0
		}
	}

	size_t dim3D = length * width * height;
	vector<int> heapIndex(dim3D);
	for (size_t n = 0; n < dim3D; n++) {
		heapIndex[n] = -1;
	}

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	//Initialize the source points
	for (k = 0; k < height; k++) 
	{
		for (i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				size_t currentIndx = x_new(i, j, length);
				if(ctImageData.imageDataPtr[k][currentIndx] == foregroundValue)
				{

					if (k > 0)
					{
						if (labelArray[k - 1][currentIndx] != 1)
						{
							updateNeighbor3D(i, j, k - 1, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}	
					}
					if (k < (height - 1))
					{
						if (labelArray[k + 1][currentIndx] != 1)
						{
							updateNeighbor3D(i, j, k + 1, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}
					}
					if (i > 0)
					{
						if (labelArray[k][x_new(i - 1, j, length)] != 1)
						{
							updateNeighbor3D(i - 1, j, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}
					}
					if (i < (length - 1))
					{
						if (labelArray[k][x_new(i + 1, j, length)] != 1)
						{
							updateNeighbor3D(i + 1, j, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}
					}
					if (j > 0)
					{
						if (labelArray[k][x_new(i, j - 1, length)] != 1)
						{
							updateNeighbor3D(i, j - 1, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}
					}
					if (j < (width - 1))
					{
						if (labelArray[k][x_new(i, j + 1, length)] != 1)
						{
							updateNeighbor3D(i, j + 1, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
						}
					}
				}
			}
		}
	}

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching3D current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		size_t currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;

		deleteRootHeap3D(inProcess, heapIndex);

		if (k > 0)
		{
			if (labelArray[k - 1][currentIndx] != 1)
			{
				updateNeighbor3D(i, j, k - 1, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
			}
		}
		if (k < (height - 1))
		{
			if (labelArray[k + 1][currentIndx] != 1)
			{
				updateNeighbor3D(i, j, k + 1, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
			}
		}
		if (i > 0)
		{
			if (labelArray[k][x_new(i - 1, j, length)] != 1)
			{
				updateNeighbor3D(i - 1, j, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
			}
		}
		if (i < (length - 1))
		{
			if (labelArray[k][x_new(i + 1, j, length)] != 1)
			{
				updateNeighbor3D(i + 1, j, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
			}
		}
		if (j > 0)
		{
			updateNeighbor3D(i, j - 1, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
		}
		if (j < (width - 1))
		{
			if (labelArray[k][x_new(i, j + 1, length)] != 1)
			{
				updateNeighbor3D(i, j + 1, k, length, width, height, distanceFuncPtr, potentialFuncPtr, labelArray, spacing, inProcess, heapIndex);
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
		delete[] potentialFuncPtr[k];
	}
	delete[] labelArray;
	delete[] potentialFuncPtr;

	return true;
}

bool fastSweepingDistanceMap(Image_Data ctImageData, dataType** distancePtr, const dataType foregroundValue)
{
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t height = ctImageData.height;
	VoxelSpacing spacing = ctImageData.spacing;
	const size_t dim2D = length * width;

	size_t length_minus = length - 1;
	size_t width_minus = width - 1;
	size_t height_minus = height - 1;

	//Initialization
	for (size_t ik = 0; ik < height; ik++)
	{
		for (size_t ij = 0; ij < dim2D; ij++) 
		{
			if(ctImageData.imageDataPtr[ik][ij] == foregroundValue)
			{
				distancePtr[ik][ij] = 0.0;
			}
			else 
			{
				distancePtr[ik][ij] = INFINITY;
			}
		}
	}

	//sweep 1
	for (int k = 0; k < height; k++) 
	{
		for (int i = 0; i < length; i++) 
		{
			for (int j = 0; j < width; j++) 
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 2
	for (int k = height_minus; k > -1; k--)
	{
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 3
	for (int k = 0; k < height; k++)
	{
		for (int i = length_minus; i > -1; i--)
		{
			for (int j = 0; j < width; j++)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 4
	for (int k = 0; k < height; k++)
	{
		for (int i = 0; i < length; i++)
		{
			for (int j = width_minus; j > -1; j--)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 5
	for (int k = height_minus; k > -1; k--)
	{
		for (int i = length_minus; i > -1; i--)
		{
			for (int j = 0; j < width; j++)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 6
	for (int k = height_minus; k > -1; k--)
	{
		for (int i = 0; i < length; i++)
		{
			for (int j = width_minus; j > -1; j--)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 7
	for (int k = 0; k < height; k++)
	{
		for (int i = length_minus; i > -1; i--)
		{
			for (int j = width_minus; j > -1; j--)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}

	//sweep 8
	for (int k = height_minus; k > -1; k--)
	{
		for (int i = length_minus; i > -1; i--)
		{
			for (int j = width_minus; j > -1; j--)
			{
				int xd = x_new((size_t)i, (size_t)j, length);
				dataType x = upwindFiniteDifferenceX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = upwindFiniteDifferenceY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = upwindFiniteDifferenceZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType coefSpeed = 1.0;
				dataType pDistance = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (pDistance < distancePtr[k][xd])
				{
					distancePtr[k][xd] = pDistance;
				}
			}
		}
	}
	
	return true;
}

bool bruteForceDistanceMap(Image_Data ctImageData, dataType** distancePtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distancePtr == NULL) {
		return false;
	}
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t height = ctImageData.height;
	VoxelSpacing spacing = ctImageData.spacing;

	double min_distance = 0.0;
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t xd = x_new(i, j, length);
				Point3D cPoint = { i, j, k };
				cPoint = getRealCoordFromImageCoord3D(cPoint, ctImageData.origin, spacing, ctImageData.orientation);

				min_distance = INFINITY;
				for (size_t tk = 0; tk < height; tk++) 
				{
					for (size_t ti = 0; ti < length; ti++) 
					{
						for (size_t tj = 0; tj < width; tj++) 
						{
							size_t txd = x_new(ti, tj, length);
							if (ctImageData.imageDataPtr[tk][txd] == foregroundValue) 
							{
								Point3D tPoint = { ti, tj, tk };
								tPoint = getRealCoordFromImageCoord3D(tPoint, ctImageData.origin, spacing, ctImageData.orientation);
								double pDistance = getPoint3DDistance(cPoint, tPoint);
								if (pDistance < min_distance) {
									min_distance = pDistance;
								}
							}
						}
					}
				}
				if (min_distance == INFINITY) {
					distancePtr[k][xd] = 0.0;
				}
				else {
					distancePtr[k][xd] = min_distance;
				}
			}
		}
	}
	return true;
}

bool rouyTourinFrontPropagation(Image_Data ctImageData, dataType** distancePtr, dataType** potential, dataType tolerance, size_t max_iteration) {

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

	dataType** previousSolution = new dataType * [height_ext];
	for (k = 0; k < height_ext; k++) {
		previousSolution[k] = new dataType[length_ext * width_ext]{ 0 };
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
	std::cout << "tau = " << tau << std::endl;

	size_t count_iteration = 0;
	dataType w = 0.0;
	size_t xd;

	while (mass > tolerance && count_iteration < max_iteration) {
		copyDataToExtendedArea(distancePtr, previousSolution, height, length, width);
		reflection3D(previousSolution, height_ext, length_ext, width_ext);
		count_iteration++;
		mass = 0.0;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					xd = x_new(i_ext, j_ext, length_ext);
					value = previousSolution[k_ext][xd];
					w = potential[k][xd];
					distancePtr[k][x_new(i, j, length)] = value + w * tau - tau * sqrt(hx_2 * max(min0(previousSolution[k_ext][x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext + 1, j_ext, length_ext)], value))
						+ hy_2 * max(min0(previousSolution[k_ext][x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext, j_ext + 1, length_ext)], value))
						+ hz_2 * max(min0(previousSolution[k_ext - 1][x_new(i_ext, j_ext, length_ext)], value), min0(previousSolution[k_ext + 1][x_new(i_ext, j_ext, length_ext)], value)));
					//Compute the mass
					mass += pow(previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)] - distancePtr[k][x_new(i, j, length)], 2);
				}
			}
		}
		mass = sqrt(mass);
	}
	std::cout << "Iteration: " << count_iteration << ", Mass: " << mass << std::endl;

	for (k = 0; k < height_ext; k++) {
		delete[] previousSolution[k];
	}
	delete[] previousSolution;

	return true;
}


