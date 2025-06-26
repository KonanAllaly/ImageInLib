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

dataType selectX(dataType* actionPtr, const size_t length, const size_t width, size_t i, size_t j) {
	
	dataType i_minus, i_plus;

	if (i == 0) {
		i_minus = INFINITY;
	}
	else {
		i_minus = actionPtr[x_new(i - 1, j, length)];
	}

	if (i == length - 1) {
		i_plus = INFINITY;
	}
	else {
		i_plus = actionPtr[x_new(i + 1, j, length)];
	}

	return min(i_minus, i_plus);
}

dataType selectY(dataType* actionPtr, const size_t length, const size_t width, size_t i, size_t j) {
	
	dataType j_minus, j_plus;

	if (j == 0) {
		j_minus = INFINITY;
	}
	else {
		j_minus = actionPtr[x_new(i, j - 1, length)];
	}

	if (j == width - 1) {
		j_plus = INFINITY;
	}
	else {
		j_plus = actionPtr[x_new(i, j + 1, length)];
	}

	return min(j_minus, j_plus);
}

dataType solve2dQuadratic(dataType X, dataType Y, dataType P, PixelSpacing h) {

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
			if (solution >= Y) {
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
			if (solution >= X)
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
			if (solution >= max(X, Y)) {
				return solution;
			}
			else {
				return (dataType)(min(X + hx * P, Y + hy * P));
			}
		}
		else {
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
	dataType seedVal2 = (dataType)(imageDataStr.imageDataPtr[x_new((size_t)seedPoints[1].x, (size_t)seedPoints[1].y, length)]);
	dataType seedVal = (dataType)(seedVal1 + seedVal2) / 2.0; //Average of the two seed points

	dataType* distanceMap = new dataType[dim2D]{ 0 };
	dataType* edgeImage = new dataType[dim2D]{ 0 };	

	dataType norm_of_gradient = 0.0, edgeValue = 0.0;
	bool isGradientComputed = false;
	Point2D grad_vector;
	PixelSpacing fVolume = imageDataStr.spacing;

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
		potentialFuncPtr[i] = fabs(imageDataStr.imageDataPtr[i] - seedVal);
	}

	////Find max difference
	//dataType maxDiff = 0.0;
	//for (i = 0; i < dim2D; i++) {
	//	if (potentialFuncPtr[i] > maxDiff) {
	//		maxDiff = potentialFuncPtr[i];
	//	}
	//}
	
	//Normalization
	dataType weight = 0.0;
	for (i = 0; i < dim2D; i++) {
		//weight = 1.0 / (1.0 + 2.0 * distanceMap[i]);
		//potentialFuncPtr[i] = (parameters.eps + potentialFuncPtr[i] / maxDiff) * weight;
		potentialFuncPtr[i] = parameters.eps + potentialFuncPtr[i];
	}

	delete[] distanceMap;
	delete[] edgeImage;

	return true;
}

void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b) {
	pointFastMarching2D temp = *a;
	*a = *b;
	*b = temp;
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
		swap2dPoints(&in_Process[pos], &in_Process[current]);
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
		swap2dPoints(&in_Process[current], &in_Process[i]);
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
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap2dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown2D(in_Process, 0);
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

	size_t i = (size_t)seedPoints[0].x;
	size_t j = (size_t)seedPoints[0].y;
	size_t currentIndx = x_new(i, j, length);
	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < length) {
		size_t jplus = j + 1;
		dataType x = selectX(distancePtr, length, width, i, jplus);
		dataType y = selectY(distancePtr, length, width, i, jplus);
		size_t indxEast = x_new(i, jplus, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		distancePtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < length) {
		size_t jminus = j - 1;
		dataType x = selectX(distancePtr, length, width, i, jminus);
		dataType y = selectY(distancePtr, length, width, i, jminus);
		size_t indxWest = x_new(i, jminus, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		distancePtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		dataType x = selectX(distancePtr, length, width, iminus, j);
		dataType y = selectY(distancePtr, length, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		distancePtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < length - 1) {
		size_t iplus = i + 1;
		dataType x = selectX(distancePtr, length, width, iplus, j);
		dataType y = selectY(distancePtr, length, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		distancePtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;
	short label = 0;
	
	while (inProcess.size() > 0) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;
		deleteRootHeap2D(inProcess);

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, iminus, j);
				dataType y = selectY(distancePtr, length, width, iminus, j);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { iminus, j, dWest };
				if (label == 3) {
					distancePtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < distancePtr[indxWest]) {
						distancePtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(inProcess, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, iplus, j);
				dataType y = selectY(distancePtr, length, width, iplus, j);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { iplus, j, dEast };
				if (label == 3) {
					distancePtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < distancePtr[indxEast]) {
						distancePtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(inProcess, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, i, jminus);
				dataType y = selectY(distancePtr, length, width, i, jminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
				if (label == 3) {
					distancePtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < distancePtr[indxNorth]) {
						distancePtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, i, jplus);
				dataType y = selectY(distancePtr, length, width, i, jplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
				if (label == 3) {
					distancePtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < distancePtr[indxSouth]) {
						distancePtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

	}

	delete[] labelArray;
}

bool partialFrontPropagation2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* endPoints, string savingPath) {

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

	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < length) {
		size_t jplus = j + 1;
		dataType x = selectX(distancePtr, length, width, i, jplus);
		dataType y = selectY(distancePtr, length, width, i, jplus);
		size_t indxEast = x_new(i, jplus, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		distancePtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < length) {
		size_t jminus = j - 1;
		dataType x = selectX(distancePtr, length, width, i, jminus);
		dataType y = selectY(distancePtr, length, width, i, jminus);
		size_t indxWest = x_new(i, jminus, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		distancePtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		dataType x = selectX(distancePtr, length, width, iminus, j);
		dataType y = selectY(distancePtr, length, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		distancePtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < length - 1) {
		size_t iplus = i + 1;
		dataType x = selectX(distancePtr, length, width, iplus, j);
		dataType y = selectY(distancePtr, length, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		distancePtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;

	//Save points for visualization
	vector<Point2D> savingList;
	size_t id_save = 0;
	size_t nb_computed_points = 0;
	size_t seedI = endPoints[1].x, seedJ = endPoints[1].y, seedIndex = x_new(seedI, seedJ, length);

	dataType max_save_action = 0.0;
	//dataType* partialDistancePtr = new dataType[length * width]{ 0 };
	size_t nb_points_processed = 0;
	while (labelArray[seedIndex] != 1) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;
		nb_computed_points++;

		if(distancePtr[currentIndx] > max_save_action) {
			max_save_action = distancePtr[currentIndx];
		}

		////Save points for visualization
		//Point2D point = { (dataType)i, (dataType)j };
		//savingList.push_back(point);
		//if (nb_computed_points % 500 == 0) {
		//	id_save++;
		//	//string saving_csv = savingPath + to_string(id_save) + ".csv";
		//	//FILE* frontPoint;
		//	//if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		//	//	printf("Enable to open");
		//	//	return false;
		//	//}
		//	//fprintf(frontPoint, "x,y\n");
		//	//for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		//	//	fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
		//	//}
		//	////Don't empty the list of points to see the whole computed points
		//	////savingList.clear();
		//	//fclose(frontPoint);
		//	string saving_file = savingPath + to_string(id_save) + ".raw";
		//	for(size_t in = 0; in < savingList.size(); in++) {
		//		size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
		//		partialDistancePtr[indx] = distancePtr[indx];
		//	}
		//	manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);
		//}
		
		nb_points_processed++;
		deleteRootHeap2D(inProcess);

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, iminus, j);
				dataType y = selectY(distancePtr, length, width, iminus, j);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { iminus, j, dWest };
				if (label == 3) {
					distancePtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < distancePtr[indxWest]) {
						distancePtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(inProcess, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, iplus, j);
				dataType y = selectY(distancePtr, length, width, iplus, j);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { iplus, j, dEast };
				if (label == 3) {
					distancePtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < distancePtr[indxEast]) {
						distancePtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(inProcess, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, i, jminus);
				dataType y = selectY(distancePtr, length, width, i, jminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
				if (label == 3) {
					distancePtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < distancePtr[indxNorth]) {
						distancePtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(distancePtr, length, width, i, jplus);
				dataType y = selectY(distancePtr, length, width, i, jplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
				if (label == 3) {
					distancePtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < distancePtr[indxSouth]) {
						distancePtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}
	}
	
	////Save points for visualization
	//if (savingList.size() != 0) {
	//	id_save++;
	//	string saving_csv = savingPath + to_string(id_save) + ".csv";
	//	FILE* frontPoint;
	//	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
	//		printf("Enable to open");
	//		return false;
	//	}
	//	fprintf(frontPoint, "x,y\n");
	//	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
	//		fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
	//	}
	//	savingList.clear();
	//	fclose(frontPoint);
	//}

	//id_save++;
	//string saving_file = savingPath + to_string(id_save) + ".raw";
	//for (size_t in = 0; in < savingList.size(); in++) {
	//	size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
	//	partialDistancePtr[indx] = distancePtr[indx];
	//}
	//manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);

	nb_points_processed += inProcess.size();
	std::cout << "Number of points processed: " << nb_points_processed << std::endl;

	//Set the distance of the end point to the maximum value
	for(i = 0; i < dim2D; i++) {
		if (labelArray[i] == 3) {
			distancePtr[i] = max_save_action + 1;
		}
	}
	inProcess.clear();
	
	//delete[] partialDistancePtr;
	delete[] labelArray;
}

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

	if(x1 > 0 && y1 >= 0 && y1 < width) {
		size_t xminus = x1 - 1;
		dataType x = selectX(actionFirstFront, length, width, xminus, y1);
		dataType y = selectY(actionFirstFront, length, width, xminus, y1);
		size_t indxWest = x_new(xminus, y1, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { xminus, y1, dWest };
		actionFirstFront[indxWest] = dWest;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[indxWest] = 2; //2 ---> in process
	}

	if(x1 < length_minus && y1 >= 0 && y1 < width) {
		size_t xplus = x1 + 1;
		dataType x = selectX(actionFirstFront, length, width, xplus, y1);
		dataType y = selectY(actionFirstFront, length, width, xplus, y1);
		size_t indxEast = x_new(xplus, y1, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { xplus, y1, dEast };
		actionFirstFront[indxEast] = dEast;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[indxEast] = 2; //2 ---> in process
	}

	if (y1 > 0 && x1 >= 0 && x1 < length) {
		size_t yminus = y1 - 1;
		dataType x = selectX(actionFirstFront, length, width, x1, yminus);
		dataType y = selectY(actionFirstFront, length, width, x1, yminus);
		size_t indxNorth = x_new(x1, yminus, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { x1, yminus, dNorth };
		actionFirstFront[indxNorth] = dNorth;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[indxNorth] = 2; //2 ---> in process
	}

	if (y1 < width_minus && x1 >= 0 && x1 < length) {
		size_t yplus = y1 + 1;
		dataType x = selectX(actionFirstFront, length, width, x1, yplus);
		dataType y = selectY(actionFirstFront, length, width, x1, yplus);
		size_t indxSouth = x_new(x1, yplus, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { x1, yplus, dSouth };
		actionFirstFront[indxSouth] = dSouth;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[indxSouth] = 2; //2 ---> in process
	}

	//Initialize neighbors for the second front

	if (x2 > 0 && y2 >= 0 && y2 < width) {
		size_t xminus = x2 - 1;
		dataType x = selectX(actionSecondFront, length, width, xminus, y2);
		dataType y = selectY(actionSecondFront, length, width, xminus, y2);
		size_t indxWest = x_new(xminus, y2, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { xminus, y2, dWest };
		actionSecondFront[indxWest] = dWest;
		narrowBandSecondFront.push_back(WestNeighbor);
		secondLabelArray[indxWest] = 2; //2 ---> in process
	}

	if(x2 < length_minus && y2 >= 0 && y2 < width) {
		size_t xplus = x2 + 1;
		dataType x = selectX(actionSecondFront, length, width, xplus, y2);
		dataType y = selectY(actionSecondFront, length, width, xplus, y2);
		size_t indxEast = x_new(xplus, y2, length);
		dataType coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { xplus, y2, dEast };
		actionSecondFront[indxEast] = dEast;
		narrowBandSecondFront.push_back(EastNeighbor);
		secondLabelArray[indxEast] = 2; //2 ---> in process
	}

	if(y2 > 0 && x2 >= 0 && x2 < length) {
		size_t yminus = y2 - 1;
		dataType x = selectX(actionSecondFront, length, width, x2, yminus);
		dataType y = selectY(actionSecondFront, length, width, x2, yminus);
		size_t indxNorth = x_new(x2, yminus, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { x2, yminus, dNorth };
		actionSecondFront[indxNorth] = dNorth;
		narrowBandSecondFront.push_back(NorthNeighbor);
		secondLabelArray[indxNorth] = 2; //2 ---> in process
	}

	if(y2 < width_minus && x2 >= 0 && x2 < length) {
		size_t yplus = y2 + 1;
		dataType x = selectX(actionSecondFront, length, width, x2, yplus);
		dataType y = selectY(actionSecondFront, length, width, x2, yplus);
		size_t indxSouth = x_new(x2, yplus, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { x2, yplus, dSouth };
		actionSecondFront[indxSouth] = dSouth;
		narrowBandSecondFront.push_back(SouthNeighbor);
		secondLabelArray[indxSouth] = 2; //2 ---> in process
	}

	//heapify 2D vector
	heapifyVector2D(narrowBandFirstFront);
	heapifyVector2D(narrowBandSecondFront);

	//Save points for visualization
	//vector<Point2D> savingList;
	//size_t id_save = 0;

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
			//narrowBandFirstFront.clear();
			//narrowBandSecondFront.clear();
			//return true; //End the propagation
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
			//narrowBandFirstFront.clear();
			//narrowBandSecondFront.clear();
			//return true; //End the propagation
		}
		else {
			secondLabelArray[indexSecond] = 1;
		}

		if (secondFrontPoint.arrival > max_save_action) {
			max_save_action = secondFrontPoint.arrival;
		}
		deleteRootHeap2D(narrowBandSecondFront);
		nb_computed_points++;

		////Save points for visualization
		//Point2D point = { (dataType)i, (dataType)j };
		//savingList.push_back(point);
		//if (nb_computed_points % 500 == 0) {
		//	id_save++;
		//	//string saving_csv = savingPath + to_string(id_save) + ".csv";
		//	//FILE* frontPoint;
		//	//if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		//	//	printf("Enable to open");
		//	//	return false;
		//	//}
		//	//fprintf(frontPoint, "x,y\n");
		//	//for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		//	//	fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
		//	//}
		//	////Don't empty the list of points to see the whole computed points
		//	////savingList.clear();
		//	//fclose(frontPoint);
		//	string saving_file = savingPath + to_string(id_save) + ".raw";
		//	for(size_t in = 0; in < savingList.size(); in++) {
		//		size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
		//		partialDistancePtr[indx] = distancePtr[indx];
		//	}
		//	manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);
		//}

		//West first front
		if (x1 > 0 && x1 < length && y1 >= 0 && y1 < width) {
			size_t xminus = x1 - 1;
			size_t indxWest = x_new(xminus, y1, length);
			short label = firstLabelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, xminus, y1);
				dataType y = selectY(actionFirstFront, length, width, xminus, y1);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { xminus, y1, dWest };
				if (label == 3) {
					actionFirstFront[indxWest] = dWest;
					firstLabelArray[indxWest] = 2;
					addPointHeap2D(narrowBandFirstFront, WestNeighbor);
				}
				else {
					if (dWest < actionFirstFront[indxWest]) {
						actionFirstFront[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, xminus, y1);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//West second front
		if(x2 > 0 && x2 < length && y2 >= 0 && y2 < width) {
			size_t xminus = x2 - 1;
			size_t indxWest = x_new(xminus, y2, length);
			short label = secondLabelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, xminus, y2);
				dataType y = selectY(actionSecondFront, length, width, xminus, y2);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { xminus, y2, dWest };
				if (label == 3) {
					actionSecondFront[indxWest] = dWest;
					secondLabelArray[indxWest] = 2;
					addPointHeap2D(narrowBandSecondFront, WestNeighbor);
				}
				else {
					if (dWest < actionSecondFront[indxWest]) {
						actionSecondFront[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, xminus, y2);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//East first front
		if (x1 < length_minus && x1 >= 0 && y1 >= 0 && y1 < width) {
			size_t xplus = x1 + 1;
			size_t indxEast = x_new(xplus, y1, length);
			short label = firstLabelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, xplus, y1);
				dataType y = selectY(actionFirstFront, length, width, xplus, y1);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { xplus, y1, dEast };
				if (label == 3) {
					actionFirstFront[indxEast] = dEast;
					firstLabelArray[indxEast] = 2;
					addPointHeap2D(narrowBandFirstFront, EastNeighbor);
				}
				else {
					if (dEast < actionFirstFront[indxEast]) {
						actionFirstFront[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, xplus, y1);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//East second front
		if (x2 < length_minus && x2 >= 0 && y2 >= 0 && y2 < width) {
			size_t xplus = x2 + 1;
			size_t indxEast = x_new(xplus, y2, length);
			short label = secondLabelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, xplus, y2);
				dataType y = selectY(actionSecondFront, length, width, xplus, y2);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { xplus, y2, dEast };
				if (label == 3) {
					actionSecondFront[indxEast] = dEast;
					secondLabelArray[indxEast] = 2;
					addPointHeap2D(narrowBandSecondFront, EastNeighbor);
				}
				else {
					if (dEast < actionSecondFront[indxEast]) {
						actionSecondFront[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, xplus, y2);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//North first front
		if(y1 > 0 && y1 < width && x1 >= 0 && x1 < length) {
			size_t yminus = y1 - 1;
			size_t indxNorth = x_new(x1, yminus, length);
			short label = firstLabelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, x1, yminus);
				dataType y = selectY(actionFirstFront, length, width, x1, yminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { x1, yminus, dNorth };
				if (label == 3) {
					actionFirstFront[indxNorth] = dNorth;
					firstLabelArray[indxNorth] = 2;
					addPointHeap2D(narrowBandFirstFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionFirstFront[indxNorth]) {
						actionFirstFront[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, x1, yminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//North second front
		if (y2 > 0 && y2 < width && x2 >= 0 && x2 < length) {
			size_t yminus = y2 - 1;
			size_t indxNorth = x_new(x2, yminus, length);
			short label = secondLabelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, x2, yminus);
				dataType y = selectY(actionSecondFront, length, width, x2, yminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { x2, yminus, dNorth };
				if (label == 3) {
					actionSecondFront[indxNorth] = dNorth;
					secondLabelArray[indxNorth] = 2;
					addPointHeap2D(narrowBandSecondFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionSecondFront[indxNorth]) {
						actionSecondFront[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, x2, yminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//South first front
		if(y1 >= 0 && y1 < width_minus && x1 >= 0 && x1 < length) {
			size_t yplus = y1 + 1;
			size_t indxSouth = x_new(x1, yplus, length);
			short label = firstLabelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionFirstFront, length, width, x1, yplus);
				dataType y = selectY(actionFirstFront, length, width, x1, yplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { x1, yplus, dSouth };
				if (label == 3) {
					actionFirstFront[indxSouth] = dSouth;
					firstLabelArray[indxSouth] = 2;
					addPointHeap2D(narrowBandFirstFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionFirstFront[indxSouth]) {
						actionFirstFront[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBandFirstFront, x1, yplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//South second front
		if(y2 >= 0 && y2 < width_minus && x2 >= 0 && x2 < length) {
			size_t yplus = y2 + 1;
			size_t indxSouth = x_new(x2, yplus, length);
			short label = secondLabelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionSecondFront, length, width, x2, yplus);
				dataType y = selectY(actionSecondFront, length, width, x2, yplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { x2, yplus, dSouth };
				if (label == 3) {
					actionSecondFront[indxSouth] = dSouth;
					secondLabelArray[indxSouth] = 2;
					addPointHeap2D(narrowBandSecondFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionSecondFront[indxSouth]) {
						actionSecondFront[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBandSecondFront, x2, yplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//if (nb_computed_points == 50000) {
		//	break;
		//}

	}

	////Save points for visualization
	//if (savingList.size() != 0) {
	//	id_save++;
	//	string saving_csv = savingPath + to_string(id_save) + ".csv";
	//	FILE* frontPoint;
	//	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
	//		printf("Enable to open");
	//		return false;
	//	}
	//	fprintf(frontPoint, "x,y\n");
	//	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
	//		fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
	//	}
	//	savingList.clear();
	//	fclose(frontPoint);
	//}

	//id_save++;
	//string saving_file = savingPath + to_string(id_save) + ".raw";
	//for (size_t in = 0; in < savingList.size(); in++) {
	//	size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
	//	partialDistancePtr[indx] = distancePtr[indx];
	//}
	//manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);

	nb_computed_points += narrowBandFirstFront.size() + narrowBandSecondFront.size();
	std::cout << "Number of points processed: " << nb_computed_points << std::endl;

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

bool frontPropagationWithKeyPointDetection(Image_Data2D actionMapStr, dataType* potentialFuncPtr, Point2D* seedPoint, const double LengthKeyPoints, vector<Point2D>& key_points, std::string path_saving) {

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

	//West
	if (i > 0 && i < length && j >= 0 && j < width) {
		size_t iminus = i - 1;
		size_t indxWest = x_new(iminus, j, length);
		dataType x = selectX(actionMapStr.imageDataPtr, length, width, iminus, j);
		dataType y = selectY(actionMapStr.imageDataPtr, length, width, iminus, j);
		dataType coefSpeed = potentialFuncPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { iminus, j, dWest };
		actionMapStr.imageDataPtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//East
	if (i < length_minus && i >= 0 && j >= 0 && j < width) {
		size_t iplus = i + 1;
		size_t indxEast = x_new(iplus, j, length);
		dataType x = selectX(actionMapStr.imageDataPtr, length, width, iplus, j);
		dataType y = selectY(actionMapStr.imageDataPtr, length, width, iplus, j);
		dataType coefSpeed = potentialFuncPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D EastNeighbor = { iplus, j, dEast };
		actionMapStr.imageDataPtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//North
	if (j > 0 && j < width && i >= 0 && i < length) {
		size_t jminus = j - 1;
		size_t indxNorth = x_new(i, jminus, length);
		dataType x = selectX(actionMapStr.imageDataPtr, length, width, i, jminus);
		dataType y = selectY(actionMapStr.imageDataPtr, length, width, i, jminus);
		dataType coefSpeed = potentialFuncPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
		actionMapStr.imageDataPtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j>= 0 && j < width_minus && i >= 0 && i < length) {
		size_t jplus = j + 1;
		size_t indxSouth = x_new(i, jplus, length);
		dataType x = selectX(actionMapStr.imageDataPtr, length, width, i, jplus);
		dataType y = selectY(actionMapStr.imageDataPtr, length, width, i, jplus);
		dataType coefSpeed = potentialFuncPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
		actionMapStr.imageDataPtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	heapifyVector2D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;
	size_t nbSourcePoint = 0;
	dataType max_weighted_distance = 0.0;

	size_t iEnd = (size_t)seedPoint[1].x;
	size_t jEnd = (size_t)seedPoint[1].y;

	//Visualize the front propagation
	size_t id_keyPoint = 1;
	std::string storing_path;
	vector<Point2D> savingList;

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
		//For vizualization
		savingList.push_back(pSourceReal);
		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			string saving_csv = path_saving + to_string(id_keyPoint) + ".csv";
			FILE* frontPoint;
			if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
				printf("Enable to open");
				return false;
			}
			fprintf(frontPoint, "x,y\n");
			for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
				fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
			}
			//savingList.clear(); //not empty the list of computed points to see the whole front
			fclose(frontPoint);

			//If the condition is true ---> new key point is found so we need to initilize it neighbors
			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionMapStr.imageDataPtr[currentIndx] = 0;

			//Initialize all the points inside the narrow band as not processed ---> label = 3
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
		
		//processed neighbors of the minimum in the narrow band

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			size_t label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionMapStr.imageDataPtr, length, width, iminus, j);
				dataType y = selectY(actionMapStr.imageDataPtr, length, width, iminus, j);
				dataType coefSpeed = potentialFuncPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxWest] = 2;
					actionMapStr.imageDataPtr[indxWest] = dWest;
					pointFastMarching2D WestNeighbor = { iminus, j, dWest };
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionMapStr.imageDataPtr[indxWest]) {
						actionMapStr.imageDataPtr[indxWest] = dWest;
						size_t pt_pos = getIndexFromHeap2D(inProcess, iminus, j);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dWest;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}
		
		//East
		if (i < length_minus && i >= 0 && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			size_t label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionMapStr.imageDataPtr, length, width, iplus, j);
				dataType y = selectY(actionMapStr.imageDataPtr, length, width, iplus, j);
				dataType coefSpeed = potentialFuncPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxEast] = 2;
					actionMapStr.imageDataPtr[indxEast] = dEast;
					pointFastMarching2D EastNeighbor = { iplus, j, dEast };
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < actionMapStr.imageDataPtr[indxEast]) {
						actionMapStr.imageDataPtr[indxEast] = dEast;
						size_t pt_pos = getIndexFromHeap2D(inProcess, iplus, j);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dEast;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			size_t label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionMapStr.imageDataPtr, length, width, i, jminus);
				dataType y = selectY(actionMapStr.imageDataPtr, length, width, i, jminus);
				dataType coefSpeed = potentialFuncPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxNorth] = 2;
					actionMapStr.imageDataPtr[indxNorth] = dNorth;
					pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionMapStr.imageDataPtr[indxNorth]) {
						actionMapStr.imageDataPtr[indxNorth] = dNorth;
						size_t pt_pos = getIndexFromHeap2D(inProcess, i, jminus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dNorth;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			size_t label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionMapStr.imageDataPtr, length, width, i, jplus);
				dataType y = selectY(actionMapStr.imageDataPtr, length, width, i, jplus);
				dataType coefSpeed = potentialFuncPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				if (label == 3) {
					labelArray[indxSouth] = 2;
					actionMapStr.imageDataPtr[indxSouth] = dSouth;
					pointFastMarching2D NorthNeighbor = { i, jplus, dSouth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dSouth < actionMapStr.imageDataPtr[indxSouth]) {
						actionMapStr.imageDataPtr[indxSouth] = dSouth;
						size_t pt_pos = getIndexFromHeap2D(inProcess, i, jplus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dSouth;
							heapifyUp2D(inProcess, pt_pos);
						}
					}
				}
			}
		}

	}

	key_points.push_back(seedPoint[1]);
	std::cout << "The source Points are: " << std::endl;
	for (int nb = 0; nb < key_points.size(); nb++) {
		std::cout << "Point " << nb + 1 << " : " << key_points[nb].x << " " << key_points[nb].y << std::endl;
	}

	//id_keyPoint++;
	string saving_csv = path_saving + to_string(id_keyPoint) + ".csv";
	FILE* frontPoint;
	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(frontPoint, "x,y\n");
	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
	}
	savingList.clear();
	fclose(frontPoint);

	delete[] labelArray;

	return true;
}

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

				//West
				if (i > 0 && i < length && j >= 0 && j < width) {
					size_t iminus = i - 1;
					size_t indxWest = x_new(iminus, j, length);
					short label = labelArray[indxWest];
					if (label != 1) {
						x = selectX(distanceFuncPtr, length, width, iminus, j);
						y = selectY(distanceFuncPtr, length, width, iminus, j);
						coefSpeed = 1.0;
						dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
						pointFastMarching2D WestNeighbor = { iminus, j, dWest};
						if (label == 3) {
							distanceFuncPtr[indxWest] = dWest;
							labelArray[indxWest] = 2;
							addPointHeap2D(inProcess, WestNeighbor);
						}
						else {
							if (dWest < distanceFuncPtr[indxWest]) {
								distanceFuncPtr[indxWest] = dWest;
								size_t pIndex = getIndexFromHeap2D(inProcess, iminus, j);
								if (pIndex != -1) {
									heapifyUp2D(inProcess, pIndex);
								}
							}
						}
					}
				}

				//East
				if (i >= 0 && i < length_minus && j >= 0 && j < width) {
					size_t iplus = i + 1;
					size_t indxEast = x_new(iplus, j, length);
					short label = labelArray[indxEast];
					if (label != 1) {
						x = selectX(distanceFuncPtr, length, width, iplus, j);
						y = selectY(distanceFuncPtr, length, width, iplus, j);
						coefSpeed = 1.0;
						dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
						pointFastMarching2D EastNeighbor = { iplus, j, dEast};
						if (label == 3) {
							distanceFuncPtr[indxEast] = dEast;
							labelArray[indxEast] = 2;
							addPointHeap2D(inProcess, EastNeighbor);
						}
						else {
							if (dEast < distanceFuncPtr[indxEast]) {
								distanceFuncPtr[indxEast] = dEast;
								size_t pIndex = getIndexFromHeap2D(inProcess, iplus, j);
								if (pIndex != -1) {
									heapifyUp2D(inProcess, pIndex);
								}
							}
						}
					}
				}

				//North
				if (j > 0 && j < width && i >= 0 && i < length) {
					size_t jminus = j - 1;
					size_t indxNorth = x_new(i, jminus, length);
					short label = labelArray[indxNorth];
					if (label != 1) {
						x = selectX(distanceFuncPtr, length, width, i, jminus);
						y = selectY(distanceFuncPtr, length, width, i, jminus);
						coefSpeed = 1.0;
						dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
						pointFastMarching2D NorthNeighbor = { i, jminus, dNorth};
						if (label == 3) {
							distanceFuncPtr[indxNorth] = dNorth;
							labelArray[indxNorth] = 2;
							addPointHeap2D(inProcess, NorthNeighbor);
						}
						else {
							if (dNorth < distanceFuncPtr[indxNorth]) {
								distanceFuncPtr[indxNorth] = dNorth;
								size_t pIndex = getIndexFromHeap2D(inProcess, i, jminus);
								if (pIndex != -1) {
									heapifyUp2D(inProcess, pIndex);
								}
							}
						}
					}
				}

				//South
				if (j >= 0 && j < width_minus && i >= 0 && i < length) {
					size_t jplus = j + 1;
					size_t indxSouth = x_new(i, jplus, length);
					short label = labelArray[indxSouth];
					if (label != 1) {
						x = selectX(distanceFuncPtr, length, width, i, jplus);
						y = selectY(distanceFuncPtr, length, width, i, jplus);
						coefSpeed = 1.0;
						dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
						pointFastMarching2D SouthNeighbor = { i, jplus, dSouth};
						if (label == 3) {
							distanceFuncPtr[indxSouth] = dSouth;
							labelArray[indxSouth] = 2;
							addPointHeap2D(inProcess, SouthNeighbor);
						}
						else {
							if (dSouth < distanceFuncPtr[indxSouth]) {
								distanceFuncPtr[indxSouth] = dSouth;
								size_t pIndex = getIndexFromHeap2D(inProcess, i, jplus);
								if (pIndex != -1) {
									heapifyUp2D(inProcess, pIndex);
								}
							}
						}
					}
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

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				x = selectX(distanceFuncPtr, length, width, iminus, j);
				y = selectY(distanceFuncPtr, length, width, iminus, j);
				coefSpeed = 1.0;
				dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { iminus, j, dWest};
				if (label == 3) {
					distanceFuncPtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < distanceFuncPtr[indxWest]) {
						distanceFuncPtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(inProcess, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				x = selectX(distanceFuncPtr, length, width, iplus, j);
				y = selectY(distanceFuncPtr, length, width, iplus, j);
				coefSpeed = 1.0;
				dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { iplus, j, dEast };
				if (label == 3) {
					distanceFuncPtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < distanceFuncPtr[indxEast]) {
						distanceFuncPtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(inProcess, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				x = selectX(distanceFuncPtr, length, width, i, jminus);
				y = selectY(distanceFuncPtr, length, width, i, jminus);
				coefSpeed = 1.0;
				dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
				if (label == 3) {
					distanceFuncPtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < distanceFuncPtr[indxNorth]) {
						distanceFuncPtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				x = selectX(distanceFuncPtr, length, width, i, jplus);
				y = selectY(distanceFuncPtr, length, width, i, jplus);
				coefSpeed = 1.0;
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
				if (label == 3) {
					distanceFuncPtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < distanceFuncPtr[indxSouth]) {
						distanceFuncPtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(inProcess, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(inProcess, pIndex);
						}
					}
				}
			}
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//=================== Functions for the 3D Fast Marching and Path Tracking ==============================
/////////////////////////////////////////////////////////////////////////////////////////////////////////

dataType select3dX(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {

	dataType x_minus, x_plus;
	if(i == 0) {
		x_minus = INFINITY;
	}
	else {
		x_minus = actionPtr[k][x_new(i - 1, j, length)];
	}
	if(i == length - 1) {
		x_plus = INFINITY;
	}
	else {
		x_plus = actionPtr[k][x_new(i + 1, j, length)];
	}
	return min(x_minus, x_plus);
}

dataType select3dY(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {

	dataType y_minus, y_plus;
	if(j == 0) {
		y_minus = INFINITY;
	}
	else {
		y_minus = actionPtr[k][x_new(i, j - 1, length)];
	}
	if(j == width - 1) {
		y_plus = INFINITY;
	}
	else {
		y_plus = actionPtr[k][x_new(i, j + 1, length)];
	}
	return min(y_minus, y_plus);
}

dataType select3dZ(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {
	
	size_t xd = x_new(i, j, length);
	dataType z_minus, z_plus;
	if(k == 0) {
		z_minus = INFINITY;
	}
	else {
		z_minus = actionPtr[k - 1][xd];
	}
	if(k == height - 1) {
		z_plus = INFINITY;
	}
	else {
		z_plus = actionPtr[k + 1][xd];
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
			if (solution >= X) {
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
			if (solution >= Y) {
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
			if (solution >= Z) {
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
			if (solution >= max(X, Y)) {
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
		if (delta > 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= max(X, Z)) {
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
			if (solution >= max(Y, Z)) {
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
			if (solution >= max(X, max(Y, Z))) {
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

void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i) {

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
		swap3dPoints(&in_Process[i], &in_Process[current]);
		heapifyDown3D(in_Process, current);
	}

}

void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i) {

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
		swap3dPoints(&in_Process[current], &in_Process[i]);
		heapifyUp3D(in_Process, current);
	}

}

void heapifyVector3D(vector<pointFastMarching3D>& in_Process) {
	int length_array = in_Process.size();
	if(length_array < 2) {
		return; //nothing to heapify
	}
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		heapifyDown3D(in_Process, ind);
	}
}

void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	if (l > 1) {
		swap3dPoints(&in_Process[0], &in_Process[l - 1]);
		in_Process.pop_back();
		heapifyDown3D(in_Process, 0);
	}
	else {
		in_Process.pop_back();
	}
}

void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	if (l < 2) 
	{
		return; //nothing to heapify
	}
	else {
		heapifyUp3D(in_Process, l - 1);
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
	//
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			xd = x_new(i, j, length);
	//			isGradientComputed = getGradient3D(ctImageData, i, j, k, &grad_vector);
	//			if (isGradientComputed == true) {
	//				norm_of_gradient = sqrt(grad_vector.x * grad_vector.x + grad_vector.y * grad_vector.y + grad_vector.z * grad_vector.z);
	//			}
	//			else {
	//				std::cout << "Error in computing gradient at point (" << i << ", " << j << ", " << k << ")" << std::endl;
	//				return false;
	//			}
	//			dataType edgeValue = gradientFunction(norm_of_gradient, parameters.K);
	//			//threshold : real image
	//			if (edgeValue <= parameters.thres) {
	//				maskThreshold[k][xd] = 1.0;
	//			}
	//			else {
	//				maskThreshold[k][xd] = 0.0;
	//			}
	//		}
	//	}
	//}

	//////////Real image
	Image_Data toDistanceMap = { height, length, width, maskThreshold, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	//fastSweepingDistanceMap(toDistanceMap, distance, 1.0);
	std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/p6/distance_fs.raw";
	manageRAWFile3D<dataType>(distance, length, width, height, storing_path.c_str(), LOAD_DATA, false);
	//storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/p6/edge_image.raw";
	//manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), STORE_DATA, false);

	////Find min and max of the distance map
	//dataType max_distance = 0.0;
	//dataType min_distance = INFINITY;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (distance[k][i] > max_distance) {
	//			max_distance = distance[k][i];
	//		}
	//		if (distance[k][i] < min_distance) {
	//			min_distance = distance[k][i];
	//		}
	//	}
	//}
	//rescaleNewRange(distance, length, width, height, 0.0, 1.0, min_distance, max_distance);

	////////Artificial image : no need to compute the edge image when empty inside
	//Image_Data toDistanceMap = { height, length, width, ctImageData.imageDataPtr, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	////fastMarching3dForDistanceMap(toDistanceMap, distance, 1.0);
	////rouyTourinDistanceMap(toDistanceMap, distance, 0.001, 1000, 1.0);
	//fastSweepingDistanceMap(toDistanceMap, distance, 1.0);
	//////bruteForceDistanceMap(toDistanceMap, distance, 1.0);
	//std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_FS.raw";
	////////std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_fstswp.raw";
	////////std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_RT.raw";
	//////std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_brtforce.raw";
	//manageRAWFile3D<dataType>(distance, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };
	seedStats = getPointNeighborhoodStats(ctImageData, seedPoint[0], parameters.radius);
	dataType value_first_pt = seedStats.mean_data;
	////seedStats = getPointNeighborhoodStats(ctImageData, seedPoint[1], parameters.radius);
	////dataType value_second_pt = seedStats.mean_data;

	//std::cout << "#############################################################" << std::endl;
	//std::cout << "Statistics on input image " << std::endl;
	//std::cout << "Seed point mean value: " << seedStats.mean_data << std::endl;
	//std::cout << "Seed point minimum value: " << seedStats.min_data << std::endl;
	//std::cout << "Seed point maximum value: " << seedStats.max_data << std::endl;
	//std::cout << "Standard deviation value: " << seedStats.sd_data << std::endl;
	//std::cout << "#############################################################" << std::endl;
	
	dataType var_epsilon = 0;//0.01;
	if (seedStats.sd_data != 0) {
		var_epsilon = seedStats.sd_data;
	}
	else {
		var_epsilon = parameters.eps;
	}
	//std::cout << "Epsilon to be used : " << var_epsilon << std::endl;
	
	dataType seedValCT = value_first_pt;
	//////dataType seedValCT = (value_first_pt + value_second_pt) / 2.0;
	//////dataType seedValCT = ctImageData.imageDataPtr[(size_t)seedPoint[0].z][x_new((size_t)seedPoint[0].x, (size_t)seedPoint[0].y, length)];

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = fabs(seedValCT - ctImageData.imageDataPtr[k][i]);
		}
	}
	
	//Image_Data toStatistics = { height, length, width, potential, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	//seedStats = getPointNeighborhoodStats(toStatistics, seedPoint[0], parameters.radius);
	//std::cout << "Statistics after difference " << std::endl;
	//std::cout << "Seed point mean value: " << seedStats.mean_data << std::endl;
	//std::cout << "Seed point minimum value: " << seedStats.min_data << std::endl;
	//std::cout << "Seed point maximum value: " << seedStats.max_data << std::endl;
	//std::cout << "Standard deviation value: " << seedStats.sd_data << std::endl;
	//std::cout << "#############################################################" << std::endl;

	//Find the max of the difference
	dataType maxImage = 0.0;
	dataType minImage = INFINITY;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (potential[k][i] > maxImage) {
				maxImage = potential[k][i];
			}
			if (potential[k][i] < minImage) {
				minImage = potential[k][i];
			}
		}
	}
	////rescaleNewRange(potential, length, width, height, 0.0, 1.0, minImage, maxImage);
	//std::cout << "#############################################################" << std::endl;
	//std::cout << "difference min : " << minImage << std::endl;
	//std::cout << "difference max : " << maxImage << std::endl;
	//std::cout << "#############################################################" << std::endl;

	//dataType maxRatio = 0.0;
	//dataType minRatio = INFINITY;
	//Normalization
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			//potential[k][i] = var_epsilon + potential[k][i];
			potential[k][i] = (var_epsilon + potential[k][i] / maxImage) * (1.0 / (1.0 + 1.0 * distance[k][i]));
		}
	}
	//std::cout << "min factor : " << minRatio << std::endl;
	//std::cout << "max factor : " << maxRatio << std::endl;

	//seedStats = getPointNeighborhoodStats(toStatistics, seedPoint[0], parameters.radius);
	//std::cout << "Statistics after dividing " << std::endl;
	//std::cout << "Seed point mean value: " << seedStats.mean_data << std::endl;
	//std::cout << "Seed point minimum value: " << seedStats.min_data << std::endl;
	//std::cout << "Seed point maximum value: " << seedStats.max_data << std::endl;
	//std::cout << "Standard deviation value: " << seedStats.sd_data << std::endl;
	
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

		if(x < 0.0 || x >= length || y < 0.0 || y >= width || z < 0.0 || z >= height) {
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

bool partialFrontPropagation(Image_Data actionPtr, dataType** potentialFuncPtr, Point3D* endPoints, std::string path_saving) {

	if (actionPtr.imageDataPtr == NULL || potentialFuncPtr == NULL || endPoints == NULL) {
		return false;
	}

	const size_t height = actionPtr.height;
	const size_t length = actionPtr.length;
	const size_t width = actionPtr.width;
	VoxelSpacing spacing = actionPtr.spacing;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

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
			actionPtr.imageDataPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the starting point
	if (endPoints[0].x < 0 || endPoints[0].x > length || endPoints[0].y < 0 || endPoints[0].y > width || endPoints[0].z < 0 || endPoints[0].z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	i = (size_t)endPoints[0].x;
	j = (size_t)endPoints[0].y;
	k = (size_t)endPoints[0].z;
	size_t currentIndx = x_new(i, j, length);
	actionPtr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	
	//Top
	if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
		size_t kminus = k - 1;
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}
	
	//Bottom
	if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
		size_t kplus = k + 1;
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	//North
	if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
		size_t jminus = j - 1;
		size_t indxNorth = x_new(i, jminus, length);
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxNorth];
		dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr.imageDataPtr[k][indxNorth] = dNorth;
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
		size_t jplus = j + 1;
		size_t indxSouth = x_new(i, jplus, length);
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxSouth];
		dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		size_t posSouth = inProcess.size();
		actionPtr.imageDataPtr[k][indxSouth] = dSouth;
		pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth };
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//East
	if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
		size_t iplus = i + 1;
		size_t indxEast = x_new(iplus, j, length);
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxEast];
		dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr.imageDataPtr[k][indxEast] = dEast;
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//West
	if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
		size_t iminus = i - 1;
		size_t indxWest = x_new(iminus, j, length);
		dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxWest];
		dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr.imageDataPtr[k][indxWest] = dWest;
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		inProcess.push_back(WestNeighbor);
	}
	
	heapifyVector3D(inProcess);

	size_t id_save = 0;
	vector <Point3D> savingList;

	if (endPoints[1].x < 0 || endPoints[1].x > length || endPoints[1].y < 0 || endPoints[1].y > width || endPoints[1].z < 0 || endPoints[1].z > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	size_t seedI = (size_t)endPoints[1].x;
	size_t seedJ = (size_t)endPoints[1].y;
	size_t seedK = (size_t)endPoints[1].z;
	size_t seedIndex = x_new(seedI, seedJ, length);

	dataType max_action = 0;
	size_t processed_point = 0;
	size_t count_mistakes = 0;
	size_t cout_laborus_update = 0;
	pointFastMarching3D current;

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		processed_point++;

		//Find the maximum action value
		if (actionPtr.imageDataPtr[k][currentIndx] > max_action) {
			max_action = actionPtr.imageDataPtr[k][currentIndx];
		}
		
		//Exit the while loop when the final point is reached/computed
		if (labelArray[seedK][seedIndex] == 1) {
			break;
		}

		////Visualize the front propagation
		//Point3D sPoint = { (dataType)i, (dataType)j, (dataType)k };
		//sPoint = getRealCoordFromImageCoord3D(sPoint, actionPtr.origin, actionPtr.spacing, actionPtr.orientation);
		//savingList.push_back(sPoint);
		////in real image save after every 10000 points
		//if (processed_point % 10000 == 0) {
		//	id_save++;
		//	string saving_csv = path_saving + to_string(id_save) + ".csv";
		//	FILE* frontPoint;
		//	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		//		printf("Enable to open");
		//		return false;
		//	}
		//	fprintf(frontPoint, "x,y,z\n");
		//	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		//		fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
		//	}
		//	savingList.clear();
		//	fclose(frontPoint);
		//}
		
		deleteRootHeap3D(inProcess);
		
		//Top
		if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kminus = k - 1;
			short label = labelArray[kminus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					labelArray[kminus][currentIndx] = 2;
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (dTop < actionPtr.imageDataPtr[kminus][currentIndx]) {
						actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kminus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dTop;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}
		
		//Bottom
		if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kplus = k + 1;
			short label = labelArray[kplus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					labelArray[kplus][currentIndx] = 2;
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (dBottom < actionPtr.imageDataPtr[kplus][currentIndx]) {
						actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kplus);
						if(pt_pos != -1) {
							inProcess[pt_pos].arrival = dBottom;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[k][indxEast];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					labelArray[k][indxEast] = 2;
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < actionPtr.imageDataPtr[k][indxEast]) {
						actionPtr.imageDataPtr[k][indxEast] = dEast;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iplus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dEast;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//West
		if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[k][indxWest];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxWest];
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					labelArray[k][indxWest] = 2;
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionPtr.imageDataPtr[k][indxWest]) {
						actionPtr.imageDataPtr[k][indxWest] = dWest;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iminus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dWest;
							heapifyUp3D(inProcess, pt_pos);
						}	
					}
				}
			}
		}
		
		//North
		if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[k][indxNorth];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				dataType coefSpeed = potentialFuncPtr[k][indxNorth];
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					labelArray[k][indxNorth] = 2;
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionPtr.imageDataPtr[k][indxNorth]) {
						actionPtr.imageDataPtr[k][indxNorth] = dNorth;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, jminus, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dNorth;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}
		
		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[k][indxSouth];
			if (label != 1) {
				dataType x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				dataType y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				dataType z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				dataType coefSpeed = potentialFuncPtr[k][indxSouth];
				dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxSouth] = dSouth;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth};
					labelArray[k][indxSouth] = 2;
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < actionPtr.imageDataPtr[k][indxSouth]) {
						actionPtr.imageDataPtr[k][indxSouth] = dSouth;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, jplus, k);
						if( pt_pos != -1) {
							inProcess[pt_pos].arrival = dSouth;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}
		
	}
	
	////Visualize the front propagation
	//id_save++;
	//string saving_csv = path_saving + to_string(id_save) + ".csv";
	//FILE* frontPoint;
	//if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(frontPoint, "x,y,z\n");
	//if (savingList.size() > 0) {
	//	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
	//		fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
	//	}
	//	savingList.clear();
	//}
	//fclose(frontPoint);
	
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

bool frontPropagation(Image_Data ctImageData, dataType** actionPtr, dataType** potentialFuncPtr, Point3D seedPoint) {

	if (actionPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	short** labelArray = new short* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D];
		if (labelArray[k] == NULL) {
			return false; // Memory allocation failed
		}
	}

	//pointer to index of inProcess vector
	size_t** indexArray = new size_t * [length * width * height] {NULL};

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
	size_t currentIndx = x_new(i, j, length);
	actionPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProces

	//Top
	if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
		size_t kminus = k - 1;
		dataType x = select3dX(actionPtr, length, width, height, i, j, kminus);
		dataType y = select3dY(actionPtr, length, width, height, i, j, kminus);
		dataType z = select3dZ(actionPtr, length, width, height, i, j, kminus);
		dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr[kminus][currentIndx] = dTop;
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//Bottom
	if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
		size_t kplus = k + 1;
		dataType x = select3dX(actionPtr, length, width, height, i, j, kplus);
		dataType y = select3dY(actionPtr, length, width, height, i, j, kplus);
		dataType z = select3dZ(actionPtr, length, width, height, i, j, kplus);
		dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr[kplus][currentIndx] = dBottom;
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	//North
	if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
		size_t jminus = j - 1;
		size_t indxNorth = x_new(i, jminus, length);
		dataType x = select3dX(actionPtr, length, width, height, i, jminus, k);
		dataType y = select3dY(actionPtr, length, width, height, i, jminus, k);
		dataType z = select3dZ(actionPtr, length, width, height, i, jminus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxNorth];
		dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr[k][indxNorth] = dNorth;
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
		size_t jplus = j + 1;
		size_t indxSouth = x_new(i, jplus, length);
		dataType x = select3dX(actionPtr, length, width, height, i, jplus, k);
		dataType y = select3dY(actionPtr, length, width, height, i, jplus, k);
		dataType z = select3dZ(actionPtr, length, width, height, i, jplus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxSouth];
		dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		size_t posSouth = inProcess.size();
		actionPtr[k][indxSouth] = dSouth;
		pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth };
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//East
	if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
		size_t iplus = i + 1;
		size_t indxEast = x_new(iplus, j, length);
		dataType x = select3dX(actionPtr, length, width, height, iplus, j, k);
		dataType y = select3dY(actionPtr, length, width, height, iplus, j, k);
		dataType z = select3dZ(actionPtr, length, width, height, iplus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxEast];
		dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr[k][indxEast] = dEast;
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//West
	if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
		size_t iminus = i - 1;
		size_t indxWest = x_new(iminus, j, length);
		dataType x = select3dX(actionPtr, length, width, height, iminus, j, k);
		dataType y = select3dY(actionPtr, length, width, height, iminus, j, k);
		dataType z = select3dZ(actionPtr, length, width, height, iminus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxWest];
		dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		actionPtr[k][indxWest] = dWest;
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		inProcess.push_back(WestNeighbor);
	}

	heapifyVector3D(inProcess);

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching3D current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;

		deleteRootHeap3D(inProcess);

		//Top
		if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kminus = k - 1;
			short label = labelArray[kminus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, i, j, kminus);
				dataType y = select3dY(actionPtr, length, width, height, i, j, kminus);
				dataType z = select3dZ(actionPtr, length, width, height, i, j, kminus);
				dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					labelArray[kminus][currentIndx] = 2;
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (dTop < actionPtr[kminus][currentIndx]) {
						actionPtr[kminus][currentIndx] = dTop;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kminus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dTop;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//Bottom
		if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kplus = k + 1;
			short label = labelArray[kplus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, i, j, kplus);
				dataType y = select3dY(actionPtr, length, width, height, i, j, kplus);
				dataType z = select3dZ(actionPtr, length, width, height, i, j, kplus);
				dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					labelArray[kplus][currentIndx] = 2;
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (dBottom < actionPtr[kplus][currentIndx]) {
						actionPtr[kplus][currentIndx] = dBottom;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kplus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dBottom;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[k][indxEast];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, iplus, j, k);
				dataType y = select3dY(actionPtr, length, width, height, iplus, j, k);
				dataType z = select3dZ(actionPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					labelArray[k][indxEast] = 2;
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < actionPtr[k][indxEast]) {
						actionPtr[k][indxEast] = dEast;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iplus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dEast;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//West
		if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[k][indxWest];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, iminus, j, k);
				dataType y = select3dY(actionPtr, length, width, height, iminus, j, k);
				dataType z = select3dZ(actionPtr, length, width, height, iminus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxWest];
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					labelArray[k][indxWest] = 2;
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionPtr[k][indxWest]) {
						actionPtr[k][indxWest] = dWest;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iminus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dWest;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[k][indxNorth];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, i, jminus, k);
				dataType y = select3dY(actionPtr, length, width, height, i, jminus, k);
				dataType z = select3dZ(actionPtr, length, width, height, i, jminus, k);
				dataType coefSpeed = potentialFuncPtr[k][indxNorth];
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					labelArray[k][indxNorth] = 2;
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionPtr[k][indxNorth]) {
						actionPtr[k][indxNorth] = dNorth;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, jminus, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dNorth;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[k][indxSouth];
			if (label != 1) {
				dataType x = select3dX(actionPtr, length, width, height, i, jplus, k);
				dataType y = select3dY(actionPtr, length, width, height, i, jplus, k);
				dataType z = select3dZ(actionPtr, length, width, height, i, jplus, k);
				dataType coefSpeed = potentialFuncPtr[k][indxSouth];
				dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr[k][indxSouth] = dSouth;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
					labelArray[k][indxSouth] = 2;
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < actionPtr[k][indxSouth]) {
						actionPtr[k][indxSouth] = dSouth;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, jplus, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dSouth;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	delete[] indexArray;

	return true;
}

bool frontPropagationWithKeyPointDetection(Image_Data actionMapStr, dataType** potentialFuncPtr, Point3D* seedPoint, const double LengthKeyPoints, vector<Point3D>& key_points, std::string path_saving) {

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

	vector <pointFastMarching3D> inProcess;
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

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		size_t kminus = k - 1;
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		size_t kplus = k + 1;
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		size_t iplus = i + 1;
		size_t indxEast = x_new(iplus, j, length);
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxEast];
		dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		actionMapStr.imageDataPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		size_t jminus = j - 1;
		size_t indxNorth = x_new(i, jminus, length);
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxNorth];
		dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		size_t indxWest = x_new(iminus, j, length);
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		dataType coefSpeed = potentialFuncPtr[k][indxWest];
		dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		actionMapStr.imageDataPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		size_t jplus = j + 1;
		size_t indxSouth = x_new(i, jplus, length);
		dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		dataType coefSpeed = potentialFuncPtr[k][indxSouth];
		dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;
	size_t nbSourcePoint = 0;
	dataType max_weighted_distance = 0.0;

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
	
	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;

		//Exit of the while loop if the end point is reached
		if (labelArray[kEnd][x_new(iEnd, jEnd, length)] == 1) {
			break;
		}

		//check if the current point is a key point
		Point3D pSource = { i, j, k };
		Point3D pSourceReal = getRealCoordFromImageCoord3D(pSource, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		Point3D pCurrent = getRealCoordFromImageCoord3D(currentSourcePoint, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		distanceToCurrentSourcePoint = getPoint3DDistance(pCurrent, pSourceReal);
		//For vizualization
		savingList.push_back(pSourceReal);
		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			string saving_csv = path_saving + to_string(id_keyPoint) + ".csv";
			FILE* frontPoint;
			if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
				printf("Enable to open");
				return false;
			}
			for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
				fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
			}
			//savingList.clear(); //not empty the list of computed points to see the whole front
			fclose(frontPoint);

			//If the condition is true ---> new key point is found so we need to initilize it neighbors
			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionMapStr.imageDataPtr[k][currentIndx] = 0;

			//Top
			if (k > 0 && j >= 0 && j < width && i >= 0 && i < length)
			{
				size_t kminus = k - 1;
				actionMapStr.imageDataPtr[kminus][currentIndx] = INFINITY;
				labelArray[kminus][currentIndx] = 3;
			}

			//Bottom
			if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width)
			{
				size_t kplus = k + 1;
				actionMapStr.imageDataPtr[kplus][currentIndx] = INFINITY;
				labelArray[kplus][currentIndx] = 3;
			}

			//East
			if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height)
			{
				size_t iplus = i + 1;
				size_t indxEast = x_new(iplus, j, length);
				actionMapStr.imageDataPtr[k][indxEast] = INFINITY;
				labelArray[k][indxEast] = 3;
			}

			//North
			if (j > 0 && i >= 0 && i < length && k >= 0 && k < height)
			{
				size_t jminus = j - 1;
				size_t indxNorth = x_new(i, jminus, length);
				actionMapStr.imageDataPtr[k][indxNorth] = INFINITY;
				labelArray[k][indxNorth] = 3;
			}

			//West
			if (i > 0 && j >= 0 && j < width && k >= 0 && k < height)
			{
				size_t iminus = i - 1;
				size_t indxWest = x_new(iminus, j, length);
				actionMapStr.imageDataPtr[k][indxWest] = INFINITY;
				labelArray[k][indxWest] = 3;
			}

			//South
			if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height)
			{
				size_t jplus = j + 1;
				size_t indxSouth = x_new(i, jplus, length);
				actionMapStr.imageDataPtr[k][indxSouth] = INFINITY;
				labelArray[k][indxSouth] = 3;
			}

			////Initialize all the points inside the narrow band as not processed ---> label = 3
			//for(size_t it = 0; it < inProcess.size(); it++) {
			//	labelArray[inProcess[it].z][x_new(inProcess[it].x, inProcess[it].y, length)] = 1;
			//}
			inProcess.clear();
			id_keyPoint++;

		}
		else {
			actionMapStr.imageDataPtr[k][currentIndx] = current.arrival;
			deleteRootHeap3D(inProcess);
		}

		//====================
		//processed neighbors of the minimum in the narrow band
		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			size_t kminus = k - 1;
			size_t label = labelArray[kminus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
				dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
				dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
				dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					labelArray[kminus][currentIndx] = 2;
					actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (dTop < actionMapStr.imageDataPtr[kminus][currentIndx]) {
						actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kminus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dTop;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kplus = k + 1;
			size_t label = labelArray[kplus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
				dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
				dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
				dataType coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					labelArray[kplus][currentIndx] = 2;
					actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (dBottom < actionMapStr.imageDataPtr[kplus][currentIndx]) {
						actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, j, kplus);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dBottom;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			size_t label = labelArray[k][indxEast];
			if (label != 1) {
				dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
				dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
				dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					labelArray[k][indxEast] = 2;
					actionMapStr.imageDataPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					addPointHeap3D(inProcess, EastNeighbor);

				}
				else {
					if (dEast < actionMapStr.imageDataPtr[k][indxEast]) {
						actionMapStr.imageDataPtr[k][indxEast] = dEast;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iplus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dEast;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			size_t label = labelArray[k][indxNorth];
			if (label != 1) {
				dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
				dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
				dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
				dataType coefSpeed = potentialFuncPtr[k][indxNorth];
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					labelArray[k][indxNorth] = 2;
					actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < actionMapStr.imageDataPtr[k][indxNorth]) {
						actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
						size_t pt_pos = getIndexFromHeap3D(inProcess, i, jminus, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dNorth;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			size_t label = labelArray[k][indxWest];
			if (label != 1) {
				dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
				dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
				dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxWest];
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					labelArray[k][indxWest] = 2;
					actionMapStr.imageDataPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < actionMapStr.imageDataPtr[k][indxWest]) {
						actionMapStr.imageDataPtr[k][indxWest] = dWest;
						size_t pt_pos = getIndexFromHeap3D(inProcess, iminus, j, k);
						if (pt_pos != -1) {
							inProcess[pt_pos].arrival = dWest;
							heapifyUp3D(inProcess, pt_pos);
						}
					}
				}
			}

			//South
			if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
				size_t jplus = j + 1;
				size_t indxSouth = x_new(i, jplus, length);
				size_t label = labelArray[k][indxSouth];
				if (label != 1) {
					dataType x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					dataType y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					dataType z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					dataType coefSpeed = potentialFuncPtr[k][indxSouth];
					dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[k][indxSouth] = 2;
						actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
						pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
						addPointHeap3D(inProcess, SouthNeighbor);
					}
					else {
						if (dSouth < actionMapStr.imageDataPtr[k][indxSouth]) {
							actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
							size_t pt_pos = getIndexFromHeap3D(inProcess, i, jplus, k);
							if (pt_pos != -1) {
								inProcess[pt_pos].arrival = dSouth;
								heapifyUp3D(inProcess, pt_pos);
							}
						}
					}
				}
			}
			//====================

		}

	}
	
	key_points.push_back(seedPoint[1]);
	std::cout << "The source Points are: " << std::endl;
	for (int nb = 0; nb < key_points.size(); nb++) {
		std::cout << "Point " << nb + 1 << " : " << key_points[nb].x << " " << key_points[nb].y << " " << key_points[nb].z << std::endl;
	}

	//id_keyPoint++;
	string saving_csv = path_saving + to_string(id_keyPoint) + ".csv";
	FILE* frontPoint;
	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
	}
	savingList.clear();
	fclose(frontPoint);

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

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

	//Top
	if (x1 >= 0 && x1 < length && y1 >= 0 && y1 < width && z1 > 0 && z1 < height) {
		size_t zminus = z1 - 1;
		dataType x = select3dX(actionFirstFront, length, width, height, x1, y1, zminus);
		dataType y = select3dY(actionFirstFront, length, width, height, x1, y1, zminus);
		dataType z = select3dZ(actionFirstFront, length, width, height, x1, y1, zminus);
		size_t indx = x_new(x1, y1, length);
		dataType coefSpeed = potentialPtr[zminus][indx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D TopNeighbor = { x1, y1, zminus, dTop };
		actionFirstFront[zminus][indx] = dTop;
		narrowBandFirstFront.push_back(TopNeighbor);
		firstLabelArray[zminus][indx] = 2; //2 ---> in process
	}

	//Bottom
	if (x1 >= 0 && x1 < length && y1 >= 0 && y1 < width && z1 >= 0 && z1 < height_minus) {
		size_t zplus = z1 + 1;
		dataType x = select3dX(actionFirstFront, length, width, height, x1, y1, zplus);
		dataType y = select3dY(actionFirstFront, length, width, height, x1, y1, zplus);
		dataType z = select3dZ(actionFirstFront, length, width, height, x1, y1, zplus);
		size_t indx = x_new(x1, y1, length);
		dataType coefSpeed = potentialPtr[zplus][indx];
		dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D BottomNeighbor = { x1, y1, zplus, dBottom };
		actionFirstFront[zplus][indx] = dBottom;
		narrowBandFirstFront.push_back(BottomNeighbor);
		firstLabelArray[zplus][indx] = 2; //2 ---> in process
	}

	//East
	if (x1 > 0 && x1 < length && y1 >= 0 && y1 < width && z1 >= 0 && z1 < height) {
		size_t xminus = x1 - 1;
		dataType x = select3dX(actionFirstFront, length, width, height, xminus, y1, z1);
		dataType y = select3dY(actionFirstFront, length, width, height, xminus, y1, z1);
		dataType z = select3dZ(actionFirstFront, length, width, height, xminus, y1, z1);
		size_t indxWest = x_new(xminus, y1, length);
		dataType coefSpeed = potentialPtr[z1][indxWest];
		dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { xminus, y1, z1, dWest };
		actionFirstFront[z1][indxWest] = dWest;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[z1][indxWest] = 2; //2 ---> in process
	}

	//West
	if (x1 >= 0 && x1 < length_minus && y1 >= 0 && y1 < width && z1 >= 0 && z1 < height) {
		size_t xplus = x1 + 1;
		dataType x = select3dX(actionFirstFront, length, width, height, xplus, y1, z1);
		dataType y = select3dY(actionFirstFront, length, width, height, xplus, y1, z1);
		dataType z = select3dZ(actionFirstFront, length, width, height, xplus, y1, z1);
		size_t indxEast = x_new(xplus, y1, length);
		dataType coefSpeed = potentialPtr[z1][indxEast];
		dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { xplus, y1, z1, dEast };
		actionFirstFront[z1][indxEast] = dEast;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[z1][indxEast] = 2; //2 ---> in process
	}

	//North
	if (y1 > 0 && y1 < length && x1 >= 0 && x1 < length && z1 >= 0 && z1 < height) {
		size_t yminus = y1 - 1;
		dataType x = select3dX(actionFirstFront, length, width, height, x1, yminus, z1);
		dataType y = select3dY(actionFirstFront, length, width, height, x1, yminus, z1);
		dataType z = select3dZ(actionFirstFront, length, width, height, x1, yminus, z1);
		size_t indxNorth = x_new(x1, yminus, length);
		dataType coefSpeed = potentialPtr[z1][indxNorth];
		dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { x1, yminus, z1, dNorth };
		actionFirstFront[z1][indxNorth] = dNorth;
		narrowBandFirstFront.push_back(WestNeighbor);
		firstLabelArray[z1][indxNorth] = 2; //2 ---> in process
	}

	//South
	if (y1 >= 0 && y1 < width_minus && x1 >= 0 && x1 < length && z1 >= 0 && z1 < height) {
		size_t yplus = y1 + 1;
		dataType x = select3dX(actionFirstFront, length, width, height, x1, yplus, z1);
		dataType y = select3dY(actionFirstFront, length, width, height, x1, yplus, z1);
		dataType z = select3dZ(actionFirstFront, length, width, height, x1, yplus, z1);
		size_t indxSouth = x_new(x1, yplus, length);
		dataType coefSpeed = potentialPtr[z1][indxSouth];
		dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { x1, yplus, z1, dSouth };
		actionFirstFront[z1][indxSouth] = dSouth;
		narrowBandFirstFront.push_back(EastNeighbor);
		firstLabelArray[z1][indxSouth] = 2; //2 ---> in process
	}

	//Initialize neighbors for the second front

	//Top
	if (x2 >= 0 && x2 < length && y2 >= 0 && y2 < width && z2 > 0 && z2 < height) {
		size_t zminus = z2 - 1;
		dataType x = select3dX(actionSecondFront, length, width, height, x2, y2, zminus);
		dataType y = select3dY(actionSecondFront, length, width, height, x2, y2, zminus);
		dataType z = select3dZ(actionSecondFront, length, width, height, x2, y2, zminus);
		size_t indx = x_new(x2, y2, length);
		dataType coefSpeed = potentialPtr[zminus][indx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D TopNeighbor = { x2, y2, zminus, dTop };
		actionSecondFront[zminus][indx] = dTop;
		narrowBandSecondFront.push_back(TopNeighbor);
		secondLabelArray[zminus][indx] = 2; //2 ---> in process
	}

	//Bottom
	if (x2 >= 0 && x2 < length && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height_minus) {
		size_t zplus = z2 + 1;
		dataType x = select3dX(actionSecondFront, length, width, height, x2, y2, zplus);
		dataType y = select3dY(actionSecondFront, length, width, height, x2, y2, zplus);
		dataType z = select3dZ(actionSecondFront, length, width, height, x2, y2, zplus);
		size_t indx = x_new(x2, y2, length);
		dataType coefSpeed = potentialPtr[zplus][indx];
		dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D BottomNeighbor = { x2, y2, zplus, dBottom };
		actionSecondFront[zplus][indx] = dBottom;
		narrowBandSecondFront.push_back(BottomNeighbor);
		secondLabelArray[zplus][indx] = 2; //2 ---> in process
	}

	//East
	if (x2 > 0 && x2 < length && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height) {
		size_t xminus = x2 - 1;
		dataType x = select3dX(actionSecondFront, length, width, height, xminus, y2, z2);
		dataType y = select3dY(actionSecondFront, length, width, height, xminus, y2, z2);
		dataType z = select3dZ(actionSecondFront, length, width, height, xminus, y2, z2);
		size_t indxWest = x_new(xminus, y2, length);
		dataType coefSpeed = potentialPtr[z2][indxWest];
		dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { xminus, y2, z2, dWest };
		actionSecondFront[z2][indxWest] = dWest;
		narrowBandSecondFront.push_back(WestNeighbor);
		secondLabelArray[z2][indxWest] = 2; //2 ---> in process
	}

	//West
	if (x2 >= 0 && x2 < length_minus && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height) {
		size_t xplus = x2 + 1;
		dataType x = select3dX(actionSecondFront, length, width, height, xplus, y2, z2);
		dataType y = select3dY(actionSecondFront, length, width, height, xplus, y2, z2);
		dataType z = select3dZ(actionSecondFront, length, width, height, xplus, y2, z2);
		size_t indxEast = x_new(xplus, y2, length);
		dataType coefSpeed = potentialPtr[z2][indxEast];
		dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { xplus, y2, z2, dEast };
		actionSecondFront[z2][indxEast] = dEast;
		narrowBandSecondFront.push_back(EastNeighbor);
		secondLabelArray[z2][indxEast] = 2; //2 ---> in process
	}

	//North
	if (y2 > 0 && y2 < width && x2 >= 0 && x2 < length && z2 >= 0 && z2 < height) {
		size_t yminus = y2 - 1;
		dataType x = select3dX(actionSecondFront, length, width, height, x2, yminus, z2);
		dataType y = select3dY(actionSecondFront, length, width, height, x2, yminus, z2);
		dataType z = select3dZ(actionSecondFront, length, width, height, x2, yminus, z2);
		size_t indxNorth = x_new(x2, yminus, length);
		dataType coefSpeed = potentialPtr[z2][indxNorth];
		dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D NorthNeighbor = { x2, yminus, z2, dNorth };
		actionSecondFront[z2][indxNorth] = dNorth;
		narrowBandSecondFront.push_back(NorthNeighbor);
		secondLabelArray[z2][indxNorth] = 2; //2 ---> in process
	}

	//South
	if (y2 >= 0 && y2 < width_minus && x2 >= 0 && x2 < length && z2 >= 0 && z2 < height) {
		size_t yplus = y2 + 1;
		dataType x = select3dX(actionSecondFront, length, width, height, x2, yplus, z2);
		dataType y = select3dY(actionSecondFront, length, width, height, x2, yplus, z2);
		dataType z = select3dZ(actionSecondFront, length, width, height, x2, yplus, z2);
		size_t indxSouth = x_new(x2, yplus, length);
		dataType coefSpeed = potentialPtr[z2][indxSouth];
		dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D SouthNeighbor = { x2, yplus, z2, dSouth };
		actionSecondFront[z2][indxSouth] = dSouth;
		narrowBandSecondFront.push_back(SouthNeighbor);
		secondLabelArray[z2][indxSouth] = 2; //2 ---> in process
	}
	

	//heapify 2D vector
	heapifyVector3D(narrowBandFirstFront);
	heapifyVector3D(narrowBandSecondFront);

	//Save points for visualization
	//vector<Point2D> savingList;
	//size_t id_save = 0;

	dataType max_save_action = 0.0;
	size_t nb_computed_points = 0;

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
		nb_computed_points++;

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
		nb_computed_points++;

		////Save points for visualization
		//Point2D point = { (dataType)i, (dataType)j };
		//savingList.push_back(point);
		//if (nb_computed_points % 500 == 0) {
		//	id_save++;
		//	//string saving_csv = savingPath + to_string(id_save) + ".csv";
		//	//FILE* frontPoint;
		//	//if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		//	//	printf("Enable to open");
		//	//	return false;
		//	//}
		//	//fprintf(frontPoint, "x,y\n");
		//	//for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
		//	//	fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
		//	//}
		//	////Don't empty the list of points to see the whole computed points
		//	////savingList.clear();
		//	//fclose(frontPoint);
		//	string saving_file = savingPath + to_string(id_save) + ".raw";
		//	for(size_t in = 0; in < savingList.size(); in++) {
		//		size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
		//		partialDistancePtr[indx] = distancePtr[indx];
		//	}
		//	manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);
		//}

		//Top first front
		if (z1 > 0 && z1 < height && x1 >= 0 && x1 < length && y1 >= 0 && y1 < width) {
			size_t zminus = z1 - 1;
			size_t indx = x_new(x1, y1, length);
			short label = firstLabelArray[zminus][indx];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, x1, y1, zminus);
				dataType y = select3dY(actionFirstFront, length, width, height, x1, y1, zminus);
				dataType z = select3dZ(actionFirstFront, length, width, height, x1, y1, zminus);
				dataType coefSpeed = potentialPtr[zminus][indx];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D TopNeighbor = { x1, y1, zminus, dTop };
				if (label == 3) {
					actionFirstFront[zminus][indx] = dTop;
					firstLabelArray[zminus][indx] = 2;
					addPointHeap3D(narrowBandFirstFront, TopNeighbor);
				}
				else {
					if (dTop < actionFirstFront[zminus][indx]) {
						actionFirstFront[zminus][indx] = dTop;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, x1, y1, zminus);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//Bottom first front
		if (z1 >= 0 && z1 < height_minus && x1 >= 0 && x1 < length && y1 >= 0 && y1 < width) {
			size_t zplus = z1 + 1;
			size_t indx = x_new(x1, y1, length);
			short label = firstLabelArray[zplus][indx];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, x1, y1, zplus);
				dataType y = select3dY(actionFirstFront, length, width, height, x1, y1, zplus);
				dataType z = select3dZ(actionFirstFront, length, width, height, x1, y1, zplus);
				dataType coefSpeed = potentialPtr[zplus][indx];
				dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D BottomNeighbor = { x1, y1, zplus, dBottom };
				if (label == 3) {
					actionFirstFront[zplus][indx] = dBottom;
					firstLabelArray[zplus][indx] = 2;
					addPointHeap3D(narrowBandFirstFront, BottomNeighbor);
				}
				else {
					if (dBottom < actionFirstFront[zplus][indx]) {
						actionFirstFront[zplus][indx] = dBottom;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, x1, y1, zplus);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//West first front
		if (x1 > 0 && x1 < length && y1 >= 0 && y1 < width && z1 >= 0 && z1 < height) {
			size_t xminus = x1 - 1;
			size_t indxWest = x_new(xminus, y1, length);
			short label = firstLabelArray[z1][indxWest];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, xminus, y1, z1);
				dataType y = select3dY(actionFirstFront, length, width, height, xminus, y1, z1);
				dataType z = select3dZ(actionFirstFront, length, width, height, xminus, y1, z1);
				dataType coefSpeed = potentialPtr[z1][indxWest];
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D WestNeighbor = { xminus, y1, z1, dWest };
				if (label == 3) {
					actionFirstFront[z1][indxWest] = dWest;
					firstLabelArray[z1][indxWest] = 2;
					addPointHeap3D(narrowBandFirstFront, WestNeighbor);
				}
				else {
					if (dWest < actionFirstFront[z1][indxWest]) {
						actionFirstFront[z1][indxWest] = dWest;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, xminus, y1, z1);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//East first front
		if (x1 < length_minus && x1 >= 0 && y1 >= 0 && y1 < width && z1 >= 0 && z1 < height) {
			size_t xplus = x1 + 1;
			size_t indxEast = x_new(xplus, y1, length);
			short label = firstLabelArray[z1][indxEast];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, xplus, y1, z1);
				dataType y = select3dY(actionFirstFront, length, width, height, xplus, y1, z1);
				dataType z = select3dZ(actionFirstFront, length, width, height, xplus, y1, z1);
				dataType coefSpeed = potentialPtr[z1][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D EastNeighbor = { xplus, y1, z1, dEast };
				if (label == 3) {
					actionFirstFront[z1][indxEast] = dEast;
					firstLabelArray[z1][indxEast] = 2;
					addPointHeap3D(narrowBandFirstFront, EastNeighbor);
				}
				else {
					if (dEast < actionFirstFront[z1][indxEast]) {
						actionFirstFront[z1][indxEast] = dEast;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, xplus, y1, z1);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//North first front
		if (y1 > 0 && y1 < width && x1 >= 0 && x1 < length && z1 >= 0 && z1 < height) {
			size_t yminus = y1 - 1;
			size_t indxNorth = x_new(x1, yminus, length);
			short label = firstLabelArray[z1][indxNorth];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, x1, yminus, z1);
				dataType y = select3dY(actionFirstFront, length, width, height, x1, yminus, z1);
				dataType z = select3dZ(actionFirstFront, length, width, height, x1, yminus, z1);
				dataType coefSpeed = potentialPtr[z1][indxNorth];
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D NorthNeighbor = { x1, yminus, z1, dNorth };
				if (label == 3) {
					actionFirstFront[z1][indxNorth] = dNorth;
					firstLabelArray[z1][indxNorth] = 2;
					addPointHeap3D(narrowBandFirstFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionFirstFront[z1][indxNorth]) {
						actionFirstFront[z1][indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, x1, yminus, z1);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//South first front
		if (y1 >= 0 && y1 < width_minus && x1 >= 0 && x1 < length && z1 >= 0 && z1 < height) {
			size_t yplus = y1 + 1;
			size_t indxSouth = x_new(x1, yplus, length);
			short label = firstLabelArray[z1][indxSouth];
			if (label != 1) {
				dataType x = select3dX(actionFirstFront, length, width, height, x1, yplus, z1);
				dataType y = select3dY(actionFirstFront, length, width, height, x1, yplus, z1);
				dataType z = select3dZ(actionFirstFront, length, width, height, x1, yplus, z1);
				dataType coefSpeed = potentialPtr[z1][indxSouth];
				dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D SouthNeighbor = { x1, yplus, z1, dSouth };
				if (label == 3) {
					actionFirstFront[z1][indxSouth] = dSouth;
					firstLabelArray[z1][indxSouth] = 2;
					addPointHeap3D(narrowBandFirstFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionFirstFront[z1][indxSouth]) {
						actionFirstFront[z1][indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap3D(narrowBandFirstFront, x1, yplus, z1);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandFirstFront, pIndex);
						}
					}
				}
			}
		}

		//==========================================

		//Top second front
		if (x2 >= 0 && x2 < length && y2 >= 0 && y2 < width && z2 > 0 && z2 < height) {
			size_t zminus = z2 - 1;
			//size_t indx = x_new(x2, y2, length);
			short label = secondLabelArray[zminus][indexSecond];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, x2, y2, zminus);
				dataType y = select3dY(actionSecondFront, length, width, height, x2, y2, zminus);
				dataType z = select3dZ(actionSecondFront, length, width, height, x2, y2, zminus);
				dataType coefSpeed = potentialPtr[zminus][indexSecond];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D TopNeighbor = { x2, y2, zminus, dTop };
				if (label == 3) {
					actionSecondFront[zminus][indexSecond] = dTop;
					secondLabelArray[zminus][indexSecond] = 2;
					addPointHeap3D(narrowBandSecondFront, TopNeighbor);
				}
				else {
					if (dTop < actionSecondFront[zminus][indexSecond]) {
						actionSecondFront[zminus][indexSecond] = dTop;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, x2, y2, zminus);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//Bottom second front
		if (x2 >= 0 && x2 < length && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height_minus) {
			size_t zplus = z2 + 1;
			//size_t indx = x_new(x2, y2, length);
			short label = secondLabelArray[zplus][indexSecond];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, x2, y2, zplus);
				dataType y = select3dY(actionSecondFront, length, width, height, x2, y2, zplus);
				dataType z = select3dZ(actionSecondFront, length, width, height, x2, y2, zplus);
				dataType coefSpeed = potentialPtr[zplus][indexSecond];
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D TopNeighbor = { x2, y2, zplus, dTop };
				if (label == 3) {
					actionSecondFront[zplus][indexSecond] = dTop;
					secondLabelArray[zplus][indexSecond] = 2;
					addPointHeap3D(narrowBandSecondFront, TopNeighbor);
				}
				else {
					if (dTop < actionSecondFront[zplus][indexSecond]) {
						actionSecondFront[zplus][indexSecond] = dTop;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, x2, y2, zplus);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//West second front
		if (x2 > 0 && x2 < length && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height) {
			size_t xminus = x2 - 1;
			size_t indxWest = x_new(xminus, y2, length);
			short label = secondLabelArray[z2][indxWest];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, xminus, y2, z2);
				dataType y = select3dY(actionSecondFront, length, width, height, xminus, y2, z2);
				dataType z = select3dZ(actionSecondFront, length, width, height, xminus, y2, z2);
				dataType coefSpeed = potentialPtr[z2][indxWest];
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D WestNeighbor = { xminus, y2, z2, dWest };
				if (label == 3) {
					actionSecondFront[z2][indxWest] = dWest;
					secondLabelArray[z2][indxWest] = 2;
					addPointHeap3D(narrowBandSecondFront, WestNeighbor);
				}
				else {
					if (dWest < actionSecondFront[z2][indxWest]) {
						actionSecondFront[z2][indxWest] = dWest;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, xminus, y2, z2);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//East second front
		if (x2 < length_minus && x2 >= 0 && y2 >= 0 && y2 < width && z2 >= 0 && z2 < height) {
			size_t xplus = x2 + 1;
			size_t indxEast = x_new(xplus, y2, length);
			short label = secondLabelArray[z2][indxEast];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, xplus, y2, z2);
				dataType y = select3dY(actionSecondFront, length, width, height, xplus, y2, z2);
				dataType z = select3dZ(actionSecondFront, length, width, height, xplus, y2, z2);
				dataType coefSpeed = potentialPtr[z2][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D EastNeighbor = { xplus, y2, z2, dEast };
				if (label == 3) {
					actionSecondFront[z2][indxEast] = dEast;
					secondLabelArray[z2][indxEast] = 2;
					addPointHeap3D(narrowBandSecondFront, EastNeighbor);
				}
				else {
					if (dEast < actionSecondFront[z2][indxEast]) {
						actionSecondFront[z2][indxEast] = dEast;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, xplus, y2, z2);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//North second front
		if (y2 > 0 && y2 < width && x2 >= 0 && x2 < length && z2 >= 0 && z2 < height) {
			size_t yminus = y2 - 1;
			size_t indxNorth = x_new(x2, yminus, length);
			short label = secondLabelArray[z2][indxNorth];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, x2, yminus, z2);
				dataType y = select3dY(actionSecondFront, length, width, height, x2, yminus, z2);
				dataType z = select3dZ(actionSecondFront, length, width, height, x2, yminus, z2);
				dataType coefSpeed = potentialPtr[z2][indxNorth];
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D NorthNeighbor = { x2, yminus, z2, dNorth };
				if (label == 3) {
					actionSecondFront[z2][indxNorth] = dNorth;
					secondLabelArray[z2][indxNorth] = 2;
					addPointHeap3D(narrowBandSecondFront, NorthNeighbor);
				}
				else {
					if (dNorth < actionSecondFront[z2][indxNorth]) {
						actionSecondFront[z2][indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, x2, yminus, z2);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//South second front
		if (y2 >= 0 && y2 < width_minus && x2 >= 0 && x2 < length && z2 >= 0 && z2 < height) {
			size_t yplus = y2 + 1;
			size_t indxSouth = x_new(x2, yplus, length);
			short label = secondLabelArray[z2][indxSouth];
			if (label != 1) {
				dataType x = select3dX(actionSecondFront, length, width, height, x2, yplus, z2);
				dataType y = select3dY(actionSecondFront, length, width, height, x2, yplus, z2);
				dataType z = select3dZ(actionSecondFront, length, width, height, x2, yplus, z2);
				dataType coefSpeed = potentialPtr[z2][indxSouth];
				dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				pointFastMarching3D SouthNeighbor = { x2, yplus, z2, dSouth };
				if (label == 3) {
					actionSecondFront[z2][indxSouth] = dSouth;
					secondLabelArray[z2][indxSouth] = 2;
					addPointHeap3D(narrowBandSecondFront, SouthNeighbor);
				}
				else {
					if (dSouth < actionSecondFront[z2][indxSouth]) {
						actionSecondFront[z2][indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap3D(narrowBandSecondFront, x2, yplus, z2);
						if (pIndex != -1) {
							heapifyUp3D(narrowBandSecondFront, pIndex);
						}
					}
				}
			}
		}

		//if (nb_computed_points == 150000) {
		//	break;
		//}

	}

	////Save points for visualization
	//if (savingList.size() != 0) {
	//	id_save++;
	//	string saving_csv = savingPath + to_string(id_save) + ".csv";
	//	FILE* frontPoint;
	//	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
	//		printf("Enable to open");
	//		return false;
	//	}
	//	fprintf(frontPoint, "x,y\n");
	//	for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
	//		fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
	//	}
	//	savingList.clear();
	//	fclose(frontPoint);
	//}

	//id_save++;
	//string saving_file = savingPath + to_string(id_save) + ".raw";
	//for (size_t in = 0; in < savingList.size(); in++) {
	//	size_t indx = x_new((size_t)savingList[in].x, (size_t)savingList[in].y, length);
	//	partialDistancePtr[indx] = distancePtr[indx];
	//}
	//manageRAWFile2D<dataType>(distancePtr, length, width, saving_file.c_str(), STORE_DATA, false);

	nb_computed_points += narrowBandFirstFront.size() + narrowBandSecondFront.size();
	std::cout << "Number of points processed: " << nb_computed_points << std::endl;
	std::cout << "Max arrival before exit: " << max_save_action << std::endl;

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
			if (ctImageData.imageDataPtr[k][i] == foregroundValue) 
			{
				distanceFuncPtr[k][i] = 0;
				labelArray[k][i] = 1;
			}
			else {
				distanceFuncPtr[k][i] = INFINITY;
				labelArray[k][i] = 3;
			}
			
		}
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
					//Top
					if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
						size_t kminus = k - 1;
						short label = labelArray[kminus][currentIndx];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
							dataType y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
							dataType coefSpeed = 1.0;
							dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[kminus][currentIndx] = dTop;
								pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
								labelArray[kminus][currentIndx] = 2;
								inProcess.push_back(TopNeighbor);
							}
							else {
								if (dTop < distanceFuncPtr[kminus][currentIndx]) {
									distanceFuncPtr[kminus][currentIndx] = dTop;
									int indexTop = getIndexFromHeap3D(inProcess, i, j, kminus);
									if (indexTop != -1) {
										inProcess[indexTop].arrival = dTop;
										//heapifyVector3D(inProcess);
									}
								}
							}
						}
					}

					//Bottom
					if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
						size_t kplus = k + 1;
						short label = labelArray[kplus][currentIndx];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
							dataType y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
							dataType coefSpeed = 1.0;
							dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[kplus][currentIndx] = dBottom;
								pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
								labelArray[kplus][currentIndx] = 2;
								inProcess.push_back(BottomNeighbor);
							}
							else {
								if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
									distanceFuncPtr[kplus][currentIndx] = dBottom;
									int indexBottom = getIndexFromHeap3D(inProcess, i, j, kplus);
									if (indexBottom != -1) {
										inProcess[indexBottom].arrival = dBottom;
									}
								}
							}
						}
					}

					//East
					if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
						size_t iplus = i + 1;
						size_t indxEast = x_new(iplus, j, length);// iplus + j * length;
						short label = labelArray[k][indxEast];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
							dataType y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
							dataType coefSpeed = 1.0;
							dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[k][indxEast] = dEast;
								pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
								labelArray[k][indxEast] = 2;
								inProcess.push_back(EastNeighbor);
							}
							else {
								if (dEast < distanceFuncPtr[k][indxEast]) {
									distanceFuncPtr[k][indxEast] = dEast;
									int indexEast = getIndexFromHeap3D(inProcess, iplus, j, k);
									if (indexEast != -1) {
										inProcess[indexEast].arrival = dEast;
									}
								}
							}
						}
					}

					//West
					if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
						size_t iminus = i - 1;
						size_t indxWest = x_new(iminus, j, length);
						short label = labelArray[k][indxWest];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
							dataType y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
							dataType coefSpeed = 1.0;
							dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[k][indxWest] = dWest;
								pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
								labelArray[k][indxWest] = 2;
								inProcess.push_back(WestNeighbor);
							}
							else {
								if (dWest < distanceFuncPtr[k][indxWest]) {
									distanceFuncPtr[k][indxWest] = dWest;
									int indexWest = getIndexFromHeap3D(inProcess, iminus, j, k);
									if (indexWest != -1) {
										inProcess[indexWest].arrival = dWest;
									}
								}
							}
						}
					}

					//North
					if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
						size_t jminus = j - 1;
						size_t indxNorth = x_new(i, jminus, length);
						short label = labelArray[k][indxNorth];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
							dataType y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
							dataType coefSpeed = 1.0;
							dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[k][indxNorth] = dNorth;
								pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
								labelArray[k][indxNorth] = 2;
								inProcess.push_back(NorthNeighbor);
							}
							else {
								if (dNorth < distanceFuncPtr[k][indxNorth]) {
									distanceFuncPtr[k][indxNorth] = dNorth;
									int indexNorth = getIndexFromHeap3D(inProcess, i, jminus, k);
									if(indexNorth != -1) {
										inProcess[indexNorth].arrival = dNorth;
									}
								}
							}
						}
					}

					//South
					if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
						size_t jplus = j + 1;
						size_t indxSouth = x_new(i, jplus, length);
						short label = labelArray[k][indxSouth];
						if (label != 1) {
							dataType x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
							dataType y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
							dataType z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
							dataType coefSpeed = 1.0;
							dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
							if (label == 3) {
								distanceFuncPtr[k][indxSouth] = dSouth;
								pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
								labelArray[k][indxSouth] = 2;
								inProcess.push_back(SouthNeighbor);
							}
							else {
								if (dSouth < distanceFuncPtr[k][indxSouth]) {
									distanceFuncPtr[k][indxSouth] = dSouth;
									int indexSouth = getIndexFromHeap3D(inProcess, i, jplus, k);
									if(indexSouth != -1) {
										inProcess[indexSouth].arrival = dSouth;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	heapifyVector3D(inProcess);

	size_t id_passage = 0;
	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		pointFastMarching3D current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		size_t currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;

		deleteRootHeap3D(inProcess);

		//Top
		if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kminus = k - 1;
			short label = labelArray[kminus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
				dataType y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
				dataType coefSpeed = 1.0;
				dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					labelArray[kminus][currentIndx] = 2;
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (dTop < distanceFuncPtr[kminus][currentIndx]) {
						distanceFuncPtr[kminus][currentIndx] = dTop;
						int indexTop = getIndexFromHeap3D(inProcess, i, j, kminus);
						if(indexTop != -1) {
							inProcess[indexTop].arrival = dTop;
							heapifyUp3D(inProcess, indexTop);
						}
					}
				}
			}
		}

		//Bottom
		if (k >= 0 && k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			size_t kplus = k + 1;
			short label = labelArray[kplus][currentIndx];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
				dataType y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
				dataType coefSpeed = 1.0;
				dataType dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					labelArray[kplus][currentIndx] = 2;
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
						distanceFuncPtr[kplus][currentIndx] = dBottom;
						int indexBottom = getIndexFromHeap3D(inProcess, i, j, kplus);
						if (indexBottom != -1) {
							inProcess[indexBottom].arrival = dBottom;
							heapifyUp3D(inProcess, indexBottom);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);// iplus + j * length;
			short label = labelArray[k][indxEast];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
				dataType y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = 1.0;
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					labelArray[k][indxEast] = 2;
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (dEast < distanceFuncPtr[k][indxEast]) {
						distanceFuncPtr[k][indxEast] = dEast;
						int indexEast = getIndexFromHeap3D(inProcess, iplus, j, k);
						if (indexEast != -1) {
							inProcess[indexEast].arrival = dEast;
							heapifyUp3D(inProcess, indexEast);
						}
					}
				}
			}
		}

		//West
		if (i > 0 && i < length && j >= 0 && j < width && k >= 0 && k < height) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[k][indxWest];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
				dataType y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
				dataType coefSpeed = 1.0;
				dataType dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					labelArray[k][indxWest] = 2;
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (dWest < distanceFuncPtr[k][indxWest]) {
						distanceFuncPtr[k][indxWest] = dWest;
						int indexWest = getIndexFromHeap3D(inProcess, iminus, j, k);
						if (indexWest != -1) {
							inProcess[indexWest].arrival = dWest;
							heapifyUp3D(inProcess, indexWest);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[k][indxNorth];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
				dataType y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
				dataType coefSpeed = 1.0;
				dataType dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					labelArray[k][indxNorth] = 2;
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (dNorth < distanceFuncPtr[k][indxNorth]) {
						distanceFuncPtr[k][indxNorth] = dNorth;
						int indexNorth = getIndexFromHeap3D(inProcess, i, jminus, k);
						if (indexNorth != -1) {
							inProcess[indexNorth].arrival = dNorth;
							heapifyUp3D(inProcess, indexNorth);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[k][indxSouth];
			if (label != 1) {
				dataType x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
				dataType y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
				dataType z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
				dataType coefSpeed = 1.0;
				dataType dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					distanceFuncPtr[k][indxSouth] = dSouth;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
					labelArray[k][indxSouth] = 2;
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (dSouth < distanceFuncPtr[k][indxSouth]) {
						distanceFuncPtr[k][indxSouth] = dSouth;
						int indexSouth = getIndexFromHeap3D(inProcess, i, jplus, k);
						if (indexSouth != -1) {
							inProcess[indexSouth].arrival = dSouth;
							heapifyUp3D(inProcess, indexSouth);
						}
					}
				}
			}
		}

	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				dataType x = select3dX(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType y = select3dY(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
				dataType z = select3dZ(distancePtr, length, width, height, (size_t)i, (size_t)j, (size_t)k);
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
				for (size_t tk = 0; tk < height; tk++) {
					for (size_t ti = 0; ti < length; ti++) {
						for (size_t tj = 0; tj < width; tj++) {
							size_t txd = x_new(ti, tj, length);
							if (ctImageData.imageDataPtr[tk][txd] == foregroundValue) {
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


