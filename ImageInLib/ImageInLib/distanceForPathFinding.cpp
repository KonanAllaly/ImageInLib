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
		a = hx_2;
		b = -2 * hx_2 * Y;
		c = (dataType)(hx_2 * Y * Y - hx_2 * hy_2 * P_2);
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
		a = hy_2;
		b = -2 * hy_2 * X;
		c = (dataType)(hy_2 * X * X - hx_2 * hy_2 * P_2);
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

	dataType seedVal = (dataType)(imageDataStr.imageDataPtr[x_new((size_t)seedPoints[0].x, (size_t)seedPoints[0].y, length)]);

	dataType* distanceMap = new dataType[dim2D]{ 0 };
	dataType* edgeImage = new dataType[dim2D]{ 0 };	

	dataType norm_of_gradient = 0.0, edgeValue = 0.0;
	bool isGradientComputed = false;
	Point2D grad_vector;
	PixelSpacing fVolume = imageDataStr.spacing;

	for (i = 0; i < length; i++) 
	{
		for (j = 0; j < width; j++) 
		{
			isGradientComputed = getGradient2D(imageDataStr.imageDataPtr, length, width, i, j, fVolume, &grad_vector);
			if (isGradientComputed == true) {
				norm_of_gradient = sqrt(grad_vector.x * grad_vector.x + grad_vector.y * grad_vector.y);
				edgeValue = gradientFunction(norm_of_gradient, parameters.K);
			}
			else {
				std::cout << "Error in computing gradient at point (" << i << ", " << j << ")" << std::endl;
				return false;
			}
			//Threshold
			if( edgeValue < parameters.thres) 
			{
				edgeValue = 0.0;
			}
			else 
			{
				edgeValue = 1.0;
			}
		}
	}

	//fastSweepingFunction_2D(distanceMap, edgeImage, length, width, 1.0, 1000000.0, 1.0);
	//string path_file = "C:/Users/Konan Allaly/Documents/Tests/output/distance_map_2d.raw";
	//manageRAWFile2D<dataType>(distanceMap, length, width, path_file.c_str(), STORE_DATA, false);

	for (i = 0; i < dim2D; i++) {
		potentialFuncPtr[i] = fabs(imageDataStr.imageDataPtr[i] - seedVal);
	}

	//Find max difference
	dataType maxDiff = 0.0;
	for (i = 0; i < dim2D; i++) {
		if (potentialFuncPtr[i] > maxDiff) {
			maxDiff = potentialFuncPtr[i];
		}
	}
	
	//Normalization
	dataType weight = 0.0;
	for (i = 0; i < dim2D; i++) {
		weight = 1.0 / (1.0 + 1.0 * distanceMap[i]);
		potentialFuncPtr[i] = (parameters.eps + potentialFuncPtr[i] / maxDiff) * weight;
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

	while (labelArray[seedIndex] != 1) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;
		nb_computed_points++;

		//Save points for visualization
		Point2D point = { (dataType)i, (dataType)j };
		savingList.push_back(point);
		if (nb_computed_points % 50 == 0) {
			id_save++;
			string saving_csv = savingPath + to_string(id_save) + ".csv";
			FILE* frontPoint;
			if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
				printf("Enable to open");
				return false;
			}
			fprintf(frontPoint, "x,y\n");
			for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
				fprintf(frontPoint, "%f,%f\n", savingList[i_n].x, savingList[i_n].y);
			}
			//Don't empty the list of points to see the whole computed points
			//savingList.clear();
			fclose(frontPoint);
		}
		
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
	
	//Save points for visualization
	if (savingList.size() != 0) {
		id_save++;
		string saving_csv = savingPath + to_string(id_save) + ".csv";
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
	}
	
	delete[] labelArray;
}

bool shortestPath2d(Image_Data2D distanceFuncPtr, Point2D* seedPoints, vector<Point2D>& path_points, Path_Parameters parameters) {

	if (distanceFuncPtr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	const size_t height = distanceFuncPtr.height;
	const size_t width = distanceFuncPtr.width;
	size_t i, j, dim2D = height * width;
	dataType dist_min = 0.0;

	dataType* gradientVectorX = new dataType[dim2D]{0};
	dataType* gradientVectorY = new dataType[dim2D]{0};

	//Normalization of the gradient
	computeImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY);
	
	dataType norm_of_gradient = 0.0, ux = 0.0, uy = 0.0;
	for (i = 0; i < dim2D; i++) {
		ux = gradientVectorX[i];
		uy = gradientVectorY[i];
		norm_of_gradient = sqrt(ux * ux + uy * uy);
		gradientVectorX[i] = gradientVectorX[i] / norm_of_gradient;
		gradientVectorY[i] = gradientVectorY[i] / norm_of_gradient;
	}

	size_t count_iter = 0;
	i = (size_t)seedPoints[1].x; 
	j = (size_t)seedPoints[1].y;
	size_t currentIndex = x_new(i, j, height);

	dataType iNew = seedPoints[1].x, jNew = seedPoints[1].y;
	
	dataType radius = 0.0, offset = 2.0;
	BoundingBox2D box = { 0, 0, 0, 0 };
	size_t n, m;

	do{

		iNew -= parameters.tau * gradientVectorX[currentIndex];
		jNew -= parameters.tau * gradientVectorY[currentIndex];
		
		Point2D current_point = { iNew, jNew };
		dist_min = getPoint2DDistance(current_point, seedPoints[0]);
		path_points.push_back(current_point);

		i = (size_t)round(iNew); 
		j = (size_t)round(jNew);
		currentIndex = x_new(i, j, height);

		count_iter++;
	}
	while(dist_min > parameters.tolerance && count_iter < parameters.max_iteration);

	delete[] gradientVectorX;
	delete[] gradientVectorY;

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
 
dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W) {

	dataType sol = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType R = sqrt(W);

	if (X == INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 2;
		b = (dataType)(- 2 * (Y + Z));
		c = (dataType)(pow(Y, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(Y, Z) + R);
		}
	}

	if (Y == INFINITY && X != INFINITY && Z != INFINITY) {
		a = 2;
		b = (dataType)(- 2 * (X + Z));
		c = (dataType)(pow(X, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, Z) + R);
		}
	}

	if (Z == INFINITY && X != INFINITY && Y != INFINITY) {
		a = 2;
		b = (dataType)(- 2 * (X + Y));
		c = (dataType)(pow(X, 2) + pow(Y, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, Y) + R);
		}
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY) {
		a = 1;
		b = (dataType)(- 2 * Z);
		c = (dataType)(pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(Z + R);
		}
	}

	if (X == INFINITY && Z == INFINITY && Y != INFINITY) {
		a = 1;
		b = (dataType)(- 2 * Y);
		c = (dataType)(pow(Y, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(Y + R);
		}
	}

	if (Y == INFINITY && Z == INFINITY && X != INFINITY) {
		a = 1;
		b = (dataType)(- 2 * X);
		c = (dataType)(pow(X, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(X + R);
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 3;
		b = (dataType)(- 2 * (X + Y + Z));
		c = (dataType)(pow(X, 2) + pow(Y, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, min(Y, Z)) + R);
		}
	}

}

dataType select3dX(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {
	//if (i == 0) {
	//	return min(actionPtr[k][x_new(i, j, length)], actionPtr[k][x_new(i + 1, j, length)]);
	//}
	//else {
	//	if (i == length - 1) {
	//		return min(actionPtr[k][x_new(i - 1, j, length)], actionPtr[k][x_new(i, j, length)]);
	//	}
	//	else {
	//		return min(actionPtr[k][x_new(i - 1, j, length)], actionPtr[k][x_new(i + 1, j, length)]);
	//	}
	//}
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
	//if (j == 0) {
	//	return min(actionPtr[k][x_new(i, j, length)], actionPtr[k][x_new(i, j + 1, length)]);
	//}
	//else {
	//	if (j == width - 1) {
	//		return min(actionPtr[k][x_new(i, j - 1, length)], actionPtr[k][x_new(i, j, length)]);
	//	}
	//	else {
	//		return min(actionPtr[k][x_new(i, j - 1, length)], actionPtr[k][x_new(i, j + 1, length)]);
	//	}
	//}
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
	//if (k == 0) {
	//	return min(actionPtr[k + 1][xd], actionPtr[k][xd]);
	//}
	//else {
	//	if (k == height - 1) {
	//		return min(actionPtr[k - 1][xd], actionPtr[k][xd]);
	//	}
	//	else {
	//		return min(actionPtr[k + 1][xd], actionPtr[k - 1][xd]);
	//	}
	//}
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
	
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
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
				//threshold : real image
				if (edgeValue <= parameters.thres) {
					maskThreshold[k][xd] = 1.0;
				}
				else {
					maskThreshold[k][xd] = 0.0;
				}
				////threshold artificial image
				//if (norm_of_gradient > 0.0) {
				//	maskThreshold[k][xd] = 0.0;
				//}
				//else {
				//	maskThreshold[k][xd] = 1.0;
				//}
			}
		}
	}
	
	//////Real image
	Image_Data toDistanceMap = { height, length, width, maskThreshold, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	std::string storing_path = "C:/Users/Konan Allaly/Documents/Tests/input/edge_image_crop_p3.raw";
	manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), LOAD_DATA, false);
	fastSweepingDistanceMap(toDistanceMap, distance, 1.0);
	//fastMarching3dForDistanceMap(toDistanceMap, distance, 1.0);
	//rouyTourinDistanceMap(toDistanceMap, distance, 0.001, 1000, 1.0);
	//////bruteForceDistanceMap(toDistanceMap, distance, 0.0);
	storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/distance_fs.raw";
	manageRAWFile3D<dataType>(distance, length, width, height, storing_path.c_str(), STORE_DATA, false);
	////storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/crop/edge_image.raw";
	////manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), STORE_DATA, false);

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
	std::cout << "Seed point mean value: " << seedStats.mean_data << std::endl;
	std::cout << "Seed point minimum value: " << seedStats.min_data << std::endl;
	std::cout << "Seed point maximum value: " << seedStats.max_data << std::endl;
	std::cout << "Standard deviation value: " << seedStats.sd_data << std::endl;
	
	dataType var_epsilon = 0;//0.01;
	if (seedStats.sd_data != 0) {
		var_epsilon = seedStats.sd_data;
	}
	else {
		var_epsilon = parameters.eps;
	}
	std::cout << "Epsilon to be used : " << var_epsilon << std::endl;
	
	dataType seedValCT = value_first_pt;
	//////dataType seedValCT = (value_first_pt + value_second_pt) / 2.0;
	//////dataType seedValCT = ctImageData.imageDataPtr[(size_t)seedPoint[0].z][x_new((size_t)seedPoint[0].x, (size_t)seedPoint[0].y, length)];

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = fabs(seedValCT - ctImageData.imageDataPtr[k][i]);
		}
	}

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
	std::cout << "difference min : " << minImage << std::endl;
	std::cout << "difference max : " << maxImage << std::endl;

	//dataType maxRatio = 0.0;
	//dataType minRatio = INFINITY;
	//Normalization
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			//potential[k][i] = (var_epsilon + potential[k][i] / maxImage) * (1.0 / (1.0 + 1.0 * distance[k][i]));
			potential[k][i] = (var_epsilon + potential[k][i]) * (1.0 / (1.0 + 1.0 * distance[k][i]));
			//potential[k][i] = (var_epsilon + potential[k][i]) * (1.0 / (1.0 + 1.0 * distance[k][i]));
			//potential[k][i] = 1.0 / (1.0 + 1.0 * distance[k][i]);
			//potential[k][i] = parameters.eps + potential[k][i] / maxImage;
			//dataType normFactor = potential[k][i] / maxImage;
			//if(normFactor > maxRatio) {
			//	maxRatio = normFactor;
			//}
			//if(normFactor < minRatio) {
			//	minRatio = normFactor;
			//}
		}
	}
	//std::cout << "min factor : " << minRatio << std::endl;
	//std::cout << "max factor : " << maxRatio << std::endl;

	for (k = 0; k < height; k++) {
		delete[] maskThreshold[k];
		delete[] distance[k];
	}
	delete[] maskThreshold;
	delete[] distance;

	return true;
}

bool fastMarching3D_N(dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D seedPoint) {

	if (distanceFuncPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	size_t** labelArray = new size_t * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new size_t[dim2D]{0};
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			distanceFuncPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the initial point
	pointFastMarching3D current;
	i = (size_t)seedPoint.x;
	j = (size_t)seedPoint.y;
	k = (size_t)seedPoint.z;
	if (i < 0 || i > length || j < 0 || j > width || k < 0 || k > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	currentIndx = x_new(i, j, length);
	distanceFuncPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dTop = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		distanceFuncPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialFuncPtr[k][indxEast];
		dEast = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		distanceFuncPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialFuncPtr[k][indxNorth];
		dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		distanceFuncPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialFuncPtr[k][indxWest];
		dWest = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		distanceFuncPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialFuncPtr[k][indxSouth];
		dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		distanceFuncPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		distanceFuncPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;

	while (inProcess.size() != 0) {

		//processed the point with minimum distance
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		distanceFuncPtr[k][currentIndx] = current.arrival;
		deleteRootHeap3D(inProcess);

		//processed neighbors of the minimum in the narrow band
		
		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			kminus = k - 1;
			label = labelArray[kminus][currentIndx];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kminus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kminus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kminus);
				coefSpeed = potentialFuncPtr[kminus][currentIndx];
				dTop = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[kminus][currentIndx] = dTop;
					labelArray[kminus][currentIndx] = 2;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (label == 2) {
						if (dTop < distanceFuncPtr[kminus][currentIndx]) {
							distanceFuncPtr[kminus][currentIndx] = dTop;
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(iplus, j, length);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialFuncPtr[k][indxEast];
				dEast = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxEast] = dEast;
					labelArray[k][indxEast] = 2;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distanceFuncPtr[k][indxEast]) {
							distanceFuncPtr[k][indxEast] = dEast;
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(i, jminus, length);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialFuncPtr[k][indxNorth];
				dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxNorth] = dNorth;
					labelArray[k][indxNorth] = 2;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distanceFuncPtr[k][indxNorth]) {
							distanceFuncPtr[k][indxNorth] = dNorth;
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(iminus, j, length);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
				y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
				z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialFuncPtr[k][indxWest];
				dWest = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxWest] = dWest;
					labelArray[k][indxWest] = 2;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distanceFuncPtr[k][indxWest]) {
							distanceFuncPtr[k][indxWest] = dWest;
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(i, jplus, length);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
				y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
				z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialFuncPtr[k][indxSouth];
				dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[k][indxSouth] = dSouth;
					labelArray[k][indxSouth] = 2;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distanceFuncPtr[k][indxSouth]) {
							distanceFuncPtr[k][indxSouth] = dSouth;
						}
					}
				}
			}
		}

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
				y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
				z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialFuncPtr[kplus][currentIndx];
				dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					distanceFuncPtr[kplus][currentIndx] = dBottom;
					labelArray[kplus][currentIndx] = 2;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (label == 2) {
						if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
							distanceFuncPtr[kplus][currentIndx] = dBottom;
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

//========================================================================================================
//====================================== Front propagation with spacing ==================================
//========================================================================================================

dataType solve3dQuadraticEikonalEquation(dataType X, dataType Y, dataType Z, dataType P, VoxelSpacing h) {

	if(h.sx <= 0 || h.sy <= 0 || h.sz <= 0) {
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

	if(X == INFINITY && Y == INFINITY && Z == INFINITY) {
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
			if(solution >= Z) {
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
			if(solution >= max(X, Y)) {
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
			if(solution >= max(X, Z)) {
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
		a = (dataType)( hy_2 + hz_2 );
		b = (dataType)( -2 * (Y * hy_2 + Z * hz_2) );
		c = (dataType)( Y * Y * hy_2 + Z * Z * hz_2 - P_2 );
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if(solution >= max(Y, Z)) {
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
		a = (dataType)( hx_2 + hy_2 + hz_2 );
		b = (dataType)( -2 * (X * hx_2 + Y * hy_2 + Z * hz_2) );
		c = (dataType)( X * X * hx_2 + Y * Y * hy_2 + Z * Z * hz_2 - P_2 );
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if(solution >= max(X, max(Y, Z))) {
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

	//pointer to index of inProcess vector
	size_t** indexArray = new size_t * [length * width * height] {NULL};

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

	pointFastMarching3D current = { i, j, k, 0 };

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

		//Visualize the front propagation
		Point3D sPoint = { (dataType)i, (dataType)j, (dataType)k };
		sPoint = getRealCoordFromImageCoord3D(sPoint, actionPtr.origin, actionPtr.spacing, actionPtr.orientation);
		savingList.push_back(sPoint);
		//in real image save after every 10000 points
		if (processed_point % 10000 == 0) {
			id_save++;
			string saving_csv = path_saving + to_string(id_save) + ".csv";
			FILE* frontPoint;
			if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
				printf("Enable to open");
				return false;
			}
			fprintf(frontPoint, "x,y,z\n");
			for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
				fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
			}
			savingList.clear();
			fclose(frontPoint);
		}
		
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
	
	//Visualize the front propagation
	id_save++;
	string saving_csv = path_saving + to_string(id_save) + ".csv";
	FILE* frontPoint;
	if (fopen_s(&frontPoint, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(frontPoint, "x,y,z\n");
	if (savingList.size() > 0) {
		for (size_t i_n = 0; i_n < savingList.size(); i_n++) {
			fprintf(frontPoint, "%f,%f,%f\n", savingList[i_n].x, savingList[i_n].y, savingList[i_n].z);
		}
		savingList.clear();
	}
	fclose(frontPoint);
	
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

	delete[] indexArray;

	return true;

}

/*
bool fastMarching3dWithSpacing(Image_Data ctImageData, dataType** actionPtr, dataType** potentialFuncPtr, Point3D seedPoint) {

	if (actionPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> inProcess;
	inProcess.reserve(length * width * height); // Reserve space to avoid frequent reallocations
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

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

	pointFastMarching3D current = { i, j, k, 0, actionPtr[k][currentIndx] };

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;

	//Top
	if (k > 0 && k < height && i >= 0 && i < length && j >= 0 && j < width) {
		size_t kminus = k - 1;
		dataType x = select3dX(actionPtr, length, width, height, i, j, kminus);
		dataType y = select3dY(actionPtr, length, width, height, i, j, kminus);
		dataType z = select3dZ(actionPtr, length, width, height, i, j, kminus);
		dataType coefSpeed = potentialFuncPtr[kminus][currentIndx];
		dataType dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		size_t posTop = inProcess.size();
		actionPtr[kminus][currentIndx] = dTop;
		pointFastMarching3D TopNeighbor = { i, j, kminus, posTop, actionPtr[kminus][currentIndx] };
		inProcess.push_back(TopNeighbor);
		indexArray[x_flat(i, j, kminus, length, width)] = &inProcess[posTop].pos;
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
		size_t posBottom = inProcess.size();
		actionPtr[kplus][currentIndx] = dBottom;
		pointFastMarching3D BottomNeighbor = { i, j, kplus, posBottom, actionPtr[kplus][currentIndx] };
		inProcess.push_back(BottomNeighbor);
		indexArray[x_flat(i, j, kplus, length, width)] = &inProcess[posBottom].pos;
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
		size_t posNorth = inProcess.size();
		actionPtr[k][indxNorth] = dNorth;
		pointFastMarching3D NorthNeighbor = { i, jminus, k, posNorth, actionPtr[k][indxNorth] };
		inProcess.push_back(NorthNeighbor);
		indexArray[x_flat(i, jminus, k, length, width)] = &inProcess[posNorth].pos;
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
		pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth, actionPtr[k][indxSouth] };
		inProcess.push_back(SouthNeighbor);
		indexArray[x_flat(i, jplus, k, length, width)] = &inProcess[posSouth].pos;
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
		size_t posEast = inProcess.size();
		actionPtr[k][indxEast] = dEast;
		pointFastMarching3D EastNeighbor = { iplus, j, k, posEast, actionPtr[k][indxEast] };
		inProcess.push_back(EastNeighbor);
		indexArray[x_flat(iplus, j, k, length, width)] = &inProcess[posEast].pos;
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
		size_t posWest = inProcess.size();
		actionPtr[k][indxWest] = dWest;
		pointFastMarching3D WestNeighbor = { iminus, j, k, posWest, actionPtr[k][indxWest] };
		inProcess.push_back(WestNeighbor);
		indexArray[x_flat(iminus, j, k, length, width)] = &inProcess[posWest].pos;
		size_t lab = labelArray[k][indxWest];
	}

	heapifyVector3D(inProcess);

	while (inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0];
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
				size_t indx = x_flat(i, j, kminus, length, width);
				if (label == 3) {
					size_t posTop = inProcess.size();
					actionPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, posTop, actionPtr[kminus][currentIndx] };
					inProcess.push_back(TopNeighbor);
					indexArray[indx] = &inProcess[posTop].pos;
					heapifyUp3D(inProcess, posTop);
					labelArray[kminus][currentIndx] = 2;
				}
				else {
					if (dTop < actionPtr[kminus][currentIndx]) {
						actionPtr[kminus][currentIndx] = dTop;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dTop;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, j, kplus, length, width);
				if (label == 3) {
					size_t posBottom = inProcess.size();
					actionPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, posBottom, actionPtr[kplus][currentIndx] };
					inProcess.push_back(BottomNeighbor);
					indexArray[indx] = &inProcess[posBottom].pos;
					heapifyUp3D(inProcess, posBottom);
					labelArray[kplus][currentIndx] = 2;
				}
				else {
					if (dBottom < actionPtr[kplus][currentIndx]) {
						actionPtr[kplus][currentIndx] = dBottom;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dBottom;
						heapifyUp3D(inProcess, pt_pos);
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
				dataType x = select3dX(actionPtr, length, width, height, iplus, j, k);
				dataType y = select3dY(actionPtr, length, width, height, iplus, j, k);
				dataType z = select3dZ(actionPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = potentialFuncPtr[k][indxEast];
				dataType dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
				size_t indx = x_flat(iplus, j, k, length, width);
				if (label == 3) {
					size_t posEast = inProcess.size();
					actionPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, posEast, actionPtr[k][indxEast] };
					inProcess.push_back(EastNeighbor);
					indexArray[indx] = &inProcess[posEast].pos;
					heapifyUp3D(inProcess, posEast);
					labelArray[k][indxEast] = 2;
				}
				else {
					if (dEast < actionPtr[k][indxEast]) {
						actionPtr[k][indxEast] = dEast;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dEast;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(iminus, j, k, length, width);
				if (label == 3) {
					size_t posWest = inProcess.size();
					actionPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, posWest, actionPtr[k][indxWest] };
					inProcess.push_back(WestNeighbor);
					size_t indx = x_flat(iminus, j, k, length, width);
					indexArray[indx] = &inProcess[posWest].pos;
					heapifyUp3D(inProcess, posWest);
					labelArray[k][indxWest] = 2;
				}
				else {
					if (dWest < actionPtr[k][indxWest]) {
						actionPtr[k][indxWest] = dWest;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dWest;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, jminus, k, length, width);
				if (label == 3) {
					size_t posNorth = inProcess.size();
					actionPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, posNorth, actionPtr[k][indxNorth] };
					inProcess.push_back(NorthNeighbor);
					indexArray[indx] = &inProcess[posNorth].pos;
					heapifyUp3D(inProcess, posNorth);
					labelArray[k][indxNorth] = 2;
				}
				else {
					if (dNorth < actionPtr[k][indxNorth]) {
						actionPtr[k][indxNorth] = dNorth;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dNorth;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, jplus, k, length, width);
				if (label == 3) {
					size_t posSouth = inProcess.size();
					actionPtr[k][indxSouth] = dSouth;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth, actionPtr[k][indxSouth] };
					inProcess.push_back(SouthNeighbor);
					indexArray[indx] = &inProcess[posSouth].pos;
					heapifyUp3D(inProcess, posSouth);
					labelArray[k][indxSouth] = 2;
				}
				else {
					if (dSouth < actionPtr[k][indxSouth]) {
						actionPtr[k][indxSouth] = dSouth;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dSouth;
						heapifyUp3D(inProcess, pt_pos);
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
*/

/*
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
		if(labelArray[kEnd][x_new(iEnd, jEnd, length)] == 1) {
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
			if(i < length_minus && j >= 0 && j < width && k >= 0 && k < height)
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
					}
				}
			}
		}
		//====================

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
*/

//========================================================================================================
//================================== Distance Map ========================================================
//========================================================================================================

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

bool fastMarching3dForDistanceMap(Image_Data ctImageData, dataType** distanceFuncPtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distanceFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> inProcess;
	inProcess.reserve(length * width * height); // Reserve space to avoid frequent reallocations
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

dataType computeHausDorffDistance(dataType** source1, dataType** source2, const size_t length, const size_t width, const size_t height) {
	
	if (source1 == NULL || source2 == NULL) {
		std::cout << "Structures wrongly defined" << std::endl;
		return -1;
	}

	dataType hausdorff_distance_a = 0.0, min_value = INFINITY;
	for (size_t k_a = 0; k_a < height; k_a++) {
		for (size_t i_a = 0; i_a < length; i_a++) {
			for (size_t j_a = 0; j_a < width; j_a++) {
				dataType value_a = source1[k_a][x_new(i_a, j_a, length)];
				min_value = INFINITY;
				for (size_t k_b = 0; k_b < height; k_b++) {
					for (size_t i_b = 0; i_b < length; i_b++) {
						for (size_t j_b = 0; j_b < width; j_b++) {
							dataType value_b = source2[k_b][x_new(i_b, j_b, length)];
							dataType h_distance = abs(value_a - value_b);
							if (h_distance < min_value) {
								min_value = h_distance;
							}
						}
					}
				}
				if (hausdorff_distance_a < min_value) {
					hausdorff_distance_a = min_value;
				}
			}
		}
	}

	dataType hausdorff_distance_b = 0.0;
	for (size_t k_b = 0; k_b < height; k_b++) {
		for (size_t i_b = 0; i_b < length; i_b++) {
			for (size_t j_b = 0; j_b < width; j_b++) {
				dataType value_b = source2[k_b][x_new(i_b, j_b, length)];
				min_value = INFINITY;
				for (size_t k_a = 0; k_a < height; k_a++) {
					for (size_t i_a = 0; i_a < length; i_a++) {
						for (size_t j_a = 0; j_a < width; j_a++) {
							dataType value_a = source1[k_a][x_new(i_a, j_a, length)];
							dataType h_distance = abs(value_a - value_b);
							if (h_distance < min_value) {
								min_value = h_distance;
							}
						}
					}
				}
				if (hausdorff_distance_b < min_value) {
					hausdorff_distance_b = min_value;
				}
			}
		}
	}
	return (dataType)(max(hausdorff_distance_a, hausdorff_distance_b));
}

//========================================================================================================
//================================= Front propagation and Path Extraction ================================
//========================================================================================================

/*
* Point3D nextKeyPointDetection(Image_Data actionMapStr, dataType** potentialFuncPtr, Point3D seedPoint, const double LengthKeyPoints) {

	if (actionMapStr.imageDataPtr == NULL || potentialFuncPtr == NULL) {
		std::cout << "Wrong input";
		return Point3D{ -1.0, -1.0, -1.0 };
	}

	const size_t length = actionMapStr.length;
	const size_t width = actionMapStr.width;
	const size_t height = actionMapStr.height;
	VoxelSpacing spacing = actionMapStr.spacing;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	size_t** labelArray = new size_t * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new size_t[dim2D]{ 0 };
	}
	if (labelArray == NULL) {
		std::cout << "Label array non defined";
		return Point3D{ -1.0, -1.0, -1.0 };
	}

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			actionMapStr.imageDataPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the initial point
	pointFastMarching3D current;
	if (seedPoint.x < 0 || seedPoint.y < 0 || seedPoint.z < 0) {
		std::cout << "Error in the input seed point" << std::endl;
		return Point3D{ -1, -1, -1 };
	}
	else {
		i = (size_t)seedPoint.x;
		j = (size_t)seedPoint.y;
		k = (size_t)seedPoint.z;
		if (i > length || j > width || k > height) {
			std::cout << "Error in the input seed point" << std::endl;
			return Point3D{ -1, -1, -1 };
		}
	}
	currentIndx = x_new(i, j, length);
	actionMapStr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	double distanceToCurrentSourcePoint = 0.0;
	Point3D currentSourcePoint = { i, j, k };

	vector <Point3D> key_points;
	//key_points.push_back(currentSourcePoint);

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialFuncPtr[kminus][currentIndx];
		//dTop = solve3dQuadratic(x, y, z, coefSpeed);
		dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialFuncPtr[k][indxEast];
		//dEast = solve3dQuadratic(x, y, z, coefSpeed);
		dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		actionMapStr.imageDataPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialFuncPtr[k][indxNorth];
		//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialFuncPtr[k][indxWest];
		//dWest = solve3dQuadratic(x, y, z, coefSpeed);
		dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		actionMapStr.imageDataPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialFuncPtr[k][indxSouth];
		//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialFuncPtr[kplus][currentIndx];
		//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;
	size_t nbSourcePoint = 0;
	dataType max_weighted_distance = 0.0;

	//size_t iEnd = (size_t)seedPoint[1].x;
	//size_t jEnd = (size_t)seedPoint[1].y;
	//size_t kEnd = (size_t)seedPoint[1].z;

	while (inProcess.size() != 0 && distanceToCurrentSourcePoint < LengthKeyPoints) {

		//processed the point with minimum distance
		current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		dataType arrival_for_max = current.arrival;

		Point3D pSource = { i, j, k };
		Point3D pSourceReal = getRealCoordFromImageCoord3D(pSource, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		Point3D pCurrent = getRealCoordFromImageCoord3D(currentSourcePoint, actionMapStr.origin, actionMapStr.spacing, actionMapStr.orientation);
		distanceToCurrentSourcePoint = getPoint3DDistance(pCurrent, pSourceReal);

		if (distanceToCurrentSourcePoint >= LengthKeyPoints) {

			currentSourcePoint = pSource;
			key_points.push_back(currentSourcePoint);
			actionMapStr.imageDataPtr[k][currentIndx] = 0;

			//Top
			if (k > 0 && j >= 0 && j < width && i >= 0 && i < length)
			{
				kminus = k - 1;
				actionMapStr.imageDataPtr[kminus][currentIndx] = INFINITY;
				labelArray[kminus][currentIndx] = 3;
			}

			//East
			if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height)
			{
				iplus = i + 1;
				indxEast = x_new(iplus, j, length);
				actionMapStr.imageDataPtr[k][indxEast] = INFINITY;
				labelArray[k][indxEast] = 3;
			}

			//North
			if (j > 0 && i >= 0 && i < length && k >= 0 && k < height)
			{
				jminus = j - 1;
				indxNorth = x_new(i, jminus, length);
				actionMapStr.imageDataPtr[k][indxNorth] = INFINITY;
				labelArray[k][indxNorth] = 3;
			}

			//West
			if (i > 0 && j >= 0 && j < width && k >= 0 && k < height)
			{
				iminus = i - 1;
				indxWest = x_new(iminus, j, length);
				actionMapStr.imageDataPtr[k][indxWest] = INFINITY;
				labelArray[k][indxWest] = 3;
			}

			//South
			if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height)
			{
				jplus = j + 1;
				indxSouth = x_new(i, jplus, length);
				actionMapStr.imageDataPtr[k][indxSouth] = INFINITY;
				labelArray[k][indxSouth] = 3;
			}

			//Bottom
			if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width)
			{
				kplus = k + 1;
				actionMapStr.imageDataPtr[kplus][currentIndx] = INFINITY;
				labelArray[kplus][currentIndx] = 3;
			}

			while (inProcess.size() != 0) {
				inProcess.pop_back();
			}
		}
		else {
			actionMapStr.imageDataPtr[k][currentIndx] = current.arrival;
			deleteRootHeap3D(inProcess);
		}

		if (arrival_for_max > max_weighted_distance) {
			max_weighted_distance = arrival_for_max;
		}

		//Continue the propagation only if key point is not detected
		if (key_points.size() == 0) {

			//processed neighbors of the minimum in the narrow band
			//Top
			if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
				kminus = k - 1;
				label = labelArray[kminus][currentIndx];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kminus);
					coefSpeed = potentialFuncPtr[kminus][currentIndx];
					//dTop = solve3dQuadratic(x, y, z, coefSpeed);
					dTop = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[kminus][currentIndx] = 2;
						actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
						pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
						addPointHeap3D(inProcess, TopNeighbor);
					}
					else {
						if (label == 2) {
							if (dTop < actionMapStr.imageDataPtr[kminus][currentIndx]) {
								actionMapStr.imageDataPtr[kminus][currentIndx] = dTop;
							}
						}
					}
				}
			}

			//East
			if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
				iplus = i + 1;
				indxEast = x_new(iplus, j, length);
				label = labelArray[k][indxEast];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iplus, j, k);
					coefSpeed = potentialFuncPtr[k][indxEast];
					//dEast = solve3dQuadratic(x, y, z, coefSpeed);
					dEast = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[k][indxEast] = 2;
						actionMapStr.imageDataPtr[k][indxEast] = dEast;
						pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
						addPointHeap3D(inProcess, EastNeighbor);

					}
					else {
						if (label == 2) {
							if (dEast < actionMapStr.imageDataPtr[k][indxEast]) {
								actionMapStr.imageDataPtr[k][indxEast] = dEast;
							}
						}
					}
				}
			}

			//North
			if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
				jminus = j - 1;
				indxNorth = x_new(i, jminus, length);
				label = labelArray[k][indxNorth];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jminus, k);
					coefSpeed = potentialFuncPtr[k][indxNorth];
					//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
					dNorth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[k][indxNorth] = 2;
						actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
						pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
						addPointHeap3D(inProcess, NorthNeighbor);
					}
					else {
						if (label == 2) {
							if (dNorth < actionMapStr.imageDataPtr[k][indxNorth]) {
								actionMapStr.imageDataPtr[k][indxNorth] = dNorth;
							}
						}
					}
				}
			}

			//West
			if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
				iminus = i - 1;
				indxWest = x_new(iminus, j, length);
				label = labelArray[k][indxWest];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, iminus, j, k);
					coefSpeed = potentialFuncPtr[k][indxWest];
					//dWest = solve3dQuadratic(x, y, z, coefSpeed);
					dWest = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[k][indxWest] = 2;
						actionMapStr.imageDataPtr[k][indxWest] = dWest;
						pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
						addPointHeap3D(inProcess, WestNeighbor);
					}
					else {
						if (label == 2) {
							if (dWest < actionMapStr.imageDataPtr[k][indxWest]) {
								actionMapStr.imageDataPtr[k][indxWest] = dWest;
							}
						}
					}
				}
			}

			//South
			if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
				jplus = j + 1;
				indxSouth = x_new(i, jplus, length);
				label = labelArray[k][indxSouth];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, jplus, k);
					coefSpeed = potentialFuncPtr[k][indxSouth];
					//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
					dSouth = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[k][indxSouth] = 2;
						actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
						pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
						addPointHeap3D(inProcess, SouthNeighbor);
					}
					else {
						if (label == 2) {
							if (dSouth < actionMapStr.imageDataPtr[k][indxSouth]) {
								actionMapStr.imageDataPtr[k][indxSouth] = dSouth;
							}
						}
					}
				}
			}

			//Bottom
			if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
				kplus = k + 1;
				label = labelArray[kplus][currentIndx];
				if (label != 1) {
					x = select3dX(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
					y = select3dY(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
					z = select3dZ(actionMapStr.imageDataPtr, length, width, height, i, j, kplus);
					coefSpeed = potentialFuncPtr[kplus][currentIndx];
					//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
					dBottom = solve3dQuadraticEikonalEquation(x, y, z, coefSpeed, spacing);
					if (label == 3) {
						labelArray[kplus][currentIndx] = 2;
						actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
						pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
						addPointHeap3D(inProcess, BottomNeighbor);
					}
					else {
						if (label == 2) {
							if (dBottom < actionMapStr.imageDataPtr[kplus][currentIndx]) {
								actionMapStr.imageDataPtr[kplus][currentIndx] = dBottom;
							}
						}
					}
				}
			}

		}

	}

	////key_points.push_back(seedPoint[1]);
	//std::cout << "The source Points are: " << std::endl;
	//for (int nb = 0; nb < key_points.size(); nb++) {
	//	std::cout << "Point " << nb + 1 << " : " << key_points[nb].x << " " << key_points[nb].y << " " << key_points[nb].z << std::endl;
	//}

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (actionMapStr.imageDataPtr[k][i] == INFINITY) {
				actionMapStr.imageDataPtr[k][i] = max_weighted_distance + 1;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return key_points[0];
}
* bool penalizedFrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialFuncPtr, Point3D* seedPoints) {

	if (inputImageData.imageDataPtr == NULL || actionMapPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
		return false;
	}

	const size_t height = inputImageData.height;
	const size_t length = inputImageData.length;
	const size_t width = inputImageData.width;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D] {0};
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			actionMapPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the initial point
	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	k = (size_t)seedPoints[0].z;
	if (i < 0 || i > length || j < 0 || j > width || k < 0 || k > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	currentIndx = x_new(i, j, length);
	actionMapPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;
	pointFastMarching3D current = { i, j, k, 0.0 };

	double radius = 3.0;
	Image_Data potentialStr = { height, length, width, potentialFuncPtr, inputImageData.origin, inputImageData.spacing, inputImageData.orientation };
	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };
	//seedStats = getStats(potentialStr, seedPoints[0], radius);
	seedStats = getPointNeighborhoodStats(potentialStr, seedPoints[0], radius);
	dataType thres_min = seedStats.mean_data - seedStats.sd_data;
	dataType thres_max = seedStats.mean_data + seedStats.sd_data;

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(actionMapPtr, length, width, height, i, j, kminus);
		y = select3dY(actionMapPtr, length, width, height, i, j, kminus);
		z = select3dZ(actionMapPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialFuncPtr[kminus][currentIndx];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dTop = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
			actionMapPtr[kminus][currentIndx] = dTop;
			inProcess.push_back(TopNeighbor);
			labelArray[kminus][currentIndx] = 2;
		}
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(actionMapPtr, length, width, height, iplus, j, k);
		y = select3dY(actionMapPtr, length, width, height, iplus, j, k);
		z = select3dZ(actionMapPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialFuncPtr[k][indxEast];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dEast = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
			actionMapPtr[k][indxEast] = dEast;
			inProcess.push_back(EastNeighbor);
			labelArray[k][indxEast] = 2;
		}
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(actionMapPtr, length, width, height, i, jminus, k);
		y = select3dY(actionMapPtr, length, width, height, i, jminus, k);
		z = select3dZ(actionMapPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialFuncPtr[k][indxNorth];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dNorth = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
			actionMapPtr[k][indxNorth] = dNorth;
			inProcess.push_back(NorthNeighbor);
			labelArray[k][indxNorth] = 2;
		}
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(actionMapPtr, length, width, height, iminus, j, k);
		y = select3dY(actionMapPtr, length, width, height, iminus, j, k);
		z = select3dZ(actionMapPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialFuncPtr[k][indxWest];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dWest = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
			actionMapPtr[k][indxWest] = dWest;
			inProcess.push_back(WestNeighbor);
			labelArray[k][indxWest] = 2;
		}
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(actionMapPtr, length, width, height, i, jplus, k);
		y = select3dY(actionMapPtr, length, width, height, i, jplus, k);
		z = select3dZ(actionMapPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialFuncPtr[k][indxSouth];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dSouth = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
			actionMapPtr[k][indxSouth] = dSouth;
			inProcess.push_back(SouthNeighbor);
			labelArray[k][indxSouth] = 2;
		}
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(actionMapPtr, length, width, height, i, j, kplus);
		y = select3dY(actionMapPtr, length, width, height, i, j, kplus);
		z = select3dZ(actionMapPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialFuncPtr[kplus][currentIndx];
		if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
			dBottom = solve3dQuadratic(x, y, z, coefSpeed);
			pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
			actionMapPtr[kplus][currentIndx] = dBottom;
			inProcess.push_back(BottomNeighbor);
			labelArray[kplus][currentIndx] = 2;
		}
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	size_t seedI = (size_t)seedPoints[1].x;
	size_t seedJ = (size_t)seedPoints[1].y;
	size_t seedK = (size_t)seedPoints[1].z;
	size_t seedIndex = x_new(seedI, seedJ, length);
	if (seedI < 0 || seedI > length || seedJ < 0 || seedJ > width || seedK < 0 || seedK > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}

	dataType max_action = 0;
	while (labelArray[seedK][seedIndex] != 1 && inProcess.size() > 0) {

		//processed the point with minimum distance
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		i = current.x;
		j = current.y;
		k = current.z;
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		actionMapPtr[k][currentIndx] = current.arrival;
		max_action = current.arrival;

		deleteRootHeap3D(inProcess);

		//processed neighbors of the minimum in the narrow band

		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			kminus = k - 1;
			label = labelArray[kminus][currentIndx];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, j, kminus);
				y = select3dY(actionMapPtr, length, width, height, i, j, kminus);
				z = select3dZ(actionMapPtr, length, width, height, i, j, kminus);
				coefSpeed = potentialFuncPtr[kminus][currentIndx];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dTop = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[kminus][currentIndx] = dTop;
						labelArray[kminus][currentIndx] = 2;
						pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
						addPointHeap3D(inProcess, TopNeighbor);
					}
					else {
						if (label == 2) {
							if (dTop < actionMapPtr[kminus][currentIndx]) {
								actionMapPtr[kminus][currentIndx] = dTop;
							}
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(iplus, j, length);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, iplus, j, k);
				y = select3dY(actionMapPtr, length, width, height, iplus, j, k);
				z = select3dZ(actionMapPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialFuncPtr[k][indxEast];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dEast = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[k][indxEast] = dEast;
						labelArray[k][indxEast] = 2;
						pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
						addPointHeap3D(inProcess, EastNeighbor);
					}
					else {
						if (label == 2) {
							if (dEast < actionMapPtr[k][indxEast]) {
								actionMapPtr[k][indxEast] = dEast;
							}
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(i, jminus, length);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, jminus, k);
				y = select3dY(actionMapPtr, length, width, height, i, jminus, k);
				z = select3dZ(actionMapPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialFuncPtr[k][indxNorth];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dNorth = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[k][indxNorth] = dNorth;
						labelArray[k][indxNorth] = 2;
						pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
						addPointHeap3D(inProcess, NorthNeighbor);
					}
					else {
						if (label == 2) {
							if (dNorth < actionMapPtr[k][indxNorth]) {
								actionMapPtr[k][indxNorth] = dNorth;
							}
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(iminus, j, length);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, iminus, j, k);
				y = select3dY(actionMapPtr, length, width, height, iminus, j, k);
				z = select3dZ(actionMapPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialFuncPtr[k][indxWest];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dWest = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[k][indxWest] = dWest;
						labelArray[k][indxWest] = 2;
						pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
						addPointHeap3D(inProcess, WestNeighbor);
					}
					else {
						if (label == 2) {
							if (dWest < actionMapPtr[k][indxWest]) {
								actionMapPtr[k][indxWest] = dWest;
							}
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(i, jplus, length);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, jplus, k);
				y = select3dY(actionMapPtr, length, width, height, i, jplus, k);
				z = select3dZ(actionMapPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialFuncPtr[k][indxSouth];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dSouth = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[k][indxSouth] = dSouth;
						labelArray[k][indxSouth] = 2;
						pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
						addPointHeap3D(inProcess, SouthNeighbor);
					}
					else {
						if (label == 2) {
							if (dSouth < actionMapPtr[k][indxSouth]) {
								actionMapPtr[k][indxSouth] = dSouth;
							}
						}
					}
				}
			}
		}

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, j, kplus);
				y = select3dY(actionMapPtr, length, width, height, i, j, kplus);
				z = select3dZ(actionMapPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialFuncPtr[kplus][currentIndx];
				if (coefSpeed >= thres_min && coefSpeed <= thres_max) {
					dBottom = solve3dQuadratic(x, y, z, coefSpeed);
					if (label == 3) {
						actionMapPtr[kplus][currentIndx] = dBottom;
						labelArray[kplus][currentIndx] = 2;
						pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
						addPointHeap3D(inProcess, BottomNeighbor);
					}
					else {
						if (label == 2) {
							if (dBottom < actionMapPtr[kplus][currentIndx]) {
								actionMapPtr[kplus][currentIndx] = dBottom;
							}
						}
					}
				}
			}
		}
	}

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (actionMapPtr[k][i] > max_action) {
				actionMapPtr[k][i] = max_action + 1;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}
* dataType solve3dQuadraticFastSweeping(dataType X, dataType Y, dataType Z, dataType P, VoxelSpacing h) {

	if (h.sx <= 0 || h.sy <= 0 || h.sz <= 0) {
		std::cout << "Error: Voxel spacing must be positive." << std::endl;
		return INFINITY; // Return a large value or handle the error appropriately
	}

	dataType sol = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P_2 = P * P;
	dataType hx_2 = 1.0 / (h.sx * h.sx);
	dataType hy_2 = 1.0 / (h.sy * h.sy);
	dataType hz_2 = 1.0 / (h.sz * h.sz);

	if (X != INFINITY && Y == INFINITY && Z == INFINITY) {
		a = hx_2;
		b = (dataType)(-2 * X * hx_2);
		c = (dataType)(X * X * hx_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
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
			return (dataType)((-b + sqrt(delta)) / (2 * a));
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
			return (dataType)((-b + sqrt(delta)) / (2 * a));
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
			sol = (dataType)((-b + sqrt(delta)) / (2 * a));
			return sol;
		}
		else {
			return (dataType)(min((X + h.sx * P), (Y + h.sy * P)));
		}
	}

	if (X != INFINITY && Z != INFINITY && Y == INFINITY) {
		a = (dataType)(hx_2 + hz_2);
		b = (dataType)(-2 * (X * hx_2 + Z * hz_2));
		c = (dataType)(X * X * hx_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol = (dataType)((-b + sqrt(delta)) / (2 * a));
			return sol;
		}
		else {
			return (dataType)(min((X + h.sx * P), (Z + h.sz * P)));
		}
	}

	if (Y != INFINITY && Z != INFINITY && X == INFINITY) {
		a = (dataType)(hy_2 + hz_2);
		b = (dataType)(-2 * (Y * hy_2 + Z * hz_2));
		c = (dataType)(Y * Y * hy_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol = (dataType)((-b + sqrt(delta)) / (2 * a));
			return sol;
		}
		else {
			return (dataType)(min((Y + h.sy * P), (Z + h.sz * P)));
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		
		a = (dataType)(hx_2 + hy_2 + hz_2);
		b = (dataType)(-2 * (X * hx_2 + Y * hy_2 + Z * hz_2));
		c = (dataType)(X * X * hx_2 + Y * Y * hy_2 + Z * Z * hz_2 - P_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol = (dataType)((-b + sqrt(delta)) / (2 * a));
			return sol;
		}
		else {
			return (dataType)(min((X + h.sx * P), min((Y + h.sy * P), (Z + h.sz * P))));
		}
		
	}

}
bool findPathFromOneGivenPointWithCircleDetection(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters) {

	//if (ctImageData.imageDataPtr == NULL)
	//	return false;
	//size_t k = 0, i = 0, j = 0, x = 0;
	//const size_t height = ctImageData.height;
	//const size_t length = ctImageData.length;
	//const size_t width = ctImageData.width;
	//const size_t dim2D = length * width;
	//double offset = 5.0;
	//BoundingBox2D box = { 0, 0, 0, 0 };
	//string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//string saving_name, extension, slice_number;

	//saving_name = path_name + "path_points.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	//dataType** actionPtr = new dataType * [height];
	//dataType** maskAction = new dataType * [height];
	//dataType** potentialPtr = new dataType * [height];
	//dataType** edgeImage = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	actionPtr[k] = new dataType[dim2D]{ 0 };
	//	maskAction[k] = new dataType[dim2D]{ 0 };
	//	potentialPtr[k] = new dataType[dim2D]{ 0 };
	//	edgeImage[k] = new dataType[dim2D]{ 0 };
	//}
	//if (actionPtr == NULL || maskAction == NULL || potentialPtr == NULL)
	//	return false;
	//compute3DPotential(ctImageData, potentialPtr, seedPoints, parameters);
	//saving_name = path_name + "potential.raw";
	//manageRAWFile3D<dataType>(potentialPtr, length, width, height, saving_name.c_str(), STORE_DATA, false);
	//const double step = 70.0;
	
	////================ find next point inside the aorta ============
	//Point3D* seeds = new Point3D[2];
	//seeds[0] = seedPoints[0];

	//Image_Data action
	//{
	//	height,
	//	length,
	//	width,
	//	actionPtr,
	//	ctImageData.origin,
	//	ctImageData.spacing,
	//	ctImageData.orientation
	//};
	//dataType tau = 0.95 * ctImageData.spacing.sx;
	//dataType tolerance = ctImageData.spacing.sx;

	////During the first propagation the mask is empty
	////partialPropagation(action, potentialPtr, maskAction, length, width, height, seeds, step);
	////fastMarching3D_N(actionPtr, potentialPtr, length, width, height, initial_point_img_coord);
	//seeds[1] = nextKeyPointDetection(action, potentialPtr, seeds[0], step);
	//partialFrontPropagation(actionPtr, potentialPtr, length, width, height, seeds);
	//saving_name = path_name + "action_initial.raw";
	//manageRAWFile3D<dataType>(actionPtr, length, width, height, saving_name.c_str(), STORE_DATA, false);
	////vector for path points
	//vector<Point3D> path_points;
	////path between initial point and second point (found point)
	//shortestPath3D(action, seeds, tau, tolerance, path_points);
	////======================= Circle detection ==========================
	////Search for circle around the found points
	//dataType* imageSlice = new dataType[dim2D]{ 0 };
	//dataType* houghSpace = new dataType[dim2D]{ 0 };

	//Point2D seed2D = { 0.0, 0.0 };
	//Point2D center2D = { 0.0, 0.0 };
	//PixelSpacing spacing = { ctImageData.spacing.sx, ctImageData.spacing.sy };
	//Image_Data2D IMAGE = {length, width, imageSlice, center2D, spacing};
	//HoughParameters hParameters = {
	//	5.0,   //minimal radius
	//	19.0,  //maximal radius
	//	0.5,   //radius step
	//	0.5,   //epsilon
	//	3,   //offset, 50 without defined seed
	//	1000,  //K edge detector coefficient
	//	0.25,  //threshold edge detector
	//};
	//string loading_path = "C:/Users/Konan Allaly/Documents/Tests/input/raw/edge image/edge_image_p1.raw";
	//string hough_path = "C:/Users/Konan Allaly/Documents/Tests/output/threshold.raw";
	//bool detected = false;
	//size_t count_detected = 0;

	//for (int n = 0; n < path_points.size() - 1; n--) {
	//	k = (size_t)path_points[n].z;
	//	Point2D seed_circle = { path_points[n].x, path_points[n].y };
	//	copyDataToAnother2dArray(edgeImage[k], imageSlice, length, width);
	//	//detected = isCircleDetected(IMAGE, houghSpace, seed_circle, hParameters, hough_path);// Add new path !!!
	//	if (detected == false) {
	//		count_detected++;
	//	}
	//}
	//std::cout << count_detected << " circles not detected" << std::endl;
	////===================================================================

	////write path points coordinates to file
	////int n = 0;
	//Point3D point_file = { 0.0, 0.0, 0.0 };
	//for (int n = path_points.size() - 1; n > -1; n--) {
	//	point_file = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	//}
	//
	////empty the vector contaning path points
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}
	//fclose(file);

	////Update the mask
	//copyDataToAnotherArray(actionPtr, maskAction, height, length, width);

	////Second segment
	//seeds[0] = seeds[1];
	////during the second propagation the mask is used to avaoid computation in already computed area
	//partialPropagation(actionPtr, potentialPtr, maskAction, length, width, height, seeds, step);

	//saving_name = path_name + "action_second.raw";
	//manageRAWFile3D<dataType>(actionPtr, length, width, height, saving_name.c_str(), STORE_DATA, false);

	////path between initial point and third point (found point)
	//shortestPath3D(actionPtr, length, width, height, ctImageData.spacing, seeds, path_points);

	//for (n = path_points.size() - 1; n > -1; n--) {
	//	point_file = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	//}

	////empty the vector contaning path points
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}

	//size_t num_slice_processed = 0;
	//int len = path_points.size();
	//std::cout << len << " point to be processed" << std::endl;
	//Point3D next_point = { 0.0, 0.0, 0.0 };
	//size_t m;
	//for (n = len - 1; n > -1; n--) {
	//	num_slice_processed++;
	//	extension = to_string(num_slice_processed);
	//	
	//	//slice number
	//	m = (size_t)path_points[n].z;
	//	slice_number = to_string(m);
	//	
	//	copyDataToAnother2dArray(ctImageData.imageDataPtr[m], imageSlice, length, width);
	//	////Slice rescaling
	//	//rescaleNewRange2D(imageSlice, length, width, 0.0, 1.0);
	//	////Slice filtering
	//	//heatImplicit2dScheme(IMAGE, filter_parameters);
	//	
	//	seed2D.x = path_points[n].x;
	//	seed2D.y = path_points[n].y;
	//	
	//	if (num_slice_processed < 10) {
	//		saving_name = path_name + "path & circle detection/threshold/threshold_000" + extension + ".raw";
	//	}
	//	else {
	//		if (num_slice_processed < 100) {
	//			saving_name = path_name + "path & circle detection/threshold/threshold_00" + extension + ".raw";
	//		}
	//		else {
	//			if (num_slice_processed < 1000) {
	//				saving_name = path_name + "path & circle detection/threshold/threshold_0" + extension + ".raw";
	//			}
	//			else {
	//				saving_name = path_name + "path & circle detection/threshold/threshold_" + extension + ".raw";
	//			}
	//		}
	//	}
	//	center2D = localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, length, width, hParameters, saving_name);
	//	
	//	max_ratio = getTheMaxValue(voteArray, length, width);
	//	if (max_ratio != 0.0) {
	//		next_point.x = center2D.x;
	//		next_point.y = center2D.y;
	//		next_point.z = m;
	//	}
	//	if (num_slice_processed < 10) {
	//		saving_name = path_name + "path & circle detection/slice/slice_000" + extension + ".raw";
	//		store2dRawData(imageSlice, length, width, saving_name.c_str());
	//		saving_name = path_name + "path & circle detection/circle/found_circle_000" + extension + ".raw";
	//		store2dRawData(houghSpace, length, width, saving_name.c_str());
	//		saving_name = path_name + "path & circle detection/ratio/ratio_000" + extension + ".raw";
	//		store2dRawData(voteArray, length, width, saving_name.c_str());
	//	}
	//	else {
	//		if (num_slice_processed < 100) {
	//			saving_name = path_name + "path & circle detection/slice/slice_00" + extension + ".raw";
	//			store2dRawData(imageSlice, length, width, saving_name.c_str());
	//			saving_name = path_name + "path & circle detection/circle/found_circle_00" + extension + ".raw";
	//			store2dRawData(houghSpace, length, width, saving_name.c_str());
	//			saving_name = path_name + "path & circle detection/ratio/ratio_00" + extension + ".raw";
	//			store2dRawData(voteArray, length, width, saving_name.c_str());
	//		}
	//		else {
	//			if (num_slice_processed < 1000) {
	//				saving_name = path_name + "path & circle detection/slice/slice_0" + extension + ".raw";
	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//				saving_name = path_name + "path & circle detection/circle/found_circle_0" + extension + ".raw";
	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//				saving_name = path_name + "path & circle detection/ratio/ratio_0" + extension + ".raw";
	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//			}
	//			else {
	//				saving_name = path_name + "path & circle detection/slice/slice_" + extension + ".raw";
	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//				saving_name = path_name + "path & circle detection/circle/found_circle_" + extension + ".raw";
	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//				saving_name = path_name + "path & circle detection/ratio/ratio_" + extension + ".raw";
	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//			}
	//		}
	//	}
	//}

	////empty the list of point before new piece of path
	//while (path_points.size() != 0) {
	//	path_points.pop_back();
	//}

	//============================================

	////size_t count_stop = 0, max_steps = 10;
	////max_ratio = 1.0;
	////size_t stop_k = (size_t)seedPoints[1].z, verif_k = (size_t)seeds[1].z;

	//
	//Point3D last_point = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//Point3D verif_point = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//double verif_dist = getPoint3DDistance(last_point, verif_point);

	//std::cout << "distance to end point : " << verif_dist << std::endl;
	//while ((stop_k != verif_k && verif_dist < 12.5) || count_step < max_steps) {
	//	//max_ratio != 0.0 && count_step < max_steps //--> old condition
	//	std::cout << "STEP " << count_step << " : " << std::endl;
	//	extension = to_string(count_step);
	//	
	//	//mask
	//	for (k = 0; k < height; k++) {
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				x = x_new(i, j, length);
	//				current_point.x = (dataType)i;
	//				current_point.y = (dataType)j;
	//				current_point.z = (dataType)k;
	//				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//				dist = getPoint3DDistance(real_world, current_point);
	//				//dist = getPoint3DDistance(initial_point, current_point);
	//				if (dist <= (step + 0.01) && actionPtr[k][x] <= value_temp) {
	//					maskDistance[k][x] = 10000.0;
	//				}
	//				else {
	//					if (maskDistance[k][x] != 10000.0) {
	//						maskDistance[k][x] = newActionPtr[k][x];
	//					}
	//				}
	//			}
	//		}
	//	}
	//	seeds[0] = seeds[1];
	//	initial_point = seeds[1];
	//	
	//	fastMarching3D_N(newActionPtr, potentialPtr, length, width, height, initial_point);

	//	//saving_name = path_name + "action_field" + extension + ".raw";
	//	//store3dRawData<dataType>(newActionPtr, length, width, height, saving_name.c_str());

	//	//Find the temporary point
	//	min_distance = BIG_VALUE;
	//	
	//	//get real world coordinate of the initial point
	//	real_world = getRealCoordFromImageCoord3D(initial_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	
	//	for (k = 0; k < height; k++) {
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				x = x_new(i, j, length);
	//				current_point.x = (dataType)i;
	//				current_point.y = (dataType)j;
	//				current_point.z = (dataType)k;
	//				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//				dist = getPoint3DDistance(real_world, current_point);
	//				//dist = getPoint3DDistance(initial_point, current_point);
	//				if (dist >= (step - 0.01) && dist <= (step + 0.01) && maskDistance[k][x] != 10000.0) {
	//					if (newActionPtr[k][x] < min_distance) {
	//						min_distance = newActionPtr[k][x];
	//						temporary_point.x = (dataType)i;
	//						temporary_point.y = (dataType)j;
	//						temporary_point.z = (dataType)k;
	//						value_temp = newActionPtr[k][x];
	//					}
	//				}
	//			}
	//		}
	//	}
	//	//path between points
	//	std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;
	//	
	//	seeds[1] = temporary_point;
	//	shortestPath3D(newActionPtr, resultedPath, length, width, height, h, seeds, path_points);

	//	//write path points coordinates to file
	//	for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
	//		point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	//	}
	//	
	//	saving_name = path_name + "path.raw";
	//	store3dRawData<dataType>(resultedPath, length, width, height, saving_name.c_str());
	//	
	//	copyDataToAnotherArray(newActionPtr, actionPtr, height, length, width);
	//	count_step++;

	//	////find circles
	//	//len = path_points.size();
	//	//std::cout << len << " point to be processed" << std::endl;
	//	//
	//	//count_stop = 0;
	//	//for (i_n = len - 1; i_n > -1; i_n--) {
	//	//	num_slice_processed++;
	//	//	extension = to_string(num_slice_processed);
	//	//	//slice number
	//	//	k_n = (size_t)path_points[i_n].z;
	//	//	slice_number = to_string(k_n);
	//	//	seed2D.x = path_points[i_n].x;
	//	//	seed2D.y = path_points[i_n].y;
	//	//	//Slice extraction
	//	//	for (i = 0; i < dim2D; i++) {
	//	//		imageSlice[i] = ctImageData.imageDataPtr[k_n][i];
	//	//	}
	//	//	//Slice rescaling
	//	//	rescaleNewRange2D(imageSlice, length, width, 0.0, 1.0);
	//	//	//Slice filtering
	//	//	heatImplicit2dScheme(IMAGE, filter_parameters);
	//	//	if (num_slice_processed < 10) {
	//	//		saving_name = path_name + "threshold_000" + extension + ".raw";
	//	//	}
	//	//	else {
	//	//		if (num_slice_processed < 100) {
	//	//			saving_name = path_name + "threshold_00" + extension + ".raw";
	//	//		}
	//	//		else {
	//	//			if (num_slice_processed < 1000) {
	//	//				saving_name = path_name + "threshold_0" + extension + ".raw";
	//	//			}
	//	//			else {
	//	//				saving_name = path_name + "threshold_" + extension + ".raw";
	//	//			}
	//	//		}
	//	//	}
	//	//	center2D = localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, length, width, params, saving_name);
	//	//	max_ratio = getTheMaxValue(voteArray, length, width);
	//	//	if (max_ratio != 0.0) {
	//	//		centeredPath[k_n][x_new((size_t)center2D.x, (size_t)center2D.y, length)] = 1.0;
	//	//		//k_center = k_n;
	//	//		next_point.x = center2D.x;
	//	//		next_point.y = center2D.y;
	//	//		next_point.z = k_n;
	//	//	}
	//	//	if (num_slice_processed < 10) {
	//	//		saving_name = path_name + "slice_000" + extension + ".raw";
	//	//		store2dRawData(imageSlice, length, width, saving_name.c_str());
	//	//		saving_name = path_name + "found_circle_000" + extension + ".raw";
	//	//		store2dRawData(houghSpace, length, width, saving_name.c_str());
	//	//		saving_name = path_name + "ratio_000" + extension + ".raw";
	//	//		store2dRawData(voteArray, length, width, saving_name.c_str());
	//	//	}
	//	//	else {
	//	//		if (num_slice_processed < 100) {
	//	//			saving_name = path_name + "slice_00" + extension + ".raw";
	//	//			store2dRawData(imageSlice, length, width, saving_name.c_str());
	//	//			saving_name = path_name + "found_circle_00" + extension + ".raw";
	//	//			store2dRawData(houghSpace, length, width, saving_name.c_str());
	//	//			saving_name = path_name + "ratio_00" + extension + ".raw";
	//	//			store2dRawData(voteArray, length, width, saving_name.c_str());
	//	//		}
	//	//		else {
	//	//			if (num_slice_processed < 1000) {
	//	//				saving_name = path_name + "slice_0" + extension + ".raw";
	//	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//	//				saving_name = path_name + "found_circle_0" + extension + ".raw";
	//	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//	//				saving_name = path_name + "ratio_0" + extension + ".raw";
	//	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//	//			}
	//	//			else {
	//	//				saving_name = path_name + "slice_" + extension + ".raw";
	//	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//	//				saving_name = path_name + "found_circle_" + extension + ".raw";
	//	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//	//				saving_name = path_name + "ratio_" + extension + ".raw";
	//	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//	//			}
	//	//		}
	//	//	}
	//	//	if (max_ratio == 0.0) {
	//	//		count_stop++;
	//	//		//break;
	//	//	}
	//	//	else {
	//	//		count_stop = 0;
	//	//	}
	//	//	if (max_ratio == 0.0 && count_stop == stop_criterium) {
	//	//		break;
	//	//	}
	//	//}

	//	//empty the list of point
	//	while (path_points.size() != 0) {
	//		path_points.pop_back();
	//	}

	//	verif_point = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	verif_dist = getPoint3DDistance(last_point, verif_point);
	//	std::cout << "distance to end point : " << verif_dist << std::endl;
	//}

	//std::cout << num_slice_processed << " slices were processed\n" << std::endl;
	//fclose(file);
	//

	//delete[] imageSlice;
	//delete[] houghSpace;
	//delete[] voteArray;

	////delete[] seeds;
	//for (k = 0; k < height; k++) {
	//	delete[] actionPtr[k];
	//	delete[] maskAction[k];
	//	delete[] potentialPtr[k];
	//	delete[] edgeImage[k];
	//}
	//delete[] actionPtr;
	//delete[] maskAction;
	//delete[] potentialPtr;
	//delete[] edgeImage;

	return true;
}

bool findPathTwoSteps(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters) {

	//size_t i, j, k, xd;
	//size_t length = ctImageData.length, width = ctImageData.width, height = ctImageData.height;
	//size_t dim2D = length * width;
	//dataType** potential = new dataType * [height];
	//dataType** action_field = new dataType * [height];
	//dataType** mask_action = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	potential[k] = new dataType[dim2D]{ 0 };
	//	action_field[k] = new dataType[dim2D]{ 0 };
	//	mask_action[k] = new dataType[dim2D]{ 0 };
	//}
	//string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//string saving_name, saving_csv = outputPath + "path_point.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//Point3D temporary_point = { 0.0, 0.0, 0.0 };
	//Point3D* seedsPath = new Point3D[2];
	//seedsPath[0] = seedPoints[0];
	//seedsPath[1] = seedPoints[1];
	//vector<Point3D> path_points;
	//dataType tau = 0.95 * ctImageData.spacing.sx;
	//dataType tolerance = ctImageData.spacing.sx;

	////======== First step ===========//
	////compute3DPotential(ctImageData, potential, seedsPath[0], parameters);
	//saving_name = outputPath + "potential.raw";
	//manageRAWFile3D<dataType>(potential, length, width, height, saving_name.c_str(), LOAD_DATA, false);
	////partialFrontPropagation(action_field, potential, length, width, height, seedsPath);
	////fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);
	//Image_Data actionDataStr = { height, length, width, action_field, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	////shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	////shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);
	////copy points to file
	//int n = 0;
	//for (n = path_points.size() - 1; n > -1; n--) {
	//	path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//}
	//while (path_points.size() != 0) {
	//	path_points.pop_back();
	//}

	//================================//
	//Point3D seed1 = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//Point3D seed2 = getRealCoordFromImageCoord3D(seedPoints[2], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//double step = 0.5 * getPoint3DDistance(seed1, seed2), dist = 0.0;
	//cout << "Step : " << step << endl;

	//////======== Second step ===========//
	////compute3DPotential(ctImageData, potential, seedsPath[1], parameters);
	//fastMarching3D_N(action_field, potential, length, width, height, seedsPath[1]);
	//dataType min_action = BIG_VALUE, value_temp = 0.0;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			xd = x_new(i, j, length);
	//			Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//			dist = getPoint3DDistance(seed1, current_point);
	//			if (dist >= (step - 0.01) && dist <= (step + 0.01) && k >= seedPoints[1].z) {
	//				if (action_field[k][xd] < min_action) {
	//					min_action = action_field[k][xd];
	//					temporary_point.x = (dataType)i;
	//					temporary_point.y = (dataType)j;
	//					temporary_point.z = (dataType)k;
	//					value_temp = action_field[k][xd];
	//				}
	//			}
	//		}
	//	}
	//}
	//cout << "found point : (" << temporary_point.x << " ," << temporary_point.y << " , " << temporary_point.z << ")" << endl;
	//seedsPath[0] = seedPoints[1];
	//seedsPath[1] = temporary_point;
	////shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	////shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);
	//for (n = path_points.size() - 1; n > -1; n--) {
	//	path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//}
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}

	////======= Last step ===============//
	//seedsPath[0] = temporary_point;
	//seedsPath[1] = seedPoints[2];
	////partialFrontPropagation(action_field, potential, length, width, height, seedsPath);
	////fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);
	////shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	////shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);
	//for (n = path_points.size() - 1; n > -1; n--) {
	//	path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//}
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}
	////=====================================

	//fclose(file);
	//delete[] seedsPath;
	//for (k = 0; k < height; k++) {
	//	delete[] potential[k];
	//	delete[] action_field[k];
	//	delete[] mask_action[k];
	//}
	//delete[] potential;
	//delete[] action_field;
	//delete[] mask_action;

	return true;
}

bool findPathFromOneGivenPoint(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters) {

	//if (ctImageData.imageDataPtr == NULL)
	//	return false;
	//size_t k = 0, i = 0, j = 0, x = 0;
	//const size_t height = ctImageData.height, length = ctImageData.length, width = ctImageData.width;
	//const size_t dim2D = length * width;
	//double epsilon = 0.01;
	//string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//string saving_name, extension;
	//saving_name = path_name + "path_points_v2.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//dataType** actionPtr = new dataType * [height];
	//dataType** newActionPtr = new dataType * [height];
	//dataType** maskDistance = new dataType * [height];
	//dataType** potentialPtr = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	actionPtr[k] = new dataType[dim2D]{ 0 };
	//	newActionPtr[k] = new dataType[dim2D]{ 0 };
	//	maskDistance[k] = new dataType[dim2D]{ 0 };
	//	potentialPtr[k] = new dataType[dim2D]{ 0 };
	//}
	//if (actionPtr == NULL || newActionPtr == NULL || maskDistance == NULL || potentialPtr == NULL)
	//	return false;

	//Point3D* seeds = new Point3D[2];
	//seeds[0] = seedPoints[0];
	//Point3D initial_point = seedPoints[0];
	//Point3D temporary_point = { 0.0, 0.0, 0.0 };
	//dataType min_distance = BIG_VALUE, value_temp = 0.0, h = 1.0;
	//Point3D seed1 = getRealCoordFromImageCoord3D(seedPoints[0], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//Point3D seed2 = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//double distance_point = getPoint3DDistance(seed1, seed2);
	//double step = 50;//0.5 * distance_point;
	//std::cout << "First step : " << step << std::endl;

	////================ find next point inside the aorta 
	//Image_Data actionDataStr = { height, length, width, actionPtr, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	//double radius_to_find_mean = 3.0;
	////compute3DPotential(ctImageData, potentialPtr, seedPoints[0], parameters);
	//saving_name = path_name + "potential.raw";
	////store3dRawData<dataType>(potentialPtr, length, width, height, saving_name.c_str());
	//manageRAWFile3D<dataType>(potentialPtr, length, width, height, saving_name.c_str(), LOAD_DATA, false);
	////fastMarching3D_N(actionPtr, potentialPtr, length, width, height, seeds[0]);
	//fastMarching3dWithSpacing(ctImageData, actionPtr, potentialPtr, seeds[0], ctImageData.spacing);
	////saving_name = path_name + "first_action.raw";
	////manageRAWFile3D<dataType>(actionPtr, length, width, height, saving_name.c_str(), LOAD_DATA, false);

	//min_distance = BIG_VALUE;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			x = x_new(i, j, length);
	//			Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//			double distanceP = getPoint3DDistance(seed1, current_point);
	//			if (distanceP >= (step - epsilon) && distanceP <= (step + epsilon)) {
	//				if (actionPtr[k][x] < min_distance) {
	//					min_distance = actionPtr[k][x];
	//					temporary_point.x = (dataType)i;
	//					temporary_point.y = (dataType)j;
	//					temporary_point.z = (dataType)k;
	//					value_temp = actionPtr[k][x];
	//				}
	//			}
	//		}
	//	}
	//}
	//seed1 = getRealCoordFromImageCoord3D(temporary_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	////distance_point = getPoint3DDistance(seed1, seed2);
	////step = 0.5 * distance_point;
	////std::cout << "Second step: " << step << std::endl;
	//std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;
	////vector for path points
	//vector<Point3D> path_points;
	////path between initial point and second point
	//seeds[1] = temporary_point;
	//dataType tau = 0.95 * ctImageData.spacing.sx;
	//dataType tolerance = ctImageData.spacing.sx;
	////shortestPath3D(actionDataStr, seeds, tau, tolerance, path_points);
	////write just path points coordinates to file
	//int i_n = 0, k_n = 0, k_center = 0;
	//Point3D point_file = { 0.0, 0.0, 0.0 };
	//for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
	//	point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	//}
	//while (path_points.size() != 0) {
	//	path_points.pop_back();
	//}

	////============================================
	//double distance_to_stop = 0.1 * distance_point;
	//size_t count_step = 0;
	//while (count_step < 3) {
	//	count_step++;
	//	std::cout << "count step = " << count_step << std::endl;
	//	//mask
	//	for (k = 0; k < height; k++) {
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				x = x_new(i, j, length);
	//				Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
	//				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//				double distanceP = getPoint3DDistance(seed1, current_point);
	//				if (distanceP <= (step + 0.01) && actionPtr[k][x] <= value_temp) {
	//					maskDistance[k][x] = 10000.0;
	//				}
	//				else {
	//					if (maskDistance[k][x] != 10000.0) {
	//						maskDistance[k][x] = newActionPtr[k][x];
	//					}
	//				}
	//			}
	//		}
	//	}
	//	seeds[0] = seeds[1];
	//	initial_point = seeds[1];
	//	//fastMarching3D_N(newActionPtr, potentialPtr, length, width, height, seeds[0]);
	//	fastMarching3dWithSpacing(ctImageData, newActionPtr, potentialPtr, seeds[0], ctImageData.spacing);
	//	//saving_name = path_name + "action_field" + extension + ".raw";
	//	//store3dRawData<dataType>(newActionPtr, length, width, height, saving_name.c_str());
	//	//Find the temporary point
	//	min_distance = BIG_VALUE;
	//	//get real world coordinate of the initial point
	//	seed1 = getRealCoordFromImageCoord3D(initial_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	for (k = 0; k < height; k++) {
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				x = x_new(i, j, length);
	//				Point3D current_point = { (dataType)i, (dataType)j, current_point.z = (dataType)k };
	//				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//				double distanceP = getPoint3DDistance(seed1, current_point);
	//				if (distanceP >= (step - 0.01) && distanceP <= (step + 0.01) && maskDistance[k][x] != 10000.0) {
	//					if (newActionPtr[k][x] < min_distance) {
	//						min_distance = newActionPtr[k][x];
	//						temporary_point.x = (dataType)i;
	//						temporary_point.y = (dataType)j;
	//						temporary_point.z = (dataType)k;
	//						value_temp = newActionPtr[k][x];
	//					}
	//				}
	//			}
	//		}
	//	}
	//	std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;
	//	seeds[1] = temporary_point;
	//	//seed2 = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//	//distance_point = getPoint3DDistance(seed1, seed2);
	//	//step = 0.5 * distance_point;
	//	//std::cout << "The step is : " << step << std::endl;
	//	actionDataStr.imageDataPtr = newActionPtr;
	//	//shortestPath3D(actionDataStr, seeds, tau, tolerance, path_points);
	//	copyDataToAnotherArray(newActionPtr, actionPtr, height, length, width);
	//	//write just path points coordinates to file
	//	for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
	//		point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	//	}
	//	while (path_points.size() != 0) {
	//		path_points.pop_back();
	//	}
	//}

	//fclose(file);
	//delete[] seeds;
	//for (k = 0; k < height; k++) {
	//	delete[] actionPtr[k];
	//	delete[] newActionPtr[k];
	//	delete[] maskDistance[k];
	//	delete[] potentialPtr[k];
	//}
	//delete[] actionPtr;
	//delete[] newActionPtr;
	//delete[] maskDistance;
	//delete[] potentialPtr;

	return true;
}

*/

/*
bool partialPropagation(Image_Data actionPtr, dataType** potentialPtr, dataType** maskPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints, const double maxLength) {

	if (actionPtr.imageDataPtr == NULL || potentialPtr == NULL || maskPtr == NULL) {
		return false;
	}

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short* [height];
	dataType** weighted_distance_saving = new dataType * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D] {0};
		weighted_distance_saving[k] = new dataType[dim2D]{ 0 };
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (maskPtr[k][i] == 0) {
				labelArray[k][i] = 3;
			}
			else {
				labelArray[k][i] = 1;
			}
			actionPtr.imageDataPtr[k][i] = INFINITY;
		}
	}

	//Processed the initial point
	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	k = (size_t)seedPoints[0].z;
	currentIndx = x_new(i, j, length);
	actionPtr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;
	pointFastMarching3D current = { i, j, k, 0.0 };

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialPtr[kminus][currentIndx];
		dTop = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialPtr[k][indxEast];
		dEast = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		actionPtr.imageDataPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialPtr[k][indxNorth];
		dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		actionPtr.imageDataPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialPtr[k][indxWest];
		dWest = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		actionPtr.imageDataPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialPtr[k][indxSouth];
		dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		actionPtr.imageDataPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialPtr[kplus][currentIndx];
		dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	//Visualize the front propagation
	size_t id_keyPoint = 1;
	std::string root_path = "C:/Users/Konan Allaly/Documents/Tests/output/";
	std::string storing_path;
	size_t saving_time = 0;

	bool status_point = false;
	double length_propagation = 0.0;
	Point3D pStart = getRealCoordFromImageCoord3D(seedPoints[0], actionPtr.origin, actionPtr.spacing, actionPtr.orientation);

	while (length_propagation < maxLength) {

		//processed the point with minimum distance
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		i = current.x;
		j = current.y;
		k = current.z;

		Point3D pEnd = { current.x, current.y, current.z };
		pEnd = getRealCoordFromImageCoord3D(pEnd, actionPtr.origin, actionPtr.spacing, actionPtr.orientation);
		length_propagation = getPoint3DDistance(pStart, pEnd);
		if (length_propagation > maxLength) {
			status_point = true;
		}

		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		actionPtr.imageDataPtr[k][currentIndx] = current.arrival;

		if (status_point == true) {
			seedPoints[1].x = inProcess[0].x;
			seedPoints[1].y = inProcess[0].y;
			seedPoints[1].z = inProcess[0].z;
		}
		else {
			seedPoints[1].x = current.x;
			seedPoints[1].y = current.y;
			seedPoints[1].z = current.z;
		}

		for (size_t ik = 0; ik < height; ik++) {
			for (size_t ii = 0; ii < dim2D; ii++) {
				if (actionPtr.imageDataPtr[ik][ii] == INFINITY) {
					weighted_distance_saving[ik][ii] = 0.0;
				}
				else {
					weighted_distance_saving[ik][ii] = actionPtr.imageDataPtr[ik][ii];
				}
			}
		}

		saving_time++;
		storing_path = root_path + "weighted distance/partial_action_" + std::to_string(id_keyPoint) + ".raw";
		manageRAWFile3D<dataType>(weighted_distance_saving, length, width, height, storing_path.c_str(), STORE_DATA, false);

		deleteRootHeap3D(inProcess);

		//processed neighbors of the minimum in the narrow band

		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			kminus = k - 1;
			label = labelArray[kminus][currentIndx];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				coefSpeed = potentialPtr[kminus][currentIndx];
				dTop = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
					labelArray[kminus][currentIndx] = 2;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (label == 2) {
						if (dTop < actionPtr.imageDataPtr[kminus][currentIndx]) {
							actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(iplus, j, length);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialPtr[k][indxEast];
				dEast = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxEast] = dEast;
					labelArray[k][indxEast] = 2;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < actionPtr.imageDataPtr[k][indxEast]) {
							actionPtr.imageDataPtr[k][indxEast] = dEast;
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(i, jminus, length);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialPtr[k][indxNorth];
				dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxNorth] = dNorth;
					labelArray[k][indxNorth] = 2;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < actionPtr.imageDataPtr[k][indxNorth]) {
							actionPtr.imageDataPtr[k][indxNorth] = dNorth;
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(iminus, j, length);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialPtr[k][indxWest];
				dWest = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxWest] = dWest;
					labelArray[k][indxWest] = 2;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < actionPtr.imageDataPtr[k][indxWest]) {
							actionPtr.imageDataPtr[k][indxWest] = dWest;
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(i, jplus, length);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialPtr[k][indxSouth];
				dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxSouth] = dSouth;
					labelArray[k][indxSouth] = 2;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < actionPtr.imageDataPtr[k][indxSouth]) {
							actionPtr.imageDataPtr[k][indxSouth] = dSouth;
						}
					}
				}
			}
		}

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialPtr[kplus][currentIndx];
				dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				if (label == 3) {
					actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
					labelArray[kplus][currentIndx] = 2;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (label == 2) {
						if (dBottom < actionPtr.imageDataPtr[kplus][currentIndx]) {
							actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
						}
					}
				}
			}
		}

	}

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (actionPtr.imageDataPtr[k][i] == INFINITY) {
				actionPtr.imageDataPtr[k][i] = 0.0;
				weighted_distance_saving[k][i] = 0.0;
			}
			else {
				weighted_distance_saving[k][i] = actionPtr.imageDataPtr[k][i];
			}
		}
	}

	saving_time++;
	storing_path = root_path + "weighted distance/partial_action_" + std::to_string(id_keyPoint) + ".raw";
	manageRAWFile3D<dataType>(weighted_distance_saving, length, width, height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
		delete[] weighted_distance_saving[k];
	}
	delete[] labelArray;
	delete[] weighted_distance_saving;

	return true;
}
*/

/*
bool partialPropagationWithSpacing(Image_Data actionPtr, dataType** potentialPtr, dataType** maskPtr, Point3D* seedPoints, const double maxLength) {

	if (actionPtr.imageDataPtr == NULL || potentialPtr == NULL || maskPtr == NULL) {
		return false;
	}

	const size_t length = actionPtr.length;
	const size_t width = actionPtr.width;
	const size_t height = actionPtr.height;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D] {0};
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;

	VoxelSpacing spacing = actionPtr.spacing;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (maskPtr[k][i] == 0) {
				labelArray[k][i] = 3;
			}
			else {
				labelArray[k][i] = 1;
			}
			actionPtr.imageDataPtr[k][i] = INFINITY;
		}
	}

	//Processed the initial point
	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	k = (size_t)seedPoints[0].z;
	currentIndx = x_new(i, j, length);
	actionPtr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;
	pointFastMarching3D current = { i, j, k, 0.0 };

	//find the neighbours of the initial point add add them to inProcess
	size_t height_minus = height - 1, length_minus = length - 1, width_minus = width - 1;
	size_t kminus = 0, kplus = 0, iminus = 0, iplus = 0, jminus = 0, jplus = 0;
	size_t indxNorth = 0, indxSouth = 0, indxWest = 0, indxEast = 0;
	//dataType dTop = 0.0, dBottom = 0.0, dNorth = 0.0, dSouth = 0.0, dWest = 0.0, dEast = 0.0, coefSpeed = 0.0;

	//Top
	if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
		kminus = k - 1;
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
		dataType coefSpeed = potentialPtr[kminus][currentIndx];
		//dTop = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
		dataType coefSpeed = potentialPtr[k][indxEast];
		//dEast = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		actionPtr.imageDataPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
		dataType coefSpeed = potentialPtr[k][indxNorth];
		//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		actionPtr.imageDataPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
		dataType coefSpeed = potentialPtr[k][indxWest];
		//dWest = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		actionPtr.imageDataPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
		dataType coefSpeed = potentialPtr[k][indxSouth];
		//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		actionPtr.imageDataPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
		dataType coefSpeed = potentialPtr[kplus][currentIndx];
		//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		dataType dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	bool status_point = false;
	double length_propagation = 0.0;
	Point3D pStart = getRealCoordFromImageCoord3D(seedPoints[0], actionPtr.origin, spacing, actionPtr.orientation);
	double max_action = 0.0;

	while (length_propagation < maxLength) {

		//processed the point with minimum distance
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		i = current.x;
		j = current.y;
		k = current.z;

		Point3D pEnd = { current.x, current.y, current.z };
		pEnd = getRealCoordFromImageCoord3D(pEnd, actionPtr.origin, spacing, actionPtr.orientation);
		length_propagation = getPoint3DDistance(pStart, pEnd);
		if (length_propagation > maxLength) {
			status_point = true;
		}

		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		actionPtr.imageDataPtr[k][currentIndx] = current.arrival;
		max_action = current.arrival;

		deleteRootHeap3D(inProcess);

		if (status_point == true) {
			seedPoints[1].x = inProcess[0].x;
			seedPoints[1].y = inProcess[0].y;
			seedPoints[1].z = inProcess[0].z;
		}
		else {
			seedPoints[1].x = current.x;
			seedPoints[1].y = current.y;
			seedPoints[1].z = current.z;
		}

		//processed neighbors of the minimum in the narrow band

		//Top
		if (k > 0 && j >= 0 && j < width && i >= 0 && i < length) {
			kminus = k - 1;
			label = labelArray[kminus][currentIndx];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kminus);
				dataType coefSpeed = potentialPtr[kminus][currentIndx];
				//dTop = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
					labelArray[kminus][currentIndx] = 2;
					pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (label == 2) {
						if (dTop < actionPtr.imageDataPtr[kminus][currentIndx]) {
							actionPtr.imageDataPtr[kminus][currentIndx] = dTop;
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(iplus, j, length);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, iplus, j, k);
				dataType coefSpeed = potentialPtr[k][indxEast];
				//dEast = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxEast] = dEast;
					labelArray[k][indxEast] = 2;
					pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < actionPtr.imageDataPtr[k][indxEast]) {
							actionPtr.imageDataPtr[k][indxEast] = dEast;
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(i, jminus, length);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jminus, k);
				dataType coefSpeed = potentialPtr[k][indxNorth];
				//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxNorth] = dNorth;
					labelArray[k][indxNorth] = 2;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < actionPtr.imageDataPtr[k][indxNorth]) {
							actionPtr.imageDataPtr[k][indxNorth] = dNorth;
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(iminus, j, length);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, iminus, j, k);
				dataType coefSpeed = potentialPtr[k][indxWest];
				//dWest = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxWest] = dWest;
					labelArray[k][indxWest] = 2;
					pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < actionPtr.imageDataPtr[k][indxWest]) {
							actionPtr.imageDataPtr[k][indxWest] = dWest;
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(i, jplus, length);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, jplus, k);
				dataType coefSpeed = potentialPtr[k][indxSouth];
				//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[k][indxSouth] = dSouth;
					labelArray[k][indxSouth] = 2;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < actionPtr.imageDataPtr[k][indxSouth]) {
							actionPtr.imageDataPtr[k][indxSouth] = dSouth;
						}
					}
				}
			}
		}

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				y = select3dY(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				z = select3dZ(actionPtr.imageDataPtr, length, width, height, i, j, kplus);
				dataType coefSpeed = potentialPtr[kplus][currentIndx];
				//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				dataType dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
				if (label == 3) {
					actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
					labelArray[kplus][currentIndx] = 2;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (label == 2) {
						if (dBottom < actionPtr.imageDataPtr[kplus][currentIndx]) {
							actionPtr.imageDataPtr[kplus][currentIndx] = dBottom;
						}
					}
				}
			}
		}

	}

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
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

bool computeDistanceToOnePoint(dataType** distancePtr, const size_t length, const size_t width, const size_t height, Point3D seed) {

	if (distancePtr == NULL)
		return false;

	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
				distancePtr[k][x_new(i, j, length)] = getPoint3DDistance(seed, current_point);
			}
		}
	}
	return true;

}
*/

/*
dataType select3dX(dataType** actionPtr, short** label, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {

	//if (i == 0) {
	//	return min(actionPtr[k][x_new(i, j, length)], actionPtr[k][x_new(i + 1, j, length)]);
	//}
	//else {
	//	if (i == length - 1) {
	//		return min(actionPtr[k][x_new(i - 1, j, length)], actionPtr[k][x_new(i, j, length)]);
	//	}
	//	else {
	//		return min(actionPtr[k][x_new(i - 1, j, length)], actionPtr[k][x_new(i + 1, j, length)]);
	//	}
	//}

	dataType x_minus, x_plus;
	if (i == 0) {
		x_minus = INFINITY;
	}
	else {
		if (label[k][x_new(i - 1, j, length)] == 1) {
			x_minus = actionPtr[k][x_new(i - 1, j, length)];
		}
		else {
			x_minus = INFINITY;
		}
	}
	if(i == length - 1) {
		x_plus = INFINITY;
	}
	else {
		if(label[k][x_new(i + 1, j, length)] == 1) {
			x_plus = actionPtr[k][x_new(i + 1, j, length)];
		}
		else {
			x_plus = INFINITY;
		}
	}
	return min(x_minus, x_plus);
}

dataType select3dY(dataType** actionPtr, short** label, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {
	//if (j == 0) {
	//	return min(actionPtr[k][x_new(i, j, length)], actionPtr[k][x_new(i, j + 1, length)]);
	//}
	//else {
	//	if (j == width - 1) {
	//		return min(actionPtr[k][x_new(i, j - 1, length)], actionPtr[k][x_new(i, j, length)]);
	//	}
	//	else {
	//		return min(actionPtr[k][x_new(i, j - 1, length)], actionPtr[k][x_new(i, j + 1, length)]);
	//	}
	//}

	dataType y_minus, y_plus;
	if(j == 0) {
		y_minus = INFINITY;
	}
	else {
		if(label[k][x_new(i, j - 1, length)] == 1) {
			y_minus = actionPtr[k][x_new(i, j - 1, length)];
		}
		else {
			y_minus = INFINITY;
		}
	}
	if(j == width - 1) {
		y_plus = INFINITY;
	}
	else {
		if(label[k][x_new(i, j + 1, length)] == 1) {
			y_plus = actionPtr[k][x_new(i, j + 1, length)];
		}
		else {
			y_plus = INFINITY;
		}
	}
	return min(y_minus, y_plus);
}

dataType select3dZ(dataType** actionPtr, short** label, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {
	size_t xd = x_new(i, j, length);
	//if (k == 0) {
	//	return min(actionPtr[k + 1][xd], actionPtr[k][xd]);
	//}
	//else {
	//	if (k == height - 1) {
	//		return min(actionPtr[k - 1][xd], actionPtr[k][xd]);
	//	}
	//	else {
	//		return min(actionPtr[k + 1][xd], actionPtr[k - 1][xd]);
	//	}
	//}

	dataType z_minus, z_plus;
	if(k == 0) {
		z_minus = INFINITY;
	}
	else {
		if(label[k - 1][xd] == 1) {
			z_minus = actionPtr[k - 1][xd];
		}
		else {
			z_minus = INFINITY;
		}
	}
	if(k == height - 1) {
		z_plus = INFINITY;
	}
	else {
		if(label[k + 1][xd] == 1) {
			z_plus = actionPtr[k + 1][xd];
		}
		else {
			z_plus = INFINITY;
		}
	}
	return min(z_minus, z_plus);
}
*/

/*
bool fastMarching3dForDistanceMap(Image_Data ctImageData, dataType** distanceFuncPtr, dataType foregroundValue) {

	if (ctImageData.imageDataPtr == NULL || distanceFuncPtr == NULL) {
		return false;
	}

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	VoxelSpacing spacing = ctImageData.spacing;

	vector <pointFastMarching3D> inProcess;
	inProcess.reserve(length * width * height); // Reserve space to avoid frequent reallocations
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

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
							size_t indx = x_flat(i, j, kminus, length, width);
							if (label == 3) {
								size_t posTop = inProcess.size();
								distanceFuncPtr[kminus][currentIndx] = dTop;
								pointFastMarching3D TopNeighbor = { i, j, kminus, posTop, distanceFuncPtr[kminus][currentIndx] };
								inProcess.push_back(TopNeighbor);
								indexArray[indx] = &inProcess[posTop].pos;
								labelArray[kminus][currentIndx] = 2;
							}
							else {
								if (dTop < distanceFuncPtr[kminus][currentIndx]) {
									distanceFuncPtr[kminus][currentIndx] = dTop;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dTop;
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
							size_t indx = x_flat(i, j, kplus, length, width);
							if (label == 3) {
								size_t posBottom = inProcess.size();
								distanceFuncPtr[kplus][currentIndx] = dBottom;
								pointFastMarching3D BottomNeighbor = { i, j, kplus, posBottom, distanceFuncPtr[kplus][currentIndx] };
								inProcess.push_back(BottomNeighbor);
								indexArray[indx] = &inProcess[posBottom].pos;
								labelArray[kplus][currentIndx] = 2;
							}
							else {
								if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
									distanceFuncPtr[kplus][currentIndx] = dBottom;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dBottom;
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
							size_t indx = x_flat(iplus, j, k, length, width);
							if (label == 3) {
								size_t posEast = inProcess.size();
								distanceFuncPtr[k][indxEast] = dEast;
								pointFastMarching3D EastNeighbor = { iplus, j, k, posEast, distanceFuncPtr[k][indxEast] };
								inProcess.push_back(EastNeighbor);
								indexArray[indx] = &inProcess[posEast].pos;
								labelArray[k][indxEast] = 2;
							}
							else {
								if (dEast < distanceFuncPtr[k][indxEast]) {
									distanceFuncPtr[k][indxEast] = dEast;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dEast;
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
							size_t indx = x_flat(iminus, j, k, length, width);
							if (label == 3) {
								size_t posWest = inProcess.size();
								distanceFuncPtr[k][indxWest] = dWest;
								pointFastMarching3D WestNeighbor = { iminus, j, k, posWest, distanceFuncPtr[k][indxWest] };
								inProcess.push_back(WestNeighbor);
								size_t indx = x_flat(iminus, j, k, length, width);
								indexArray[indx] = &inProcess[posWest].pos;
								labelArray[k][indxWest] = 2;
							}
							else {
								if (dWest < distanceFuncPtr[k][indxWest]) {
									distanceFuncPtr[k][indxWest] = dWest;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dWest;
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
							size_t indx = x_flat(i, jminus, k, length, width);
							if (label == 3) {
								size_t posNorth = inProcess.size();
								distanceFuncPtr[k][indxNorth] = dNorth;
								pointFastMarching3D NorthNeighbor = { i, jminus, k, posNorth, distanceFuncPtr[k][indxNorth] };
								inProcess.push_back(NorthNeighbor);
								indexArray[indx] = &inProcess[posNorth].pos;
								labelArray[k][indxNorth] = 2;
							}
							else {
								if (dNorth < distanceFuncPtr[k][indxNorth]) {
									distanceFuncPtr[k][indxNorth] = dNorth;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dNorth;
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
							size_t indx = x_flat(i, jplus, k, length, width);
							if (label == 3) {
								size_t posSouth = inProcess.size();
								distanceFuncPtr[k][indxSouth] = dSouth;
								pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth, distanceFuncPtr[k][indxSouth] };
								inProcess.push_back(SouthNeighbor);
								indexArray[indx] = &inProcess[posSouth].pos;
								labelArray[k][indxSouth] = 2;
							}
							else {
								if (dSouth < distanceFuncPtr[k][indxSouth]) {
									distanceFuncPtr[k][indxSouth] = dSouth;
									size_t pt_pos = *indexArray[indx];
									inProcess[pt_pos].arrival = dSouth;
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
		indexArray[x_flat(i, j, k, length, width)] = NULL; // Clear the index pointer for this point

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
				size_t indx = x_flat(i, j, kminus, length, width);
				if (label == 3) {
					size_t posTop = inProcess.size();
					distanceFuncPtr[kminus][currentIndx] = dTop;
					pointFastMarching3D TopNeighbor = { i, j, kminus, posTop, distanceFuncPtr[kminus][currentIndx] };
					inProcess.push_back(TopNeighbor);
					indexArray[indx] = &inProcess[posTop].pos;
					heapifyUp3D(inProcess, posTop);
					labelArray[kminus][currentIndx] = 2;
				}
				else {
					if (dTop < distanceFuncPtr[kminus][currentIndx]) {
						distanceFuncPtr[kminus][currentIndx] = dTop;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dTop;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, j, kplus, length, width);
				if (label == 3) {
					size_t posBottom = inProcess.size();
					distanceFuncPtr[kplus][currentIndx] = dBottom;
					pointFastMarching3D BottomNeighbor = { i, j, kplus, posBottom, distanceFuncPtr[kplus][currentIndx] };
					inProcess.push_back(BottomNeighbor);
					indexArray[indx] = &inProcess[posBottom].pos;
					heapifyUp3D(inProcess, posBottom);
					labelArray[kplus][currentIndx] = 2;
				}
				else {
					if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
						distanceFuncPtr[kplus][currentIndx] = dBottom;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dBottom;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(iplus, j, k, length, width);
				if (label == 3) {
					size_t posEast = inProcess.size();
					distanceFuncPtr[k][indxEast] = dEast;
					pointFastMarching3D EastNeighbor = { iplus, j, k, posEast, distanceFuncPtr[k][indxEast] };
					inProcess.push_back(EastNeighbor);
					indexArray[indx] = &inProcess[posEast].pos;
					heapifyUp3D(inProcess, posEast);
					labelArray[k][indxEast] = 2;
				}
				else {
					if (dEast < distanceFuncPtr[k][indxEast]) {
						distanceFuncPtr[k][indxEast] = dEast;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dEast;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(iminus, j, k, length, width);
				if (label == 3) {
					size_t posWest = inProcess.size();
					distanceFuncPtr[k][indxWest] = dWest;
					pointFastMarching3D WestNeighbor = { iminus, j, k, posWest, distanceFuncPtr[k][indxWest] };
					inProcess.push_back(WestNeighbor);
					size_t indx = x_flat(iminus, j, k, length, width);
					indexArray[indx] = &inProcess[posWest].pos;
					heapifyUp3D(inProcess, posWest);
					labelArray[k][indxWest] = 2;
				}
				else {
					if (dWest < distanceFuncPtr[k][indxWest]) {
						distanceFuncPtr[k][indxWest] = dWest;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dWest;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, jminus, k, length, width);
				if (label == 3) {
					size_t posNorth = inProcess.size();
					distanceFuncPtr[k][indxNorth] = dNorth;
					pointFastMarching3D NorthNeighbor = { i, jminus, k, posNorth, distanceFuncPtr[k][indxNorth] };
					inProcess.push_back(NorthNeighbor);
					indexArray[indx] = &inProcess[posNorth].pos;
					heapifyUp3D(inProcess, posNorth);
					labelArray[k][indxNorth] = 2;
				}
				else {
					if (dNorth < distanceFuncPtr[k][indxNorth]) {
						distanceFuncPtr[k][indxNorth] = dNorth;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dNorth;
						heapifyUp3D(inProcess, pt_pos);
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
				size_t indx = x_flat(i, jplus, k, length, width);
				if (label == 3) {
					size_t posSouth = inProcess.size();
					distanceFuncPtr[k][indxSouth] = dSouth;
					pointFastMarching3D SouthNeighbor = { i, jplus, k, posSouth, distanceFuncPtr[k][indxSouth] };
					inProcess.push_back(SouthNeighbor);
					indexArray[indx] = &inProcess[posSouth].pos;
					heapifyUp3D(inProcess, posSouth);
					labelArray[k][indxSouth] = 2;
				}
				else {
					if (dSouth < distanceFuncPtr[k][indxSouth]) {
						distanceFuncPtr[k][indxSouth] = dSouth;
						size_t pt_pos = *indexArray[indx];
						inProcess[pt_pos].arrival = dSouth;
						heapifyUp3D(inProcess, pt_pos);
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
*/

/*
void computePotentialMeanVariance(Image_Data2D ctImageData, dataType* potentialPtr, Point2D seedPoint, dataType radiusInitial, dataType radiusMax, dataType radiusStep) {

	size_t height = ctImageData.height;
	size_t width = ctImageData.width;
	size_t i = 0, j = 0;
	dataType lambda1 = 1.0, lambda2 = 1.0, omega = 0.00001;

	Statistics stats = get2DPointNeighborhoodStats(ctImageData, seedPoint, radiusInitial);
	dataType meanInitial = stats.mean_data / radiusInitial, mean = 0.0;
	dataType varianceInitial = stats.variance / radiusInitial, variance = 0.0;

	dataType radius = 0.0, value = 0.0, minValue = 0.0;
	dataType meanVal = 0.0, varianceVal = 0.0, optimalRadius = 0.0;

	dataType* radiusArray = new dataType[height * width]{ 0 };

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			radius = radiusInitial;
			meanVal = 0.0;
			varianceVal = 0.0;
			optimalRadius = 0.0;
			Point2D current_point = { i, j };
			while (radius <= radiusMax) {
				stats = get2DPointNeighborhoodStats(ctImageData, current_point, radius);
				mean = stats.mean_data / radius;
				variance = stats.variance / radius;
				if (mean >= meanVal) {
					meanVal = mean;
					//optimalRadius = radius;
				}
				if (variance >= varianceVal) {
					varianceVal = variance;
					optimalRadius = radius;
				}
				radius = radius + radiusStep;
			}
			value = omega + lambda1 * pow(meanVal - meanInitial, 2) + lambda2 * pow(varianceVal - varianceInitial, 2);
			potentialPtr[x_new(i, j, height)] = value;
			radiusArray[x_new(i, j, height)] = optimalRadius;
		}
	}

	delete[] radiusArray;
}
*/
