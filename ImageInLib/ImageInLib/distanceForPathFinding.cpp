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

dataType solve2dQuadratic(dataType X, dataType Y, dataType W) {

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P = sqrt(W);

	a = 2.0; 
	if (X == INFINITY) {
		X = 0; 
		a = a - 1;;
	}
	if (Y == INFINITY) {
		Y = 0; 
		a = a - 1;
	}

	//a can not be equal to 0

	b = -2 * (X + Y); 
	c = (dataType)(pow(X, 2) + pow(Y, 2) - W);
	delta = (dataType)(pow(b, 2) - 4 * a * c);

	if (delta >= 0) {
		solution = (dataType)((-b + sqrt(delta)) / (2 * a));
	}
	else {
		solution = (dataType)(min(X, Y) + P);
	}

	return solution;
}

dataType selectX(dataType* actionPtr, const size_t height, const size_t width, const size_t i, const size_t j) {
	
	// this function return the minimum in the upwind principle
	// x--->j and y--->i
	dataType i_minus, i_plus;

	if (i == 0) {
		i_minus = INFINITY;
	}
	else {
		i_minus = actionPtr[x_new(i - 1, j, height)];
	}

	if (i == height - 1) {
		i_plus = INFINITY;
	}
	else {
		i_plus = actionPtr[x_new(i + 1, j, height)];
	}

	return min(i_minus, i_plus);
}

dataType selectY(dataType* actionPtr, const size_t height, const size_t width, const size_t i, const size_t j) {
	
	// this function return the minimum in the upwind principle
	// x--->j and y--->i
	dataType j_minus, j_plus;

	if (j == 0) {
		j_minus = BIG_VALUE;
	}
	else {
		j_minus = actionPtr[x_new(i, j - 1, height)];
	}

	if (j == width - 1) {
		j_plus = BIG_VALUE;
	}
	else {
		j_plus = actionPtr[x_new(i, j + 1, height)];
	}

	return min(j_minus, j_plus);
}

bool computePotential(Image_Data2D imageDataStr, dataType* potentialFuncPtr, Point2D * seedPoints)
{
	//This function is used to compute the potential function
	//were epsilon can be used also as parameter
	//epislon = 0.01 was perfect for our experimentations

	if (imageDataStr.imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;


	const size_t height = imageDataStr.height;
	const size_t width = imageDataStr.width;
	size_t i, dim2D = height * width;
	size_t i1 = (size_t)seedPoints[0].x, j1 = (size_t)seedPoints[0].y;
	size_t i2 = (size_t)seedPoints[1].x, j2 = (size_t)seedPoints[1].y;

	dataType seedVal = (dataType)((imageDataStr.imageDataPtr[x_new(i1, j1, height)] + imageDataStr.imageDataPtr[x_new(i2, j2, height)]) / 2.0);
	dataType epsilon = 0.01;
	dataType K = 0.00005;
	dataType coef_dist = 1000, coef = 0.0;

	dataType* gradientVectorX = new dataType[dim2D]{ 0 };
	dataType* gradientVectorY = new dataType[dim2D]{ 0 };
	dataType* edgeDetector = new dataType[dim2D]{ 0 };

	computeImageGradient(imageDataStr, gradientVectorX, gradientVectorY);

	dataType ux = 0.0, uy = 0.0, norm_of_gradient_square = 0.0;
	for (i = 0; i < dim2D; i++) {
		ux = gradientVectorX[i];
		uy = gradientVectorY[i];
		norm_of_gradient_square = ux * ux + uy * uy;
		edgeDetector[i] = 0.01 + sqrt(norm_of_gradient_square);
	}
	
	//Normalization
	dataType weight = 0.0, edgeValue = 0.0, norm_of_gradient = 0.0;
	for (i = 0; i < dim2D; i++) {
		potentialFuncPtr[i] = edgeDetector[i];
	}

	delete[] gradientVectorX;
	delete[] gradientVectorY;

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
		if (val_current < val_parent) {
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

bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, Point2D* seedPoints) {

	size_t dim2D = height * width;

	short* labelArray = new short[dim2D] {0};

	if (imageDataPtr == NULL || distancePtr == NULL || potentialPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0;
	vector<pointFastMarching2D> inProcess;

	dataType x = 0.0, y = 0.0, coefSpeed = 0.0;

	//Compute the potential function
	//computePotential(imageDataPtr, potentialPtr, height, width, seedPoints);

	//STEP 1
	//In labelAray we have : 
	//1 ---> already processed, 
	//2 ---> in process and 
	//3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distancePtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	size_t currentIndx = x_new(i, j, height);
	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < height) {
		size_t jplus = j + 1;
		x = selectX(distancePtr, height, width, i, jplus);
		y = selectY(distancePtr, height, width, i, jplus);
		size_t indxEast = x_new(i, jplus, height);
		coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		distancePtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < height) {
		size_t jminus = j - 1;
		x = selectX(distancePtr, height, width, i, jminus);
		y = selectY(distancePtr, height, width, i, jminus);
		size_t indxWest = x_new(i, jminus, height);
		coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		distancePtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		x = selectX(distancePtr, height, width, iminus, j);
		y = selectY(distancePtr, height, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, height);
		coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		distancePtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < height - 1) {
		size_t iplus = i + 1;
		x = selectX(distancePtr, height, width, iplus, j);
		y = selectY(distancePtr, height, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, height);
		coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		distancePtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;
	short label = 0;
	
	while (inProcess.size() != 0) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, height);
		labelArray[currentIndx] = 1;
		distancePtr[currentIndx] = current.arrival;
		deleteRootHeap2D(inProcess);

		//East
		if (j < width - 1 && i >= 0 && i < height) {
			size_t jplus = j + 1;
			size_t indxEast = x_new(i, jplus, height);
			label = labelArray[indxEast];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jplus);
				y = selectY(distancePtr, height, width, i, jplus);
				coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					pointFastMarching2D EastNeighbor = { i, jplus, dEast };
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distancePtr[indxEast]) {
							distancePtr[indxEast] = dEast;
						}
					}
				}
			}
		}

		//West
		if (j > 0 && i >= 0 && i < height) {
			size_t jminus = j - 1;
			size_t indxWest = x_new(i, jminus, height);
			label = labelArray[indxWest];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jminus);
				y = selectY(distancePtr, height, width, i, jminus);
				coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					pointFastMarching2D WestNeighbor = { i, jminus, dWest };
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distancePtr[indxWest]) {
							distancePtr[indxWest] = dWest;
						}
					}
				}
			}
		}

		//North
		if (j >= 0 && j < width && i > 0) {
			size_t iminus = i - 1;
			size_t indxNorth = x_new(iminus, j, height);
			label = labelArray[indxNorth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iminus, j);
				y = selectY(distancePtr, height, width, iminus, j);
				coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distancePtr[indxNorth]) {
							distancePtr[indxNorth] = dNorth;
						}
					}
				}
			}


		}

		//South
		if (j >= 0 && j < width && i < height - 1) {
			size_t iplus = i + 1;
			size_t indxSouth = x_new(iplus, j, height);
			label = labelArray[indxSouth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iplus, j);
				y = selectY(distancePtr, height, width, iplus, j);
				coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distancePtr[indxSouth]) {
							distancePtr[indxSouth] = dSouth;
						}
					}
				}
			}
		}

	}

	for (i = 0; i < dim2D; i++) {
		if (distancePtr[i] == BIG_VALUE) {
			distancePtr[i] = 0.0;
		}
	}

	delete[] labelArray;
}

bool partialFrontPropagation2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, Point2D* seedPoints) {

	size_t dim2D = height * width;

	short* labelArray = new short[dim2D] {0};

	if (imageDataPtr == NULL || distancePtr == NULL || potentialPtr == NULL || seedPoints == NULL || labelArray == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0;
	vector<pointFastMarching2D> inProcess;

	dataType x = 0.0, y = 0.0, coefSpeed = 0.0;

	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	size_t currentIndx = x_new(i, j, height);

	////Compute the potential function
	//computePotential(imageDataPtr, potentialPtr, height, width, seedPoints);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
	for (k = 0; k < dim2D; k++) {
		distancePtr[k] = INFINITY;
		labelArray[k] = 3;
	}

	distancePtr[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < height) {
		size_t jplus = j + 1;
		x = selectX(distancePtr, height, width, i, jplus);
		y = selectY(distancePtr, height, width, i, jplus);
		size_t indxEast = x_new(i, jplus, height);
		coefSpeed = potentialPtr[indxEast];
		dataType dEast = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		distancePtr[indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0 && i >= 0 && i < height) {
		size_t jminus = j - 1;
		x = selectX(distancePtr, height, width, i, jminus);
		y = selectY(distancePtr, height, width, i, jminus);
		size_t indxWest = x_new(i, jminus, height);
		coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		distancePtr[indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		x = selectX(distancePtr, height, width, iminus, j);
		y = selectY(distancePtr, height, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, height);
		coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		distancePtr[indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < height - 1) {
		size_t iplus = i + 1;
		x = selectX(distancePtr, height, width, iplus, j);
		y = selectY(distancePtr, height, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, height);
		coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		distancePtr[indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(inProcess);

	pointFastMarching2D current;
	short label = 0;

	size_t seedI = seedPoints[1].x, seedJ = seedPoints[1].y, seedIndex = x_new(seedI, seedJ, height);
	while (labelArray[seedIndex] != 1) {

		current = inProcess[0];
		i = current.x;
		j = current.y;
		currentIndx = x_new(i, j, height);
		labelArray[currentIndx] = 1;
		distancePtr[currentIndx] = current.arrival;
		deleteRootHeap2D(inProcess);

		//East
		if (j < width - 1 && i >= 0 && i < height) {
			size_t jplus = j + 1;
			size_t indxEast = x_new(i, jplus, height);
			label = labelArray[indxEast];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jplus);
				y = selectY(distancePtr, height, width, i, jplus);
				coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					pointFastMarching2D EastNeighbor = { i, jplus, dEast };
					addPointHeap2D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distancePtr[indxEast]) {
							distancePtr[indxEast] = dEast;
						}
					}
				}
			}
		}

		//West
		if (j > 0 && i >= 0 && i < height) {
			size_t jminus = j - 1;
			size_t indxWest = x_new(i, jminus, height);
			label = labelArray[indxWest];
			if (label != 1) {
				x = selectX(distancePtr, height, width, i, jminus);
				y = selectY(distancePtr, height, width, i, jminus);
				coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					pointFastMarching2D WestNeighbor = { i, jminus, dWest };
					addPointHeap2D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distancePtr[indxWest]) {
							distancePtr[indxWest] = dWest;
						}
					}
				}
			}
		}

		//North
		if (j >= 0 && j < width && i > 0) {
			size_t iminus = i - 1;
			size_t indxNorth = x_new(iminus, j, height);
			label = labelArray[indxNorth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iminus, j);
				y = selectY(distancePtr, height, width, iminus, j);
				coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
					addPointHeap2D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distancePtr[indxNorth]) {
							distancePtr[indxNorth] = dNorth;
						}
					}
				}
			}


		}

		//South
		if (j >= 0 && j < width && i < height - 1) {
			size_t iplus = i + 1;
			size_t indxSouth = x_new(iplus, j, height);
			label = labelArray[indxSouth];
			if (label != 1) {
				x = selectX(distancePtr, height, width, iplus, j);
				y = selectY(distancePtr, height, width, iplus, j);
				coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed);
				if (label == 3) {
					distancePtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
					addPointHeap2D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distancePtr[indxSouth]) {
							distancePtr[indxSouth] = dSouth;
						}
					}
				}
			}
		}

	}

	for (i = 0; i < dim2D; i++) {
		if (distancePtr[i] == BIG_VALUE) {
			distancePtr[i] = 0.0;
		}
	}

	delete[] labelArray;
}

bool shortestPath2d(Image_Data2D distanceFuncPtr, dataType* resultedPath, Point2D* seedPoints) {

	if (distanceFuncPtr.imageDataPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	const size_t height = distanceFuncPtr.height;
	const size_t width = distanceFuncPtr.width;
	size_t i, j, dim2D = height * width, maxIter = 10000;
	dataType tau = 0.8, dist_min = 0.0, tol = 1.0;

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

		iNew = iNew - tau * gradientVectorX[currentIndex];
		jNew = jNew - tau * gradientVectorY[currentIndex];
		
		Point2D current_point = { iNew, jNew };
		dist_min = getPoint2DDistance(current_point, seedPoints[0]);

		i = (size_t)round(iNew); 
		j = (size_t)round(jNew);
		currentIndex = x_new(i, j, height);

		count_iter++;
	}
	while(dist_min > tol && count_iter < maxIter);

	delete[] gradientVectorX;
	delete[] gradientVectorY;

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

	dataType i_minus, i_plus;

	if (i == 0) {
		i_minus = INFINITY;
	}
	else {
		i_minus = actionPtr[k][x_new(i - 1, j, length)];
	}

	if (i == length - 1) {
		i_plus = INFINITY;
	}
	else {
		i_plus = actionPtr[k][x_new(i + 1, j, length)];
	}

	return min(i_minus, i_plus);
}

dataType select3dY(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {

	dataType j_minus, j_plus;

	if (j == 0) {
		j_minus = INFINITY;
	}
	else {
		j_minus = actionPtr[k][x_new(i, j - 1, length)];
	}

	if (j == width - 1) {
		j_plus = INFINITY;
	}
	else {
		j_plus = actionPtr[k][x_new(i, j + 1, length)];
	}

	return min(j_minus, j_plus);
}

dataType select3dZ(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k) {

	dataType k_minus, k_plus;

	size_t xd = x_new(i, j, length);

	if (k == 0) {
		k_minus = INFINITY;
	}
	else {
		k_minus = actionPtr[k - 1][xd];
	}

	if (k == height - 1) {
		k_plus = INFINITY;
	}
	else {
		k_plus = actionPtr[k + 1][xd];
	}

	return min(k_minus, k_plus);
}

bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height, VoxelSpacing spacing) {

	if (imageDataPtr == NULL || gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, currentInd = 0;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;

	dataType sx = spacing.sx, sy = spacing.sy, sz = spacing.sz;

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {

				currentInd = x_new(i, j, length);

				if (k == 0) {
					uz = (imageDataPtr[k + 1][currentInd] - imageDataPtr[k][currentInd]) / sz;
				}
				else {
					if (k == (height - 1)) {
						uz = (imageDataPtr[k][currentInd] - imageDataPtr[k - 1][currentInd]) / sz;
					}
					else {
						uz = (imageDataPtr[k + 1][currentInd] - imageDataPtr[k - 1][currentInd]) / (2 * sz);
					}
				}

				if (i == 0) {
					ux = (imageDataPtr[k][x_new(i + 1, j, length)] - imageDataPtr[k][x_new(i, j, length)]) / sx;
				}
				else {
					if (i == (length - 1)) {
						ux = (imageDataPtr[k][x_new(i, j, length)] - imageDataPtr[k][x_new(i - 1, j, length)]) / sx;
					}
					else {
						ux = (imageDataPtr[k][x_new(i + 1, j, length)] - imageDataPtr[k][x_new(i - 1, j, length)]) / (2 * sx);
					}
				}

				if (j == 0) {
					uy = (imageDataPtr[k][x_new(i, j + 1, length)] - imageDataPtr[k][x_new(i, j, length)]) / sy;
				}
				else {
					if (j == (width - 1)) {
						uy = (imageDataPtr[k][x_new(i, j, length)] - imageDataPtr[k][x_new(i, j - 1, length)]) / sy;
					}
					else {
						uy = (imageDataPtr[k][x_new(i, j + 1, length)] - imageDataPtr[k][x_new(i, j - 1, length)]) / (2 * sy);
					}
				}

				gradientVectorX[k][currentInd] = (dataType)ux;
				gradientVectorY[k][currentInd] = (dataType)uy;
				gradientVectorZ[k][currentInd] = (dataType)uz;
			}
		}
	}
	return true;
}

void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b) {
	pointFastMarching3D temp = *a;
	*a = *b;
	*b = temp;
}

void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i) {

	int length_array = in_Process.size();
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

void heapifyVector3D(vector<pointFastMarching3D>& in_Process) {
	int length_array = in_Process.size();
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		heapifyDown3D(in_Process, ind);
	}
}

void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap3dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown3D(in_Process, 0);
}

void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp3D(in_Process, l - 1);
}

void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i) {

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

bool compute3DPotential(Image_Data ctImageData, dataType** potential, Point3D* seedPoint, Potential_Parameters parameters) {

	if (ctImageData.imageDataPtr == NULL || potential == NULL) {
		return false;
	}
		
	size_t i = 0, k = 0;
	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t dim2D = length * width;

	dataType** gradientVectorX = new dataType * [height];
	dataType** gradientVectorY = new dataType * [height];
	dataType** gradientVectorZ = new dataType * [height];
	dataType** distance = new dataType * [height];
	dataType** edgeDetector = new dataType * [height];
	dataType** maskThreshold = new dataType * [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType[dim2D]{0};
		gradientVectorY[k] = new dataType[dim2D]{0};
		gradientVectorZ[k] = new dataType[dim2D]{0};
		distance[k] = new dataType[dim2D]{ 0 };
		edgeDetector[k] = new dataType[dim2D]{ 0 };
		maskThreshold[k] = new dataType[dim2D]{ 0 };
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL || distance == NULL || edgeDetector == NULL || maskThreshold == NULL)
		return false;
	
	dataType norm_of_gradient = 0.0;;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;
	
	
	compute3dImageGradient(ctImageData.imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, ctImageData.spacing);
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			
			ux = gradientVectorX[k][i];
			uy = gradientVectorY[k][i];
			uz = gradientVectorZ[k][i];
			norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz);
			edgeDetector[k][i] = gradientFunction(norm_of_gradient, parameters.K);
			
			//threshold
			if (edgeDetector[k][i] < parameters.thres) {
				maskThreshold[k][i] = 0.0;
			}
			else {
				maskThreshold[k][i] = 1.0;
			}
		}
	}
	fastSweepingFunction_3D(distance, maskThreshold, length, width, height, ctImageData.spacing.sx, 10000000.0, 0.0);

	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };
	seedStats = getPointNeighborhoodStats(ctImageData, seedPoint[0], parameters.radius);
	dataType value_first_pt = seedStats.mean_data;
	seedStats = getPointNeighborhoodStats(ctImageData, seedPoint[1], parameters.radius);
	dataType value_second_pt = seedStats.mean_data;

	dataType seedValCT = (value_first_pt + value_second_pt) / 2.0;
	//dataType seedValCT = value_first_pt;

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = fabs(seedValCT - ctImageData.imageDataPtr[k][i]);
		}
	}

	//Find the max of each potential for normalization
	dataType maxImage = 0.0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (potential[k][i] > maxImage) {
				maxImage = potential[k][i];
			}
		}
	}

	//Normalization
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			dataType weight_dist = 1.0 / (1.0 + distance[k][i]);
			potential[k][i] = parameters.eps + weight_dist * potential[k][i] / maxImage;
		}
	}

	for (k = 0; k < height; k++) {
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
		delete[] edgeDetector[k];
		delete[] maskThreshold[k];
		delete[] distance[k];
	}
	delete[] gradientVectorX;
	delete[] gradientVectorY;
	delete[] gradientVectorZ;
	delete[] edgeDetector;
	delete[] maskThreshold;
	delete[] distance;

	return true;
}

bool shortestPath3D(Image_Data actionMapStr, Point3D* seedPoints, dataType tau, dataType tolerance, vector<Point3D>& path_points) {

	if (actionMapStr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	const size_t length = actionMapStr.length;
	const size_t width = actionMapStr.width;
	const size_t height = actionMapStr.height;
	VoxelSpacing spacing = actionMapStr.spacing;

	size_t i = 0, j = 0, k = 0, xd = 0, dim2D = length * width;
	//size_t max_iter = 1000;
	//dataType tau = 0.8, tol = 1.0;

	//estimate the max number of iterations
	Point3D pEnd = getRealCoordFromImageCoord3D(seedPoints[0], actionMapStr.origin, spacing, actionMapStr.orientation);
	Point3D pStart = getRealCoordFromImageCoord3D(seedPoints[1], actionMapStr.origin, spacing, actionMapStr.orientation);
	double totalLenght = getPoint3DDistance(pEnd, pStart);
	size_t max_iter = 1000;//(size_t)(1.5 * totalLenght);
	//std::cout << "max number of iterations : " << max_iter << std::endl;

	dataType** gradientVectorX = new dataType * [height];
	dataType** gradientVectorY = new dataType * [height];
	dataType** gradientVectorZ = new dataType * [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType[dim2D];
		gradientVectorY[k] = new dataType[dim2D];
		gradientVectorZ[k] = new dataType[dim2D];
		if (gradientVectorX[k] == NULL || gradientVectorY[k] == NULL || gradientVectorZ[k] == NULL)
			return false;
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;

	//Normalization of the gradient
	compute3dImageGradient(actionMapStr.imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, spacing);

	dataType ux = 0.0, uy = 0.0, uz = 0.0, norm_of_gradient = 0.0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			ux = gradientVectorX[k][i];
			uy = gradientVectorY[k][i];
			uz = gradientVectorZ[k][i];
			norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz);
			if (norm_of_gradient != 0) {
				gradientVectorX[k][i] = (dataType)(ux / norm_of_gradient);
				gradientVectorY[k][i] = (dataType)(uy / norm_of_gradient);
				gradientVectorZ[k][i] = (dataType)(uz / norm_of_gradient);
			}
			else {
				gradientVectorX[k][i] = 0;
				gradientVectorY[k][i] = 0;
				gradientVectorZ[k][i] = 0;
			}
		}
	}

	//Find the closest point till the last point
	i = (size_t)seedPoints[1].x;
	j = (size_t)seedPoints[1].y;
	k = (size_t)seedPoints[1].z;
	size_t currentIndx = x_new(i, j, length);

	dataType iNew = seedPoints[1].x;
	dataType jNew = seedPoints[1].y;
	dataType kNew = seedPoints[1].z;
	double dist_to_end = 0.0;

	size_t count_iter = 1;

	Point3D final_point = getRealCoordFromImageCoord3D(seedPoints[0], actionMapStr.origin, spacing, actionMapStr.orientation);

	do {

		currentIndx = x_new(i, j, length);

		iNew = iNew - tau * gradientVectorX[k][currentIndx];
		jNew = jNew - tau * gradientVectorY[k][currentIndx];
		kNew = kNew - tau * gradientVectorZ[k][currentIndx];
		
		Point3D point_current = { iNew, jNew, kNew };
		Point3D pDistance = getRealCoordFromImageCoord3D(point_current, actionMapStr.origin, spacing, actionMapStr.orientation);
		path_points.push_back(point_current);
		
		//compute distance current Point - last point
		dist_to_end = getPoint3DDistance(pDistance, final_point);
		
		i = (size_t)(round(iNew));
		j = (size_t)(round(jNew));
		k = (size_t)(round(kNew));
		
		count_iter++;

	} while (dist_to_end > tolerance && count_iter < max_iter);

	for (k = 0; k < height; k++) {
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
	}
	delete[] gradientVectorX;
	delete[] gradientVectorY;
	delete[] gradientVectorZ;

	return true;
}

bool findPathFromOneGivenPointWithCircleDetection(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters, size_t stop_criterium) {

	if (ctImageData.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	size_t k = 0, i = 0, j = 0, x = 0;

	const size_t height = ctImageData.height;
	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t dim2D = length * width;

	double offset = 5.0;
	BoundingBox2D box = { 0, 0, 0, 0 };

	string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, extension, slice_number;

	saving_name = path_name + "path_points.csv";
	FILE* file;
	if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType** actionPtr = new dataType * [height];
	dataType** maskAction = new dataType * [height];
	dataType** potentialPtr = new dataType * [height];
	for (k = 0; k < height; k++) {
		actionPtr[k] = new dataType[dim2D]{ 0 };
		maskAction[k] = new dataType[dim2D]{ 0 };
		potentialPtr[k] = new dataType[dim2D]{ 0 };
	}
	if (actionPtr == NULL || maskAction == NULL || potentialPtr == NULL)
		return false;

	//compute3DPotential(ctImageData, potentialPtr, seedPoints[0], parameters);
	saving_name = path_name + "potential.raw";
	manageRAWFile3D<dataType>(potentialPtr, length, width, height, saving_name.c_str(), LOAD_DATA, false);

	const double step = 70.0;

	//================ find next point inside the aorta ============

	Point3D* seeds = new Point3D[2];
	seeds[0] = seedPoints[0];

	Image_Data action
	{
		height,
		length,
		width,
		actionPtr,
		ctImageData.origin,
		ctImageData.spacing,
		ctImageData.orientation
	};

	dataType tau = 0.95 * ctImageData.spacing.sx;
	dataType tolerance = ctImageData.spacing.sx;

	//During the first propagation the mask is empty
	partialPropagation(action, potentialPtr, maskAction, length, width, height, seeds, step);
	//fastMarching3D_N(actionPtr, potentialPtr, length, width, height, initial_point_img_coord);

	saving_name = path_name + "action_initial.raw";
	manageRAWFile3D<dataType>(actionPtr, length, width, height, saving_name.c_str(), STORE_DATA, false);

	//vector for path points
	vector<Point3D> path_points;

	////path between initial point and second point (found point)
	//shortestPath3D(actionPtr, length, width, height, ctImageData.spacing, seeds, path_points);
	shortestPath3D(action, seeds, tau, tolerance, path_points);

	//write path points coordinates to file
	int n = 0;
	Point3D point_file = { 0.0, 0.0, 0.0 };
	for (n = path_points.size() - 1; n > -1; n--) {
		point_file = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	}

	////empty the vector contaning path points
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}

	fclose(file);

	
	//Update the mask
	copyDataToAnotherArray(actionPtr, maskAction, height, length, width);

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
	

	//=========== Hough transform parameters ==================
	
	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* houghSpace = new dataType[dim2D]{ 0 };
	dataType* voteArray = new dataType[dim2D]{ 0 };

	//Slice extraction
	
	////filtering parameters
	//Filter_Parameters filter_parameters;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	Point2D seed2D = { 0.0, 0.0 };
	Point2D center2D = { 0.0, 0.0 };
	PixelSpacing spacing = { ctImageData.spacing.sx, ctImageData.spacing.sy };
	Image_Data2D IMAGE
	{
		length,
		width,
		imageSlice,
		center2D,
		spacing
	};

	//HoughParameters hParameters = { 6.0, 19.0, 0.5, 5.0, 1000, 1.0, 0.2, 0.5, spacing };
	//dataType max_ratio = 1.0;

	//================ Use circle detection as path finding stopping criterium =====================
	
	
	//size_t num_slice_processed = 0;
	//int len = path_points.size();
	//std::cout << len << " point to be processed" << std::endl;
	//
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

	delete[] seeds;
	for (k = 0; k < height; k++) {
		delete[] actionPtr[k];
		delete[] maskAction[k];
		delete[] potentialPtr[k];
	}
	delete[] actionPtr;
	delete[] maskAction;
	delete[] potentialPtr;

	return true;
}

bool findPathTwoSteps(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters) {

	size_t i, j, k, xd;
	size_t length = ctImageData.length, width = ctImageData.width, height = ctImageData.height;
	size_t dim2D = length * width;

	dataType** potential = new dataType*[height];
	dataType** action_field = new dataType *[height];
	dataType** mask_action = new dataType *[height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		action_field[k] = new dataType[dim2D]{ 0 };
		mask_action[k] = new dataType[dim2D]{ 0 };
	}

	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, saving_csv = outputPath + "path_point.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	Point3D temporary_point = { 0.0, 0.0, 0.0 };
	Point3D* seedsPath = new Point3D[2];
	seedsPath[0] = seedPoints[0];
	seedsPath[1] = seedPoints[1];
	vector<Point3D> path_points;

	dataType tau = 0.95 * ctImageData.spacing.sx;
	dataType tolerance = ctImageData.spacing.sx;

	//======== First step ===========//
	//compute3DPotential(ctImageData, potential, seedsPath[0], parameters);

	saving_name = outputPath + "potential.raw";
	manageRAWFile3D<dataType>(potential, length, width, height, saving_name.c_str(), LOAD_DATA, false);

	partialFrontPropagation(action_field, potential, length, width, height, seedsPath);
	//fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);

	Image_Data actionDataStr = { height, length, width, action_field, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };
	//shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);
	
	//copy points to file
	int n = 0;
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}
	
	//================================//

	Point3D seed1 = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	Point3D seed2 = getRealCoordFromImageCoord3D(seedPoints[2], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	double step = 0.5 * getPoint3DDistance(seed1, seed2), dist = 0.0;
	cout << "Step : " << step << endl;

	////======== Second step ===========//
	//compute3DPotential(ctImageData, potential, seedsPath[1], parameters);
	fastMarching3D_N(action_field, potential, length, width, height, seedsPath[1]);

	dataType min_action = BIG_VALUE, value_temp = 0.0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(i, j, length);
				Point3D current_point = {(dataType)i, (dataType)j, (dataType)k};
				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
				dist = getPoint3DDistance(seed1, current_point);
				if (dist >= (step - 0.01) && dist <= (step + 0.01) && k >= seedPoints[1].z) {
					if (action_field[k][xd] < min_action) {
						min_action = action_field[k][xd];
						temporary_point.x = (dataType)i;
						temporary_point.y = (dataType)j;
						temporary_point.z = (dataType)k;
						value_temp = action_field[k][xd];
					}
				}
			}
		}
	}
	cout << "found point : (" << temporary_point.x << " ," << temporary_point.y << " , " << temporary_point.z << ")" << endl;
	
	seedsPath[0] = seedPoints[1];
	seedsPath[1] = temporary_point;
	//shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);

	for (n = path_points.size() - 1; n > - 1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	//======= Last step ===============//
	seedsPath[0] = temporary_point;
	seedsPath[1] = seedPoints[2];
	
	partialFrontPropagation(action_field, potential, length, width, height, seedsPath);
	//fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);
	
	//shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	shortestPath3D(actionDataStr, seedsPath, tau, tolerance, path_points);
	
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}
	//=====================================
	
	fclose(file);
	delete[] seedsPath;
	for (k = 0; k < height; k++) {
		delete[] potential[k];
		delete[] action_field[k];
		delete[] mask_action[k];
	}
	delete[] potential;
	delete[] action_field;
	delete[] mask_action;

	return true;
}

bool findPathFromOneGivenPoint(Image_Data ctImageData, Point3D * seedPoints, Potential_Parameters parameters) {

	if (ctImageData.imageDataPtr == NULL)
		return false;

	size_t k = 0, i = 0, j = 0, x = 0;
	const size_t height = ctImageData.height, length = ctImageData.length, width = ctImageData.width;
	const size_t dim2D = length * width;

	double epsilon = 0.01;

	string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, extension;

	saving_name = path_name + "path_points_v2.csv";
	FILE* file;
	if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType** actionPtr = new dataType * [height];
	dataType** newActionPtr = new dataType * [height];
	dataType** maskDistance = new dataType * [height];
	dataType** potentialPtr = new dataType * [height];
	for (k = 0; k < height; k++) {
		actionPtr[k] = new dataType[dim2D]{ 0 };
		newActionPtr[k] = new dataType[dim2D]{ 0 };
		maskDistance[k] = new dataType[dim2D]{ 0 };
		potentialPtr[k] = new dataType[dim2D]{ 0 };
	}
	if (actionPtr == NULL || newActionPtr == NULL || maskDistance == NULL || potentialPtr == NULL)
		return false;

	Point3D* seeds = new Point3D[2];
	seeds[0] = seedPoints[0];
	Point3D initial_point = seedPoints[0];
	Point3D temporary_point = { 0.0, 0.0, 0.0 };

	dataType min_distance = BIG_VALUE, value_temp = 0.0, h = 1.0;

	Point3D seed1 = getRealCoordFromImageCoord3D(seedPoints[0], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	Point3D seed2 = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	
	double distance_point = getPoint3DDistance(seed1, seed2);
	double step = 50;//0.5 * distance_point;
	std::cout << "First step : " << step << std::endl;

	//================ find next point inside the aorta 

	Image_Data actionDataStr = { height, length, width, actionPtr, ctImageData.origin, ctImageData.spacing, ctImageData.orientation };

	double radius_to_find_mean = 3.0;
	//compute3DPotential(ctImageData, potentialPtr, seedPoints[0], parameters);

	saving_name = path_name + "potential.raw";
	//store3dRawData<dataType>(potentialPtr, length, width, height, saving_name.c_str());
	manageRAWFile3D<dataType>(potentialPtr, length, width, height, saving_name.c_str(), LOAD_DATA, false);

	//fastMarching3D_N(actionPtr, potentialPtr, length, width, height, seeds[0]);
	fastMarching3dWithSpacing(ctImageData, actionPtr, potentialPtr, seeds[0], ctImageData.spacing);
	//saving_name = path_name + "first_action.raw";
	//manageRAWFile3D<dataType>(actionPtr, length, width, height, saving_name.c_str(), LOAD_DATA, false);

	min_distance = BIG_VALUE;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				x = x_new(i, j, length);
				Point3D current_point = { (dataType)i, (dataType)j, (dataType)k};
				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
				double distanceP = getPoint3DDistance(seed1, current_point);
				if (distanceP >= (step - epsilon) && distanceP <= (step + epsilon)) {
					if (actionPtr[k][x] < min_distance) {
						min_distance = actionPtr[k][x];
						temporary_point.x = (dataType)i;
						temporary_point.y = (dataType)j;
						temporary_point.z = (dataType)k;
						value_temp = actionPtr[k][x];
					}
				}
			}
		}
	}
	
	seed1 = getRealCoordFromImageCoord3D(temporary_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//distance_point = getPoint3DDistance(seed1, seed2);
	//step = 0.5 * distance_point;
	//std::cout << "Second step: " << step << std::endl;
	std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;

	//vector for path points
	vector<Point3D> path_points;

	//path between initial point and second point
	seeds[1] = temporary_point;
	
	dataType tau = 0.95 * ctImageData.spacing.sx;
	dataType tolerance = ctImageData.spacing.sx;
	shortestPath3D(actionDataStr, seeds, tau, tolerance, path_points);

	//write just path points coordinates to file
	int i_n = 0, k_n = 0, k_center = 0;
	Point3D point_file = { 0.0, 0.0, 0.0 };
	
	for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
		point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}

	//============================================

	double distance_to_stop = 0.1 * distance_point;

	size_t count_step = 0;

	
	while (count_step < 3) {
		
		count_step++;
		std::cout << "count step = " << count_step << std::endl;

		//mask
		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {
					x = x_new(i, j, length);
					Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
					current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
					double distanceP = getPoint3DDistance(seed1, current_point);
					if (distanceP <= (step + 0.01) && actionPtr[k][x] <= value_temp) {
						maskDistance[k][x] = 10000.0;
					}
					else {
						if (maskDistance[k][x] != 10000.0) {
							maskDistance[k][x] = newActionPtr[k][x];
						}
					}
				}
			}
		}
		seeds[0] = seeds[1];
		initial_point = seeds[1];

		//fastMarching3D_N(newActionPtr, potentialPtr, length, width, height, seeds[0]);
		fastMarching3dWithSpacing(ctImageData, newActionPtr, potentialPtr, seeds[0], ctImageData.spacing);

		//saving_name = path_name + "action_field" + extension + ".raw";
		//store3dRawData<dataType>(newActionPtr, length, width, height, saving_name.c_str());

		//Find the temporary point
		min_distance = BIG_VALUE;

		//get real world coordinate of the initial point
		seed1 = getRealCoordFromImageCoord3D(initial_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);

		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {
					x = x_new(i, j, length);
					Point3D current_point = { (dataType)i, (dataType)j, current_point.z = (dataType)k };
					current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
					double distanceP = getPoint3DDistance(seed1, current_point);
					if (distanceP >= (step - 0.01) && distanceP <= (step + 0.01) && maskDistance[k][x] != 10000.0) {
						if (newActionPtr[k][x] < min_distance) {
							min_distance = newActionPtr[k][x];
							temporary_point.x = (dataType)i;
							temporary_point.y = (dataType)j;
							temporary_point.z = (dataType)k;
							value_temp = newActionPtr[k][x];
						}
					}
				}
			}
		}
		
		std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;

		seeds[1] = temporary_point;
		//seed2 = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		//distance_point = getPoint3DDistance(seed1, seed2);
		//step = 0.5 * distance_point;
		//std::cout << "The step is : " << step << std::endl;

		actionDataStr.imageDataPtr = newActionPtr;
		shortestPath3D(actionDataStr, seeds, tau, tolerance, path_points);

		copyDataToAnotherArray(newActionPtr, actionPtr, height, length, width);

		//write just path points coordinates to file
		for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
			point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
			fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
		}
		while (path_points.size() != 0) {
			path_points.pop_back();
		}

	}
	

	fclose(file);

	delete[] seeds;

	for (k = 0; k < height; k++) {
		delete[] actionPtr[k];
		delete[] newActionPtr[k];
		delete[] maskDistance[k];
		delete[] potentialPtr[k];
	}
	delete[] actionPtr;
	delete[] newActionPtr;
	delete[] maskDistance;
	delete[] potentialPtr;

	return true;
}

bool partialFrontPropagation(dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints) {

	if (distanceFuncPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	short** labelArray = new short *[height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new short[dim2D]{0};
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
	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	k = (size_t)seedPoints[0].z;
	if (i < 0 || i > length || j < 0 || j > width || k < 0 || k > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	currentIndx = x_new(i, j, length);
	distanceFuncPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;
	pointFastMarching3D current = {i, j, k, 0.0};

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

	size_t seedI = (size_t)seedPoints[1].x;
	size_t seedJ = (size_t)seedPoints[1].y;
	size_t seedK = (size_t)seedPoints[1].z;
	size_t seedIndex = x_new(seedI, seedJ, length);
	if (seedI < 0 || seedI > length || seedJ < 0 || seedJ > width || seedK < 0 || seedK > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}

	while (labelArray[seedK][seedIndex] != 1) {

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
		for (i = 0; i < dim2D; i++) {
			if (distanceFuncPtr[k][i] == INFINITY) {
				distanceFuncPtr[k][i] = 0.0;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

bool partialPropagation(Image_Data actionPtr, dataType** potentialPtr, dataType** maskPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints, const double maxLength) {

	if (actionPtr.imageDataPtr == NULL || potentialPtr == NULL || maskPtr == NULL) {
		return false;
	}

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
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

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

//=====================================

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

	//string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/radius.raw";
	//store2dRawData<dataType>(radiusArray, height, width, outputPath.c_str());

	delete[] radiusArray;
}

dataType min0(dataType x, dataType y) {
	if (y - x > 0)
		return pow(x - y, 2);
	else
		return 0;
}

bool rouyTourinDistanceMap(Image_Data ctImageData, dataType** distancePtr, dataType foregroundValue, dataType tolerance, dataType tau) {
	
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
	}

	double mass = 10.0;
	dataType hx = ctImageData.spacing.sx;
	dataType hy = ctImageData.spacing.sy;
	dataType hz = ctImageData.spacing.sz;
	dataType value = 0.0;

	size_t count_iteration = 0;
	size_t count_point_to_processed = 0;
	size_t max_iteration = 1000;

	while (mass > tolerance && count_iteration < max_iteration) {

		copyDataToExtendedArea(distancePtr, previousSolution, height, length, width);
		reflection3D(previousSolution, height_ext, length_ext, width_ext);

		count_iteration++;
		mass = 0.0;
		for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
			for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
				for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
					if (ctImageData.imageDataPtr[k][x_new(i, j, length)] == foregroundValue) {
						if (count_iteration == 1) {
							count_point_to_processed++;
						}
						value = previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)];
						distancePtr[k][x_new(i, j, length)] = value + tau - tau * sqrt((1.0 / pow(hx, 2)) * max(min0(previousSolution[k_ext][x_new(i_ext - 1, j_ext, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext + 1, j_ext, length_ext)], value))
							+ (1.0 / pow(hy, 2)) * max(min0(previousSolution[k_ext][x_new(i_ext, j_ext - 1, length_ext)], value), min0(previousSolution[k_ext][x_new(i_ext, j_ext + 1, length_ext)], value))
							+ (1.0 / pow(hz, 2)) * max(min0(previousSolution[k_ext - 1][x_new(i_ext, j_ext, length_ext)], value), min0(previousSolution[k_ext + 1][x_new(i_ext, j_ext, length_ext)], value)));

						//Compute the mass
						mass += pow(previousSolution[k_ext][x_new(i_ext, j_ext, length_ext)] - distancePtr[k][x_new(i, j, length)], 2);
					}
				}
			}
		}
		mass = sqrt(mass);
		
		//if (count_iteration % 10 == 0) {
		//	std::cout << "Mass = " << mass << std::endl;
		//	std::cout << count_iteration << " iterations" << std::endl;
		//}

	}
	//std::cout << count_point_to_processed << " points were processed" << std::endl;

	for (k = 0; k < height_ext; k++) {
		delete[] previousSolution[k];
	}
	delete[] previousSolution;

	return true;
}

//======================================

dataType solve3dQuadraticFastMarching(dataType X, dataType Y, dataType Z, dataType W, VoxelSpacing h) {

	dataType sol = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType R = sqrt(W);

	if (X == INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 1.0 / pow(h.sy,2) + 1.0 / pow(h.sz,2);
		b = (dataType)(-2 * (Y / pow(h.sy,2) + Z / pow(h.sz,2)));
		c = (dataType)(pow(Y / h.sy, 2) + pow(Z / h.sz, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(Y, Z) + R);
		}
	}

	if (Y == INFINITY && X != INFINITY && Z != INFINITY) {
		a = 1.0 / pow(h.sx, 2) + 1.0 / pow(h.sz, 2);
		b = (dataType)(-2 * (X / pow(h.sx,2) + Z / pow(h.sz,2)));
		c = (dataType)(pow(X / h.sx, 2) + pow(Z / h.sz, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, Z) + R);
		}
	}

	if (Z == INFINITY && X != INFINITY && Y != INFINITY) {
		a = 1.0 / pow(h.sx, 2) + 1.0 / pow(h.sy, 2);
		b = (dataType)(-2 * (X / pow(h.sx, 2) + Y / pow(h.sy, 2)));
		c = (dataType)(pow(X / h.sx, 2) + pow(Y / h.sy, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, Y) + R);
		}
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY) {
		a = 1.0 / pow(h.sz, 2);
		b = (dataType)(-2 * Z / pow(h.sz,2));
		c = (dataType)(pow(Z / h.sz, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(Z + R);
		}
	}

	if (X == INFINITY && Z == INFINITY && Y != INFINITY) {
		a = 1.0 / pow(h.sy, 2);
		b = (dataType)(-2 * Y / pow(h.sy, 2));
		c = (dataType)(pow(Y / h.sy, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(Y + R);
		}
	}

	if (Y == INFINITY && Z == INFINITY && X != INFINITY) {
		a = 1.0 / pow(h.sx, 2);
		b = (dataType)(-2 * X / pow(h.sx, 2));
		c = (dataType)(pow(X / h.sx, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(X + R);
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 1.0 / pow(h.sx, 2) + 1.0 / pow(h.sy, 2) + 1.0 / pow(h.sz, 2);
		b = (dataType)( -2 * (X / pow(h.sx, 2) + Y / pow(h.sy, 2) + Z / pow(h.sz, 2)) );
		c = (dataType)(pow(X / h.sx, 2) + pow(Y / h.sy, 2) + pow(Z / h.sz, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			return (dataType)((-b + sqrt(delta)) / (2 * a));
		}
		else {
			return (dataType)(min(X, min(Y, Z)) + R);
		}
	}

}

bool fastMarching3dWithSpacing(Image_Data ctImageData, dataType** distanceFuncPtr, dataType** potentialFuncPtr, Point3D seedPoint, VoxelSpacing spacing) {

	if (distanceFuncPtr == NULL || potentialFuncPtr == NULL) {
		return false;
	}

	const size_t length = ctImageData.length;
	const size_t width = ctImageData.width;
	const size_t height = ctImageData.height;

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0, dim2D = length * width;

	size_t** labelArray = new size_t * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new size_t[dim2D]{ 0 };
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
		//dTop = solve3dQuadratic(x, y, z, coefSpeed);
		dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		//dEast = solve3dQuadratic(x, y, z, coefSpeed);
		dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		//dWest = solve3dQuadratic(x, y, z, coefSpeed);
		dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dTop = solve3dQuadratic(x, y, z, coefSpeed);
				dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dEast = solve3dQuadratic(x, y, z, coefSpeed);
				dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
				dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dWest = solve3dQuadratic(x, y, z, coefSpeed);
				dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
				dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
				dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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

dataType FastSweepingUpdate(dataType** distancePtr, size_t length, size_t width, size_t height, int indx, int indy, int indz, VoxelSpacing h) {

	dataType x = select3dX(distancePtr, length, width, height, indx, indy, indz);
	dataType y = select3dY(distancePtr, length, width, height, indx, indy, indz);
	dataType z = select3dZ(distancePtr, length, width, height, indx, indy, indz);
	
	dataType dist = 0.0;
	dataType dist_new = 0.0;

	dataType a = -max(-x, max(-y, -z));
	dataType c = max(x, max(y, z));
	dataType b = (x + y + z) - (a + c);

	dist_new = a + h.sx;

	return dist;
}

bool penalizedFrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialFuncPtr, Point3D* seedPoints) {

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
		//dTop = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		//actionMapPtr[kminus][currentIndx] = dTop;
		//inProcess.push_back(TopNeighbor);
		//labelArray[kminus][currentIndx] = 2;
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
		//dEast = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		//actionMapPtr[k][indxEast] = dEast;
		//inProcess.push_back(EastNeighbor);
		//labelArray[k][indxEast] = 2;
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
		//dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		//actionMapPtr[k][indxNorth] = dNorth;
		//inProcess.push_back(NorthNeighbor);
		//labelArray[k][indxNorth] = 2;
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
		//dWest = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		//actionMapPtr[k][indxWest] = dWest;
		//inProcess.push_back(WestNeighbor);
		//labelArray[k][indxWest] = 2;
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
		//dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		//actionMapPtr[k][indxSouth] = dSouth;
		//inProcess.push_back(SouthNeighbor);
		//labelArray[k][indxSouth] = 2;
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
		//dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		//pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		//actionMapPtr[kplus][currentIndx] = dBottom;
		//inProcess.push_back(BottomNeighbor);
		//labelArray[kplus][currentIndx] = 2;
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

//=======================================

bool frontPropagationWithKeyPointDetection(Image_Data actionMapStr, dataType** potentialFuncPtr, Point3D* seedPoint, const double LengthKeyPoints, vector<Point3D>& key_points) {

	if (actionMapStr.imageDataPtr == NULL || potentialFuncPtr == NULL) {
		return false;
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
	if (labelArray == NULL)
		return false;

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
	i = (size_t)seedPoint[0].x;
	j = (size_t)seedPoint[0].y;
	k = (size_t)seedPoint[0].z;
	if (i < 0 || i > length || j < 0 || j > width || k < 0 || k > height) {
		std::cout << "Error in the input seed point" << std::endl;
		return false;
	}
	currentIndx = x_new(i, j, length);
	actionMapStr.imageDataPtr[k][currentIndx] = 0.0;
	labelArray[k][currentIndx] = 1;

	double distanceToCurrentSourcePoint = 0.0;
	Point3D currentSourcePoint = { i, j, k };
	key_points.push_back(currentSourcePoint);

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
		dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
		dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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

	size_t iEnd = (size_t)seedPoint[1].x;
	size_t jEnd = (size_t)seedPoint[1].y;
	size_t kEnd = (size_t)seedPoint[1].z;
	
	while (inProcess.size() != 0 && labelArray[kEnd][x_new(iEnd, jEnd, length)] != 1) {

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
			if(i < length_minus && j >= 0 && j < width && k >= 0 && k < height)
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
				dTop = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				dEast = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				dNorth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				dWest = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				dSouth = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
				dBottom = solve3dQuadraticFastMarching(x, y, z, coefSpeed, spacing);
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
	
	key_points.push_back(seedPoint[1]);
	std::cout << "The source Points are: " << std::endl;
	for (int nb = 0; nb < key_points.size(); nb++) {
		std::cout << "Point " << nb + 1 << " : " << key_points[nb].x << " " << key_points[nb].y << " " << key_points[nb].z << std::endl;
	}

	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (actionMapStr.imageDataPtr[k][i] == INFINITY) {
	//			actionMapStr.imageDataPtr[k][i] = max_weighted_distance + 1;
	//		}
	//	}
	//}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}