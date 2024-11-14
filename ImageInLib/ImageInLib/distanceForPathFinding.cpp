#include <iostream>
#include <sstream>  
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include<tuple>
#include "distanceForPathFinding.h"
#include<template_functions.h>
#include"filtering.h"
#include"../src/heat_equation.h"
#include"hough_transform.h"
#include"../src/distance_function.h"

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

bool computePotential(dataType* imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D * seedPoints)
{
	//This function is used to compute the potential function
	//were epsilon can be used also as parameter
	//epislon = 0.01 was perfect for our experimentations

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, dim2D = height * width;
	size_t i1 = (size_t)seedPoints[0].x, j1 = (size_t)seedPoints[0].y;
	size_t i2 = (size_t)seedPoints[1].x, j2 = (size_t)seedPoints[1].y;

	dataType seedVal = (dataType)((imageDataPtr[x_new(i1, j1, height)] + imageDataPtr[x_new(i2, j2, height)]) / 2.0);
	dataType epsilon = 0.01;
	dataType K = 0.00005;
	dataType coef_dist = 1000, coef = 0.0;

	dataType* gradientVectorX = new dataType[dim2D]{ 0 };
	dataType* gradientVectorY = new dataType[dim2D]{ 0 };
	dataType* edgeDetector = new dataType[dim2D]{ 0 };

	computeImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, height, width, 1.0);

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
	computePotential(imageDataPtr, potentialPtr, height, width, seedPoints);

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> in process and 3 ---> not processed
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
	computePotential(imageDataPtr, potentialPtr, height, width, seedPoints);

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

bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, saving_csv = outputPath + "path2D.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	size_t i, j, dim2D = height * width, maxIter = 10000;
	dataType tau = 0.8, dist_min = 0.0, tol = 1.0;

	dataType* gradientVectorX = new dataType[dim2D]{0};
	dataType* gradientVectorY = new dataType[dim2D]{0};

	//Normalization of the gradient
	computeImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY, height, width, 1.0);
	
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
	//resultedPath[currentIndex] = 1.0;

	dataType iNew = seedPoints[1].x, jNew = seedPoints[1].y;
	
	dataType radius = 0.0, offset = 2.0;
	BoundingBox2D box = { 0, 0, 0, 0 };
	size_t n, m;

	fprintf(file, "x,y\n");

	do{

		iNew = iNew - tau * gradientVectorX[currentIndex];
		jNew = jNew - tau * gradientVectorY[currentIndex];
		fprintf(file, "%f,%f\n", iNew, jNew);
		Point2D current_point = { iNew, jNew };
		dist_min = getPoint2DDistance(current_point, seedPoints[0]);

		i = (size_t)round(iNew); 
		j = (size_t)round(jNew);
		currentIndex = x_new(i, j, height);

		count_iter++;
	}
	while(dist_min > tol && count_iter < maxIter);

	fclose(file);

	//delete[] radiusArray;
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

	//dataType sx = 1.0, sy = 1.0, sz = 1.0;
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

bool compute3DPotential(Image_Data ctImageData, dataType** potential, Point3D seedPoint, double radius, Potential_Parameters parameters) {

	if (ctImageData.imageDataPtr == NULL || potential == NULL) {
		return false;
	}
		
	size_t i = 0, k = 0;
	const size_t height = ctImageData.height, length = ctImageData.length, width = ctImageData.width;
	const size_t dim2D = length * width;

	dataType** gradientVectorX = new dataType * [height];
	dataType** gradientVectorY = new dataType * [height];
	dataType** gradientVectorZ = new dataType * [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType[dim2D]{0};
		gradientVectorY[k] = new dataType[dim2D]{0};
		gradientVectorZ[k] = new dataType[dim2D]{0};
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;

	compute3dImageGradient(ctImageData.imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, ctImageData.spacing);

	//use distance map in potential
	dataType** distance = new dataType * [height];
	dataType** edgeDetector = new dataType * [height];
	dataType** maskThreshold = new dataType * [height];
	for (k = 0; k < height; k++) {
		distance[k] = new dataType[dim2D]{0};
		edgeDetector[k] = new dataType[dim2D]{0};
		maskThreshold[k] = new dataType[dim2D]{0};
	}
	
	dataType norm_of_gradient = 0.0, edge_coef = 1000.0;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			
			ux = gradientVectorX[k][i];
			uy = gradientVectorY[k][i];
			uz = gradientVectorZ[k][i];
			norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz);
			edgeDetector[k][i] = gradientFunction(norm_of_gradient, edge_coef);
			
			//threshold
			if (edgeDetector[k][i] < 0.15) {
				maskThreshold[k][i] = 0.0;
			}
			else {
				maskThreshold[k][i] = 1.0;
			}

		}
	}
	
	fastSweepingFunction_3D(distance, maskThreshold, length, width, height, 1.0, 10000000.0, 0.0);

	//get real world coordinates of the seed point
	Point3D initial_point = getRealCoordFromImageCoord3D(seedPoint, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);

	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };
	seedStats = getStats(ctImageData, initial_point, radius);
	dataType seedValCT = seedStats.mean_data;

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = abs(seedValCT - ctImageData.imageDataPtr[k][i]);
		}
	}

	//Find the max of each potential for normalization
	dataType maxImage = 0.0, maxMean = 0.0;
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

bool shortestPath3D(dataType** distanceFuncPtr, const size_t length, const size_t width, const size_t height, VoxelSpacing spacing, Point3D* seedPoints, vector<Point3D>& path_points) {

	if (distanceFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i = 0, j = 0, k = 0, xd = 0, dim2D = length * width, max_iter = 1000;
	dataType tau = 0.8, tol = 1.0;

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
	compute3dImageGradient(distanceFuncPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, spacing);

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

	do {

		currentIndx = x_new(i, j, length);

		iNew = iNew - tau * gradientVectorX[k][currentIndx];
		jNew = jNew - tau * gradientVectorY[k][currentIndx];
		kNew = kNew - tau * gradientVectorZ[k][currentIndx];
		
		Point3D point_current = { iNew, jNew, kNew };
		path_points.push_back(point_current);
		
		//compute distance current Point - last point
		dist_to_end = getPoint3DDistance(point_current, seedPoints[0]);
		
		i = (size_t)(round(iNew));
		j = (size_t)(round(jNew));
		k = (size_t)(round(kNew));
		
		count_iter++;

	} while (dist_to_end > tol && count_iter < max_iter);

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

/*
bool findPathFromOneGivenPointWithCircleDetection(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters, size_t stop_criterium) {

	if (ctImageData.imageDataPtr == NULL || meanImagePtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t k = 0, i = 0, j = 0, x = 0;
	const size_t height = ctImageData.height, length = ctImageData.length, width = ctImageData.width;
	const size_t dim2D = length * width;
	double offset = 5.0;
	BoundingBox2D box = { 0, 0, 0, 0 };

	string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, extension, slice_number;

	saving_name = path_name + "path_points_test.csv";
	FILE* file;
	if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	std::cout << "Initial point : (" << seedPoints[0].x << ", " << seedPoints[0].y << " , " << seedPoints[0].z << ")" << std::endl;

	dataType** actionPtr = new dataType * [height];
	dataType** newActionPtr = new dataType * [height];
	dataType** maskDistance = new dataType * [height];
	dataType** potentialPtr = new dataType * [height];
	dataType** centeredPath = new dataType * [height];
	for (k = 0; k < height; k++) {
		actionPtr[k] = new dataType[dim2D]{ 0 };
		newActionPtr[k] = new dataType[dim2D]{ 0 };
		maskDistance[k] = new dataType[dim2D]{ 0 };
		potentialPtr[k] = new dataType[dim2D]{ 0 };
		centeredPath[k] = new dataType[dim2D]{ 0 };
	}
	if (actionPtr == NULL || newActionPtr == NULL || maskDistance == NULL || potentialPtr == NULL || centeredPath == NULL)
		return false;

	Point3D* seeds = new Point3D[2]; // dynamic array needed for path finding.
	seeds[0] = seedPoints[0];
	Point3D initial_point = seedPoints[0];
	Point3D temporary_point = { 0.0, 0.0, 0.0 };
	Point3D current_point = { 0.0, 0.0, 0.0 };
	Point3D real_world = { 0.0, 0.0, 0.0 };

	double step = 50.0, dist = 0.0;
	dataType min_distance = BIG_VALUE, value_temp = 0.0, h = 1.0;

	//================ find next point inside the aorta ============
	
	double radius_to_find_mean = 3.0;
	compute3DPotential(ctImageData, potentialPtr, seedPoints[0], radius_to_find_mean, parameters);
	
	//saving_name = "C:/Users/Konan Allaly/Documents/Tests/input/potential.raw";
	//load3dArrayRAW<dataType>(potentialPtr, length, width, height, saving_name.c_str(), false);
	//saving_name = path_name + "potential.raw";
	//store3dRawData<dataType>(potentialPtr, length, width, height, saving_name.c_str());
	//saving_name = "C:/Users/Konan Allaly/Documents/Tests/input/action_field.raw";
	//load3dArrayRAW<dataType>(actionPtr, length, width, height, saving_name.c_str(), false);
	
	fastMarching3D_N(actionPtr, potentialPtr, length, width, height, initial_point);
	
	//saving_name = path_name + "action_field_initial.raw";
	//store3dRawData<dataType>(actionPtr, length, width, height, saving_name.c_str());
	
	//get real world coordinate of the initial point
	real_world = getRealCoordFromImageCoord3D(initial_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);

	int count_step = 1;
	min_distance = BIG_VALUE;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				x = x_new(i, j, length);
				current_point.x = (dataType)i;
				current_point.y = (dataType)j;
				current_point.z = (dataType)k;
				current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
				dist = getPoint3DDistance(real_world, current_point);
				//dist = getPoint3DDistance(initial_point, current_point);
				if (dist >= (step - 0.01) && dist <= (step + 0.01)) {
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
	count_step++;
	std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;

	//vector for path points
	vector<Point3D> path_points;

	//path between initial point and second point
	seeds[1] = temporary_point;
	shortestPath3D(actionPtr, resultedPath, length, width, height, h, seeds, path_points);

	//write path points coordinates to file
	int i_n = 0, k_n = 0, k_center = 0;
	Point3D point_file = { 0.0, 0.0, 0.0 };
	for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
		point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	}

	exit(0);
	
	//saving_name = path_name + "path.raw";
	//store3dRawData<dataType>(resultedPath, length, width, height, saving_name.c_str());
	
	//=========== Hough transform parameters ==================
	
	//dataType* imageSlice = new dataType[dim2D]{ 0 };
	//dataType* houghSpace = new dataType[dim2D]{ 0 };
	//dataType* voteArray = new dataType[dim2D]{ 0 };

	////Slice extraction
	//Point2D seed2D = { 0.0, 0.0 };

	////filtering parameters
	//Filter_Parameters filter_parameters;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	//Image_Data2D IMAGE;
	//IMAGE.imageDataPtr = imageSlice;
	//IMAGE.height = length; IMAGE.width = width;

	size_t m = 0, n = 0, p = 0;

	//Point2D center2D = { 0.0, 0.0 };
	//PixelSpacing spacing = { ctImageData.spacing.sx, ctImageData.spacing.sy };
	//HoughParameters params = { 6.0, 19.0, 0.5, 5.0, 1000, 1.0, 0.2, 0.5, spacing };
	//dataType max_ratio = 1.0;

	//================ Use circle detection as path finding stopping criterium =====================
	
	size_t num_slice_processed = 0;
	//int len = path_points.size();
	//std::cout << len << " point to be processed" << std::endl;
	//
	//Point3D next_point = { 0.0, 0.0, 0.0 };
	//for (i_n = len - 1; i_n > -1; i_n--) {
	//	num_slice_processed++;
	//	extension = to_string(num_slice_processed);
	//	//slice number
	//	k_n = (size_t)path_points[i_n].z;
	//	slice_number = to_string(k_n);
	//	//Slice extraction
	//	for (i = 0; i < dim2D; i++) {
	//		imageSlice[i] = ctImageData.imageDataPtr[k_n][i];
	//	}
	//	//Slice rescaling
	//	rescaleNewRange2D(imageSlice, length, width, 0.0, 1.0);
	//	//Slice filtering
	//	heatImplicit2dScheme(IMAGE, filter_parameters);
	//	seed2D.x = path_points[i_n].x;
	//	seed2D.y = path_points[i_n].y;
	//	if (num_slice_processed < 10) {
	//		saving_name = path_name + "threshold_000" + extension + ".raw";
	//	}
	//	else {
	//		if (num_slice_processed < 100) {
	//			saving_name = path_name + "threshold_00" + extension + ".raw";
	//		}
	//		else {
	//			if (num_slice_processed < 1000) {
	//				saving_name = path_name + "threshold_0" + extension + ".raw";
	//			}
	//			else {
	//				saving_name = path_name + "threshold_" + extension + ".raw";
	//			}
	//		}
	//	}
	//	center2D = localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, length, width, params, saving_name);
	//	
	//	max_ratio = getTheMaxValue(voteArray, length, width);
	//	if (max_ratio != 0.0) {
	//		centeredPath[k_n][x_new((size_t)center2D.x, (size_t)center2D.y, length)] = 1.0;
	//		//k_center = k_n;
	//		next_point.x = center2D.x;
	//		next_point.y = center2D.y;
	//		next_point.z = k_n;
	//	}
	//	if (num_slice_processed < 10) {
	//		saving_name = path_name + "slice_000" + extension + ".raw";
	//		store2dRawData(imageSlice, length, width, saving_name.c_str());
	//		saving_name = path_name + "found_circle_000" + extension + ".raw";
	//		store2dRawData(houghSpace, length, width, saving_name.c_str());
	//		saving_name = path_name + "ratio_000" + extension + ".raw";
	//		store2dRawData(voteArray, length, width, saving_name.c_str());
	//	}
	//	else {
	//		if (num_slice_processed < 100) {
	//			saving_name = path_name + "slice_00" + extension + ".raw";
	//			store2dRawData(imageSlice, length, width, saving_name.c_str());
	//			saving_name = path_name + "found_circle_00" + extension + ".raw";
	//			store2dRawData(houghSpace, length, width, saving_name.c_str());
	//			saving_name = path_name + "ratio_00" + extension + ".raw";
	//			store2dRawData(voteArray, length, width, saving_name.c_str());
	//		}
	//		else {
	//			if (num_slice_processed < 1000) {
	//				saving_name = path_name + "slice_0" + extension + ".raw";
	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//				saving_name = path_name + "found_circle_0" + extension + ".raw";
	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//				saving_name = path_name + "ratio_0" + extension + ".raw";
	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//			}
	//			else {
	//				saving_name = path_name + "slice_" + extension + ".raw";
	//				store2dRawData(imageSlice, length, width, saving_name.c_str());
	//				saving_name = path_name + "found_circle_" + extension + ".raw";
	//				store2dRawData(houghSpace, length, width, saving_name.c_str());
	//				saving_name = path_name + "ratio_" + extension + ".raw";
	//				store2dRawData(voteArray, length, width, saving_name.c_str());
	//			}
	//		}
	//	}
	//}

	//empty the list of point before new piece of path
	while (path_points.size() != 0) {
		path_points.pop_back();
	}

	//============================================
	size_t count_stop = 0, max_steps = 10;
	
	//max_ratio = 1.0;

	size_t stop_k = (size_t)seedPoints[1].z, verif_k = (size_t)seeds[1].z;

	Point3D last_point = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	Point3D verif_point = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	double verif_dist = getPoint3DDistance(last_point, verif_point);

	std::cout << "distance to end point : " << verif_dist << std::endl;

	while ((stop_k != verif_k && verif_dist < 12.5) || count_step < max_steps) {
		//max_ratio != 0.0 && count_step < max_steps //--> old condition
		std::cout << "STEP " << count_step << " : " << std::endl;
		extension = to_string(count_step);
		
		//mask
		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {
					x = x_new(i, j, length);
					current_point.x = (dataType)i;
					current_point.y = (dataType)j;
					current_point.z = (dataType)k;
					current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
					dist = getPoint3DDistance(real_world, current_point);
					//dist = getPoint3DDistance(initial_point, current_point);
					if (dist <= (step + 0.01) && actionPtr[k][x] <= value_temp) {
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
		
		fastMarching3D_N(newActionPtr, potentialPtr, length, width, height, initial_point);

		//saving_name = path_name + "action_field" + extension + ".raw";
		//store3dRawData<dataType>(newActionPtr, length, width, height, saving_name.c_str());

		//Find the temporary point
		min_distance = BIG_VALUE;
		
		//get real world coordinate of the initial point
		real_world = getRealCoordFromImageCoord3D(initial_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		
		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {
					x = x_new(i, j, length);
					current_point.x = (dataType)i;
					current_point.y = (dataType)j;
					current_point.z = (dataType)k;
					current_point = getRealCoordFromImageCoord3D(current_point, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
					dist = getPoint3DDistance(real_world, current_point);
					//dist = getPoint3DDistance(initial_point, current_point);
					if (dist >= (step - 0.01) && dist <= (step + 0.01) && maskDistance[k][x] != 10000.0) {
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
		//path between points
		std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;
		
		seeds[1] = temporary_point;
		shortestPath3D(newActionPtr, resultedPath, length, width, height, h, seeds, path_points);

		//write path points coordinates to file
		for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
			point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
			fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
		}
		
		saving_name = path_name + "path.raw";
		store3dRawData<dataType>(resultedPath, length, width, height, saving_name.c_str());
		
		copyDataToAnotherArray(newActionPtr, actionPtr, height, length, width);
		count_step++;

		////find circles
		//len = path_points.size();
		//std::cout << len << " point to be processed" << std::endl;
		//
		//count_stop = 0;
		//for (i_n = len - 1; i_n > -1; i_n--) {
		//	num_slice_processed++;
		//	extension = to_string(num_slice_processed);
		//	//slice number
		//	k_n = (size_t)path_points[i_n].z;
		//	slice_number = to_string(k_n);
		//	seed2D.x = path_points[i_n].x;
		//	seed2D.y = path_points[i_n].y;
		//	//Slice extraction
		//	for (i = 0; i < dim2D; i++) {
		//		imageSlice[i] = ctImageData.imageDataPtr[k_n][i];
		//	}
		//	//Slice rescaling
		//	rescaleNewRange2D(imageSlice, length, width, 0.0, 1.0);
		//	//Slice filtering
		//	heatImplicit2dScheme(IMAGE, filter_parameters);
		//	if (num_slice_processed < 10) {
		//		saving_name = path_name + "threshold_000" + extension + ".raw";
		//	}
		//	else {
		//		if (num_slice_processed < 100) {
		//			saving_name = path_name + "threshold_00" + extension + ".raw";
		//		}
		//		else {
		//			if (num_slice_processed < 1000) {
		//				saving_name = path_name + "threshold_0" + extension + ".raw";
		//			}
		//			else {
		//				saving_name = path_name + "threshold_" + extension + ".raw";
		//			}
		//		}
		//	}
		//	center2D = localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, length, width, params, saving_name);
		//	max_ratio = getTheMaxValue(voteArray, length, width);
		//	if (max_ratio != 0.0) {
		//		centeredPath[k_n][x_new((size_t)center2D.x, (size_t)center2D.y, length)] = 1.0;
		//		//k_center = k_n;
		//		next_point.x = center2D.x;
		//		next_point.y = center2D.y;
		//		next_point.z = k_n;
		//	}
		//	if (num_slice_processed < 10) {
		//		saving_name = path_name + "slice_000" + extension + ".raw";
		//		store2dRawData(imageSlice, length, width, saving_name.c_str());
		//		saving_name = path_name + "found_circle_000" + extension + ".raw";
		//		store2dRawData(houghSpace, length, width, saving_name.c_str());
		//		saving_name = path_name + "ratio_000" + extension + ".raw";
		//		store2dRawData(voteArray, length, width, saving_name.c_str());
		//	}
		//	else {
		//		if (num_slice_processed < 100) {
		//			saving_name = path_name + "slice_00" + extension + ".raw";
		//			store2dRawData(imageSlice, length, width, saving_name.c_str());
		//			saving_name = path_name + "found_circle_00" + extension + ".raw";
		//			store2dRawData(houghSpace, length, width, saving_name.c_str());
		//			saving_name = path_name + "ratio_00" + extension + ".raw";
		//			store2dRawData(voteArray, length, width, saving_name.c_str());
		//		}
		//		else {
		//			if (num_slice_processed < 1000) {
		//				saving_name = path_name + "slice_0" + extension + ".raw";
		//				store2dRawData(imageSlice, length, width, saving_name.c_str());
		//				saving_name = path_name + "found_circle_0" + extension + ".raw";
		//				store2dRawData(houghSpace, length, width, saving_name.c_str());
		//				saving_name = path_name + "ratio_0" + extension + ".raw";
		//				store2dRawData(voteArray, length, width, saving_name.c_str());
		//			}
		//			else {
		//				saving_name = path_name + "slice_" + extension + ".raw";
		//				store2dRawData(imageSlice, length, width, saving_name.c_str());
		//				saving_name = path_name + "found_circle_" + extension + ".raw";
		//				store2dRawData(houghSpace, length, width, saving_name.c_str());
		//				saving_name = path_name + "ratio_" + extension + ".raw";
		//				store2dRawData(voteArray, length, width, saving_name.c_str());
		//			}
		//		}
		//	}
		//	if (max_ratio == 0.0) {
		//		count_stop++;
		//		//break;
		//	}
		//	else {
		//		count_stop = 0;
		//	}
		//	if (max_ratio == 0.0 && count_stop == stop_criterium) {
		//		break;
		//	}
		//}

		//empty the list of point
		while (path_points.size() != 0) {
			path_points.pop_back();
		}

		verif_point = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		verif_dist = getPoint3DDistance(last_point, verif_point);

		std::cout << "distance to end point : " << verif_dist << std::endl;

	}

	std::cout << num_slice_processed << " slices were processed\n" << std::endl;
	fclose(file);

	//delete[] imageSlice;
	//delete[] houghSpace;
	//delete[] voteArray;

	for (k = 0; k < height; k++) {
		delete[] actionPtr[k];
		delete[] newActionPtr[k];
		delete[] maskDistance[k];
		delete[] potentialPtr[k];
		delete[] centeredPath[k];
	}
	delete[] actionPtr;
	delete[] newActionPtr;
	delete[] maskDistance;
	delete[] potentialPtr;
	delete[] centeredPath;

	delete[] seeds;

	return true;
}
*/

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

	double radius = 3.0;

	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, saving_csv = outputPath + "path_p2.csv";
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

	//======== First step ===========//
	compute3DPotential(ctImageData, potential, seedsPath[0], radius, parameters);

	partialFrontPropagation(action_field, potential, length, width, height, seedsPath);
	//fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);
	
	//saving_name = outputPath + "action3_1.raw";
	//store3dRawData<dataType>(action_field, length, width, height, saving_name.c_str());

	shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	
	fprintf(file, "x,y,z\n");
	for (i = 0; i < path_points.size(); i++) {
		path_points[i] = getRealCoordFromImageCoord3D(path_points[i], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[i].x, path_points[i].y, path_points[i].z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}
	
	//================================//

	Point3D seed1 = getRealCoordFromImageCoord3D(seedPoints[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	Point3D seed2 = getRealCoordFromImageCoord3D(seedPoints[2], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	double step = 0.75 * getPoint3DDistance(seed1, seed2), dist = 0.0;
	cout << "Step : " << step << endl;

	//======== Second step ===========//
	fastMarching3D_N(action_field, potential, length, width, height, seedsPath[1]);
	
	//saving_name = outputPath + "action3_2.raw";
	//store3dRawData<dataType>(action_field, length, width, height, saving_name.c_str());

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
	shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);

	//copy points to file
	for (i = 0; i < path_points.size(); i++) {
		path_points[i] = getRealCoordFromImageCoord3D(path_points[i], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[i].x, path_points[i].y, path_points[i].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	//======= Last step ===============//
	seedsPath[0] = temporary_point;
	seedsPath[1] = seedPoints[2];
	
	fastMarching3D_N(action_field, potential, length, width, height, seedsPath[0]);
	
	//saving_name = outputPath + "action3_3.raw";
	//store3dRawData<dataType>(action_field, length, width, height, saving_name.c_str());
	
	shortestPath3D(action_field, length, width, height, ctImageData.spacing, seedsPath, path_points);
	
	for (i = 0; i < path_points.size(); i++) {
		path_points[i] = getRealCoordFromImageCoord3D(path_points[i], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", path_points[i].x, path_points[i].y, path_points[i].z);
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

/*
bool findPathFromOneGivenPoint(Image_Data ctImageData, Point3D * seedPoints, Potential_Parameters parameters) {

	if (ctImageData.imageDataPtr == NULL)
		return false;

	size_t k = 0, i = 0, j = 0, x = 0;
	const size_t height = ctImageData.height, length = ctImageData.length, width = ctImageData.width;
	const size_t dim2D = length * width;

	double epsilon = 0.01;

	string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, extension;

	saving_name = path_name + "path_points.csv";
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
	double step = 0.5 * distance_point;

	//================ find next point inside the aorta 

	// saving_name = path_name + "potential.raw";
	//store3dRawData<dataType>(potentialPtr, length, width, height, saving_name.c_str());

	double radius_to_find_mean = 3.0;
	compute3DPotential(ctImageData, potentialPtr, seedPoints[0], radius_to_find_mean, parameters);

	fastMarching3D_N(actionPtr, potentialPtr, length, width, height, seeds[0]);

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
	distance_point = getPoint3DDistance(seed1, seed2);
	std::cout << "found point : (" << temporary_point.x << ", " << temporary_point.y << ", " << temporary_point.z << ")" << std::endl;

	//vector for path points
	vector<Point3D> path_points;

	//path between initial point and second point
	seeds[1] = temporary_point;
	shortestPath3D(actionPtr, length, width, height, ctImageData.spacing, seeds, path_points);

	//write just path points coordinates to file
	int i_n = 0, k_n = 0, k_center = 0;
	Point3D point_file = { 0.0, 0.0, 0.0 };
	
	for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
		point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
	}

	//============================================

	double distance_to_stop = 0.1 * distance_point;

	size_t count_step = 0;

	while (distance_to_stop < distance_point) {
		
		count_step++;
		std::cout << "count step = " << count_step << std::endl;

		//count_step++;
		//std::cout << "STEP " << count_step << " : " << std::endl;
		//extension = to_string(count_step);

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

		fastMarching3D_N(newActionPtr, potentialPtr, length, width, height, seeds[0]);

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
		
		//slice_stop = (size_t)temporary_point.z;
		//cout << "current slice stop : " << slice_stop << endl;

		seeds[1] = temporary_point;
		//path between points
		//seed1 = getRealCoordFromImageCoord3D(seeds[0], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		seed2 = getRealCoordFromImageCoord3D(seeds[1], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
		distance_point = getPoint3DDistance(seed1, seed2);
		shortestPath3D(newActionPtr, length, width, height, ctImageData.spacing, seeds, path_points);

		copyDataToAnotherArray(newActionPtr, actionPtr, height, length, width);

		//write just path points coordinates to file
		for (i_n = path_points.size() - 1; i_n > -1; i_n--) {
			point_file = getRealCoordFromImageCoord3D(path_points[i_n], ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
			fprintf(file, "%f,%f,%f\n", point_file.x, point_file.y, point_file.z);
		}
		
		if (count_step == 2) {
			exit(0);
		}
		
		//saving_name = path_name + "path.raw";
		//store3dRawData<dataType>(resultedPath, length, width, height, saving_name.c_str());
		
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
*/

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

bool partialPropagation(dataType** actionMapPtr, dataType** potentialPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints, const double maxLength) {

	if (actionMapPtr == NULL || potentialPtr == NULL) {
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
			actionMapPtr[k][i] = INFINITY;
			labelArray[k][i] = 3;
		}
	}

	//Processed the initial point
	i = (size_t)seedPoints[0].x;
	j = (size_t)seedPoints[0].y;
	k = (size_t)seedPoints[0].z;
	currentIndx = x_new(i, j, length);
	actionMapPtr[k][currentIndx] = 0.0;
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
		x = select3dX(actionMapPtr, length, width, height, i, j, kminus);
		y = select3dY(actionMapPtr, length, width, height, i, j, kminus);
		z = select3dZ(actionMapPtr, length, width, height, i, j, kminus);
		coefSpeed = potentialPtr[kminus][currentIndx];
		dTop = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D TopNeighbor = { i, j, kminus, dTop };
		actionMapPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(iplus, j, length);
		x = select3dX(actionMapPtr, length, width, height, iplus, j, k);
		y = select3dY(actionMapPtr, length, width, height, iplus, j, k);
		z = select3dZ(actionMapPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialPtr[k][indxEast];
		dEast = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D EastNeighbor = { iplus, j, k, dEast };
		actionMapPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(i, jminus, length);
		x = select3dX(actionMapPtr, length, width, height, i, jminus, k);
		y = select3dY(actionMapPtr, length, width, height, i, jminus, k);
		z = select3dZ(actionMapPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialPtr[k][indxNorth];
		dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D NorthNeighbor = { i, jminus, k, dNorth };
		actionMapPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(iminus, j, length);
		x = select3dX(actionMapPtr, length, width, height, iminus, j, k);
		y = select3dY(actionMapPtr, length, width, height, iminus, j, k);
		z = select3dZ(actionMapPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialPtr[k][indxWest];
		dWest = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D WestNeighbor = { iminus, j, k, dWest };
		actionMapPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(i, jplus, length);
		x = select3dX(actionMapPtr, length, width, height, i, jplus, k);
		y = select3dY(actionMapPtr, length, width, height, i, jplus, k);
		z = select3dZ(actionMapPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialPtr[k][indxSouth];
		dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D SouthNeighbor = { i, jplus, k, dSouth };
		actionMapPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(actionMapPtr, length, width, height, i, j, kplus);
		y = select3dY(actionMapPtr, length, width, height, i, j, kplus);
		z = select3dZ(actionMapPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialPtr[kplus][currentIndx];
		dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D BottomNeighbor = { i, j, kplus, dBottom };
		actionMapPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	heapifyVector3D(inProcess);
	size_t label = 0;

	double length_propagation = 0.0;

	while (length_propagation < maxLength) {

		//processed the point with minimum distance
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		i = current.x;
		j = current.y;
		k = current.z;
		
		seedPoints[1].x = current.x;
		seedPoints[1].y = current.y;
		seedPoints[1].z = current.z;
		length_propagation = getPoint3DDistance(seedPoints[0], seedPoints[1]);
		
		currentIndx = x_new(i, j, length);
		labelArray[k][currentIndx] = 1;
		actionMapPtr[k][currentIndx] = current.arrival;

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
				coefSpeed = potentialPtr[kminus][currentIndx];
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

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(iplus, j, length);
			label = labelArray[k][indxEast];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, iplus, j, k);
				y = select3dY(actionMapPtr, length, width, height, iplus, j, k);
				z = select3dZ(actionMapPtr, length, width, height, iplus, j, k);
				coefSpeed = potentialPtr[k][indxEast];
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

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(i, jminus, length);
			label = labelArray[k][indxNorth];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, jminus, k);
				y = select3dY(actionMapPtr, length, width, height, i, jminus, k);
				z = select3dZ(actionMapPtr, length, width, height, i, jminus, k);
				coefSpeed = potentialPtr[k][indxNorth];
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

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(iminus, j, length);
			label = labelArray[k][indxWest];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, iminus, j, k);
				y = select3dY(actionMapPtr, length, width, height, iminus, j, k);
				z = select3dZ(actionMapPtr, length, width, height, iminus, j, k);
				coefSpeed = potentialPtr[k][indxWest];
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

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(i, jplus, length);
			label = labelArray[k][indxSouth];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, jplus, k);
				y = select3dY(actionMapPtr, length, width, height, i, jplus, k);
				z = select3dZ(actionMapPtr, length, width, height, i, jplus, k);
				coefSpeed = potentialPtr[k][indxSouth];
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

		//Bottom
		if (k < height_minus && i >= 0 && i < length && j >= 0 && j < width) {
			kplus = k + 1;
			label = labelArray[kplus][currentIndx];
			if (label != 1) {
				x = select3dX(actionMapPtr, length, width, height, i, j, kplus);
				y = select3dY(actionMapPtr, length, width, height, i, j, kplus);
				z = select3dZ(actionMapPtr, length, width, height, i, j, kplus);
				coefSpeed = potentialPtr[kplus][currentIndx];
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

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (actionMapPtr[k][i] == INFINITY) {
				actionMapPtr[k][i] = 0.0;
			}
		}
	}

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////
//=================== Functions for the 4D fast marching ==============================
///////////////////////////////////////////////////////////////////////////////////////

double getDistancePoints4D(Point4D p1, Point4D p2) {
	return (dataType)(sqrt(pow(p1.r - p2.r, 2) + pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2)));
}

dataType select4DX(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rlength, const size_t i, const size_t j, const size_t k, const size_t l) {

	dataType i_minus, i_plus;
	size_t zr = x_new(k, l, height);

	if (i == 0) {
		i_minus = INFINITY;
	}
	else {
		if (i > 0) {
			i_minus = actionMapPtr[zr][x_new(i - 1, j, length)];
		}
		//the case i < 0 is impossible 
	}

	if (i == length - 1) {
		i_plus = INFINITY;
	}
	else {
		if (i < length - 1) {
			i_plus = actionMapPtr[zr][x_new(i + 1, j, length)];
		}
		//the case i > length - 1 is impossible 
	}

	return min(i_minus, i_plus);
}

dataType select4DY(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rlength, const size_t i, const size_t j, const size_t k, const size_t l) {

	dataType j_minus, j_plus;
	size_t zr = x_new(k, l, height);

	if (j == 0) {
		j_minus = INFINITY;
	}
	else {
		if (j > 0) {
			j_minus = actionMapPtr[zr][x_new(i, j - 1, length)];
		}
		//the case j < 0 is impossible 
	}

	if (j == width - 1) {
		j_plus = INFINITY;
	}
	else {
		if (j < width - 1) {
			j_plus = actionMapPtr[zr][x_new(i, j + 1, length)];
		}
		//the case j > width - 1 is impossible 
	}

	return min(j_minus, j_plus);
}

dataType select4DZ(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rlength, const size_t i, const size_t j, const size_t k, const size_t l) {

	dataType k_minus, k_plus;

	size_t xy = x_new(i, j, length);

	if (k == 0) {
		k_minus = INFINITY;
	}
	else {
		if (k > 0) {
			k_minus = actionMapPtr[x_new(k - 1, l, height)][xy];
		}
		//the case k < 0 is impossible 
	}

	if (k == height - 1) {
		k_plus = INFINITY;
	}
	else {
		if (k < height - 1) {
			k_plus = actionMapPtr[x_new(k + 1, l, height)][xy];
		}
		//the case k > height - 1 is impossible 
	}

	return min(k_minus, k_plus);
}

dataType select4DR(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rlength, const size_t i, const size_t j, const size_t k, const size_t l) {

	dataType r_minus, r_plus;
	size_t xy = x_new(i, j, length);

	if (l == 0) {
		r_minus = INFINITY;
	}
	else {
		if (l > 0) {
			r_minus = actionMapPtr[x_new(k, l - 1, height)][xy];
		}
		//the case r < 0 is impossible 
	}

	if (l == rlength - 1) {
		r_plus = INFINITY;
	}
	else {
		if (l < rlength - 1) {
			r_plus = actionMapPtr[x_new(k, l + 1, height)][xy];
		}
		//the case l > rlength - 1 is impossible 
	}

	return min(r_minus, r_plus);
}

dataType solve4DQuadratic(dataType X, dataType Y, dataType Z, dataType R, dataType W) {

	dataType sol1 = 0.0, sol2 = 0.0, solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;

	if (X != INFINITY && Y != INFINITY && Z != INFINITY && R != INFINITY) {
		a = 4;
		b = (dataType)(-2 * (X + Y + Z + R));
		c = (dataType)(pow(X, 2) + pow(Y, 2) + pow(Z, 2) + pow(R,2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(min(min(X, Y),Z),R) + W);
		}
		//cout << "case 1 happened" << endl;
	}


	if (X != INFINITY && Y != INFINITY && Z != INFINITY && R == INFINITY) {
		a = 3;
		b = (dataType)(-2 * (X + Y + Z));
		c = (dataType)(pow(X, 2) + pow(Y, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(min(X, Y),Z) + W);
		}
		//cout << "case 2 happened" << endl;
	}

	if (X != INFINITY && Y != INFINITY && Z == INFINITY && R != INFINITY) {
		a = 3;
		b = (dataType)(-2 * (X + Y + R));
		c = (dataType)(pow(X, 2) + pow(Y, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(min(X, Y), R) + W);
		}
		//cout << "case 3 happened" << endl;
	}

	if (X != INFINITY && Y == INFINITY && Z != INFINITY && R != INFINITY) {
		a = 3;
		b = (dataType)(-2 * (X + Z + R));
		c = (dataType)(pow(X, 2) + pow(Z, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(min(X, Z), R) + W);
		}
		//cout << "case 4 happened" << endl;
	}

	if (X == INFINITY && Y != INFINITY && Z != INFINITY && R != INFINITY) {
		a = 3;
		b = (dataType)(-2 * (Y + Z + R));
		c = (dataType)(pow(Y, 2) + pow(Z, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(min(Y, Z), R) + W);
		}
		//cout << "case 5 happened" << endl;
	}


	if (X != INFINITY && Y != INFINITY && Z == INFINITY && R == INFINITY) {
		a = 2;
		b = (dataType)(-2 * (X + Y));
		c = (dataType)(pow(X, 2) + pow(Y, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(X, Y) + W);
		}
		//cout << "case 6 happened" << endl;
	}

	if (X != INFINITY && Y == INFINITY && Z != INFINITY && R == INFINITY) {
		a = 2;
		b = (dataType)(-2 * (X + Z));
		c = (dataType)(pow(X, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(X, Z) + W);
		}
		//cout << "case 7 happened" << endl;
	}

	if (X != INFINITY && Y == INFINITY && Z == INFINITY && R != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (X + R));
		c = (dataType)(pow(X, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(X, R) + W);
		}
		//cout << "case 8 happened" << endl;
	}

	if (X == INFINITY && Y != INFINITY && Z != INFINITY && R == INFINITY) {
		a = 2;
		b = (dataType)(-2 * (Y + Z));
		c = (dataType)(pow(Y, 2) + pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(Y, Z) + W);
		}
		//cout << "case 9 happened" << endl;
	}

	if (X == INFINITY && Y != INFINITY && Z == INFINITY && R != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (Y + R));
		c = (dataType)(pow(Y, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(Y, R) + W);
		}
		//cout << "case 10 happened" << endl;
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY && R != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (Z + R));
		c = (dataType)(pow(Z, 2) + pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(min(Z, R) + W);
		}
		//cout << "case 11 happened" << endl;
	}


	if (X != INFINITY && Y == INFINITY && Z == INFINITY && R == INFINITY) {
		a = 1;
		b = (dataType)(-2 * X);
		c = (dataType)(pow(X, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(X + W);
		}
		//cout << "case 12 happened" << endl;
	}

	if (X == INFINITY && Y != INFINITY && Z == INFINITY && R == INFINITY) {
		a = 1;
		b = (dataType)(-2 * Y);
		c = (dataType)(pow(Y, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(Y + W);
		}
		//cout << "case 13 happened" << endl;
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY && R == INFINITY) {
		a = 1;
		b = (dataType)(-2 * Z);
		c = (dataType)(pow(Z, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(Z + W);
		}
		//cout << "case 14 happened" << endl;
	}

	if (X == INFINITY && Y == INFINITY && Z == INFINITY && R != INFINITY) {
		a = 1;
		b = (dataType)(-2 * R);
		c = (dataType)(pow(R, 2) - W);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			sol1 = (dataType)((-b + sqrt(delta)) / (2 * a));
			sol2 = (dataType)((-b - sqrt(delta)) / (2 * a));
			solution = max(sol1, sol2);
		}
		else {
			solution = (dataType)(R + W);
		}
		//cout << "case 15 happened" << endl;
	}

	return solution;
}

void swap4DPoints(pointFastMarching4D* a, pointFastMarching4D* b) {
	pointFastMarching4D temp = *a;
	*a = *b;
	*b = temp;
}

void heapifyDown4D(vector<pointFastMarching4D>& in_Process, int i) {

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
		swap4DPoints(&in_Process[i], &in_Process[current]);
		heapifyDown4D(in_Process, current);
	}

}

void heapifyUp4D(vector<pointFastMarching4D>& in_Process, int i) {

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
		swap4DPoints(&in_Process[current], &in_Process[i]);
		heapifyUp4D(in_Process, current);
	}

}

void heapifyVector4D(vector<pointFastMarching4D>& in_Process) {
	int length_array = in_Process.size();
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		heapifyDown4D(in_Process, ind);
	}
}

void deleteRootHeap4D(vector<pointFastMarching4D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap4DPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown4D(in_Process, 0);
}

void addPointHeap4D(vector<pointFastMarching4D>& in_Process, pointFastMarching4D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp4D(in_Process, l - 1);
}

void computePotential4D(Image_Data ctImageData, dataType** potentialPtr, Point3D seedPoint, const size_t radiusLength,  dataType radiusInitial, dataType radiusScale) {
	
	size_t height = ctImageData.height;
	size_t length = ctImageData.length;
	size_t width = ctImageData.width;
	size_t i, j, k, l, ij, kl;
	dataType lamba1 = 1.0, lambda2 = 1.0;

	dataType radius = 0.0;
	Statistics stats = getPointNeighborhoodStats(ctImageData, seedPoint, radiusInitial);
	dataType coefInitial = stats.mean_data / radiusInitial;
	dataType varianceInitial = sqrt(stats.sd_data) / radiusInitial;
	dataType mean = 0.0, variance = 0.0;

	for (l = 0; l < radiusLength; l++) {
		radius = (l + 1) * radiusScale;
		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {
					kl = x_new(k, l, height);
					ij = x_new(i, j, length);
					Point3D current_point = { i, j, k };
					stats = getPointNeighborhoodStats(ctImageData, current_point, radius);
					mean = stats.mean_data / radius;
					variance = sqrt(stats.sd_data) / radius;
					potentialPtr[kl][ij] = lamba1 * pow(abs(mean - coefInitial), 2) + lambda2 * pow(abs(variance - varianceInitial), 2);
				}
			}
		}
	}

}

bool fastMarching4D(dataType** actionMapPtr, dataType** potentialPtr, const size_t length, const size_t width, const size_t height, const size_t radiusLength, Point3D seedPoint, dataType raduisScale, dataType radiusInitial) {

	if (actionMapPtr == NULL || potentialPtr == NULL) {
		return false;
	}

	size_t kl = 0, ij = 0;
	const size_t newHeight = radiusLength * height;
	size_t dim2D = width * length;
	size_t** labelArray = new size_t * [newHeight];
	for (kl = 0; kl < newHeight; kl++) {
		labelArray[kl] = new size_t[dim2D]{0};
	}
	if (labelArray == NULL)
		return false;

	//size_t current_kl = 0, current_ij = 0;

	//Initialization
	//All the points are notProcessed ---> label = 3
	for (kl = 0; kl < newHeight; kl++) {
		for (ij = 0; ij < dim2D; ij++) {
			actionMapPtr[kl][ij] = INFINITY;
			labelArray[kl][ij] = 3;
		}
	}

	//Processed the initial point
	pointFastMarching4D current;
	size_t i = (size_t)seedPoint.x;
	size_t j = (size_t)seedPoint.y;
	size_t k = (size_t)seedPoint.z;
	size_t l = (size_t)(radiusInitial / raduisScale);
	//cout << "l = " << l << endl;
	actionMapPtr[x_new(k, l, height)][x_new(i, j, length)] = 0.0;
	labelArray[x_new(k, l, height)][x_new(i, j, length)] = 1;

	//find the neighbors of the initial point add them to inProcess
	vector <pointFastMarching4D> inProcess;
	 
	//Top : k - 1
	if (i >= 0 && i < length && j >= 0 && j < width && k > 0 && l >= 0 && l < radiusLength) {
		size_t kminus = k - 1;
		size_t current_ij = x_new(i, j, length);
		size_t current_kl = x_new(kminus, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, j, kminus, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//Bottom : k + 1
	if (i >= 0 && i < length && j >= 0 && j < width && k < (height - 1) && l >= 0 && l < radiusLength) {
		size_t kplus = k + 1;
		size_t current_ij = x_new(i, j, length);
		size_t current_kl = x_new(kplus, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, j, kplus, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//East : i + 1
	if (i < (length - 1) && j >= 0 && j < width && k >= 0 && k < height && l >= 0 && l < radiusLength) {
		size_t iplus = i + 1;
		size_t current_ij = x_new(iplus, j, length);
		size_t current_kl = x_new(k, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { iplus, j, k, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//West : i - 1
	if (i > 0 && j >= 0 && j < width && k >= 0 && k < height && l >= 0 && l < radiusLength) {
		size_t iminus = i - 1;
		size_t current_ij = x_new(iminus, j, length);
		size_t current_kl = x_new(k, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { iminus, j, k, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//North : j - 1
	if (i >= 0 && i < length && j > 0 && k >= 0 && k < height && l >= 0 && l < radiusLength) {
		size_t jminus = j - 1;
		size_t current_ij = x_new(i, jminus, length);
		size_t current_kl = x_new(k, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, jminus, k, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//South : j + 1
	if (i >= 0 && i < length && j < (width - 1) && k >= 0 && k < height && l >= 0 && l < radiusLength) {
		size_t jplus = j + 1;
		size_t current_ij = x_new(i, jplus, length);
		size_t current_kl = x_new(k, l, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, jplus, k, l, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//max radius : l + 1
	if (i >= 0 && i < length && j >= 0 && j < width && k >= 0 && k < height && l < (radiusLength - 1)) {
		size_t lplus = l + 1;
		size_t current_ij = x_new(i, j, length);
		size_t current_kl = x_new(k, lplus, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, j, k, lplus, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	//min radius : l - 1
	if (i >= 0 && i < length && j >= 0 && j < width && k >= 0 && k < height && l > 0) {
		size_t lminus = l - 1;
		size_t current_ij = x_new(i, j, length);
		size_t current_kl = x_new(k, lminus, height);
		dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
		dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
		dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
		dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
		dataType coefSpeed = potentialPtr[current_kl][current_ij];
		dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
		pointFastMarching4D neighbor = { i, j, k, lminus, action };
		actionMapPtr[current_kl][current_ij] = action;
		inProcess.push_back(neighbor);
		labelArray[current_kl][current_ij] = 2;
	}

	heapifyVector4D(inProcess);
	size_t label = 0;

	while (inProcess.size() != 0) {

		//processed the point with minimum action
		current = inProcess[0];
		i = current.x;
		j = current.y;
		k = current.z;
		l = current.r;
		labelArray[x_new(k, l, height)][x_new(i, j, length)] = 1;
		actionMapPtr[x_new(k, l, height)][x_new(i, j, length)] = current.arrival;
		deleteRootHeap4D(inProcess);

		//processed neighbors of the minimum in the narrow band

		//Top : k - 1
		if (i >= 0 && i < length && j >= 0 && j < width && k > 0 && l >= 0 && l < radiusLength) {
			size_t kminus = k - 1;
			size_t current_ij = x_new(i, j, length);
			size_t current_kl = x_new(kminus, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, kminus, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { i, j, kminus, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//Bottom : k + 1
		if (i >= 0 && i < length && j >= 0 && j < width && k < (height - 1) && l >= 0 && l < radiusLength) {
			size_t kplus = k + 1;
			size_t current_ij = x_new(i, j, length);
			size_t current_kl = x_new(kplus, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, kplus, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { i, j, kplus, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//East : i + 1
		if (i < (length - 1) && j >= 0 && j < width && k >= 0 && k < height && l >= 0 && l < radiusLength) {
			size_t iplus = i + 1;
			size_t current_ij = x_new(iplus, j, length);
			size_t current_kl = x_new(k, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, iplus, j, k, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { iplus, j, k, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//West : i - 1
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height && l >= 0 && l < radiusLength) {
			size_t iminus = i - 1;
			size_t current_ij = x_new(iminus, j, length);
			size_t current_kl = x_new(k, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, iminus, j, k, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { iminus, j, k, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//North : j - 1
		if (i >= 0 && i < length && j > 0 && k >= 0 && k < height && l >= 0 && l < radiusLength) {
			size_t jminus = j - 1;
			size_t current_ij = x_new(i, jminus, length);
			size_t current_kl = x_new(k, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, jminus, k, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { i, jminus, k, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//South : j + 1
		if (i >= 0 && i < length && j < (width - 1) && k >= 0 && k < height && l >= 0 && l < radiusLength) {
			size_t jplus = j + 1;
			size_t current_ij = x_new(i, jplus, length);
			size_t current_kl = x_new(k, l, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, jplus, k, l);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { i, jplus, k, l, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//Radius min : l - 1
		if (i >= 0 && i < length && j >= 0 && j < width && k >= 0 && k < height && l > 0) {
			size_t lminus = l - 1;
			size_t current_ij = x_new(i, j, length);
			size_t current_kl = x_new(k, lminus, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, k, lminus);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D minNeighbor = { i, j, k, lminus, action };
					addPointHeap4D(inProcess, minNeighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

		//Radius max : l + 1
		if (i >= 0 && i < length && j >= 0 && j < width && k >= 0 && k < height && l < (radiusLength - 1)) {
			size_t lplus = l + 1;
			size_t current_ij = x_new(i, j, length);
			size_t current_kl = x_new(k, lplus, height);
			label = labelArray[current_kl][current_ij];
			if (label != 1) {
				dataType x = select4DX(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
				dataType y = select4DY(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
				dataType z = select4DZ(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
				dataType r = select4DR(actionMapPtr, length, width, height, radiusLength, i, j, k, lplus);
				dataType coefSpeed = potentialPtr[current_kl][current_ij];
				dataType action = solve4DQuadratic(x, y, z, r, coefSpeed);
				if (label == 3) {
					actionMapPtr[current_kl][current_ij] = action;
					labelArray[current_kl][current_ij] = 2;
					pointFastMarching4D neighbor = { i, j, k, lplus, action };
					addPointHeap4D(inProcess, neighbor);
				}
				else {
					if (label == 2) {
						if (action < actionMapPtr[current_kl][current_ij]) {
							actionMapPtr[current_kl][current_ij] = action;
						}
					}
				}
			}
		}

	}

	//====================
	//Saving
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, extension;
	dataType** actionField = new dataType * [height];
	for (k = 0; k < height; k++) {
		actionField[k] = new dataType[dim2D]{ 0 };
	}
	for (l = 0; l < radiusLength; l++) {
		extension = to_string(l);
		for (k = 0; k < height; k++) {
			for (i = 0; i < dim2D; i++) {
				actionField[k][i] = actionMapPtr[x_new(k, l, height)][i];
			}
		}
		saving_name = outputPath + "action4D_" + extension + ".raw";
		store3dRawData<dataType>(actionField, length, width, height, saving_name.c_str());
	}
	for (k = 0; k < height; k++) {
		delete[] actionField[k];
	}
	delete[] actionField;
	saving_name = outputPath + "action4D.raw";
	store3dRawData<dataType>(actionMapPtr, length, width, newHeight, saving_name.c_str());
	//=====================

	for (kl = 0; kl < newHeight; kl++) {
		delete[] labelArray[kl];
	}
	delete[] labelArray;

	return true;
}

bool compute4DimageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, dataType** gradientVectorR, const size_t length, const size_t width, const size_t height, const size_t radiusLength, dataType h) {

	if (imageDataPtr == NULL || gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL || gradientVectorR == NULL) {
		return false;
	}

	size_t i = 0, j = 0, k = 0, l = 0;
	dataType ux = 0.0, uy = 0.0, uz = 0.0, ur = 0.0;
	size_t index_ij = 0, index_kl = 0;

	for (l = 0; l < radiusLength; l++) {
		for (k = 0; k < height; k++) {
			for (i = 0; i < length; i++) {
				for (j = 0; j < width; j++) {

					index_kl = x_new(k, l, height);
					index_ij = x_new(i, j, length);

					if (l == 0) {
						ur = (imageDataPtr[x_new(k, l + 1, height)][index_ij] - imageDataPtr[x_new(k, l, height)][index_ij]) / h;
					}
					else {
						if (l == (radiusLength - 1)) {
							ur = (imageDataPtr[x_new(k, l, height)][index_ij] - imageDataPtr[x_new(k, l - 1, height)][index_ij]) / h;
						}
						else {
							ur = (imageDataPtr[x_new(k, l + 1, height)][index_ij] - imageDataPtr[x_new(k, l - 1, height)][index_ij]) / (2 * h);
						}
					}

					if (k == 0) {
						uz = (imageDataPtr[x_new(k + 1, l, height)][index_ij] - imageDataPtr[x_new(k, l, height)][index_ij]) / h;
					}
					else {
						if (k == (height - 1)) {
							uz = (imageDataPtr[x_new(k, l, height)][index_ij] - imageDataPtr[x_new(k - 1, l, height)][index_ij]) / h;
						}
						else {
							uz = (imageDataPtr[x_new(k + 1, l, height)][index_ij] - imageDataPtr[x_new(k - 1, l, height)][index_ij]) / (2 * h);
						}
					}

					if (i == 0) {
						ux = (imageDataPtr[index_kl][x_new(i + 1, j, length)] - imageDataPtr[index_kl][x_new(i, j, length)]) / h;
					}
					else {
						if (i == (length - 1)) {
							ux = (imageDataPtr[index_kl][x_new(i, j, length)] - imageDataPtr[index_kl][x_new(i - 1, j, length)]) / h;
						}
						else {
							ux = (imageDataPtr[index_kl][x_new(i + 1, j, length)] - imageDataPtr[index_kl][x_new(i - 1, j, length)]) / (2 * h);
						}
					}

					if (j == 0) {
						uy = (imageDataPtr[index_kl][x_new(i, j + 1, length)] - imageDataPtr[index_kl][x_new(i, j, length)]) / h;
					}
					else {
						if (j == (width - 1)) {
							uy = (imageDataPtr[index_kl][x_new(i, j, length)] - imageDataPtr[index_kl][x_new(i, j - 1, length)]) / h;
						}
						else {
							uy = (imageDataPtr[index_kl][x_new(i, j + 1, length)] - imageDataPtr[index_kl][x_new(i, j - 1, length)]) / (2 * h);
						}
					}

					gradientVectorX[index_kl][index_ij] = (dataType)ux;
					gradientVectorY[index_kl][index_ij] = (dataType)uy;
					gradientVectorZ[index_kl][index_ij] = (dataType)uz;
					gradientVectorR[index_kl][index_ij] = (dataType)ur;
				}
			}
		}
	}
	return true;
}

bool shortestPath4D(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, const size_t radiusLength, dataType raduisScale, dataType radiusInitial, dataType radiusFinal) {

	if (actionMapPtr == NULL || seedPoints == NULL)
		return false;	
	
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string saving_name, saving_csv = outputPath + "path4D_fm.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	size_t k, i, j, kl, ij, l, dim2D = length * width, max_iter = 100;
	double tau = 0.5, tol = 1.0;
	
	dataType i_current = seedPoints[1].x;
	dataType j_current = seedPoints[1].y;
	dataType k_current = seedPoints[1].z;
	dataType l_current = radiusFinal / raduisScale;

	dataType i_seed = seedPoints[0].x;
	dataType j_seed = seedPoints[0].y;
	dataType k_seed = seedPoints[0].z;
	dataType l_seed = radiusInitial / raduisScale;

	Point4D final_point = { i_seed, j_seed, k_seed, l_seed };

	const size_t newHeight = radiusLength * height;
	dataType** gradientVectorX = new dataType * [newHeight];
	dataType** gradientVectorY = new dataType * [newHeight];
	dataType** gradientVectorZ = new dataType * [newHeight];
	dataType** gradientVectorR = new dataType * [newHeight];
	for (k = 0; k < newHeight; k++) {
		gradientVectorX[k] = new dataType[dim2D]{0};
		gradientVectorY[k] = new dataType[dim2D]{0};
		gradientVectorZ[k] = new dataType[dim2D]{0};
		gradientVectorR[k] = new dataType[dim2D]{0};
		if (gradientVectorX[k] == NULL || gradientVectorY[k] == NULL || gradientVectorZ[k] == NULL || gradientVectorR[k] == NULL)
			return false;
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL || gradientVectorR == NULL)
		return false;

	//Normalization of the gradient
	compute4DimageGradient(actionMapPtr, gradientVectorX, gradientVectorY, gradientVectorZ, gradientVectorR, length, width, height, radiusLength, 1.0);

	dataType ux = 0.0, uy = 0.0, uz = 0.0, ur = 0.0, norm_of_gradient = 0.0;
	for (k = 0; k < newHeight; k++) {
		for (i = 0; i < dim2D; i++) {
			ux = gradientVectorX[k][i];
			uy = gradientVectorY[k][i];
			uz = gradientVectorZ[k][i];
			ur = gradientVectorR[k][i];
			norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz + ur * ur);
			gradientVectorX[k][i] = ux / norm_of_gradient;
			gradientVectorY[k][i] = uy / norm_of_gradient;
			gradientVectorZ[k][i] = uz / norm_of_gradient;
			gradientVectorR[k][i] = ur / norm_of_gradient;
		}
	}

	size_t count = 0;
	double dist_to_end = 0.0;

	fprintf(file,"x,y,z,radius\n");

	do {

		i = (size_t)i_current; //round(i_current);
		j = (size_t)j_current; //round(j_current);
		k = (size_t)k_current; //round(k_current);
		l = (size_t)l_current; //round(l_current);
		kl = x_new(k, l, height);
		ij = x_new(i, j, length);

		l_current = l_current - tau * gradientVectorR[kl][ij];
		k_current = k_current - tau * gradientVectorZ[kl][ij];
		i_current = i_current - tau * gradientVectorX[kl][ij];
		j_current = j_current - tau * gradientVectorY[kl][ij];

		fprintf(file,"%f,%f,%f,%f\n", i_current, j_current, k_current, l_current * raduisScale);
		
		Point4D current_point = { i_current, j_current, k_current, l_current };
		dist_to_end = getDistancePoints4D(current_point, final_point);
		count++;

	} while (dist_to_end > tol && count < max_iter);

	cout << "distance to end point : " << dist_to_end << endl;

	fclose(file);
	for (k = 0; k < newHeight; k++) {
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
		delete[] gradientVectorR[k];
	}
	delete[] gradientVectorX;
	delete[] gradientVectorY;
	delete[] gradientVectorZ;
	delete[] gradientVectorR;

	return true;
}

bool shortestPath3DPlus(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, const size_t radiusLength, dataType raduisScale, dataType radiusInitial, dataType radiusFinal) {

	if (actionMapPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j, k, l, dim2D = length * width;
	size_t kl, ij, max_iter = 100;
	double tau = 0.5, tol = 1.0;

	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	string extension, saving_name, saving_csv;

	dataType i_seed = seedPoints[0].x;
	dataType j_seed = seedPoints[0].y;
	dataType k_seed = seedPoints[0].z;
	Point3D final_point = { i_seed, j_seed, k_seed };

	dataType** actionField = new dataType * [height];
	dataType** gradientVectorX = new dataType * [height];
	dataType** gradientVectorY = new dataType * [height];
	dataType** gradientVectorZ = new dataType * [height];
	for (k = 0; k < height; k++) {
		actionField[k] = new dataType[dim2D]{ 0 };
		gradientVectorX[k] = new dataType[dim2D]{ 0 };
		gradientVectorY[k] = new dataType[dim2D]{ 0 };
		gradientVectorZ[k] = new dataType[dim2D]{ 0 };
		if (gradientVectorX[k] == NULL || gradientVectorY[k] == NULL || gradientVectorZ[k] == NULL)
			return false;
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;

	dataType ux = 0.0, uy = 0.0, uz = 0.0, norm_of_gradient = 0.0;

	size_t count = 0;
	double dist_to_end = 0.0;

	VoxelSpacing spacing = { 1.0, 1.0, 1.0 };

	for (l = 0; l < radiusLength; l++) {

		dataType i_current = seedPoints[1].x;
		dataType j_current = seedPoints[1].y;
		dataType k_current = seedPoints[1].z;

		extension = to_string(l);
		saving_csv = outputPath + "path_point_" + extension + ".csv";
		FILE* file;
		if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
			printf("Enable to open");
			return false;
		}

		//extract action
		for (k = 0; k < height; k++) {
			for (i = 0; i < dim2D; i++) {
				actionField[k][i] = actionMapPtr[x_new(k, l, height)][i];
			}
		}

		compute3dImageGradient(actionField, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, spacing);

		for (k = 0; k < height; k++) {
			for (i = 0; i < dim2D; i++) {
				ux = gradientVectorX[k][i];
				uy = gradientVectorY[k][i];
				uz = gradientVectorZ[k][i];
				norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz);
				gradientVectorX[k][i] = ux / norm_of_gradient;
				gradientVectorY[k][i] = uy / norm_of_gradient;
				gradientVectorZ[k][i] = uz / norm_of_gradient;
			}
		}
		
		dist_to_end = 0.0;
		count = 0;

		fprintf(file, "x,y,z\n");
		do {

			i = (size_t)i_current;
			j = (size_t)j_current;
			k = (size_t)k_current;

			ij = x_new(i, j, length);
			k_current = k_current - tau * gradientVectorZ[k][ij];
			i_current = i_current - tau * gradientVectorX[k][ij];
			j_current = j_current - tau * gradientVectorY[k][ij];

			fprintf(file, "%f,%f,%f\n", i_current, j_current, k_current);

			Point3D current_point = { i_current, j_current, k_current};
			dist_to_end = getPoint3DDistance(current_point, final_point);
			count++;

		} while (dist_to_end > tol && count < max_iter);
		fclose(file);
	}

	cout << "distance to end point : " << dist_to_end << endl;
	for (k = 0; k < height; k++) {
		delete[] actionField[k];
		delete[] gradientVectorX[k];
		delete[] gradientVectorY[k];
		delete[] gradientVectorZ[k];
	}
	delete[] actionField;
	delete[] gradientVectorX;
	delete[] gradientVectorY;
	delete[] gradientVectorZ;

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

void compute3DPotentialMeanVariance(Image_Data ctImageData, dataType** potentialPtr, Point3D seedPoint, dataType radiusInitial, dataType radiusMin, dataType radiusMax, dataType radiusStep) {

	size_t width = ctImageData.width;
	size_t length = ctImageData.length;
	size_t height = ctImageData.height;
	size_t i = 0, j = 0, k = 0;
	dataType lambda1 = 1.0, lambda2 = 1.0, omega = 0.00001;

	Statistics stats = getPointNeighborhoodStats(ctImageData, seedPoint, radiusInitial);
	dataType meanInitial = stats.mean_data / radiusInitial, mean = 0.0;
	dataType varianceInitial = stats.variance / radiusInitial, variance = 0.0;

	dataType radius = 0.0, value = 0.0, minValue = 0.0;
	dataType meanVal = 0.0, varianceVal = 0.0;
	
	dataType offset = 2.0;
	Statistics result = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	BoundingBox3D box = { 0, 0, 0, 0, 0, 0 };

	//===========================================

	/*
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				Point3D current_point = { i, j, k };
				box = findBoundingBox3D(current_point, length, width, height, radiusMax, offset);
				for (size_t kn = box.k_min; kn <= box.k_max; kn++) {
					for (size_t in = box.i_min; in <= box.i_max; in++) {
						for (size_t jn = box.j_min; jn <= box.j_max; jn++) {
							radius = radiusMin;
							while (radius <= radiusMax) {

							}
						}

					}
				}
			}
		}
	}

	Statistics result = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	BoundingBox3D box = { 0, 0, 0, 0, 0, 0 };
	dataType offset = 2.0;
	//box = findBoundingBox3D(point_of_interest, imageData.length, imageData.width, imageData.height, radius, offset);
	Point3D seed = getRealCoordFromImageCoord3D(point_of_interest, imageData.origin, imageData.spacing, imageData.orientation);

	dataType sum = 0.0;
	result.max_data = 0.0;
	result.min_data = 1000000000000000;
	int nb_point = 0;

	for (k = box.k_min; k <= box.k_max; k++) {
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, imageData.origin, imageData.spacing, imageData.orientation);
				double dist = getPoint3DDistance(seed, current_point);
				if (dist <= radius) {
					dataType voxel_value = imageData.imageDataPtr[k][x_new(i, j, length)];
					if (voxel_value < result.min_data) {
						result.min_data = voxel_value;
					}
					if (voxel_value > result.max_data) {
						result.max_data = voxel_value;
					}
					sum = sum + voxel_value;
					nb_point = nb_point + 1;
				}
			}
		}
	}

	if (nb_point != 0) {
		result.mean_data = sum / (dataType)nb_point;
		dataType sum_diff = 0.0;
		for (k = box.k_min; k <= box.k_max; k++) {
			for (i = box.i_min; i <= box.i_max; i++) {
				for (j = box.j_min; j <= box.j_max; j++) {
					Point3D current_point = { i, j, k };
					current_point = getRealCoordFromImageCoord3D(current_point, imageData.origin, imageData.spacing, imageData.orientation);
					double dist = getPoint3DDistance(seed, current_point);
					if (dist <= radius) {
						sum_diff = sum_diff + pow(imageData.imageDataPtr[k][x_new(i, j, length)] - result.mean_data, 2);
					}
				}
			}
		}
		result.variance = sum_diff / (dataType)nb_point;
		result.sd_data = sqrt(sum_diff / (dataType)nb_point);
	}
	else {
		result.mean_data = 0;
		result.sd_data = 0.0;
		result.max_data = 0.0;
		result.min_data = 0.0;
		result.variance = 0.0;
	}

	*/

	//============================================

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				radius = radiusMin;
				meanVal = 0.0;
				varianceVal = 0.0;
				Point3D current_point = { i, j, k };
				while (radius <= radiusMax) {
					stats = getPointNeighborhoodStats(ctImageData, current_point, radius);
					mean = stats.mean_data / radius;
					variance = stats.variance / radius;
					if (mean > meanVal) {
						meanVal = mean;
					}
					if (variance > varianceVal) {
						varianceVal = variance;
					}
					radius = radius + radiusStep;
				}
				value = omega + lambda1 * pow(meanVal - meanInitial, 2) + lambda2 * pow(varianceVal - varianceInitial, 2);
				potentialPtr[k][x_new(i, j, length)] = value;
			}
		}
	}

}