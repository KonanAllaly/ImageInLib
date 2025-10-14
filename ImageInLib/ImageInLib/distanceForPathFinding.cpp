#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include<cmath>
#include <omp.h>
#include<vector>
#include<tuple>
#include "distanceForPathFinding.h"
#include<template_functions.h>

#define BIG_VALUE INFINITY

using namespace std;

//J.A Sethian, A Fast Marching Level Set method for Monotonically advancing fronts, 1995, page 8 and 10.
//link to article ---> http://ugweb.cs.ualberta.ca/~vis/courses/CompVis/readings/modelrec/sethian95fastlev.pdf

//Functions for 3D images

dataType upwindFiniteDifferenceX(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z) {

	dataType x_minus, x_plus;

	if (x == 0) 
	{
		x_minus = BIG_VALUE;
	}
	else 
	{
		x_minus = actionMapPtr[z][x_new(x - 1, y, dimX)];
	}

	if (x == dimX - 1) 
	{
		x_plus = BIG_VALUE;
	}
	else 
	{
		x_plus = actionMapPtr[z][x_new(x + 1, y, dimX)];
	}

	return min(x_minus, x_plus);
}

dataType upwindFiniteDifferenceY(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z) {

	dataType y_minus, y_plus;

	if (y == 0) 
	{
		y_minus = BIG_VALUE;
	}
	else 
	{
		y_minus = actionMapPtr[z][x_new(x, y - 1, dimX)];
	}

	if (y == dimY - 1) 
	{
		y_plus = BIG_VALUE;
	}
	else {
		y_plus = actionMapPtr[z][x_new(x, y + 1, dimX)];
	}

	return min(y_minus, y_plus);
}

dataType upwindFiniteDifferenceZ(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z) {

	dataType z_minus, z_plus;

	size_t xd = x_new(x, y, dimX);

	if (z == 0) 
	{
		z_minus = BIG_VALUE;
	}
	else 
	{
		z_minus = actionMapPtr[z - 1][xd];
	}

	if (z == dimZ - 1) 
	{
		z_plus = BIG_VALUE;
	}
	else 
	{
		z_plus = actionMapPtr[z + 1][xd];
	}

	return min(z_minus, z_plus);
}

// 3U^2 - 2U(X+Y+Z) + (X^2 + Y^2 + Z^2 - W) = 0 ---> aU + 2bU + c = 0
dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType P) {

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType P_2 = P * P;

	if (X == INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (Y + Z));
		c = (dataType)(pow(Y, 2) + pow(Z, 2) - P_2);
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(Y, Z) + P;
		}
	}

	if (Y == INFINITY && X != INFINITY && Z != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (X + Z));
		c = (dataType)(pow(X, 2) + pow(Z, 2) - P_2);
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, Z) + P;
		}
	}

	if (Z == INFINITY && X != INFINITY && Y != INFINITY) {
		a = 2;
		b = (dataType)(-2 * (X + Y));
		c = (dataType)(pow(X, 2) + pow(Y, 2) - P_2);
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, Y) + P;
		}
	}

	if (X == INFINITY && Y == INFINITY && Z != INFINITY) {
		a = 1;
		b = -2 * Z;
		c = pow(Z, 2) - P_2;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return Z + P;
		}
	}

	if (X == INFINITY && Z == INFINITY && Y != INFINITY) {
		a = 1;
		b = -2 * Y;
		c = pow(Y, 2) - P_2;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return Y + P;
		}
	}

	if (Y == INFINITY && Z == INFINITY && X != INFINITY) {
		a = 1;
		b = -2 * X;
		c = pow(X, 2) - P_2;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return X + P;
		}
	}

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) {
		a = 3;
		b = -2 * (X + Y + Z);
		c = pow(X, 2) + pow(Y, 2) + pow(Z, 2) - P_2;
		delta = b * b - 4 * a * c;
		if (delta >= 0) {
			return (-b + sqrt(delta)) / (2 * a);
		}
		else {
			return min(X, min(Y, Z)) + P;
		}
	}

}

bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints) {

	if (imageDataPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL)
		return false;

	size_t i, j, k;
	const size_t dim2D = length * width;
	size_t i0 = seedPoints[0].y, j0 = seedPoints[0].x, k0 = seedPoints[0].z;
	size_t i1 = seedPoints[1].y, j1 = seedPoints[1].x, k1 = seedPoints[1].z;

	dataType** gradientVectorX = new dataType* [height];
	dataType** gradientVectorY = new dataType* [height];
	dataType** gradientVectorZ = new dataType* [height];
	for (k = 0; k < height; k++) {
		gradientVectorX[k] = new dataType [dim2D];
		gradientVectorY[k] = new dataType [dim2D];
		gradientVectorZ[k] = new dataType [dim2D];
	}
	if (gradientVectorX == NULL || gradientVectorY == NULL || gradientVectorZ == NULL)
		return false;
	
	//compute3dImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, length, width, height, 1.0);

	size_t seedIndice = x_new(j0, i0, width), currentIndx = 0;
	dataType seedVal = (imageDataPtr[k0][x_new(j0, i0, width)] + imageDataPtr[k1][x_new(j1, i1, width)]) / 2;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;
	dataType epsilon = 0.01, K = 0.00005;

	//Computation of potential function
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				potentialFuncPtr[k][currentIndx] = abs(seedVal - imageDataPtr[k][currentIndx]);
			}
		}
	}

	//Find max
	dataType max_potential = -1 * INFINITY;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				if (potentialFuncPtr[k][currentIndx] > max_potential) {
					max_potential = potentialFuncPtr[k][currentIndx];
				}
			}
		}
	}

	//Normalization
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				ux = gradientVectorX[k][currentIndx];
				uy = gradientVectorY[k][currentIndx];
				uz = gradientVectorZ[k][currentIndx];
				potentialFuncPtr[k][currentIndx] = epsilon + (potentialFuncPtr[k][currentIndx] / max_potential) * (1 + K * (ux * ux + uy * uy + uz * uz));
			}
		}
	}

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

//function to swap 3D points in the fast marching contest
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
	if(l == 0) {
		return; //nothing to delete
	}
	else {
		if(l == 1) {
			in_Process.pop_back();
			return; //only one element
		}
		else {
			swap3dPoints(&in_Process[0], &in_Process[l - 1]);
			in_Process.pop_back();
			heapifyDown3D(in_Process, 0);
		}
	}
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

int getIndexFromHeap3D(vector<pointFastMarching3D>& in_Process, size_t i, size_t j, size_t k) {
	//we use type int for indexes because we do operations like pos--
	for (int ind = 0; ind < in_Process.size(); ind++) {
		if (in_Process[ind].x == i && in_Process[ind].y == j && in_Process[ind].z == k) {
			return ind;
		}
	}
	return -1; //not found
}

bool fastMarching3D_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints) {

	if (imageDataPtr == NULL || distanceFuncPtr == NULL || potentialFuncPtr == NULL || seedPoints == NULL) {
		return false;
	}

	vector <pointFastMarching3D> inProcess;
	size_t i = 0, j = 0, k = 0;

	size_t** labelArray = new size_t * [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new size_t[width * length];
	}
	if (labelArray == NULL)
		return false;

	dataType x = 0.0, y = 0.0, z = 0.0;
	size_t currentIndx = 0;


	//Initialization
	//All the points are notProcessed ---> label = 3
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				currentIndx = x_new(j, i, width);
				distanceFuncPtr[k][currentIndx] = INFINITY;
				labelArray[k][currentIndx] = 3;
			}
		}
	}

	/*
	compute3dPotential(imageDataPtr, potentialFuncPtr, length, width, height, seedPoints);

	//Processed the initial point
	pointFastMarching3D current;
	j = seedPoints[0].x;
	i = seedPoints[0].y;
	k = seedPoints[0].z;
	currentIndx = x_new(j, i, width);
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
		pointFastMarching3D TopNeighbor = { j, i, kminus, dTop };
		distanceFuncPtr[kminus][currentIndx] = dTop;
		inProcess.push_back(TopNeighbor);
		labelArray[kminus][currentIndx] = 2;
	}

	//Bottom
	if (k < height_minus && j >= 0 && j < width && i >= 0 && i < length) {
		kplus = k + 1;
		x = select3dX(distanceFuncPtr, length, width, height, i, j, kplus);
		y = select3dY(distanceFuncPtr, length, width, height, i, j, kplus);
		z = select3dZ(distanceFuncPtr, length, width, height, i, j, kplus);
		coefSpeed = potentialFuncPtr[kplus][currentIndx];
		dBottom = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D BottomNeighbor = { j, i, kplus, dBottom };
		distanceFuncPtr[kplus][currentIndx] = dBottom;
		inProcess.push_back(BottomNeighbor);
		labelArray[kplus][currentIndx] = 2;
	}

	//North
	if (k >= 0 && k < height && j > 0 && i >= 0 && i < length) {
		jminus = j - 1;
		indxNorth = x_new(jminus, i, width);
		x = select3dX(distanceFuncPtr, length, width, height, i, jminus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jminus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jminus, k);
		coefSpeed = potentialFuncPtr[k][indxNorth];
		dNorth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D NorthNeighbor = { jminus, i, k, dNorth };
		distanceFuncPtr[k][indxNorth] = dNorth;
		inProcess.push_back(NorthNeighbor);
		labelArray[k][indxNorth] = 2;
	}

	//South
	if (k >= 0 && k < height && j < width_minus && i >= 0 && i < length) {
		jplus = j + 1;
		indxSouth = x_new(jplus, i, width);
		x = select3dX(distanceFuncPtr, length, width, height, i, jplus, k);
		y = select3dY(distanceFuncPtr, length, width, height, i, jplus, k);
		z = select3dZ(distanceFuncPtr, length, width, height, i, jplus, k);
		coefSpeed = potentialFuncPtr[k][indxSouth];
		dSouth = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D SouthNeighbor = { jplus, i, k, dSouth };
		distanceFuncPtr[k][indxSouth] = dSouth;
		inProcess.push_back(SouthNeighbor);
		labelArray[k][indxSouth] = 2;
	}

	//East
	if (k >= 0 && k < height && j >= 0 && j < width && i < length_minus) {
		iplus = i + 1;
		indxEast = x_new(j, iplus, width);
		x = select3dX(distanceFuncPtr, length, width, height, iplus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iplus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iplus, j, k);
		coefSpeed = potentialFuncPtr[k][indxEast];
		dEast = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D EastNeighbor = { j, iplus, k, dEast };
		distanceFuncPtr[k][indxEast] = dEast;
		inProcess.push_back(EastNeighbor);
		labelArray[k][indxEast] = 2;
	}

	//West
	if (k >= 0 && k < height && j >= 0 && j < width && i > 0) {
		iminus = i - 1;
		indxWest = x_new(j, iminus, width);
		x = select3dX(distanceFuncPtr, length, width, height, iminus, j, k);
		y = select3dY(distanceFuncPtr, length, width, height, iminus, j, k);
		z = select3dZ(distanceFuncPtr, length, width, height, iminus, j, k);
		coefSpeed = potentialFuncPtr[k][indxWest];
		dWest = solve3dQuadratic(x, y, z, coefSpeed);
		pointFastMarching3D WestNeighbor = { j, iminus, k, dWest };
		distanceFuncPtr[k][indxWest] = dWest;
		inProcess.push_back(WestNeighbor);
		labelArray[k][indxWest] = 2;
	}

	//int n = inProcess.size();

	heapifyVector3D(inProcess);
	size_t label = 0;

	int l = 0, m = 0;

	while (inProcess.size() != 0) {

		//processed the minimum
		current = inProcess[0]; // index 0 exist because, if not we will be out of the while loop
		j = current.x;
		i = current.y;
		k = current.z;
		currentIndx = x_new(j, i, width);
		labelArray[k][currentIndx] = 1;
		distanceFuncPtr[k][currentIndx] = current.arrival;
		
		deleteRootHeap3D(inProcess);

		//Treat neighbors of the minimum in the narrow band
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
					pointFastMarching3D BottomNeighbor = { j, i, kplus, dBottom };
					addPointHeap3D(inProcess, BottomNeighbor);
				}
				else {
					if (label == 2) {
						if (dBottom < distanceFuncPtr[kplus][currentIndx]) {
							distanceFuncPtr[kplus][currentIndx] = dBottom;
							int index = getIndexFromHeap3D(inProcess, j, i, kplus);
							if (index != -1) {
								inProcess[index].arrival = dBottom;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}

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
					pointFastMarching3D TopNeighbor = { j, i, kminus, dTop };
					addPointHeap3D(inProcess, TopNeighbor);
				}
				else {
					if (label == 2) {
						if (dTop < distanceFuncPtr[kminus][currentIndx]) {
							distanceFuncPtr[kminus][currentIndx] = dTop;
							int index = getIndexFromHeap3D(inProcess, j, i, kminus);
							if (index != -1) {
								inProcess[index].arrival = dTop;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}

		//East
		if (i < length_minus && j >= 0 && j < width && k >= 0 && k < height) {
			iplus = i + 1;
			indxEast = x_new(j, iplus, width);
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
					pointFastMarching3D EastNeighbor = { j, iplus, k, dEast };
					addPointHeap3D(inProcess, EastNeighbor);
				}
				else {
					if (label == 2) {
						if (dEast < distanceFuncPtr[k][indxEast]) {
							distanceFuncPtr[k][indxEast] = dEast;
							int index = getIndexFromHeap3D(inProcess, j, iplus, k);
							if (index != -1) {
								inProcess[index].arrival = dEast;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}

		//West
		if (i > 0 && j >= 0 && j < width && k >= 0 && k < height) {
			iminus = i - 1;
			indxWest = x_new(j, iminus, width);
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
					pointFastMarching3D WestNeighbor = { j, iminus, k, dWest };
					addPointHeap3D(inProcess, WestNeighbor);
				}
				else {
					if (label == 2) {
						if (dWest < distanceFuncPtr[k][indxWest]) {
							distanceFuncPtr[k][indxWest] = dWest;
							int index = getIndexFromHeap3D(inProcess, j, iminus, k);
							if (index != -1) {
								inProcess[index].arrival = dWest;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}

		//North
		if (j > 0 && i >= 0 && i < length && k >= 0 && k < height) {
			jminus = j - 1;
			indxNorth = x_new(jminus, i, width);
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
					pointFastMarching3D NorthNeighbor = { jminus, i, k, dNorth };
					addPointHeap3D(inProcess, NorthNeighbor);
				}
				else {
					if (label == 2) {
						if (dNorth < distanceFuncPtr[k][indxNorth]) {
							distanceFuncPtr[k][indxNorth] = dNorth;
							int index = getIndexFromHeap3D(inProcess, jminus, i, k);
							if (index != -1) {
								inProcess[index].arrival = dNorth;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}

		//South
		if (j < width_minus && i >= 0 && i < length && k >= 0 && k < height) {
			jplus = j + 1;
			indxSouth = x_new(jplus, i, width);
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
					pointFastMarching3D SouthNeighbor = { jplus, i, k, dSouth };
					addPointHeap3D(inProcess, SouthNeighbor);
				}
				else {
					if (label == 2) {
						if (dSouth < distanceFuncPtr[k][indxSouth]) {
							distanceFuncPtr[k][indxSouth] = dSouth;
							int index = getIndexFromHeap3D(inProcess, jplus, i, k);
							if (index != -1) {
								inProcess[index].arrival = dSouth;
								heapifyUp3D(inProcess, index);
							}
						}
					}
				}
			}
		}
	}
	*/

	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
	}
	delete[] labelArray;

	return true;
}

bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, point3d* seedPoints) {

	if (distanceFuncPtr == NULL || resultedPath == NULL || seedPoints == NULL)
		return false;

	size_t dim2d = length * width, max_iter = 1000000000; //width* length* height;
	dataType tau = 0.8, tol = 1.0;
	size_t i_init = seedPoints[0].y, j_init = seedPoints[0].x, k_init = seedPoints[0].z;
	size_t i_end = seedPoints[1].y, j_end = seedPoints[1].x, k_end = seedPoints[1].z;
	
	//Find the closest point till the last point
	size_t cpt = 0;
	size_t i_current = i_end;
	size_t j_current = j_end;
	size_t k_current = k_end;
	size_t currentIndx = x_new(j_current, i_current, width);
	resultedPath[k_current][currentIndx] = 1;

	dataType iNew = i_current;
	dataType jNew = j_current;
	dataType kNew = k_current;
	dataType currentDist = 0.0;
	dataType dist_min = 0.0;

	const FiniteVolumeSize3D spacing = { 1.0, 1.0, 1.0 };
	Point3D v_gradient;
	dataType norm_of_gradient = 0;

	do {

		getGradient3D(distanceFuncPtr, length, width, height, i_current, j_current, k_current, spacing, &v_gradient);
		norm_of_gradient = sqrt(v_gradient.x * v_gradient.x + v_gradient.y * v_gradient.y + v_gradient.z * v_gradient.z);	

		currentIndx = x_new(j_current, i_current, width);
		iNew = iNew - tau * (v_gradient.x / norm_of_gradient);
		jNew = jNew - tau * (v_gradient.y / norm_of_gradient);
		kNew = kNew - tau * (v_gradient.z / norm_of_gradient);

		dist_min = sqrt((iNew - i_init) * (iNew - i_init) + (jNew - j_init) * (jNew - j_init) + (kNew - k_init) * (kNew - k_init));

		i_current = (size_t)(round(iNew)); 
		j_current = (size_t)(round(jNew));
		k_current = (size_t)(round(kNew));
		resultedPath[k_current][x_new(j_current, i_current, width)] = 1;
		cpt++;

	} while (dist_min > tol && cpt < max_iter);

	cout << "\nDistance to the end point : " << dist_min << endl;
	cout << "\nNumber of iterations : " << cpt << endl;

	return true;
}

/*
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

		if (x < 0.0 || x >= length || y < 0.0 || y >= width || z < 0.0 || z >= height)
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
*/
