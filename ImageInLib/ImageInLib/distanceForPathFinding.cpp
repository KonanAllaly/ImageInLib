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

#include "imageInterpolation.h"

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

	if (X == INFINITY && Y == INFINITY && Z == INFINITY) 
	{
		return INFINITY; // No solution if all coordinates are infinite
		std::cout << "Error: All coordinates are infinite." << std::endl;
	}

	if (X != INFINITY && Y == INFINITY && Z == INFINITY) 
	{
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

	if (Y != INFINITY && X == INFINITY && Z == INFINITY) 
	{
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

	if (Z != INFINITY && X == INFINITY && Y == INFINITY) 
	{
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

	if (X != INFINITY && Y != INFINITY && Z == INFINITY) 
	{
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

	if (X != INFINITY && Z != INFINITY && Y == INFINITY) 
	{
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

	if (Y != INFINITY && Z != INFINITY && X == INFINITY) 
	{
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

	if (X != INFINITY && Y != INFINITY && Z != INFINITY) 
	{
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

void heapifyDown3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex, int i) {

	int length_array = in_Process.size();
	if (length_array == 0) {
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
		size_t idx1 = in_Process[i].index;
		size_t idx2 = in_Process[current].index;

		swap_elts(&heapIndex[idx1], &heapIndex[idx2], sizeof(int));
		swap_elts(&in_Process[i], &in_Process[current], sizeof(pointFastMarching3D));

		heapifyDown3D(in_Process, heapIndex, current);
	}

}

void heapifyUp3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex, int i) {

	if (i <= 0 || i >= in_Process.size()) {
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
		size_t ind1 = in_Process[i].index;
		size_t ind2 = in_Process[current].index;

		swap_elts(&heapIndex[ind1], &heapIndex[ind2], sizeof(int));
		swap_elts(&in_Process[current], &in_Process[i], sizeof(pointFastMarching3D));

		heapifyUp3D(in_Process, heapIndex, current);
	}

}

void heapifyVector3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex) {
	int length_array = in_Process.size();
	if (length_array < 2) {
		return; //nothing to heapify
	}
	int ind, start = length_array / 2 - 1;
	for (ind = start; ind >= 0; ind--) {
		heapifyDown3D(in_Process, heapIndex, ind);
	}
}

void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process, vector<int>& heapIndex) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	if (l > 1) 
	{
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
	else 
	{
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
		addPointHeap3D(narrowBand, heapIndex, neighbor);
		action[ind_z][xd] = solution;
		labelArray[ind_z][xd] = 2;
	}
	else if (labelArray[ind_z][xd] == 2 && solution < action[ind_z][xd])
	{
		action[ind_z][xd] = solution;
		int pIndex = heapIndex[pos];
		if (pIndex != -1)
		{
			heapifyUp3D(narrowBand, heapIndex, pIndex);
		}
	}
}

bool fastMarching3D_N(Image_Data ctImageData, dataType** actionPtr, dataType** potentialFuncPtr, Point3D seedPoint) {

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

	const FiniteVolumeSize3D fVolume = { actionMapStr.spacing.sx, actionMapStr.spacing.sy, actionMapStr.spacing.sz };

	bool isGradientComputed = false;
	Point3D grad_vector;
	dataType norm_of_gradient = 0.0;
	do {

		isGradientComputed = getGradient3D(actionMapStr.imageDataPtr, length, width, height, i, j, k, fVolume, &grad_vector);
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
	if (labelArray == NULL || potentialFuncPtr == NULL) {
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
				if (ctImageData.imageDataPtr[k][currentIndx] == foregroundValue)
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
