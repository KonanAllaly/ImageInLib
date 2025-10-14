#pragma once

#include<iostream>
#include<vector>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>

using namespace std;

	//==============================================================
	 
	//3D functions
	typedef struct {
		size_t x, y, z;
	} point3d;

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching3D;

	typedef struct {
		dataType K; //edge detection coef
		dataType thres;//edge detector threshold
		dataType eps; //path smothing parameter
		double radius;
	} Potential_Parameters;

	dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W);

	dataType upwindFiniteDifferenceX(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z);

	dataType upwindFiniteDifferenceY(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z);

	dataType upwindFiniteDifferenceZ(dataType** actionMapPtr, const size_t dimX, const size_t dimY, const size_t dimZ, const size_t x, const size_t y, const size_t z);

	bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints);

	//heap functions
	void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b);

	void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i);

	void heapifyVector3D(vector<pointFastMarching3D>& in_Process);

	void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i);

	void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process);

	void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point);

	int getIndexFromHeap3D(vector<pointFastMarching3D>& in_Process, size_t i, size_t j, size_t k);

	bool fastMarching3D_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, point3d* seedPoints);

	bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, point3d* seedPoints);

