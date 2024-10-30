
//#ifdef __cplusplus
//extern "C" {
//#endif

#pragma once

#include<iostream>
#include<vector>
#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>
#include<../src/imageInterpolation.h>

using namespace std;

	typedef struct {
		size_t x, y;
	} point2D;

	typedef struct {
		size_t x, y;
		dataType arrival;
	}pointFastMarching2D;

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching3D;

	typedef struct {
		dataType K, eps, h;
		double radius;
	} Potential_Parameters;

	// 2D functions

	/// <summary>
	/// swap two given points 
	/// </summary>
	/// <param name="a">first point</param>
	/// <param name="b">second point</param>
	void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b);

	/// <summary>
	/// heap operation from top to down
	/// </summary>
	/// <param name="in_Process">list of 2D point</param>
	/// <param name="pos">index of pixel 2D --> 1D</param>
	void heapifyDown2D(vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// heap operation from bottom to top
	/// </summary>
	/// <param name="in_Process">list of 2D point</param>
	/// <param name="pos">index of pixel 2D --> 1D</param>
	void heapifyUp2D(vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// create heap structure
	/// </summary>
	/// <param name="in_Process">list of 2D point</param>
	void heapifyVector2D(vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// remove the smallest element of the heap structure
	/// </summary>
	/// <param name="in_Process">list of 2D point</param>
	void deleteRootHeap2D(vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// add new element in the heap structure
	/// </summary>
	/// <param name="in_Process">list of 2D point</param>
	/// <param name="point"></param>
	void addPointHeap2D(vector<pointFastMarching2D>& in_Process, pointFastMarching2D point);

	/// <summary>
	/// compute arrival time for each point in given image, given starting point and given local speed
	/// </summary>
	/// <param name="imageDataPtr">input image</param>
	/// <param name="distancePtr">action map</param>
	/// <param name="potentialPtr">speed map</param>
	/// <param name="height">Height of the input image</param>
	/// <param name="width">Width of the input image</param>
	/// <param name="seedPoints">starting point</param>
	/// <returns></returns>
	bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, Point2D* seedPoints);

	bool partialFrontPropagation2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, Point2D* seedPoints);

	//fast marching functions

	/// <summary>
	/// solve quadratic equation : 	This fuction is used the solve the following quadratic coming the discretization 
	/// by upwind principle in the implementation of the fast marching method 
	/// aU^2 -2U(X+Y) + (X^2 + Y^2 - W) = 0
	/// </summary>
	/// <param name="X">coefficient related to the first variable</param>
	/// <param name="Y">coefficient related to the second variable</param>
	/// <param name="W">additional coefficient</param>
	/// <returns></returns>
	dataType solve2dQuadratic(dataType X, dataType Y, dataType W);

	/// <summary>
	/// selection of neighboring pixels in x direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="I">y coordinate</param>
	/// <param name="J">x coordinate</param>
	/// <returns></returns>
	dataType selectX(dataType* actionPtr, const size_t height, const size_t width, const size_t i, const size_t j);

	/// <summary>
	/// selection of neighboring pixels in y direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="I">y coordinate</param>
	/// <param name="J">x coordinate</param>
	/// <returns></returns>
	dataType selectY(dataType* actionPtr, const size_t height, const size_t width, const size_t i, const size_t j);

	/// <summary>
	/// compute speed for each pixel for the fast marching
	/// </summary>
	/// <param name="imageDataPtr">input image data</param>
	/// <param name="potentialFuncPtr">speed map</param>
	/// <param name="height">image length</param>
	/// <param name="width">image width</param>
	/// <param name="seedPoints">starting point</param>
	/// <returns></returns>
	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, Point2D* seedPoints);

	/// <summary>
	/// find the shortest path between two given points
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="resultedPath">found points</param>
	/// <param name="height">image length</param>
	/// <param name="width">image width</param>
	/// <param name="h">space step</param>
	/// <param name="seedPoints">initial and final points</param>
	/// <returns></returns>
	bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, Point2D* seedPoints);

	//================================================================

	//3D functions

	/// <summary>
	/// solve quadratic equation
	/// </summary>
	/// <param name="X">first quadratic coefficient</param>
	/// <param name="Y">second quadratic coefficient</param>
	/// <param name="Z">third quadratic coefficient</param>
	/// <param name="W">second member</param>
	/// <returns></returns>
	dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W);

	/// <summary>
	/// discretization in x-direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="dimK">image height</param>
	/// <param name="I">coordinate in y direction</param>
	/// <param name="J">coordinate in x direction</param>
	/// <param name="K">coordinate in z direction</param>
	/// <returns></returns>
	dataType select3dX(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k);

	/// <summary>
	/// discretization in y-direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="dimK">image height</param>
	/// <param name="I">coordinate in y direction</param>
	/// <param name="J">coordinate in x direction</param>
	/// <param name="K">coordinate in z direction</param>
	/// <returns></returns>
	dataType select3dY(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k);

	/// <summary>
	/// discretization in z-direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="dimK">image height</param>
	/// <param name="I">coordinate in y direction</param>
	/// <param name="J">coordinate in x direction</param>
	/// <param name="K">coordinate in z direction</param>
	/// <returns></returns>
	dataType select3dZ(dataType** actionPtr, const size_t length, const size_t width, const size_t height, const size_t i, const size_t j, const size_t k);

	/// <summary>
	/// Compute 3D image gradient
	/// </summary>
	/// <param name="imageDataPtr">input image</param>
	/// <param name="gradientVectorX">gradient X component</param>
	/// <param name="gradientVectorY">gradient Y component</param>
	/// <param name="gradientVectorZ">gradient Z component</param>
	/// <param name="lenght">image length</param>
	/// <param name="width">image width</param>
	/// <param name="height">image height</param>
	/// <param name="h">space discretization</param>
	/// <returns>return true after succes</returns>
	bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t lenght, const size_t width, const size_t height, VoxelSpacing spacing);

	/// <summary>
	/// swap 3D variables
	/// </summary>
	/// <param name="a">first variable</param>
	/// <param name="b">second variable</param>
	void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b);

	/// <summary>
	/// create heap sorting from last element 
	/// </summary>
	/// <param name="in_Process">vector containing all the elements</param>
	/// <param name="i">index where to start heapifying</param>
	void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i);

	/// <summary>
	/// create heap structure from vector
	/// </summary>
	/// <param name="in_Process">vector containing the elements</param>
	void heapifyVector3D(vector<pointFastMarching3D>& in_Process);

	/// <summary>
	/// heapify from the root element
	/// </summary>
	/// <param name="in_Process">vector of all the elements</param>
	/// <param name="i">index of element where to start the heapifying</param>
	void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i);

	/// <summary>
	/// delete root element from the heap structure
	/// </summary>
	/// <param name="in_Process">vector containing all the elements</param>
	void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process);

	/// <summary>
	/// add new point in the heap structure and heapifung
	/// </summary>
	/// <param name="in_Process">vector containing all the elements</param>
	/// <param name="point">element to be added</param>
	void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point);

	/// <summary>
	/// 3D fast marching (front propagation), starting by one seed point
	/// </summary>
	/// <param name="distanceFuncPtr">pointer to caintaining computed action field</param>
	/// <param name="potentialFuncPtr">pointer containing the potential (speed) in each point</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="height">image height</param>
	/// <param name="seedPoint">ssed point</param>
	/// <returns>return true after succes</returns>
	bool fastMarching3D_N(dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D seedPoint);

	/// <summary>
	/// compute speed (potential) in each point
	/// </summary>
	/// <param name="ctImageData">input image data</param>
	/// <param name="meanImagePtr">image generate by averaging around each point</param>
	/// <param name="potential">potential (speed) function</param>
	/// <param name="seedPoint">seed point</param>
	/// <param name="radius">radius for averaging</param>
	/// <param name="parameters">parameters for potential computation</param>
	/// <returns></returns>
	bool compute3DPotential(Image_Data ctImageData, dataType** potential, Point3D seedPoint, double radius, Potential_Parameters parameters);

	/// <summary>
	/// find the shortest path between two given points
	/// </summary>
	/// <param name="distanceFuncPtr">action field computed by fast marching from given seed point</param>
	/// <param name="resultedPath">shortest path between two given points</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="height">image height</param>
	/// <param name="h">space discretization</param>
	/// <param name="seedPoints">seed points</param>
	/// <param name="path_points">vector to keep all the path points coordinates</param>
	/// <returns></returns>
	bool shortestPath3D(dataType** distanceFuncPtr, const size_t length, const size_t width, const size_t height, VoxelSpacing spacing, Point3D* seedPoints, vector<Point3D>& path_points);

	/// <summary>
	/// find a path inside the aorta from given one point
	/// </summary>
	/// <param name="ctImageData">input image data</param>
	/// <param name="meanImagePtr">mean image obtained by averaging pixel value around each pixel</param>
	/// <param name="resultedPath">shortest path inside the aorta</param>
	/// <param name="seedPoints">seed points</param>
	/// <param name="parameters">parameters for potential function computation</param>
	/// <param name="stop_criterium">stopping criterium</param>
	/// <returns></returns>
	bool findPathTwoSteps(Image_Data ctImageData, Point3D* seedPoints, Potential_Parameters parameters);

	bool partialFrontPropagation(dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints);

//#ifdef __cplusplus
//}
//#endif

	//=============== 4D functions =========================
	
	typedef struct {
		size_t x, y, z, r;
		dataType arrival;
	}pointFastMarching4D;

	typedef struct {
		size_t x, y, z, r;
	}Point4D;

	double getDistancePoints4D(Point4D p1, Point4D p2);

	dataType select4DX(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rMax, const size_t i, const size_t j, const size_t k, const size_t r);

	dataType select4DY(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rMax, const size_t i, const size_t j, const size_t k, const size_t r);

	dataType select4DZ(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rMax, const size_t i, const size_t j, const size_t k, const size_t r);

	dataType select4DR(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, const size_t rMax, const size_t i, const size_t j, const size_t k, const size_t r);

	dataType solve4DQuadratic(dataType X, dataType Y, dataType Z, dataType R, dataType W);

	void swap4DPoints(pointFastMarching4D* a, pointFastMarching4D* b);

	void heapifyDown4D(vector<pointFastMarching4D>& in_Process, int i);

	void heapifyUp4D(vector<pointFastMarching4D>& in_Process, int i);

	void heapifyVector4D(vector<pointFastMarching4D>& in_Process);

	void deleteRootHeap4D(vector<pointFastMarching4D>& in_Process);

	void addPointHeap4D(vector<pointFastMarching4D>& in_Process, pointFastMarching4D point);

	void computePotential4D(Image_Data ctImageData, dataType** potentialPtr, Point3D seedPoint, const size_t radiusLength, dataType radiusInitial, dataType radiusScale);
	
	bool fastMarching4D(dataType** actionMapPtr, dataType** potentialPtr, const size_t length, const size_t width, const size_t height, const size_t radiusLength, Point3D seedPoint, dataType raduisScale, dataType radiusInitial);

	bool shortestPath4D(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, const size_t radiusLength, dataType raduisScale, dataType radiusInitial, dataType radiusFinal);

	bool computeDistanceToOnePoint(dataType** distancePtr, const size_t length, const size_t width, const size_t height, Point3D seed);

	bool shortestPath3DPlus(dataType** actionMapPtr, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, const size_t radiusLength, dataType raduisScale, dataType radiusInitial, dataType radiusFinal);

	//==================================
	void computePotentialMeanVariance(Image_Data2D ctImageData, dataType* potentialPtr, Point2D seedPoint, dataType radiusInitial, dataType radiusMax, dataType radiusStep);

	void compute3DPotentialMeanVariance(Image_Data ctImageData, dataType** potentialPtr, Point3D seedPoint, dataType radiusInitial, dataType radiusMin, dataType radiusMax, dataType radiusStep);