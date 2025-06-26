
#pragma once

#include<iostream>
#include<vector>

//#include<unordered_map>
//#include<tuple>
//#include<algorithm>

#include "common_functions.h"
#include "../src/data_load.h"
#include "../src/endianity_bl.h"
#include <stdio.h>
#include <string.h>
#include<../src/imageInterpolation.h>

using namespace std;

	typedef struct {
		size_t x, y;
		dataType arrival;
	} pointFastMarching2D;

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	} pointFastMarching3D;

	typedef struct {
		dataType K; //edge detection coef
		dataType thres;//edge detector threshold
		dataType eps; //path smothing parameter
		double radius;
	} Potential_Parameters;

	typedef struct {
		dataType tau; // used in the descent gradient
		size_t max_iteration; // maximal iteration to stop descent gradient
		double tolerance; //minimal distance to stop
	} Path_Parameters;

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
	/// solve quadratic equation : 	This fuction is used the solve the following quadratic coming the discretization 
	/// by upwind principle in the implementation of the fast marching method 
	/// aU^2 -2U(X+Y) + (X^2 + Y^2 - W) = 0
	/// </summary>
	/// <param name="X">coefficient related to the first variable</param>
	/// <param name="Y">coefficient related to the second variable</param>
	/// <param name="W">additional coefficient</param>
	/// <returns></returns>
	dataType solve2dQuadratic(dataType X, dataType Y, dataType P, PixelSpacing h);

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

	int getIndexFromHeap2D(vector<pointFastMarching2D>& in_Process, size_t i, size_t j);

	/// <summary>
	/// 
	/// </summary>
	/// <param name="imageData"></param>
	/// <param name="distancePtr"></param>
	/// <param name="potentialPtr"></param>
	/// <param name="endPoints"></param>
	/// <param name="savingPath"></param>
	/// <returns></returns>
	bool partialFrontPropagation2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* endPoints, string savingPath);

	bool doubleFrontPropagation2D(Image_Data2D imageData, dataType* actionFirstFront, dataType* actionSecondFront, dataType* potentialPtr, Point2D* endPoints, string savingPath);

	bool frontPropagationWithKeyPointDetection(Image_Data2D actionMapStr, dataType* potentialFuncPtr, Point2D* seedPoint, const double LengthKeyPoints, vector<Point2D>& key_points, std::string path_saving);

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
	bool fastMarching2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* seedPoints);
	
	/// <summary>
	/// Computes the potential function for a given 2D image using specified seed points and parameters.
	/// </summary>
	/// <param name="imageDataStr">The 2D image data to process.</param>
	/// <param name="potentialFuncPtr">Pointer to the output array where the computed potential function values will be stored.</param>
	/// <param name="seedPoints">Pointer to an array of seed points used in the potential computation.</param>
	/// <param name="parameters">The parameters controlling the potential computation.</param>
	/// <returns>True if the potential computation was successful; otherwise, false.</returns>
	bool computePotential(Image_Data2D imageDataStr, dataType* potentialFuncPtr, Point2D* seedPoints, Potential_Parameters parameters);

	/// <summary>
	/// Finds the shortest path in a 2D image using a action map and a set of seed points.
	/// </summary>
	/// <param name="distanceFuncPtr">A 2D image or function representing the distance or cost map.</param>
	/// <param name="seedPoints">Pointer to an array of seed points from which the path search begins.</param>
	/// <param name="path_points">Reference to a vector where the computed shortest path points will be stored.</param>
	/// <param name="parameters">Additional parameters controlling the pathfinding algorithm.</param>
	/// <returns>True if a valid shortest path is found; otherwise, false.</returns>
	bool shortestPath2d(Image_Data2D distanceFuncPtr, Point2D* seedPoints, vector<Point2D>& path_points, Path_Parameters parameters);

	bool bruteForceDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType foregroundValue);

	bool fastMarchingForDistanceMap(Image_Data2D ctImageData, dataType* distanceFuncPtr, dataType foregroundValue);

	bool fastSweepingDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType foregroundValue);

	bool rouyTourinDistanceMap2D(Image_Data2D ctImageData, dataType* distancePtr, dataType tolerance, size_t max_iteration, dataType foregroundValue);

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

	int getIndexFromHeap3D(vector<pointFastMarching3D>& in_Process, size_t i, size_t j, size_t k);

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
	bool compute3DPotential(Image_Data ctImageData, dataType** potential, Point3D* seedPoint, Potential_Parameters parameters);

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
	bool shortestPath3D(Image_Data actionMapStr, Point3D* seedPoints, vector<Point3D>& path_points, Path_Parameters parameters);

	//=================================

	/// <summary>
	/// Solves a 3D quadratic equation using the fast marching method.
	/// </summary>
	/// <param name="X">The value at the first neighboring voxel.</param>
	/// <param name="Y">The value at the second neighboring voxel.</param>
	/// <param name="Z">The value at the third neighboring voxel.</param>
	/// <param name="W">The right-hand side or source term for the equation.</param>
	/// <param name="h">The spacing between voxels in each dimension.</param>
	/// <returns>The computed solution for the 3D quadratic equation at the current voxel.</returns>
	dataType solve3dQuadraticEikonalEquation(dataType X, dataType Y, dataType Z, dataType P, VoxelSpacing h);

	/// <summary>
	/// Performs partial front propagation on image data and optionally saves the resulting path.
	/// </summary>
	/// <param name="actionPtr">The image data on which to perform the front propagation.</param>
	/// <param name="potentialFuncPtr">A pointer to a 2D array representing the potential function used during propagation.</param>
	/// <param name="endPoints">A pointer to an array of 3D points specifying the endpoints for the propagation.</param>
	/// <param name="path_saving">A string specifying the file path where the resulting path should be saved.</param>
	/// <returns>Returns true if the partial front propagation and path saving were successful; otherwise, returns false.</returns>
	bool partialFrontPropagation(Image_Data actionPtr, dataType** potentialFuncPtr, Point3D* endPoints, std::string path_saving);

	/// <summary>
	/// Performs the 3D Fast Marching Method on a CT image with specified voxel spacing.
	/// </summary>
	/// <param name="ctImageData">The input CT image data to process.</param>
	/// <param name="distanceFuncPtr">Pointer to the output distance function array.</param>
	/// <param name="potentialFuncPtr">Pointer to the output potential function array.</param>
	/// <param name="seedPoint">The starting point (seed) for the fast marching algorithm.</param>
	/// <param name="spacing">The spacing between voxels in the 3D image.</param>
	/// <returns>True if the fast marching computation was successful; otherwise, false.</returns>
	bool frontPropagation(Image_Data ctImageData, dataType** distanceFuncPtr, dataType** potentialFuncPtr, Point3D seedPoint);

	/// <summary>
	/// Performs front propagation with key point detection on an action map, starting from a seed point, and saves the resulting path points.
	/// </summary>
	/// <param name="actionMapStr">The action map data on which front propagation is performed.</param>
	/// <param name="potentialFuncPtr">Pointer to a 2D array representing the potential function used during propagation.</param>
	/// <param name="seedPoint">Pointer to the starting point (seed) for the propagation.</param>
	/// <param name="LengthKeyPoints">The length parameter used to determine key points along the path.</param>
	/// <param name="path_points">Reference to a vector where the detected path points will be stored.</param>
	/// <param name="path_saving">The file path where the resulting path will be saved.</param>
	/// <returns>True if the front propagation and key point detection succeed; otherwise, false.</returns>
	bool frontPropagationWithKeyPointDetection(Image_Data actionMapStr, dataType** potentialFuncPtr, Point3D* seedPoint, const double LengthKeyPoints, vector<Point3D>& key_points, std::string path_saving);

	bool doubleFrontPropagation(Image_Data imageData, dataType** actionFirstFront, dataType** actionSecondFront, dataType** potentialPtr, Point3D* endPoints, string savingPath);

	/// <summary>
	/// Returns the smaller of two values.
	/// </summary>
	/// <param name="x">The first value to compare.</param>
	/// <param name="y">The second value to compare.</param>
	/// <returns>The lesser of x and y.</returns>
	dataType min0(dataType x, dataType y);

	/// <summary>
	/// Computes a distance map using the Rouy-Tourin algorithm on the given image data.
	/// </summary>
	/// <param name="ctImageData">The input image data to process.</param>
	/// <param name="distancePtr">A pointer to a 2D array where the computed distance map will be stored.</param>
	/// <param name="foregroundValue">The value in the image that represents the foreground or object of interest.</param>
	/// <param name="tolerance">The convergence tolerance for the distance computation.</param>
	/// <param name="tau">The time step parameter for the Rouy-Tourin algorithm.</param>
	/// <returns>Returns true if the distance map was successfully computed; otherwise, returns false.</returns>
	bool rouyTourinDistanceMap(Image_Data ctImageData, dataType** distancePtr, dataType tolerance, size_t max_iteration, dataType foregroundValue);

	/// <summary>
	/// Performs the 3D Fast Marching Method to compute a distance map from a given image.
	/// </summary>
	/// <param name="ctImageData">The input image data to process.</param>
	/// <param name="distanceFuncPtr">A pointer to a 2D array where the computed distance map will be stored.</param>
	/// <param name="foregroundValue">The value in the image that represents the foreground or object of interest.</param>
	/// <returns>Returns true if the distance map was successfully computed; otherwise, returns false.</returns>
	bool fastMarchingDistanceMap(Image_Data ctImageData, dataType** distanceFuncPtr, dataType foregroundValue);

	/// <summary>
	/// Computes a distance map from the given image data using the fast sweeping method.
	/// </summary>
	/// <param name="ctImageData">The input image data to process.</param>
	/// <param name="distancePtr">A pointer to a 2D array where the computed distance map will be stored.</param>
	/// <param name="backgroundValue">The value in the image that represents the background.</param>
	/// <returns>True if the distance map was successfully computed; otherwise, false.</returns>
	bool fastSweepingDistanceMap(Image_Data ctImageData, dataType** distancePtr, const dataType backgroundValue);

	/// <summary>
	/// Computes a distance map for the given image using a brute-force approach.
	/// </summary>
	/// <param name="ctImageData">The image data to process.</param>
	/// <param name="distancePtr">A pointer to a 2D array where the computed distance map will be stored.</param>
	/// <param name="foregroundValue">The value in the image that represents the foreground.</param>
	/// <returns>True if the distance map was successfully computed; otherwise, false.</returns>
	bool bruteForceDistanceMap(Image_Data ctImageData, dataType** distancePtr, dataType foregroundValue);
