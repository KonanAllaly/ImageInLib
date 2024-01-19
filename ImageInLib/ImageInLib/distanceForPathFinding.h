
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

	//2D functions

	typedef struct {
		size_t x, y;
	} point2D;

	typedef struct {
		size_t x, y;
		dataType arrival;
	}pointFastMarching2D;

	// heap functions

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
	bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, point2D* seedPoints);

	//fast marching functions

	/// <summary>
	/// solve quadratic equation
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
	dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	/// <summary>
	/// selection of neighboring pixels in y direction
	/// </summary>
	/// <param name="distanceFuncPtr">action map</param>
	/// <param name="dimI">image length</param>
	/// <param name="dimJ">image width</param>
	/// <param name="I">y coordinate</param>
	/// <param name="J">x coordinate</param>
	/// <returns></returns>
	dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	//bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width, dataType h);

	/// <summary>
	/// compute norm of 2D gradient
	/// </summary>
	/// <param name="gradientVectorX">gradient coordinate in x direction</param>
	/// <param name="gradientVectorY">gradient coordinate in y direction</param>
	/// <param name="height">image length</param>
	/// <param name="width">image width</param>
	/// <returns></returns>
	dataType computeGradientNorm2d(dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width);

	/// <summary>
	/// compute speed for each pixel for the fast marching
	/// </summary>
	/// <param name="imageDataPtr">input image data</param>
	/// <param name="potentialFuncPtr">speed map</param>
	/// <param name="height">image length</param>
	/// <param name="width">image width</param>
	/// <param name="seedPoints">starting point</param>
	/// <returns></returns>
	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, point2D* seedPoints);

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
	bool shortestPath2d(dataType* distanceFuncPtr, dataType* resultedPath, const size_t height, const size_t width, dataType h, point2D* seedPoints);

	//==============================================================

	typedef struct
	{
		dataType** minimum;
		dataType** maximum;
		dataType** mean;
		dataType** sd;
	} statictics_Pointers;

	typedef struct {
		dataType K, epsilon;
		dataType c_min, c_max, c_mean, c_sd, c_ct, c_pet;
	} Potential_Parameters;

	/// <summary>
	/// : This function create statistics images from statistics computed in each voxel
	/// in a small ball defined by the radius
	/// </summary>
	/// <param name="imageData"> : Structure to hold the input image</param>
	/// <param name="statsImage"> : Structure to hold the statistics images</param>
	/// <param name="radius"> : radius to be considered for the statistics</param>
	/// <returns></returns>
	bool generateStatisticsImages(Image_Data imageData, statictics_Pointers statsImage, double radius);

	/// <summary>
	/// compute potential function for fast marching
	/// </summary>
	/// <param name="ctImageData">structure for input CT image data</param>
	/// <param name="petImageData">structure for input PET image data</param>
	/// <param name="statsImage">point for images created from statistics</param>
	/// <param name="potential">point to hold computed potential</param>
	/// <param name="seedPoints">initial and final points</param>
	/// <param name="radius">radius to find statistics</param>
	/// <param name="params">potential function parameters</param>
	/// <returns></returns>
	bool computePotential_N(Image_Data ctImageData, Image_Data petImageData, statictics_Pointers statsImage, dataType** potential, Point3D* seedPoints, double radius, Potential_Parameters params);

	//================================================================

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching3D;

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
	dataType select3dX(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

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
	dataType select3dY(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

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
	dataType select3dZ(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	//double computeGradientNorm3d(dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height);

	bool compute3dImageGradient(dataType** imageDataPtr, dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t lenght, const size_t width, const size_t height, double h);

	bool compute3dPotential(dataType** imageDataPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D* seedPoints);

	//heap functions
	void swap3dPoints(pointFastMarching3D* a, pointFastMarching3D* b);

	void heapifyDown3D(vector<pointFastMarching3D>& in_Process, int i);

	void heapifyVector3D(vector<pointFastMarching3D>& in_Process);

	void heapifyUp3D(vector<pointFastMarching3D>& in_Process, int i);

	void deleteRootHeap3D(vector<pointFastMarching3D>& in_Process);

	void addPointHeap3D(vector<pointFastMarching3D>& in_Process, pointFastMarching3D point);

	//================================================================

	bool fastMarching3D_N(dataType** imageDataPtr, dataType** distanceFuncPtr, dataType** potentialFuncPtr, const size_t length, const size_t width, const size_t height, Point3D seedPoint);

	bool shortestPath3d(dataType** distanceFuncPtr, dataType** resultedPath, dataType** ctImageData, dataType** potentialPtr, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, string path_curvature, FILE* file_curve);

	//================================================================

	bool computePotentialNew(Image_Data ctImageData, dataType** meanImagePtr, dataType** potential, Point3D* seedPoints, double radius, Potential_Parameters parameters);

	bool findPathFromOneGivenPoint(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters);

	bool findPathFromOneGivenPointWithStopingCriterium(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters);

	bool shortestPath3dWithoutCurvature(dataType** distanceFuncPtr, dataType** resultedPath, const size_t length, const size_t width, const size_t height, dataType h, Point3D* seedPoints, vector<Point3D>& path_points);

	//=================================================================

	typedef struct {
		Point3D origin;
		VoxelSpacing spacing;
		size_t dimension[3];
	}imageMetaData;

	/// <summary>
	/// This function find automatically the bounding box of input binary image
	/// </summary>
	/// <param name="ctImageData"> : structure holding the image to be cropped</param>
	/// <param name="offset"> : offset</param>
	/// <returns> : the meta data of the image cropped image</returns>
	imageMetaData croppImage3D(Image_Data ctImageData, const size_t offset);

	//=================================================================
	bool findPathFromOneGivenPointWithCircleDetection(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters, size_t stop_criterium);

	//=================================================================
	bool findPath(Image_Data ctImageData, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters);

	bool computePotentialFromOnePoint(Image_Data ctImageData, dataType** potential, Point3D seed, double radius, Potential_Parameters parameters);

//#ifdef __cplusplus
//}
//#endif