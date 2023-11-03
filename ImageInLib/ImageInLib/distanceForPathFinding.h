
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

	void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b);

	void heapifyDown2D(vector<pointFastMarching2D>& in_Process, int pos);

	void heapifyUp2D(vector<pointFastMarching2D>& in_Process, int pos);

	void heapifyVector2D(vector<pointFastMarching2D>& in_Process);

	void deleteRootHeap2D(vector<pointFastMarching2D>& in_Process);

	void addPointHeap2D(vector<pointFastMarching2D>& in_Process, pointFastMarching2D point);

	bool fastMarching2D(dataType* imageDataPtr, dataType* distancePtr, dataType* potentialPtr, const size_t height, const size_t width, point2D* seedPoints);

	//fast marching functions

	dataType solve2dQuadratic(dataType X, dataType Y, dataType W);

	dataType selectX(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	dataType selectY(dataType* distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t I, const size_t J);

	bool computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width, dataType h);

	dataType computeGradientNorm2d(dataType* gradientVectorX, dataType* gradientVectorY, const size_t height, const size_t width);

	bool computePotential(dataType * imageDataPtr, dataType* potentialFuncPtr, const size_t height, const size_t width, point2D* seedPoints);

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

	bool computePotential_N(Image_Data ctImageData, Image_Data petImageData, statictics_Pointers statsImage, dataType** potential, Point3D* seedPoints, double radius, Potential_Parameters params);

	bool potentialOnEdgeImage(dataType** imageData, dataType** potential, Point3D* seeds, const size_t length, const size_t width, const size_t height);

	//================================================================

	typedef struct {
		size_t x, y, z;
		dataType arrival;
	}pointFastMarching3D;

	dataType solve3dQuadratic(dataType X, dataType Y, dataType Z, dataType W);

	dataType select3dX(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dY(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	dataType select3dZ(dataType** distanceFuncPtr, const size_t dimI, const size_t dimJ, const size_t dimK, const size_t I, const size_t J, const size_t K);

	double computeGradientNorm3d(dataType** gradientVectorX, dataType** gradientVectorY, dataType** gradientVectorZ, const size_t length, const size_t width, const size_t height);

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

	//bool findPathBetweenTwoGivenPoints(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters);

	bool findPathFromOneGivenPoint(Image_Data ctImageData, dataType** meanImagePtr, dataType** resultedPath, Point3D* seedPoints, Potential_Parameters parameters);

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

//#ifdef __cplusplus
//}
//#endif